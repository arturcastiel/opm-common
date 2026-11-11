// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 NORCE.
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::PTFlash
 */
#ifndef OPM_CHI_FLASH_HPP
#define OPM_CHI_FLASH_HPP

#include <opm/material/constraintsolvers/FlashStatistics.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/constraintsolvers/RachfordRice.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/Constants.hpp>
#include <opm/material/eos/CubicEOS.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/classname.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace Opm {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components for the chiwoms problem.
 *
 */
template <class Scalar, class FluidSystem, bool isThermal = false>
class PTFlash
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx};
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx};
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx};
    static constexpr int numMiscibleComponents = FluidSystem::numMiscibleComponents;
    static constexpr int numMisciblePhases = FluidSystem::numMisciblePhases; //oil, gas
    static constexpr int numEq = numMisciblePhases + numMisciblePhases * numMiscibleComponents;

    using EOSType = CompositionalConfig::EOSType;

    // Tangent-plane distance below which a trial phase does not split off.
    // Close to a critical point an unstable feed lowers it by as little as
    // 1e-7 well inside the two-phase region.
    static constexpr Scalar stabilityTolerance = 1e-8;
    // Successive substitution has converged, or has collapsed onto the
    // trivial solution with all equilibrium ratios at one.
    static constexpr Scalar substitutionTolerance = 1e-10;
    static constexpr Scalar trivialSolutionTolerance = 1e-5;
    // Plain substitution steps of the two-phase split between two
    // extrapolations, see successiveSubstitutionComposition_().
    static constexpr int substitutionAccelerationInterval = 5;
    // The hand-over from successive substitution to Newton in the hybrid
    // method. Substitution continues while its observed convergence rate
    // predicts convergence within maxPredictedSubstitutionSteps more steps,
    // up to maxHybridSubstitutionSteps: a substitution step costs a fraction
    // of a Newton step. Otherwise Newton takes over after
    // substitutionStepsBeforeNewton steps.
    static constexpr int substitutionStepsBeforeNewton = 5;
    static constexpr int maxHybridSubstitutionSteps = 20;
    static constexpr Scalar maxPredictedSubstitutionSteps = 8;
    // Plain substitution steps of a stability trial between two
    // extrapolations, and the largest change of a ln K that one extrapolation
    // may make, see checkStability_().
    static constexpr int stabilityAccelerationInterval = 5;
    static constexpr Scalar maxStabilityExtrapolation = 1.0;
    // Substitution steps of a stability trial before second-order steps take
    // over, and the largest change of a ln W that one such step may make,
    // see stabilityNewtonStep_().
    static constexpr int stabilityNewtonStart = 3;
    static constexpr Scalar maxStabilityNewtonStep = 1.0;
    // The convergence test leaves the mole number sum S of a trial uncertain
    // by about the root of substitutionTolerance. A converged trial with S
    // that close to one takes up to maxResolvingSteps more second-order steps,
    // until its residual is below resolvedStabilityTolerance.
    static constexpr Scalar unresolvedStabilityBand = 1e-4;
    static constexpr Scalar resolvedStabilityTolerance = 1e-20;
    static constexpr int maxResolvingSteps = 4;
    // Substitution steps a cold start takes on the split before the
    // stability test, see flashFirst_(): without a phase below the tangent
    // plane of the feed, and in all.
    static constexpr int flashFirstIterations = 2;
    static constexpr int flashFirstMaxIterations = 3;
    // The clamps of CubicEOS::computeFugacityCoefficient(); a clamped value
    // does not give a tangent-plane distance.
    static constexpr Scalar maxFugacityCoefficient = 1e10;
    static constexpr Scalar minFugacityCoefficient = 1e-10;

public:
    /*!
     * \brief Calculates the fluid state from the global mole fractions of the components and the phase pressures
     *
     */
    template <class FluidState>
    static bool solve(FluidState& fluid_state,
                      PTFlashMethod twoPhaseMethod,
                      Scalar flash_tolerance,
                      const EOSType& eos_type,
                      int verbosity = 0)
    {
        using ScalarFluidState = CompositionalFluidState<Scalar, FluidSystem>;
        // Solve on a value-only copy; derivatives are reconstructed on the
        // caller's state after the scalar flash has converged.
        ScalarFluidState scalar_state;

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            scalar_state.setKvalue(compIdx, Opm::getValue(fluid_state.K(compIdx)));
            scalar_state.setMoleFraction(compIdx, Opm::getValue(fluid_state.moleFraction(compIdx)));
        }

        scalar_state.setLvalue(Opm::getValue(fluid_state.L()));
        scalar_state.setPressure(FluidSystem::oilPhaseIdx,
                                 Opm::getValue(fluid_state.pressure(FluidSystem::oilPhaseIdx)));
        scalar_state.setPressure(FluidSystem::gasPhaseIdx,
                                 Opm::getValue(fluid_state.pressure(FluidSystem::gasPhaseIdx)));

        scalar_state.setTemperature(Opm::getValue(fluid_state.temperature(0)));

        const auto [is_single_phase, newton_jacobian]
            = flashScalar_(scalar_state, twoPhaseMethod, flash_tolerance, eos_type, verbosity);

        // Transfer the converged values before reconstructing their derivatives.
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            const auto liquid_mole_fraction = scalar_state.moleFraction(oilPhaseIdx, compIdx);
            fluid_state.setMoleFraction(oilPhaseIdx, compIdx, liquid_mole_fraction);
            const auto vapour_mole_fraction = scalar_state.moleFraction(gasPhaseIdx, compIdx);
            fluid_state.setMoleFraction(gasPhaseIdx, compIdx, vapour_mole_fraction);
        }

        updateDerivatives_(scalar_state, fluid_state, eos_type, is_single_phase, newton_jacobian);

        // The equilibrium ratios of a two-phase answer are y / x: a warm
        // start of the next flash from this state starts from them.
        if (!is_single_phase) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                const auto& x = fluid_state.moleFraction(oilPhaseIdx, compIdx);
                if (Opm::getValue(x) > 0.) {
                    fluid_state.setKvalue(compIdx, fluid_state.moleFraction(gasPhaseIdx, compIdx) / x);
                }
            }
        }

        return is_single_phase;
    } //end solve

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * This is a convenience method which assumes that the capillary pressure is
     * zero...
     */
    template <class FluidState, class ComponentVector>
    static void solve(FluidState& fluid_state,
                      const ComponentVector& globalMolarities,
                      Scalar tolerance = 0.0)
    {
        using MaterialTraits = NullMaterialTraits<Scalar, numPhases>;
        using MaterialLaw = NullMaterial<MaterialTraits>;
        using MaterialLawParams = typename MaterialLaw::Params;

        MaterialLawParams matParams;
        solve<MaterialLaw>(fluid_state, matParams, globalMolarities, tolerance);
    }

    /*!
     * \brief The isothermal flash on a scalar fluid state, without derivatives.
     *
     * The stages are:
     * -# a stability test when the state does not already carry a two-phase
     *    split (\f$L \le 0\f$ or \f$L \ge 1\f$), which either declares the
     *    mixture single-phase or returns starting equilibrium ratios;
     * -# for a two-phase mixture, the Rachford-Rice equation for the
     *    starting liquid fraction, then the phase compositions from
     *    equality of fugacities by successive substitution and/or Newton;
     * -# a finite negative-flash liquid fraction is labelled single-phase;
     * -# a rejected split, or coincident phases from a warm start, triggers a
     *    new stability test from Wilson estimates and at most one retry; and
     * -# for a single-phase mixture, the phase label.
     *
     * The state also receives the fugacity coefficients of the phase
     * compositions it returns, which are the feed for a single phase.
     *
     * \return Whether the mixture is single-phase.
     */
    template <typename FluidState>
    static bool flash_solve_scalar(FluidState& fluid_state,
                                   PTFlashMethod twoPhaseMethod,
                                   const Scalar flash_tolerance,
                                   const EOSType& eos_type,
                                   const int verbosity = 0)
    {
        const bool is_single_phase
            = flashScalar_(fluid_state, twoPhaseMethod, flash_tolerance, eos_type, verbosity)
                  .single_phase;
        // The solve leaves those of its last evaluation, which need not be the result.
        using ParamCache =
            typename FluidSystem::template ParameterCache<typename FluidState::ValueType>;
        ParamCache param_cache(eos_type);
        for (unsigned phaseIdx = 0; phaseIdx < numMisciblePhases; ++phaseIdx) {
            param_cache.updatePhase(fluid_state, phaseIdx);
            updateFugacityCoefficients_(fluid_state, param_cache, phaseIdx);
        }
        return is_single_phase;
    }

protected:
    /*!
     * \brief Stores the fugacity coefficients of all components of a phase in
     *        the fluid state, from a parameter cache updated for the phase.
     *
     * A fluid system with fugacityCoefficients() evaluates the terms the
     * components share once; otherwise each component is evaluated alone.
     */
    template <class FlashFluidState, class ParamCache>
    static void updateFugacityCoefficients_(FlashFluidState& fluid_state,
                                            const ParamCache& param_cache,
                                            const unsigned phase_idx)
    {
        std::array<typename FlashFluidState::ValueType, numComponents> coefficients;
        if constexpr (requires {
                          FluidSystem::fugacityCoefficients(
                              fluid_state, param_cache, phase_idx, coefficients);
                      }) {
            FluidSystem::fugacityCoefficients(fluid_state, param_cache, phase_idx, coefficients);
        }
        else {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                coefficients[compIdx] = FluidSystem::fugacityCoefficient(
                    fluid_state, param_cache, phase_idx, compIdx);
            }
        }
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluid_state.setFugacityCoefficient(phase_idx, compIdx, coefficients[compIdx]);
        }
    }

    static constexpr std::size_t numFlashUnknowns = numMisciblePhases * numMiscibleComponents + 1;
    // x, y and L, then both phase pressures and the temperature.
    using JacobianState = std::array<Scalar, numFlashUnknowns + 3>;

    //! The composition Jacobian dF/d(x, y, L) where newtonComposition_()
    //! converged, and the state it was assembled at.
    struct NewtonJacobian
    {
        Dune::FieldMatrix<Scalar, numFlashUnknowns, numFlashUnknowns> jacobian;
        JacobianState state;
    };

    //! A scalar flash result, with the composition Jacobian when Newton
    //! finished the split.
    struct ScalarFlashResult
    {
        bool single_phase;
        std::optional<NewtonJacobian> newton_jacobian;
    };

    //! flash_solve_scalar(), also returning the composition Jacobian when
    //! Newton finished the split.
    template <typename FluidState>
    static ScalarFlashResult flashScalar_(FluidState& fluid_state,
                                          PTFlashMethod twoPhaseMethod,
                                          const Scalar flash_tolerance,
                                          const EOSType& eos_type,
                                          const int verbosity)
    {
        std::optional<NewtonJacobian> newton_jacobian;
        // The stable feed root of the last stability test.
        std::optional<FeedRoot> stable_feed_root;
        // States without a physical two-phase fraction require stability
        // analysis. Every completed split is classified below.
        bool is_single_phase = false;
        bool is_negative_flash = false;
        auto L_scalar = fluid_state.L();
        using ScalarVector = Dune::FieldVector<Scalar, numComponents>;
        ScalarVector K_scalar, z_scalar;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            K_scalar[compIdx] = fluid_state.K(compIdx);
            z_scalar[compIdx] = fluid_state.moleFraction(compIdx);
        }

        // Here, a cold start means there is no usable previous two-phase split
        // (L is outside (0, 1)); a warm start reuses previous K and L values.
        const bool is_cold_start = L_scalar <= 0 || L_scalar >= 1;
        // A cold start first tries the split itself; one that lowers the
        // Gibbs energy below the feed's makes the stability test unnecessary.
        const bool split_from_flash_first = is_cold_start
            && flashFirst_(K_scalar, L_scalar, newton_jacobian, fluid_state, z_scalar,
                           twoPhaseMethod, flash_tolerance, eos_type, verbosity);
        if (is_cold_start && !split_from_flash_first) {
            if (verbosity >= 1) {
                OpmLog::debug("Perform stability test (L <= 0 or L >= 1)!");
            }
            stable_feed_root = phaseStabilityTest_(
                is_single_phase, K_scalar, fluid_state, z_scalar, eos_type, verbosity);
        }
        if (verbosity >= 1) {
            OpmLog::debug(fmt::format("Inputs after stability test are K = [{}], L = [{}], z = [{}], P = {}, and T = {}",
                                     fmt::join(K_scalar, " "), L_scalar, fmt::join(z_scalar, " "),
                                     fluid_state.pressure(0), fluid_state.temperature(0)));
        }
        // Update the composition if cell is two-phase
        if (!is_single_phase && !split_from_flash_first) {
            // Only the first split of a warm start skips the stability test.
            bool is_warm_split = !is_cold_start;
            const auto try_solve_split = [&](const auto& on_failure) {
                newton_jacobian.reset();
                try {
                    OPM_FLASH_COUNT(rachfordRiceSolves);
                    L_scalar = RachfordRice::solve(K_scalar, z_scalar, L_scalar, verbosity);
                    // Stored ratios whose split lies outside (0, 1) predict a
                    // single phase from an earlier state. Substitution from
                    // them can settle on a negative-flash fixed point of the
                    // new state that is not its answer, so only the stability
                    // test decides.
                    if (is_warm_split && !(L_scalar > 0. && L_scalar < 1.)) {
                        OPM_THROW_NOLOG(NumericalProblem,
                                        fmt::format("The stored equilibrium ratios imply no "
                                                    "split (L = {})", L_scalar));
                    }
                    newton_jacobian = flash_2ph(z_scalar, twoPhaseMethod, K_scalar, L_scalar,
                                                fluid_state, flash_tolerance, eos_type, verbosity);
                    return true;
                }
                catch (const NumericalProblem& error) {
                    on_failure(error);
                }
                catch (const Dune::FMatrixError& error) {
                    on_failure(error);
                }
                return false;
            };
            const auto classify_negative_flash = [&]() {
                const Scalar liquid_fraction = Opm::getValue(L_scalar);
                is_negative_flash = liquid_fraction <= 0. || liquid_fraction >= 1.;
                is_single_phase = is_negative_flash;
                if (is_negative_flash) {
                    L_scalar = liquid_fraction <= 0. ? 0. : 1.;
                }
            };

            // Rachford Rice equation to get initial L for composition solver
            const bool split_succeeded = try_solve_split([&](const auto& error) {
                if (verbosity >= 1) {
                    OpmLog::debug(fmt::format(
                        "The two-phase flash was rejected ({}); reassessing phase stability.",
                        error.what()));
                }
            });

            // A finite negative-flash root outside [0, 1] is the single-phase
            // criterion, not a failed composition solve.
            if (split_succeeded) {
                classify_negative_flash();
            }

            // A warm start must be reassessed from Wilson estimates because its
            // stored K may be a trivial fixed point. For a cold start, a
            // coincident split contradicts the preceding instability result.
            // Coincident phases solve the equations for any L, so a warm
            // start that ends on them is reassessed whatever its L, also one
            // the negative-flash test would call single-phase.
            const bool phases_coincide = split_succeeded && (!is_single_phase || !is_cold_start)
                && phasesCoincide_(fluid_state);
            if (phases_coincide && is_cold_start) {
                OPM_THROW_NOLOG(NumericalProblem,
                                "The two-phase flash converged to coincident phases after "
                                "stability analysis identified an unstable mixture");
            }
            if (!split_succeeded || (phases_coincide && !is_cold_start)) {
                if (phases_coincide && verbosity >= 1) {
                    OpmLog::debug("The two-phase flash converged to coincident phases; "
                                  "reassessing phase stability.");
                }

                is_single_phase = false;
                is_negative_flash = false;
                // A uniform K vector also represents coincident normalised
                // compositions and can be a trivial fixed point of the
                // stability trials. Restart from Wilson estimates to obtain an
                // independent test.
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    K_scalar[compIdx] = fluid_state.wilsonK_(compIdx);
                }
                stable_feed_root = phaseStabilityTest_(
                    is_single_phase, K_scalar, fluid_state, z_scalar, eos_type, verbosity);
                if (!is_single_phase) {
                    OPM_FLASH_COUNT(splitRetries);
                    is_warm_split = false;
                    // Stability found a distinct trial phase and returned a
                    // fresh K estimate. Retry the split once from that state.
                    try_solve_split([](const auto& error) {
                        OPM_THROW_NOLOG(
                            NumericalProblem,
                            fmt::format("Two-phase flash retry failed after stability "
                                        "analysis: {}",
                                        error.what()));
                    });
                    classify_negative_flash();
                    if (!is_single_phase && phasesCoincide_(fluid_state)) {
                        OPM_THROW_NOLOG(NumericalProblem,
                                        "The two-phase flash retry converged to coincident phases "
                                        "after stability analysis identified an unstable mixture");
                    }
                }
            }
        }
        if (is_single_phase) {
            newton_jacobian.reset();
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                fluid_state.setMoleFraction(gasPhaseIdx, compIdx, z_scalar[compIdx]);
                fluid_state.setMoleFraction(oilPhaseIdx, compIdx, z_scalar[compIdx]);
            }
            if (!is_negative_flash) {
                // A feed whose cubic has a liquid and a vapour root is labelled
                // by the root with the lower Gibbs energy. A single root larger
                // than the pseudo-critical volume of the feed is a vapour; a
                // smaller one is labelled by Li's correlation.
                // The stability test that found the feed stable has its root.
                const auto [stable_phase_idx, has_two_roots, molar_volume] = stable_feed_root
                    ? *stable_feed_root : feedStableRoot_(fluid_state, z_scalar, eos_type);
                if (has_two_roots) {
                    L_scalar = stable_phase_idx == static_cast<int>(oilPhaseIdx) ? 1.0 : 0.0;
                    if (verbosity >= 1) {
                        OpmLog::debug(fmt::format("Cell is single-phase, labelled by its stable "
                                                  "cubic root (L = {})", L_scalar));
                    }
                }
                else {
                    const Scalar critical_volume = pseudoCriticalVolume_(z_scalar);
                    if (molar_volume > critical_volume) {
                        L_scalar = 0.0;
                        if (verbosity >= 1) {
                            OpmLog::debug(fmt::format("Cell is single-phase, vapor (L = 0.0): its cubic "
                                                      "root is larger than the pseudo-critical volume "
                                                      "({} > {} m3/mol)", molar_volume, critical_volume));
                        }
                    }
                    else {
                        L_scalar = li_single_phase_label_(fluid_state, z_scalar, verbosity);
                    }
                }
            }
        }
        fluid_state.setLvalue(L_scalar);
        return {is_single_phase, newton_jacobian};
    }

    //! The stable cubic root of a feed, see feedStableRoot_().
    struct FeedRoot
    {
        int phase_idx; //!< phase index of the stable root
        bool has_two_roots; //!< whether the cubic has a separate liquid and vapour root
        Scalar molar_volume; //!< molar volume of the stable root [m^3/mol]
    };

    /*!
     * \brief The cubic root of the feed with the lower Gibbs energy.
     *
     * The feed composition is evaluated on the root of each phase label: the
     * smallest root for oil, the largest for gas. Where the cubic has a
     * physical liquid and vapour root, the one with the lower
     * \f$\sum_i z_i \ln\varphi_i\f$ is the stable state of the feed. Where
     * the smallest root is negative, the liquid volume is floored, and the
     * fugacity coefficients on it are clamped to their maximum, so it is never
     * the stable one. With a single root both labels give the same root.
     *
     * \return the phase index of the stable root, whether the cubic has a
     *         separate liquid and vapour root, and the molar volume of the
     *         stable root
     */
    template <class FlashFluidState, class ComponentVector>
    static FeedRoot feedStableRoot_(const FlashFluidState& fluid_state,
                                    const ComponentVector& z,
                                    const EOSType& eos_type)
    {
        using FlashEval = typename FlashFluidState::ValueType;
        using ParamCache = typename FluidSystem::template ParameterCache<FlashEval>;
        ParamCache feed_cache(eos_type);
        std::array<FlashEval, numComponents> stable_coefficients;
        return feedStableRoot_(fluid_state, z, feed_cache, stable_coefficients);
    }

    /*!
     * \brief feedStableRoot_() on a given parameter cache, also returning the
     *        fugacity coefficients of the feed on its stable root, the
     *        reference of the stability trials.
     *
     * The cache is left updated for the feed on both roots; a later update
     * for another composition at the same temperature and pressure keeps
     * its pure-component parameters.
     */
    template <class FlashFluidState, class ComponentVector, class ParamCache>
    static FeedRoot feedStableRoot_(const FlashFluidState& fluid_state,
                                    const ComponentVector& z,
                                    ParamCache& feed_cache,
                                    std::array<typename FlashFluidState::ValueType, numComponents>&
                                        stable_coefficients)
    {
        using FlashEval = typename FlashFluidState::ValueType;
        using CubicEOS = typename Opm::CubicEOS<Scalar, FluidSystem>;

        constexpr int oil = static_cast<int>(oilPhaseIdx);
        constexpr int gas = static_cast<int>(gasPhaseIdx);
        FlashFluidState feed_state = fluid_state;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            feed_state.setMoleFraction(oil, compIdx, z[compIdx]);
            feed_state.setMoleFraction(gas, compIdx, z[compIdx]);
        }
        // One cache holds both phases, which share the pure-component parameters.
        feed_cache.updatePhase(feed_state, oil);
        feed_cache.updatePhase(feed_state, gas);

        // A separate liquid root is far smaller than the vapour root; a single
        // root gives two equal volumes up to rounding.
        const Scalar liquid_volume = Opm::getValue(feed_cache.molarVolume(oil));
        const Scalar vapour_volume = Opm::getValue(feed_cache.molarVolume(gas));
        if (!(liquid_volume < (1.0 - 1e-6) * vapour_volume)) {
            CubicEOS::computeFugacityCoefficients(feed_state, feed_cache, oil, stable_coefficients);
            return {oil, false, vapour_volume};
        }

        const auto sum_z_ln_phi = [&](const std::array<FlashEval, numComponents>& phis) {
            Scalar sum = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (!(Opm::getValue(z[compIdx]) > 0.)) {
                    continue;
                }
                const auto& phi = phis[compIdx];
                sum += Opm::getValue(z[compIdx]) * std::log(Opm::getValue(phi));
            }
            return sum;
        };
        std::array<FlashEval, numComponents> liquid_coefficients;
        CubicEOS::computeFugacityCoefficients(feed_state, feed_cache, gas, stable_coefficients);
        CubicEOS::computeFugacityCoefficients(feed_state, feed_cache, oil, liquid_coefficients);
        const bool vapour_is_stable
            = sum_z_ln_phi(stable_coefficients) < sum_z_ln_phi(liquid_coefficients);
        if (!vapour_is_stable) {
            stable_coefficients = liquid_coefficients;
        }
        return {vapour_is_stable ? gas : oil, true, vapour_is_stable ? vapour_volume : liquid_volume};
    }

    /*!
     * \brief The pseudo-critical molar volume of the feed by Kay's rule,
     *        \f$\sum_i z_i V_{c,i}\f$ [m^3/mol].
     *
     * A single cubic root larger than it is less dense than the critical
     * state of the feed: a vapour, also below Li's temperature.
     */
    template <class ComponentVector>
    static Scalar pseudoCriticalVolume_(const ComponentVector& z)
    {
        Scalar critical_volume = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // The fluid systems give the critical volume in m^3/kmol.
            critical_volume += Opm::getValue(z[compIdx])
                * FluidSystem::criticalVolume(compIdx) / prefix::kilo;
        }
        return critical_volume;
    }

    template <class Vector, class FlashFluidState>
    static typename Vector::field_type li_single_phase_label_(const FlashFluidState& fluid_state, const Vector& z, int verbosity)
    {
        // Calculate intermediate sum
        typename Vector::field_type sumVz = 0.0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            // Get component information
            const auto& V_crit = FluidSystem::criticalVolume(compIdx);

            // Sum calculation
            sumVz += (V_crit * z[compIdx]);
        }

        // Calculate approximate (pseudo) critical temperature using Li's method
        typename Vector::field_type Tc_est = 0.0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            // Get component information
            const auto& V_crit = FluidSystem::criticalVolume(compIdx);
            const auto& T_crit = FluidSystem::criticalTemperature(compIdx);

            // Sum calculation
            Tc_est += (V_crit * T_crit * z[compIdx] / sumVz);
        }

        // Get temperature
        const auto& T = fluid_state.temperature(0);

        // If temperature is below estimated critical temperature --> phase = liquid; else vapor
        typename Vector::field_type L;
        if (T < Tc_est) {
            // Liquid
            L = 1.0;

            // Print
            if (verbosity >= 1) {
                OpmLog::debug(fmt::format("Cell is single-phase, liquid (L = 1.0) due to Li's phase labeling method giving T < Tc_est ({} < {})!",
                                         T, Tc_est));
            }
        }
        else {
            // Vapor
            L = 0.0;

            // Print
            if (verbosity >= 1) {
                OpmLog::debug(fmt::format("Cell is single-phase, vapor (L = 0.0) due to Li's phase labeling method giving T >= Tc_est ({} >= {})!",
                                         T, Tc_est));
            }
        }

        return L;
    }

    /*!
     * \brief The tangent plane of the Gibbs energy at the feed,
     *        \f$d_i = \ln z_i + \ln\varphi_i(\mathbf z)\f$, on the stable cubic
     *        root of the feed (the lower \f$\sum_i z_i d_i\f$).
     *
     * \return nothing where a fugacity coefficient is at its clamp
     */
    template <class FlashFluidState, class ComponentVector>
    static std::optional<Dune::FieldVector<Scalar, numComponents>>
    feedTangentPlane_(const FlashFluidState& fluid_state,
                      const ComponentVector& z,
                      const EOSType& eos_type)
    {
        using ParamCache =
            typename FluidSystem::template ParameterCache<typename FlashFluidState::ValueType>;
        using PlaneVector = Dune::FieldVector<Scalar, numComponents>;
        FlashFluidState feed_state = fluid_state;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            feed_state.setMoleFraction(oilPhaseIdx, compIdx, z[compIdx]);
            feed_state.setMoleFraction(gasPhaseIdx, compIdx, z[compIdx]);
        }
        std::optional<PlaneVector> plane;
        Scalar plane_gibbs = 0.0;
        // One cache for both roots: temperature and pressure are fixed.
        ParamCache param_cache(eos_type);
        for (unsigned phaseIdx = 0; phaseIdx < numMisciblePhases; ++phaseIdx) {
            param_cache.updatePhase(feed_state, phaseIdx);
            updateFugacityCoefficients_(feed_state, param_cache, phaseIdx);
            PlaneVector phase_plane(0.0);
            Scalar phase_gibbs = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                const Scalar z_i = Opm::getValue(z[compIdx]);
                if (!(z_i > 0.)) {
                    continue;
                }
                const Scalar phi = Opm::getValue(feed_state.fugacityCoefficient(phaseIdx, compIdx));
                if (!(phi > minFugacityCoefficient && phi < maxFugacityCoefficient)) {
                    return std::nullopt;
                }
                phase_plane[compIdx] = std::log(z_i * phi);
                phase_gibbs += z_i * phase_plane[compIdx];
            }
            if (!plane || phase_gibbs < plane_gibbs) {
                plane = phase_plane;
                plane_gibbs = phase_gibbs;
            }
        }
        return plane;
    }

    /*!
     * \brief Flash first, stability on demand.
     *
     * A cold start without a phase state (the incoming \f$L\f$ is neither 0
     * nor 1) first takes substitution steps of the two-phase split from the
     * incoming ratios, the order Michelsen recommends and thermopack
     * (tp_solver.f90, twoPhaseTPflash) and neqsim (TPflash.java) use. After
     * each step both phases are measured against the tangent plane
     * \f$\mathbf d\f$ of the feed,
     *
     * \f[ \mathrm{tpd}(\mathbf w) = \sum_i w_i \left(\ln w_i + \ln\varphi_i(\mathbf w) - d_i\right). \f]
     *
     * The Gibbs energy of the split exceeds that of the feed by
     * \f$L\,\mathrm{tpd}(\mathbf x) + (1 - L)\,\mathrm{tpd}(\mathbf y)\f$.
     * When this is below \f$-\f$stabilityTolerance, a phase lies below the
     * tangent plane: the feed is unstable, which is all the stability test
     * could establish, so the test is skipped and the split, already below
     * the feed, is finished by \p twoPhaseMethod.
     *
     * Each step needs an interior Rachford-Rice root and fugacity
     * coefficients off their clamps. After flashFirstIterations steps the
     * attempt continues, up to flashFirstMaxIterations, only while a phase is
     * below the plane. A split whose phase is below the plane but whose
     * Gibbs energy is not yet below the feed is not finished from here:
     * Newton from such a start can take hundreds of steps where the
     * stability test's ratios converge in a few.
     *
     * Without acceptance, or when the finished split is not a distinct
     * interior one, the fluid state is restored and the caller runs the
     * stability test from the incoming ratios, as without this step. An
     * incoming \f$L\f$ of 0 or 1 is a single-phase result of the previous
     * call; such a state is most likely still single-phase and goes to the
     * stability test directly.
     *
     * \return whether the split was accepted; K, L, the fluid state and the
     *         Newton Jacobian then hold it
     */
    template <class FlashFluidState, class ComponentVector>
    static bool flashFirst_(ComponentVector& K,
                            typename FlashFluidState::ValueType& L,
                            std::optional<NewtonJacobian>& newton_jacobian,
                            FlashFluidState& fluid_state,
                            const ComponentVector& z,
                            const PTFlashMethod twoPhaseMethod,
                            const Scalar flash_tolerance,
                            const EOSType& eos_type,
                            const int verbosity)
    {
        using ValueType = typename FlashFluidState::ValueType;
        using ParamCache = typename FluidSystem::template ParameterCache<ValueType>;

        if (L == 0. || L == 1.) {
            return false;
        }
        const FlashFluidState incoming_state = fluid_state;
        ComponentVector trial_K = K;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!(Opm::getValue(z[compIdx]) > 0.)) {
                trial_K[compIdx] = 1.;
            }
        }
        ValueType trial_L = 0.;
        std::optional<Dune::FieldVector<Scalar, numComponents>> feed_plane;
        bool feed_is_unstable = false;
        // One cache for all trial steps: temperature and pressure are fixed.
        ParamCache param_cache(eos_type);
        for (int iteration = 0; iteration < flashFirstMaxIterations; ++iteration) {
            try {
                OPM_FLASH_COUNT(rachfordRiceSolves);
                trial_L = RachfordRice::solve(trial_K, z, 0);
            }
            catch (const NumericalProblem&) {
                break;
            }
            if (!(trial_L > 0. && trial_L < 1.)) {
                break;
            }
            OPM_FLASH_COUNT(ssiIterations);
            computeLiquidVapor_(fluid_state, trial_L, trial_K, z);
            bool is_clamped = false;
            for (unsigned phaseIdx = 0; phaseIdx < numMisciblePhases; ++phaseIdx) {
                param_cache.updatePhase(fluid_state, phaseIdx);
                updateFugacityCoefficients_(fluid_state, param_cache, phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = Opm::getValue(fluid_state.fugacityCoefficient(phaseIdx, compIdx));
                    is_clamped = is_clamped
                        || !(phi > minFugacityCoefficient && phi < maxFugacityCoefficient);
                }
            }
            if (is_clamped) {
                break;
            }
            if (!feed_plane) {
                feed_plane = feedTangentPlane_(incoming_state, z, eos_type);
                if (!feed_plane) {
                    break;
                }
            }
            const auto tangent_plane_distance = [&](const unsigned phaseIdx) {
                Scalar distance = 0.0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar w = Opm::getValue(fluid_state.moleFraction(phaseIdx, compIdx));
                    if (Opm::getValue(z[compIdx]) > 0. && w > 0.) {
                        const Scalar phi
                            = Opm::getValue(fluid_state.fugacityCoefficient(phaseIdx, compIdx));
                        distance += w * (std::log(w * phi) - (*feed_plane)[compIdx]);
                    }
                }
                return distance;
            };
            const Scalar liquid_distance = tangent_plane_distance(oilPhaseIdx);
            const Scalar vapour_distance = tangent_plane_distance(gasPhaseIdx);
            const Scalar liquid_fraction = Opm::getValue(trial_L);
            if (liquid_fraction * liquid_distance + (1. - liquid_fraction) * vapour_distance
                < -stabilityTolerance) {
                feed_is_unstable = true;
                break;
            }
            const bool phase_below_plane = liquid_distance < -stabilityTolerance
                || vapour_distance < -stabilityTolerance;
            if (!phase_below_plane && iteration + 1 >= flashFirstIterations) {
                break;
            }
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (Opm::getValue(z[compIdx]) > 0.) {
                    trial_K[compIdx] *= fluid_state.fugacity(oilPhaseIdx, compIdx)
                        / fluid_state.fugacity(gasPhaseIdx, compIdx);
                }
            }
        }

        if (feed_is_unstable) {
            if (verbosity >= 1) {
                OpmLog::debug("A phase of the split below the tangent plane of the feed: "
                              "unstable, no stability test.");
            }
            try {
                auto jacobian = flash_2ph(z, twoPhaseMethod, trial_K, trial_L, fluid_state,
                                          flash_tolerance, eos_type, verbosity);
                const Scalar liquid_fraction = Opm::getValue(trial_L);
                if (liquid_fraction > 0. && liquid_fraction < 1.
                    && !phasesCoincide_(fluid_state)) {
                    K = trial_K;
                    L = trial_L;
                    newton_jacobian = std::move(jacobian);
                    return true;
                }
            }
            catch (const NumericalProblem&) {
            }
            catch (const Dune::FMatrixError&) {
            }
        }
        fluid_state = incoming_state;
        return false;
    }

    /*!
     * \brief Michelsen's stability analysis of the mixture.
     *
     * A vapour-like and a liquid-like trial phase are each grown from the
     * mixture with checkStability_(). The mixture is stable, and so single
     * phase, when neither trial ends with a mole number sum above one. When
     * it is not, the trial compositions give the starting equilibrium ratios
     * \f$K_i = y_i / x_i\f$ for the two-phase flash.
     *
     * \param[out] isStable Whether the mixture stays in one phase.
     * \param[in,out] K Starting ratios in, refined ratios out for a two-phase mixture.
     * \param[in,out] fluid_state The mixture; a stable one gets \f$z\f$ as both phase compositions.
     * \param z The mixture composition.
     * \param eos_type The equation of state.
     * \param verbosity Level of debug logging.
     * \return The stable cubic root of the feed, which labels a stable one.
     */
    template <class FlashFluidState, class ComponentVector>
    static FeedRoot phaseStabilityTest_(bool& isStable, ComponentVector& K, FlashFluidState& fluid_state, const ComponentVector& z, const EOSType& eos_type, int verbosity)
    {
        OPM_FLASH_COUNT(stabilityTests);
        bool liquid_trial_is_trivial, vapour_trial_is_trivial;
        ComponentVector liquid_trial_composition, vapour_trial_composition;
        typename FlashFluidState::ValueType liquid_trial_sum, vapour_trial_sum;
        ComponentVector vapour_trial_k = K;
        ComponentVector liquid_trial_k = K;

        // Both trials are measured against the feed on its stable root,
        // evaluated once. One parameter cache serves the feed and both
        // trials: temperature and pressure are fixed.
        using ParamCache =
            typename FluidSystem::template ParameterCache<typename FlashFluidState::ValueType>;
        ParamCache param_cache(eos_type);
        std::array<typename FlashFluidState::ValueType, numComponents> reference_coefficients;
        const FeedRoot feed_root = feedStableRoot_(fluid_state, z, param_cache, reference_coefficients);
        const int reference_phase_idx = feed_root.phase_idx;

        // Grow a vapour-like trial phase.
        if (verbosity == 3 || verbosity == 4) {
            OpmLog::debug("Stability test for vapor phase:");
        }
        checkStability_(fluid_state,
                        vapour_trial_is_trivial,
                        vapour_trial_k,
                        vapour_trial_composition,
                        vapour_trial_sum,
                        z,
                        /*is_vapour_trial=*/true,
                        reference_phase_idx,
                        reference_coefficients,
                        param_cache,
                        verbosity);
        const bool vapour_trial_finds_no_instability
            = (vapour_trial_sum < (1.0 + stabilityTolerance)) || vapour_trial_is_trivial;

        // Grow a liquid-like trial phase.
        if (verbosity == 3 || verbosity == 4) {
            OpmLog::debug("Stability test for liquid phase:");
        }
        checkStability_(fluid_state,
                        liquid_trial_is_trivial,
                        liquid_trial_k,
                        liquid_trial_composition,
                        liquid_trial_sum,
                        z,
                        /*is_vapour_trial=*/false,
                        reference_phase_idx,
                        reference_coefficients,
                        param_cache,
                        verbosity);
        const bool liquid_trial_finds_no_instability
            = (liquid_trial_sum < (1.0 + stabilityTolerance)) || liquid_trial_is_trivial;

        // Neither trial phase found a negative tangent-plane distance.
        isStable = liquid_trial_finds_no_instability && vapour_trial_finds_no_instability;
        if (isStable) {
            // A single phase has the feed composition under either phase label.
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                fluid_state.setMoleFraction(gasPhaseIdx, compIdx, z[compIdx]);
                fluid_state.setMoleFraction(oilPhaseIdx, compIdx, z[compIdx]);
            }
        }
        // The two stationary trial compositions seed the equilibrium ratios.
        else {
            ComponentVector starting_k = K;
            for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
                if (Opm::getValue(z[compIdx]) > 0.) {
                    K[compIdx]
                        = vapour_trial_composition[compIdx] / liquid_trial_composition[compIdx];
                } else {
                    K[compIdx] = 1.;
                }
            }
            // Both trials can reach the same stationary point, as where the
            // feed splits into two liquids, and their ratios of one would
            // collapse the split. The feed is then the other phase. The trial
            // is the vapour if it is richer than the feed in the components
            // with the larger starting ratios, else the liquid.
            if (trivialSolutionMeasure_(K) < trivialSolutionTolerance) {
                const auto& trial_composition = vapour_trial_sum > liquid_trial_sum
                    ? vapour_trial_composition : liquid_trial_composition;
                Scalar lightness = 0.;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar ratio = Opm::getValue(starting_k[compIdx]);
                    if (Opm::getValue(z[compIdx]) > 0. && ratio > 0.) {
                        lightness += Opm::getValue(trial_composition[compIdx] - z[compIdx])
                            * std::log(ratio);
                    }
                }
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (Opm::getValue(z[compIdx]) > 0.) {
                        K[compIdx] = lightness > 0.
                            ? trial_composition[compIdx] / z[compIdx]
                            : z[compIdx] / trial_composition[compIdx];
                    }
                }
            }
        }
        return feed_root;
    }

    template <class ComponentVector>
    static Scalar trivialSolutionMeasure_(const ComponentVector& equilibrium_ratios)
    {
        Scalar measure = 0.;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const Scalar equilibrium_ratio = Opm::getValue(equilibrium_ratios[compIdx]);
            if (!(equilibrium_ratio > 0.) || !std::isfinite(equilibrium_ratio)) {
                return std::numeric_limits<Scalar>::infinity();
            }
            const Scalar log_ratio = std::log(equilibrium_ratio);
            measure += log_ratio * log_ratio;
        }
        return measure;
    }

    /*!
     * \brief One trial phase of Michelsen's stability analysis.
     *
     * A trial phase with mole numbers \f$W_i = K_i z_i\f$ (vapour-like) or
     * \f$W_i = z_i / K_i\f$ (liquid-like) is driven to a stationary point of
     * the tangent-plane-distance function by successive substitution,
     *
     * \f[ K_i \leftarrow K_i R_i, \qquad
     *     R_i = \frac{f_i(\mathbf z)}{S f_i(\mathbf W / S)} \text{ (vapour)}, \quad
     *     R_i = \frac{S f_i(\mathbf W / S)}{f_i(\mathbf z)} \text{ (liquid)}, \f]
     *
     * with \f$S = \sum_i W_i\f$. The iteration has converged when
     * \f$\sum_i (R_i - 1)^2\f$ is small. It has collapsed onto the mixture
     * itself, the trivial solution, when \f$\sum_i \ln^2 K_i\f$ is small. At a
     * stationary point the tangent-plane distance of the normalised trial
     * composition is \f$-\ln S\f$, so a trial that ends with \f$S > 1\f$
     * establishes that the mixture is unstable.
     *
     * Both trials iterate \f$\ln W_i \leftarrow d_i - \ln\varphi_i(\mathbf W / S)\f$
     * with \f$d_i = \ln z_i + \ln\varphi_i(\mathbf z)\f$, a descent method for
     * the modified tangent-plane distance
     *
     * \f[ tm(\mathbf W) = 1 + \sum_i W_i \big(\ln W_i + \ln\varphi_i - d_i - 1\big)
     *                   = 1 - S \mp \sum_i W_i \ln R_i \f]
     *
     * (upper sign vapour), whose stationary points are those of the trial
     * (Michelsen, Fluid Phase Equilib. 9 (1982) 1-19). Near a critical point
     * or the stability limit the steps in \f$\ln K\f$ shrink by a factor close
     * to one and plain substitution takes hundreds of steps. Every
     * stabilityAccelerationInterval plain steps the remaining steps are
     * therefore summed from the last three, as in extrapolateSubstitution_(),
     * and the extrapolated ratios are kept only if they lower \f$tm\f$ below
     * its value at the last evaluated ratios; otherwise the iteration
     * continues from the plain step. An extrapolation is shortened to change
     * no \f$\ln K_i\f$ by more than maxStabilityExtrapolation: where a fugacity
     * coefficient is clamped the identity above no longer holds, and a far
     * extrapolation could then lower \f$tm\f$ without approaching a
     * stationary point. The stationary points, the convergence
     * and trivial-solution tests are those of plain substitution.
     *
     * \param fluid_state The mixture at the pressure and temperature of the flash.
     * \param[out] isTrivial Whether the trial collapsed onto the mixture.
     * \param[in,out] K Starting ratios in, the trial's converged ratios out.
     * \param[out] trial_composition The normalised trial composition.
     * \param[out] trial_sum The mole number sum \f$S\f$ of the trial phase.
     * \param z The mixture composition.
     * \param is_vapour_trial Whether to grow a vapour-like or liquid-like trial.
     * \param reference_phase_idx The phase index of the feed's stable root,
     *        which supplies the reference fugacities.
     * \param eos_type The equation of state.
     * \param verbosity Level of debug logging.
     */
    template <class FlashFluidState, class ComponentVector>
    static void checkStability_(const FlashFluidState& fluid_state,
                                bool& isTrivial,
                                ComponentVector& K,
                                ComponentVector& trial_composition,
                                typename FlashFluidState::ValueType& trial_sum,
                                const ComponentVector& z,
                                bool is_vapour_trial,
                                const int reference_phase_idx,
                                const EOSType& eos_type,
                                int verbosity)
    {
        using FlashEval = typename FlashFluidState::ValueType;
        using CubicEOS = typename Opm::CubicEOS<Scalar, FluidSystem>;
        using ParamCache = typename FluidSystem::template ParameterCache<FlashEval>;
        FlashFluidState reference_state = fluid_state;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            reference_state.setMoleFraction(reference_phase_idx, compIdx, z[compIdx]);
        }
        ParamCache param_cache(eos_type);
        param_cache.updatePhase(reference_state, reference_phase_idx);
        std::array<FlashEval, numComponents> reference_coefficients;
        CubicEOS::computeFugacityCoefficients(
            reference_state, param_cache, reference_phase_idx, reference_coefficients);
        checkStability_(fluid_state, isTrivial, K, trial_composition, trial_sum, z,
                        is_vapour_trial, reference_phase_idx, reference_coefficients,
                        param_cache, verbosity);
    }

    /*!
     * \brief checkStability_() with the fugacity coefficients of the feed on
     *        its stable root given, and a parameter cache at the temperature
     *        and pressure of the flash, which the trial updates.
     */
    template <class FlashFluidState, class ComponentVector, class ParamCache>
    static void checkStability_(const FlashFluidState& fluid_state,
                                bool& isTrivial,
                                ComponentVector& K,
                                ComponentVector& trial_composition,
                                typename FlashFluidState::ValueType& trial_sum,
                                const ComponentVector& z,
                                bool is_vapour_trial,
                                const int reference_phase_idx,
                                const std::array<typename FlashFluidState::ValueType, numComponents>&
                                    reference_coefficients,
                                ParamCache& trial_param_cache,
                                int verbosity)
    {
        using FlashEval = typename FlashFluidState::ValueType;
        using CubicEOS = typename Opm::CubicEOS<Scalar, FluidSystem>;

        // The trial composition is evaluated on its candidate phase root. The
        // feed composition supplies the reference fugacities on its stable
        // root, the state whose tangent plane the trial is measured against.
        FlashFluidState trial_state = fluid_state;
        FlashFluidState reference_state = fluid_state;
        const int trial_phase_idx
            = is_vapour_trial ? static_cast<int>(gasPhaseIdx) : static_cast<int>(oilPhaseIdx);

        // An absent component has no fugacity ratio. Setting its undefined K
        // to one excludes it from both convergence measures below.
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!(Opm::getValue(z[compIdx]) > 0.)) {
                K[compIdx] = 1.;
            }
        }

        // Setup output
        if (verbosity >= 3) {
            OpmLog::debug(fmt::format("{:>10}{:>16}{:>16}", "Iteration", "K-Norm", "R-Norm"));
        }

        // The reference phase is the feed on its stable root throughout.
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            reference_state.setMoleFraction(reference_phase_idx, compIdx, z[compIdx]);
            reference_state.setFugacityCoefficient(
                reference_phase_idx, compIdx, reference_coefficients[compIdx]);
        }
        // The last three plain steps in ln K, newest first, and the plain
        // point to return to when an extrapolation does not lower tm.
        using FieldType = typename ComponentVector::field_type;
        std::array<ComponentVector, 3> steps{};
        int plain_steps = 0;
        bool extrapolated = false;
        ComponentVector plain_K = K;
        FieldType tm_before = 0.;
        int resolving_steps = 0;
        bool newton_step = false;
        bool newton_rejected = false;

        // Grow the trial phase out of the mixture and find a stationary point
        // of the tangent-plane-distance function.
        constexpr int max_iterations = 20000;
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            OPM_FLASH_COUNT(stabilityIterations);
            // Mole numbers of the trial phase, then its normalised composition.
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                trial_composition[compIdx]
                    = is_vapour_trial ? K[compIdx] * z[compIdx] : z[compIdx] / K[compIdx];
            }
            trial_sum = std::accumulate(
                trial_composition.begin(), trial_composition.end(), FlashEval(0.0));
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                trial_composition[compIdx] /= trial_sum;
                trial_state.setMoleFraction(trial_phase_idx, compIdx, trial_composition[compIdx]);
            }

            trial_param_cache.updatePhase(trial_state, trial_phase_idx);

            // Fugacity coefficients of the trial phase.
            std::array<FlashEval, numComponents> trial_fugacity_coefficients;
            CubicEOS::computeFugacityCoefficients(
                trial_state, trial_param_cache, trial_phase_idx, trial_fugacity_coefficients);
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                trial_state.setFugacityCoefficient(
                    trial_phase_idx, compIdx, trial_fugacity_coefficients[compIdx]);
            }

            // The substitution factors R_i of the brief.
            ComponentVector substitution_factor;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (!(Opm::getValue(z[compIdx]) > 0.)) {
                    substitution_factor[compIdx] = 1.;
                    continue;
                }
                const auto trial_fugacity = trial_state.fugacity(trial_phase_idx, compIdx);
                const auto reference_fugacity
                    = reference_state.fugacity(reference_phase_idx, compIdx);
                substitution_factor[compIdx] = is_vapour_trial
                    ? reference_fugacity / trial_fugacity / trial_sum
                    : trial_fugacity / reference_fugacity * trial_sum;
            }

            // The plain steps in ln K and the modified tangent-plane distance.
            ComponentVector step(0.);
            FieldType tm = 1. - trial_sum;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (!(Opm::getValue(z[compIdx]) > 0.)) {
                    continue;
                }
                step[compIdx] = log(substitution_factor[compIdx]);
                const FieldType mole_number = trial_sum * trial_composition[compIdx];
                tm += is_vapour_trial ? -mole_number * step[compIdx] : mole_number * step[compIdx];
            }
            if (extrapolated) {
                extrapolated = false;
                if (!(tm < tm_before) || !std::isfinite(Opm::getValue(tm))) {
                    // Continue from the plain step instead, and by
                    // substitution after a rejected second-order step.
                    K = plain_K;
                    plain_steps = 0;
                    newton_rejected = newton_rejected || newton_step;
                    continue;
                }
            }

            const ComponentVector evaluated_K = K;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                K[compIdx] *= substitution_factor[compIdx];
            }
            Scalar substitution_residual = 0.0;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (!(Opm::getValue(z[compIdx]) > 0.)) {
                    continue;
                }
                const auto substitution_error = Opm::getValue(substitution_factor[compIdx]) - 1.0;
                substitution_residual += substitution_error * substitution_error;
            }
            // A converged trial may have converged onto the mixture itself.
            // Where its mole number sum is within the resolution of the
            // convergence test of one, the second-order steps continue until
            // the sum decides.
            const bool resolves_decision = substitution_residual < resolvedStabilityTolerance
                || std::abs(Opm::getValue(trial_sum) - 1.) > unresolvedStabilityBand
                || iteration + 1 < stabilityNewtonStart || newton_rejected
                || resolving_steps >= maxResolvingSteps;
            if (substitution_residual < substitutionTolerance && !resolves_decision) {
                ++resolving_steps;
            }
            else if (substitution_residual < substitutionTolerance) {
                isTrivial = trivialSolutionMeasure_(K) < trivialSolutionTolerance;
                if (verbosity >= 3) {
                    OpmLog::debug(fmt::format("{:>10}{:>16}", iteration, substitution_residual));
                }
                return;
            }

            steps = {step, steps[0], steps[1]};
            ++plain_steps;

            // A second-order step from the evaluated point, once substitution
            // has had a few steps to settle, kept on the same terms as an
            // extrapolation.
            newton_step = false;
            if (!newton_rejected && iteration + 1 >= stabilityNewtonStart) {
                bool is_clamped = false;
                std::array<Scalar, numComponents> composition, gradient;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = Opm::getValue(trial_fugacity_coefficients[compIdx]);
                    is_clamped = is_clamped
                        || !(phi > minFugacityCoefficient && phi < maxFugacityCoefficient);
                    const bool is_present = Opm::getValue(z[compIdx]) > 0.;
                    composition[compIdx]
                        = is_present ? Opm::getValue(trial_composition[compIdx]) : 0.;
                    // g_i = ln W_i + ln phi_i - d_i, the negative plain step
                    // of a vapour trial and the plain step of a liquid one.
                    gradient[compIdx] = is_vapour_trial ? -Opm::getValue(step[compIdx])
                                                        : Opm::getValue(step[compIdx]);
                }
                std::array<Scalar, numComponents> growth;
                if (!is_clamped
                    && stabilityNewtonStep_(trial_param_cache, trial_phase_idx,
                                            Opm::getValue(trial_state.pressure(trial_phase_idx)),
                                            Opm::getValue(trial_state.temperature(trial_phase_idx)),
                                            composition, gradient, growth)) {
                    plain_K = K;
                    tm_before = tm;
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                        K[compIdx] = is_vapour_trial ? evaluated_K[compIdx] * growth[compIdx]
                                                     : evaluated_K[compIdx] / growth[compIdx];
                    }
                    newton_step = true;
                }
            }

            // The trial has collapsed onto the mixture when the point it
            // evaluates next is the trivial solution.
            const Scalar trivial_residual = trivialSolutionMeasure_(K);
            if (verbosity >= 3) {
                OpmLog::debug(fmt::format(
                    "{:>10}{:>16}{:>16}", iteration, trivial_residual, substitution_residual));
            }
            isTrivial = trivial_residual < trivialSolutionTolerance;
            if (isTrivial) {
                return;
            }
            if (newton_step) {
                extrapolated = true;
                plain_steps = 0;
                continue;
            }

            ComponentVector extrapolation;
            if (plain_steps >= stabilityAccelerationInterval
                && extrapolateSubstitution_(steps, extrapolation)) {
                plain_K = K;
                tm_before = tm;
                const FieldType largest = extrapolation.infinity_norm();
                const FieldType shortening = largest > maxStabilityExtrapolation
                    ? maxStabilityExtrapolation / largest : FieldType{1.};
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    K[compIdx] *= exp(shortening * extrapolation[compIdx]);
                }
                extrapolated = true;
                plain_steps = 0;
            }
        }
        OPM_THROW_NOLOG(NumericalProblem, " Stability test did not converge");
    }//end checkStability

    /*!
     * \brief The composition derivatives of the fugacity coefficients of a
     *        phase, \f$n\,\partial \ln\varphi_i / \partial n_j\f$ at constant
     *        temperature and pressure, with \f$n\f$ the total mole number.
     *
     * For the general two-parameter cubic of CubicEOS,
     *
     * \f[ \ln\varphi_i = \frac{B_i}{B}(Z - 1) - \ln(Z - B)
     *   + \beta \left(\frac{2 a_i}{A} - \frac{B_i}{B}\right), \quad
     *   \beta = \frac{A}{(m_1 - m_2) B} \ln\frac{Z + m_2 B}{Z + m_1 B}, \quad
     *   a_i = \sum_k x_k A_{ik}, \f]
     *
     * is differentiated with the mole fractions as independent variables,
     * \f$\partial A / \partial x_j = 2 a_j\f$, \f$\partial B / \partial x_j = B_j\f$
     * and \f$\partial Z / \partial x_j = -(F_A\,2 a_j + F_B B_j) / F_Z\f$ from the
     * cubic \f$F(Z; A, B) = 0\f$. Since \f$\ln\varphi_i\f$ is homogeneous of degree
     * zero in the mole numbers,
     * \f$n\,\partial \ln\varphi_i / \partial n_j = D_{ij} - \sum_k x_k D_{ik}\f$
     * with \f$D_{ij} = \partial \ln\varphi_i / \partial x_j\f$.
     *
     * \param cache The parameter cache, updated for the phase.
     * \param phase_idx The phase.
     * \param x The normalised phase composition.
     * \param[out] derivatives The matrix, row i for \f$\ln\varphi_i\f$.
     * \return false where the cubic is degenerate at the cached root
     */
    template <class ParamCache>
    static bool fugacityCompositionDerivatives_(
        const ParamCache& cache,
        const unsigned phase_idx,
        const Scalar pressure,
        const Scalar temperature,
        const std::array<Scalar, numComponents>& x,
        Dune::FieldMatrix<Scalar, numComponents, numComponents>& derivatives)
    {
        constexpr int n = numComponents;
        const Scalar A = Opm::getValue(cache.A(phase_idx));
        const Scalar B = Opm::getValue(cache.B(phase_idx));
        const Scalar m1 = Opm::getValue(cache.m1(phase_idx));
        const Scalar m2 = Opm::getValue(cache.m2(phase_idx));
        const Scalar Z = pressure * Opm::getValue(cache.molarVolume(phase_idx))
            / (Constants<Scalar>::R * temperature);
        if (!(A > 0. && B > 0. && Z > B)) {
            return false;
        }

        // The cubic of CubicEOS::computeMolarVolume() and its partial derivatives.
        const Scalar c2 = (m1 + m2 - 1.) * B - 1.;
        const Scalar c1 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1.);
        const Scalar F_Z = (3. * Z + 2. * c2) * Z + c1;
        const Scalar F_A = Z - B;
        const Scalar F_B = (m1 + m2 - 1.) * Z * Z
            + (2. * m1 * m2 * B - (m1 + m2) * (2. * B + 1.)) * Z
            - A - m1 * m2 * (3. * B * B + 2. * B);
        if (!(std::abs(F_Z) > 1e-12 * (1. + std::abs(c1)))) {
            return false;
        }

        // Per component j: dZ/dx_j, and with it the x_j-derivatives of the
        // terms of ln phi_i that do not depend on i.
        const Scalar inv_A = 1. / A;
        const Scalar inv_B = 1. / B;
        const Scalar Zp2 = Z + m2 * B;
        const Scalar Zp1 = Z + m1 * B;
        const Scalar inv_Zp2 = 1. / Zp2;
        const Scalar inv_Zp1 = 1. / Zp1;
        const Scalar inv_ZmB = 1. / (Z - B);
        const Scalar inv_F_Z = 1. / F_Z;
        const Scalar coefficient = A * inv_B / (m1 - m2);
        const Scalar beta = coefficient * std::log(Zp2 * inv_Zp1);
        std::array<std::array<Scalar, n>, n> A_ij;
        std::array<Scalar, n> a, Bi, gamma, u, v, dbeta;
        Scalar mean_u = 0., mean_v = 0., mean_dbeta = 0.;
        for (int j = 0; j < n; ++j) {
            a[j] = 0.;
            for (int k = 0; k < n; ++k) {
                A_ij[j][k] = Opm::getValue(cache.aCache(phase_idx, j, k));
                a[j] += A_ij[j][k] * x[k];
            }
            Bi[j] = Opm::getValue(cache.Bi(phase_idx, j));
            gamma[j] = 2. * a[j] * inv_A - Bi[j] * inv_B;
            const Scalar dZ = -(F_A * 2. * a[j] + F_B * Bi[j]) * inv_F_Z;
            u[j] = (dZ - (Z - 1.) * Bi[j] * inv_B) * inv_B;
            v[j] = (Bi[j] - dZ) * inv_ZmB;
            dbeta[j] = beta * gamma[j]
                + coefficient * ((dZ + m2 * Bi[j]) * inv_Zp2 - (dZ + m1 * Bi[j]) * inv_Zp1);
            mean_u += x[j] * u[j];
            mean_v += x[j] * v[j];
            mean_dbeta += x[j] * dbeta[j];
        }
        // D_ij = d ln phi_i / d x_j
        //      = B_i u_j + v_j + gamma_i dbeta_j
        //        + beta (2 A_ij / A - 4 a_i a_j / A^2 + B_i B_j / B^2),
        // and with sum_j x_j A_ij = a_i, sum_j x_j a_j = A, sum_j x_j B_j = B
        // its x-weighted row sum is B_i mean_u + mean_v + gamma_i (mean_dbeta - beta).
        // The difference is symmetric, so one triangle is evaluated.
        const Scalar c_A = 2. * beta * inv_A;
        const Scalar c_aa = 4. * beta * inv_A * inv_A;
        const Scalar c_BB = beta * inv_B * inv_B;
        for (int j = 0; j < n; ++j) {
            u[j] -= mean_u;
            v[j] -= mean_v;
            dbeta[j] += beta - mean_dbeta;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                const Scalar value = Bi[i] * u[j] + v[j] + gamma[i] * dbeta[j]
                    + c_A * A_ij[i][j]
                    - c_aa * a[i] * a[j] + c_BB * Bi[i] * Bi[j];
                derivatives[i][j] = value;
                derivatives[j][i] = value;
            }
        }
        return true;
    }

    /*!
     * \brief A second-order step of a stability trial.
     *
     * In the variables \f$\alpha_i = 2\sqrt{W_i}\f$ the modified tangent-plane
     * distance \f$tm\f$ of checkStability_() has the gradient
     * \f$\sqrt{W_i}\,g_i\f$, \f$g_i = \ln W_i + \ln\varphi_i - d_i\f$, and, up to
     * a term \f$\delta_{ij} g_i / 2\f$ that vanishes at a stationary point, the
     * Hessian
     *
     * \f[ H_{ij} = \delta_{ij} + \sqrt{W_i W_j}\,
     *              \frac{\partial \ln\varphi_i}{\partial n_j}
     *            = \delta_{ij} + \sqrt{w_i w_j}\, S
     *              \frac{\partial \ln\varphi_i}{\partial n_j}, \f]
     *
     * the identity for an ideal mixture, where the step is a substitution
     * step (Michelsen, Fluid Phase Equilib. 9 (1982) 1-19; Michelsen and
     * Mollerup, Thermodynamic Models, 2nd ed. (2007), ch. 10). The step
     * solves \f$H\,\Delta\alpha = -\sqrt{W}\,g\f$ and is shortened to change no
     * \f$\ln W_i\f$ by more than maxStabilityNewtonStep.
     *
     * \param cache The parameter cache, updated for the trial phase.
     * \param phase_idx The trial phase.
     * \param w The normalised trial composition, zero for an absent component.
     * \param g The gradient \f$g_i\f$, zero for an absent component.
     * \param[out] growth The factors \f$W_i^{new} / W_i\f$.
     * \return false where \f$H\f$ is not positive definite or the cubic is
     *         degenerate; the trial then continues by substitution
     */
    template <class ParamCache>
    static bool stabilityNewtonStep_(const ParamCache& cache,
                                     const unsigned phase_idx,
                                     const Scalar pressure,
                                     const Scalar temperature,
                                     const std::array<Scalar, numComponents>& w,
                                     const std::array<Scalar, numComponents>& g,
                                     std::array<Scalar, numComponents>& growth)
    {
        constexpr int n = numComponents;
        Dune::FieldMatrix<Scalar, n, n> H;
        if (!fugacityCompositionDerivatives_(cache, phase_idx, pressure, temperature, w, H)) {
            return false;
        }
        std::array<Scalar, n> root_w, rhs;
        for (int i = 0; i < n; ++i) {
            root_w[i] = std::sqrt(w[i]);
            rhs[i] = -root_w[i] * g[i];
        }

        // H = L D L^T of the lower triangle of H, then the solve.
        std::array<Scalar, n> inv_diagonal;
        for (int j = 0; j < n; ++j) {
            Scalar diagonal = 1. + w[j] * H[j][j];
            for (int k = 0; k < j; ++k) {
                diagonal -= H[j][k] * H[j][k] * inv_diagonal[k];
            }
            if (!(diagonal > 1e-12)) {
                return false;
            }
            inv_diagonal[j] = 1. / diagonal;
            for (int i = j + 1; i < n; ++i) {
                // H[i][j] holds L_ij D_j while the column is built.
                Scalar entry = root_w[i] * root_w[j] * H[i][j];
                for (int k = 0; k < j; ++k) {
                    entry -= H[i][k] * H[j][k] * inv_diagonal[k];
                }
                H[i][j] = entry;
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < i; ++k) {
                rhs[i] -= H[i][k] * inv_diagonal[k] * rhs[k];
            }
        }
        for (int i = n - 1; i >= 0; --i) {
            rhs[i] *= inv_diagonal[i];
            for (int k = i + 1; k < n; ++k) {
                rhs[i] -= H[k][i] * inv_diagonal[i] * rhs[k];
            }
        }

        // sqrt(W_i^new / W_i) = 1 + t_i with t_i = Delta alpha_i / alpha_i,
        // shortened so that |ln W_i^new - ln W_i| <= maxStabilityNewtonStep.
        static const Scalar max_root_growth = std::exp(0.5 * maxStabilityNewtonStep);
        static const Scalar min_root_growth = std::exp(-0.5 * maxStabilityNewtonStep);
        Scalar shortening = 1.;
        for (int i = 0; i < n; ++i) {
            growth[i] = w[i] > 0. ? 0.5 * rhs[i] / root_w[i] : 0.;
            if (!std::isfinite(growth[i])) {
                return false;
            }
            if (growth[i] > max_root_growth - 1.) {
                shortening = std::min(shortening, (max_root_growth - 1.) / growth[i]);
            }
            else if (growth[i] < min_root_growth - 1.) {
                shortening = std::min(shortening, (min_root_growth - 1.) / growth[i]);
            }
        }
        for (int i = 0; i < n; ++i) {
            const Scalar root_growth = 1. + shortening * growth[i];
            growth[i] = root_growth * root_growth;
        }
        return true;
    }

    /*!
     * \brief The phase compositions from the equilibrium ratios and the
     *        liquid fraction,
     *
     * \f[ x_i = \frac{z_i}{L + (1 - L) K_i}, \qquad y_i = K_i x_i, \f]
     *
     * each renormalised to sum to one.
     */
    template <class FlashFluidState, class ComponentVector>
    static void computeLiquidVapor_(FlashFluidState& fluid_state, typename FlashFluidState::ValueType& L, ComponentVector& K, const ComponentVector& z)
    {
        // Calculate x and y, and normalize
        ComponentVector x;
        ComponentVector y;
        typename FlashFluidState::ValueType sumx=0;
        typename FlashFluidState::ValueType sumy=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            x[compIdx] = z[compIdx]/(L + (1-L)*K[compIdx]);
            sumx += x[compIdx];
            y[compIdx] = (K[compIdx]*z[compIdx])/(L + (1-L)*K[compIdx]);
            sumy += y[compIdx];
        }
        x /= sumx;
        y /= sumy;

        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            fluid_state.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluid_state.setMoleFraction(gasPhaseIdx, compIdx, y[compIdx]);
        }
    }

    /*!
     * \brief The two-phase compositions by the requested method.
     *
     * \return The composition Jacobian when Newton finished the solve.
     */
    template <class FluidState, class ComponentVector>
    static auto flash_2ph(const ComponentVector& z_scalar,
                          PTFlashMethod flash_2p_method,
                          ComponentVector& K_scalar,
                          typename FluidState::ValueType& L_scalar,
                          FluidState& fluid_state_scalar,
                          const Scalar flash_tolerance,
                          const EOSType& eos_type,
                          int verbosity = 0) -> std::optional<NewtonJacobian> {
        if (verbosity >= 1) {
            OpmLog::debug(fmt::format("Cell is two-phase! Solve Rachford-Rice with initial K = [{}]",
                                     fmt::join(K_scalar, " ")));
        }

        // Calculate composition using nonlinear solver
        bool converged = false;
        std::optional<NewtonJacobian> newton_jacobian;
        switch (flash_2p_method) {
        case PTFlashMethod::Newton:
            if (verbosity >= 1) {
                OpmLog::debug("Calculate composition using Newton.");
            }
            newton_jacobian = newtonComposition_(K_scalar, L_scalar, fluid_state_scalar, z_scalar,
                                                 flash_tolerance, eos_type, verbosity);
            converged = true;
            break;

        case PTFlashMethod::Ssi:
            if (verbosity >= 1) {
                OpmLog::debug("Calculate composition using Successive Substitution.");
            }
            converged = successiveSubstitutionComposition_(K_scalar, L_scalar, fluid_state_scalar, z_scalar, false, flash_tolerance, eos_type, verbosity);
            break;

        case PTFlashMethod::SsiAccelerated:
            if (verbosity >= 1) {
                OpmLog::debug("Calculate composition using extrapolated Successive Substitution.");
            }
            converged = successiveSubstitutionComposition_(K_scalar, L_scalar, fluid_state_scalar, z_scalar, false, flash_tolerance, eos_type, verbosity,
                                                           /*extrapolate=*/true);
            break;

        case PTFlashMethod::SsiNewton:
            converged = successiveSubstitutionComposition_(K_scalar, L_scalar, fluid_state_scalar, z_scalar, true, flash_tolerance, eos_type, verbosity);
            if (!converged) {
                // In this method Newton is an accelerator, not the
                // guarantee: when its composition update fails to converge —
                // or a diverged iterate makes the EoS report non-finite
                // mixing parameters — switch back to plain successive
                // substitution and finish the solve with the method that is
                // robust on these states, instead of failing the whole
                // solve. Newton throws rather than returning a bad result;
                // it may already have mutated some fluid state (e.g. via
                // computeLiquidVapor_()) before throwing, but the successive-
                // substitution fallback re-derives the composition from K and
                // L and overwrites that state, so restarting it here is safe.
                const auto fall_back_to_ssi = [&](const auto& error) {
                    if (verbosity >= 1) {
                        OpmLog::debug(fmt::format(
                            "Newton did not finish the composition update ({}); "
                            "switching back to successive substitution.", error.what()));
                    }
                    converged = successiveSubstitutionComposition_(K_scalar, L_scalar, fluid_state_scalar, z_scalar, false, flash_tolerance, eos_type, verbosity,
                                                                   /*extrapolate=*/true);
                };
                try {
                    newton_jacobian = newtonComposition_(K_scalar, L_scalar, fluid_state_scalar,
                                                         z_scalar, flash_tolerance, eos_type,
                                                         verbosity);
                    converged = true;
                }
                catch (const NumericalProblem& error) {
                    fall_back_to_ssi(error);
                }
                catch (const Dune::FMatrixError& error) {
                    // FieldMatrix::solve reports a singular Jacobian with
                    // FMatrixError, which is outside the std::runtime_error hierarchy.
                    fall_back_to_ssi(error);
                }
            }
            break;
        }

        if (!converged) {
            OPM_THROW_NOLOG(NumericalProblem,
                            fmt::format("flash calculation did not get converged with {}",
                                        ptFlashMethodToString(flash_2p_method)));
        }

        // A finite L outside [0, 1] is a valid negative-flash result and is
        // classified as single-phase by flashScalar_().
        if (!std::isfinite(Opm::getValue(L_scalar))
            || !phaseCompositionsArePhysical_(fluid_state_scalar)) {
            OPM_THROW_NOLOG(
                NumericalProblem,
                fmt::format("Two-phase flash converged to a non-physical composition "
                            "with L = {} using {}",
                            Opm::getValue(L_scalar),
                            ptFlashMethodToString(flash_2p_method)));
        }
        return newton_jacobian;
    }

    template <class FlashFluidState>
    static bool phaseCompositionsArePhysical_(const FlashFluidState& fluid_state)
    {
        // Keep the variable bound independent of the fugacity-residual tolerance.
        const Scalar bound_slack
            = std::max(Scalar{1.e-8}, std::numeric_limits<Scalar>::epsilon());
        const auto is_out_of_range = [bound_slack](const Scalar value) {
            return !std::isfinite(value) || value < -bound_slack
                || value > Scalar{1} + bound_slack;
        };

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const auto liquid_mole_fraction
                = Opm::getValue(fluid_state.moleFraction(oilPhaseIdx, compIdx));
            const auto vapour_mole_fraction
                = Opm::getValue(fluid_state.moleFraction(gasPhaseIdx, compIdx));
            if (is_out_of_range(liquid_mole_fraction)
                || is_out_of_range(vapour_mole_fraction)) {
                return false;
            }
        }
        return true;
    }

    /*!
     * \brief Whether a candidate two-phase state has coincident phases.
     *
     * \f$\mathbf x = \mathbf y\f$ solves the fugacity and closure equations for
     * any \f$L\f$, so the component balance \f$z_i - L x_i - (1 - L) y_i\f$ loses
     * its dependence on \f$L\f$ and makes the Newton Jacobian singular in the
     * liquid-fraction direction. Round-off can then produce a small residual at
     * an unbounded, non-physical value of \f$L\f$.
     *
     * For each active component, the measure
     * \f$\sum_i (\ln(y_i/x_i))^2\f$ approaches zero as the phases coincide.
     * Computing the ratio from the phase compositions is important because
     * SSI may retain a uniformly scaled \f$K\f$ vector after normalising both
     * phases. A component whose ratio is undefined is omitted from this measure;
     * final-state bounds validation independently rejects non-physical mole
     * fractions. At least two defined ratios are required because normalising a
     * single active component makes its phase compositions equal for every K.
     */
    template <class FlashFluidState>
    static bool phasesCoincide_(const FlashFluidState& fluid_state)
    {
        Dune::FieldVector<Scalar, numComponents> equilibrium_ratios(1.);
        unsigned active_components = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const Scalar liquid_mole_fraction
                = Opm::getValue(fluid_state.moleFraction(oilPhaseIdx, compIdx));
            const Scalar vapour_mole_fraction
                = Opm::getValue(fluid_state.moleFraction(gasPhaseIdx, compIdx));
            if (!(liquid_mole_fraction > 0.) || !(vapour_mole_fraction > 0.)) {
                continue;
            }
            ++active_components;
            equilibrium_ratios[compIdx] = vapour_mole_fraction / liquid_mole_fraction;
        }
        return active_components >= 2
            && trivialSolutionMeasure_(equilibrium_ratios) < trivialSolutionTolerance;
    }

    /*!
     * \brief Newton's method for the two-phase compositions.
     *
     * The unknowns are \f$(\mathbf x, \mathbf y, L)\f$ and the residuals,
     * assembled by assembleNewton_(), are the component balances, equal
     * fugacities and the closure
     *
     * \f[ z_i - L x_i - (1 - L) y_i = 0, \qquad
     *     f_i^{L} - f_i^{V} = 0, \qquad
     *     \sum_i x_i - \sum_i y_i = 0. \f]
     *
     * The Jacobian comes from automatic differentiation of the residuals in
     * the unknowns. The convergence test follows its assembly, so the
     * Jacobian returned belongs to the converged state.
     *
     * \throws NumericalProblem when Newton does not converge.
     */
    template <class FlashFluidState, class ComponentVector>
    static auto newtonComposition_(ComponentVector& K,
                                   typename FlashFluidState::ValueType& L,
                                   FlashFluidState& fluid_state,
                                   const ComponentVector& z,
                                   const Scalar tolerance,
                                   const EOSType& eos_type,
                                   int verbosity) -> NewtonJacobian
    {
        constexpr std::size_t num_equations = numMisciblePhases * numMiscibleComponents + 1;
        constexpr std::size_t num_unknowns = numMisciblePhases * numMiscibleComponents + 1;
        using NewtonVector = Dune::FieldVector<Scalar, num_equations>;
        using NewtonMatrix = Dune::FieldMatrix<Scalar, num_equations, num_unknowns>;

        NewtonVector newton_step(0.);
        NewtonVector residual(0.);
        NewtonMatrix jacobian(0.);

        // Compute x and y from K, L and Z
        computeLiquidVapor_(fluid_state, L, K, z);

        // Print initial condition
        if (verbosity >= 1) {
            OpmLog::debug(fmt::format(" the current L is {}", Opm::getValue(L)));
            std::vector<Scalar> liquid_composition_values(numComponents);
            std::vector<Scalar> vapour_composition_values(numComponents);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                liquid_composition_values[compIdx]
                    = Opm::getValue(fluid_state.moleFraction(oilPhaseIdx, compIdx));
                vapour_composition_values[compIdx]
                    = Opm::getValue(fluid_state.moleFraction(gasPhaseIdx, compIdx));
            }
            OpmLog::debug(fmt::format("Initial guess: x = [{}], y = [{}], and L = {}",
                                      fmt::join(liquid_composition_values, " "),
                                      fmt::join(vapour_composition_values, " "),
                                      Opm::getValue(L)));
        }

        // Print header
        if (verbosity == 2 || verbosity == 4) {
            OpmLog::debug(fmt::format("{:>10}{:>16}{:>16}", "Iteration", "Norm2(step)", "Norm2(Residual)"));
        }

        // Each flash unknown is an independent automatic-differentiation variable.
        using NewtonEval = DenseAd::Evaluation<Scalar, num_unknowns>;
        std::vector<NewtonEval> liquid_composition(numComponents);
        std::vector<NewtonEval> vapour_composition(numComponents);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            liquid_composition[compIdx]
                = NewtonEval(fluid_state.moleFraction(oilPhaseIdx, compIdx), compIdx);
            const unsigned idx = compIdx + numComponents;
            vapour_composition[compIdx]
                = NewtonEval(fluid_state.moleFraction(gasPhaseIdx, compIdx), idx);
        }
        auto liquid_fraction = NewtonEval(L, num_unknowns - 1);

        // This state carries derivatives with respect to the Newton unknowns.
        CompositionalFluidState<NewtonEval, FluidSystem> newton_state;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            newton_state.setMoleFraction(
                FluidSystem::oilPhaseIdx, compIdx, liquid_composition[compIdx]);
            newton_state.setMoleFraction(
                FluidSystem::gasPhaseIdx, compIdx, vapour_composition[compIdx]);
            newton_state.setKvalue(compIdx,
                                   vapour_composition[compIdx] / liquid_composition[compIdx]);
        }
        newton_state.setLvalue(liquid_fraction);
        newton_state.setPressure(FluidSystem::oilPhaseIdx,
                                 fluid_state.pressure(FluidSystem::oilPhaseIdx));
        newton_state.setPressure(FluidSystem::gasPhaseIdx,
                                 fluid_state.pressure(FluidSystem::gasPhaseIdx));

        // Preserve the caller's saturations, although they do not enter the
        // current flash residual.
        newton_state.setSaturation(FluidSystem::gasPhaseIdx,
                                   fluid_state.saturation(FluidSystem::gasPhaseIdx));
        newton_state.setSaturation(FluidSystem::oilPhaseIdx,
                                   fluid_state.saturation(FluidSystem::oilPhaseIdx));
        newton_state.setTemperature(fluid_state.temperature(0));

        // Reject the degenerate root before assembling its singular Jacobian.
        // The hybrid method can recover by returning to successive substitution.
        if (phasesCoincide_(newton_state)) {
            OPM_THROW_NOLOG(NumericalProblem,
                            " Newton was started from the trivial solution, where the "
                            "composition Jacobian is singular in the liquid fraction");
        }

        using ParamCache = typename FluidSystem::template ParameterCache<NewtonEval>;
        ParamCache paramCache(eos_type);

        // The state at the current unknowns, with its fugacity coefficients,
        // residual and Jacobian; returns the residual norm, infinite where it
        // is not finite.
        const auto evaluate = [&]() {
            OPM_FLASH_COUNT(newtonIterations);
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                newton_state.setMoleFraction(
                    FluidSystem::oilPhaseIdx, compIdx, liquid_composition[compIdx]);
                newton_state.setMoleFraction(
                    FluidSystem::gasPhaseIdx, compIdx, vapour_composition[compIdx]);
                newton_state.setKvalue(compIdx,
                                       vapour_composition[compIdx] / liquid_composition[compIdx]);
            }
            newton_state.setLvalue(liquid_fraction);
            for (unsigned phaseIdx = 0; phaseIdx < numMisciblePhases; ++phaseIdx) {
                paramCache.updatePhase(newton_state, phaseIdx);
                updateFugacityCoefficients_(newton_state, paramCache, phaseIdx);
            }
            assembleNewton_<CompositionalFluidState<NewtonEval, FluidSystem>,
                            ComponentVector,
                            num_unknowns,
                            num_equations>(newton_state, z, jacobian, residual);
            const Scalar norm = residual.two_norm();
            return std::isfinite(norm) ? norm : std::numeric_limits<Scalar>::infinity();
        };

        bool converged = false;
        unsigned iteration = 0;
        constexpr unsigned max_iterations = 1000;
        // Newton that has not halved its residual for this many steps is not
        // converging to the root it started near: it wanders, or stalls on a
        // local minimum of the residual norm where the Jacobian is nearly
        // singular, and would only stop at the iteration limit. Stop it, so
        // that the hybrid method falls back to substitution from the same
        // ratios at once. Converging Newton halves the residual every step.
        constexpr unsigned max_steps_without_progress = 10;
        constexpr Scalar progress_factor = 0.5;
        unsigned steps_without_progress = 0;

        Scalar residual_norm = evaluate();
        Scalar progress_residual = residual_norm;
        while (iteration < max_iterations) {
            if (verbosity >= 1) {
                OpmLog::debug(fmt::format(" newton residual is {}", residual_norm));
            }
            converged = residual_norm < tolerance;
            if (converged) {
                break;
            }
            if (residual_norm < progress_factor * progress_residual) {
                progress_residual = residual_norm;
                steps_without_progress = 0;
            }
            else if (++steps_without_progress >= max_steps_without_progress
                     || !std::isfinite(residual_norm)) {
                OPM_THROW_NOLOG(NumericalProblem,
                                fmt::format(" Newton composition update did not halve its "
                                            "residual in {} iterations (residual {})",
                                            max_steps_without_progress, residual_norm));
            }

            jacobian.solve(newton_step, residual);
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                liquid_composition[compIdx] -= newton_step[compIdx];
                vapour_composition[compIdx] -= newton_step[compIdx + numComponents];
            }
            liquid_fraction -= newton_step[num_equations - 1];
            residual_norm = evaluate();
            ++iteration;
        }
        if (verbosity >= 1) {
            fmt::memory_buffer buf;
            for (unsigned equation_idx = 0; equation_idx < num_equations; ++equation_idx) {
                for (unsigned unknown_idx = 0; unknown_idx < num_unknowns; ++unknown_idx) {
                    fmt::format_to(
                        std::back_inserter(buf), " {}", jacobian[equation_idx][unknown_idx]);
                }
                fmt::format_to(std::back_inserter(buf), "\n");
            }
            OpmLog::debug(fmt::to_string(buf));
        }
        if (!converged) {
            OPM_THROW_NOLOG(NumericalProblem,
                            fmt::format(" Newton composition update did not converge "
                                        "within {} iterations",
                                        max_iterations));
        }

        // A small equation residual does not enforce composition bounds. A
        // finite L outside [0, 1] remains a valid negative-flash result.
        if (!std::isfinite(Opm::getValue(liquid_fraction))
            || !phaseCompositionsArePhysical_(newton_state)) {
            OPM_THROW_NOLOG(
                NumericalProblem,
                fmt::format(" Newton composition update converged to a non-physical "
                            "composition with L = {}",
                            Opm::getValue(liquid_fraction)));
        }

        // Copy the converged unknowns to the caller and the explicit K/L outputs.
        for (unsigned idx = 0; idx < numComponents; ++idx) {
            const auto liquid_mole_fraction
                = Opm::getValue(newton_state.moleFraction(oilPhaseIdx, idx));
            fluid_state.setMoleFraction(FluidSystem::oilPhaseIdx, idx, liquid_mole_fraction);
            const auto vapour_mole_fraction
                = Opm::getValue(newton_state.moleFraction(gasPhaseIdx, idx));
            fluid_state.setMoleFraction(FluidSystem::gasPhaseIdx, idx, vapour_mole_fraction);
            const auto equilibrium_ratio = Opm::getValue(newton_state.K(idx));
            fluid_state.setKvalue(idx, equilibrium_ratio);
            K[idx] = equilibrium_ratio;
        }
        L = Opm::getValue(liquid_fraction);
        fluid_state.setLvalue(L);
        return {jacobian, jacobianState_(fluid_state)};
    }

    /*!
     * \brief The residual and Jacobian of newtonComposition_() at the current
     *        state, whose values carry derivatives in the unknowns.
     */
    template <typename FlashFluidState,
              typename ComponentVector,
              std::size_t num_unknowns,
              std::size_t num_equations>
    static void assembleNewton_(const FlashFluidState& fluid_state,
                                const ComponentVector& overall_composition,
                                Dune::FieldMatrix<double, num_equations, num_unknowns>& jacobian,
                                Dune::FieldVector<double, num_equations>& residual)
    {
        using Eval = DenseAd::Evaluation<double, num_unknowns>;
        std::vector<Eval> liquid_composition(numComponents);
        std::vector<Eval> vapour_composition(numComponents);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            liquid_composition[compIdx] = fluid_state.moleFraction(oilPhaseIdx, compIdx);
            vapour_composition[compIdx] = fluid_state.moleFraction(gasPhaseIdx, compIdx);
        }
        const Eval& liquid_fraction = fluid_state.L();
        jacobian = 0.;
        residual = 0.;
        // Fugacities are of the order of the pressure, so their rows are
        // scaled by it: an absolute tolerance on raw differences sits at the
        // round-off floor of the EoS at reservoir pressure.
        const double pref = Opm::getValue(fluid_state.pressure(oilPhaseIdx));
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            {
                // z - L*x - (1-L) * y
                const auto component_balance = -overall_composition[compIdx]
                    + liquid_fraction * liquid_composition[compIdx]
                    + (1 - liquid_fraction) * vapour_composition[compIdx];
                residual[compIdx] = Opm::getValue(component_balance);
                for (unsigned unknown_idx = 0; unknown_idx < num_unknowns; ++unknown_idx) {
                    jacobian[compIdx][unknown_idx] = component_balance.derivative(unknown_idx);
                }
            }

            {
                // (f_liquid - f_vapor) / p = 0
                const auto fugacity_difference = (fluid_state.fugacity(oilPhaseIdx, compIdx)
                    - fluid_state.fugacity(gasPhaseIdx, compIdx)) / pref;
                residual[compIdx + numComponents] = Opm::getValue(fugacity_difference);
                for (unsigned unknown_idx = 0; unknown_idx < num_unknowns; ++unknown_idx) {
                    jacobian[compIdx + numComponents][unknown_idx]
                        = fugacity_difference.derivative(unknown_idx);
                }
            }
        }
        // sum(x) - sum(y) = 0
        Eval liquid_composition_sum = 0.;
        Eval vapour_composition_sum = 0.;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            liquid_composition_sum += liquid_composition[compIdx];
            vapour_composition_sum += vapour_composition[compIdx];
        }
        const auto composition_closure = liquid_composition_sum - vapour_composition_sum;
        residual[num_equations - 1] = Opm::getValue(composition_closure);
        for (unsigned unknown_idx = 0; unknown_idx < num_unknowns; ++unknown_idx) {
            jacobian[num_equations - 1][unknown_idx] = composition_closure.derivative(unknown_idx);
        }
    }

    //! The values the composition Jacobian depends on.
    template <class FlashFluidState>
    static JacobianState jacobianState_(const FlashFluidState& fluid_state)
    {
        JacobianState state{};
        for (unsigned compIdx = 0; compIdx < numMiscibleComponents; ++compIdx) {
            state[compIdx] = Opm::getValue(fluid_state.moleFraction(oilPhaseIdx, compIdx));
            state[compIdx + numMiscibleComponents]
                = Opm::getValue(fluid_state.moleFraction(gasPhaseIdx, compIdx));
        }
        // L is the last flash unknown, as in newtonComposition_().
        constexpr std::size_t liquid_fraction_idx = numFlashUnknowns - 1;
        state[liquid_fraction_idx] = Opm::getValue(fluid_state.L());
        state[liquid_fraction_idx + 1] = Opm::getValue(fluid_state.pressure(oilPhaseIdx));
        state[liquid_fraction_idx + 2] = Opm::getValue(fluid_state.pressure(gasPhaseIdx));
        state[liquid_fraction_idx + 3] = Opm::getValue(fluid_state.temperature(0));
        return state;
    }

    /*!
     * \brief The Jacobian of the newtonComposition_() residuals in
     *        \f$\mathbf u = (\mathbf x, \mathbf y, L)\f$ at a scalar state.
     */
    template <typename FlashFluidStateScalar>
    static Dune::FieldMatrix<Scalar, numFlashUnknowns, numFlashUnknowns>
    compositionJacobian_(const FlashFluidStateScalar& fluid_state_scalar, const EOSType& eos_type)
    {
        constexpr std::size_t num_equations = numFlashUnknowns;
        constexpr std::size_t num_flash_unknowns = numFlashUnknowns;
        using UnknownEval = Opm::DenseAd::Evaluation<double, num_flash_unknowns>;
        using ScalarComponentVector = Dune::FieldVector<double, numComponents>;
        using UnknownFluidState = Opm::CompositionalFluidState<UnknownEval, FluidSystem>;

        UnknownFluidState unknown_derivative_state;
        // The overall composition enters the residual but not its Jacobian.
        ScalarComponentVector overall_composition;
        for (unsigned  comp_idx = 0; comp_idx < numComponents; ++comp_idx) {
            overall_composition[comp_idx]
                = Opm::getValue(fluid_state_scalar.moleFraction(comp_idx));
        }
        for (unsigned comp_idx = 0; comp_idx < numComponents; ++comp_idx) {
            const auto liquid_mole_fraction
                = UnknownEval(fluid_state_scalar.moleFraction(oilPhaseIdx, comp_idx), comp_idx);
            unknown_derivative_state.setMoleFraction(oilPhaseIdx, comp_idx, liquid_mole_fraction);
            const unsigned idx = comp_idx + numComponents;
            const auto vapour_mole_fraction
                = UnknownEval(fluid_state_scalar.moleFraction(gasPhaseIdx, comp_idx), idx);
            unknown_derivative_state.setMoleFraction(gasPhaseIdx, comp_idx, vapour_mole_fraction);
            unknown_derivative_state.setKvalue(comp_idx,
                                               vapour_mole_fraction / liquid_mole_fraction);
        }
        const auto liquid_fraction = UnknownEval(fluid_state_scalar.L(), num_flash_unknowns - 1);
        unknown_derivative_state.setLvalue(liquid_fraction);
        unknown_derivative_state.setPressure(oilPhaseIdx, fluid_state_scalar.pressure(oilPhaseIdx));
        unknown_derivative_state.setPressure(gasPhaseIdx, fluid_state_scalar.pressure(gasPhaseIdx));
        unknown_derivative_state.setTemperature(fluid_state_scalar.temperature(0));

        using UnknownParamCache = typename FluidSystem::template ParameterCache<UnknownEval>;
        UnknownParamCache unknown_param_cache(eos_type);
        for (unsigned phase_idx = 0; phase_idx < numMisciblePhases; ++phase_idx) {
            unknown_param_cache.updatePhase(unknown_derivative_state, phase_idx);
            updateFugacityCoefficients_(unknown_derivative_state, unknown_param_cache, phase_idx);
        }

        using UnknownResidual = Dune::FieldVector<Scalar, num_equations>;
        using UnknownJacobian = Dune::FieldMatrix<Scalar, num_equations, num_flash_unknowns>;
        UnknownResidual unknown_residual;
        UnknownJacobian unknown_jacobian;

        assembleNewton_<UnknownFluidState,
                        ScalarComponentVector,
                        num_flash_unknowns,
                        num_equations>(
            unknown_derivative_state, overall_composition, unknown_jacobian, unknown_residual);
        return unknown_jacobian;
    }

    /*!
     * \brief Gives the converged flash result the derivatives of the caller.
     *
     * The flash is solved on scalars; the compositions it returns are then
     * functions of the caller's variables through the equations they
     * satisfy. Single-phase results need no work, as \f$x = y = z\f$.
     */
    template <typename FlashFluidStateScalar, typename FluidState>
    static void updateDerivatives_(const FlashFluidStateScalar& fluid_state_scalar,
                                   FluidState& fluid_state,
                                   const EOSType& eos_type,
                                   bool is_single_phase,
                                   const std::optional<NewtonJacobian>& newton_jacobian)
    {
        if (!is_single_phase) {
            updateDerivativesTwoPhase_(fluid_state_scalar, fluid_state, eos_type, newton_jacobian);
        } else {
            updateDerivativesSinglePhase_(fluid_state_scalar, fluid_state);
        }
    }

    /*!
     * \brief Derivatives of a two-phase result by the implicit function theorem.
     *
     * With the flash residuals \f$F(\mathbf u, \mathbf s) = 0\f$ in the
     * unknowns \f$\mathbf u = (\mathbf x, \mathbf y, L)\f$ and the state
     * variables \f$\mathbf s = (p, [T,] \mathbf z)\f$,
     *
     * \f[ \frac{\partial \mathbf u}{\partial \mathbf s}
     *   = - \left( \frac{\partial F}{\partial \mathbf u} \right)^{-1}
     *       \frac{\partial F}{\partial \mathbf s}, \f]
     *
     * and the chain rule through the caller's derivatives of \f$\mathbf s\f$
     * gives those of \f$\mathbf u\f$. A Newton solve that converged at exactly
     * this state already assembled \f$\partial F / \partial \mathbf u\f$.
     */
    template <typename FlashFluidStateScalar, typename FluidState>
    static void updateDerivativesTwoPhase_(const FlashFluidStateScalar& fluid_state_scalar,
                                           FluidState& fluid_state,
                                           const EOSType& eos_type,
                                           const std::optional<NewtonJacobian>& newton_jacobian)
    {
        using InputEval = typename FluidState::ValueType;
        using ComponentVector = Dune::FieldVector<InputEval, numComponents>;
        ComponentVector z;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            z[compIdx] = fluid_state.moleFraction(compIdx);
        }

        constexpr std::size_t num_equations = numMisciblePhases * numMiscibleComponents + 1;
        constexpr std::size_t num_state_variables
            = isThermal ? numComponents + 2 : numComponents + 1;
        using StateEval = Opm::DenseAd::Evaluation<double, num_state_variables>;
        using StateComponentVector = Dune::FieldVector<StateEval, numComponents>;
        using StateFluidState = Opm::CompositionalFluidState<StateEval, FluidSystem>;

        // Evaluate F with s = (p, [T], z) as independent AD variables while
        // holding the converged flash unknowns u = (x, y, L) fixed.
        StateFluidState state_derivative_state;
        StateComponentVector state_composition;
        const StateEval state_pressure
            = StateEval(fluid_state_scalar.pressure(FluidSystem::oilPhaseIdx), 0);
        state_derivative_state.setPressure(FluidSystem::oilPhaseIdx, state_pressure);
        state_derivative_state.setPressure(FluidSystem::gasPhaseIdx, state_pressure);

        if constexpr (isThermal) {
            // set the temperature with derivatives
            const StateEval state_temperature = StateEval(fluid_state_scalar.temperature(0), 1);
            state_derivative_state.setTemperature(state_temperature);
        } else {
            // set the temperature as a scalar
            state_derivative_state.setTemperature(Opm::getValue(fluid_state_scalar.temperature(0)));
        }

        // composition variable offset: 2 for thermal (pressure=0, temperature=1, z starts at 2)
        //                               1 for non-thermal (pressure=0, z starts at 1)
        constexpr unsigned composition_offset = isThermal ? 2 : 1;
        for (unsigned idx = 0; idx < numComponents; ++idx) {
            state_composition[idx] = StateEval(Opm::getValue(z[idx]), idx + composition_offset);
        }
        // The converged flash unknowns are constants in this derivative pass.
        for (unsigned idx = 0; idx < numComponents; ++idx) {
            const auto liquid_mole_fraction = fluid_state_scalar.moleFraction(oilPhaseIdx, idx);
            state_derivative_state.setMoleFraction(
                FluidSystem::oilPhaseIdx, idx, liquid_mole_fraction);
            const auto vapour_mole_fraction = fluid_state_scalar.moleFraction(gasPhaseIdx, idx);
            state_derivative_state.setMoleFraction(
                FluidSystem::gasPhaseIdx, idx, vapour_mole_fraction);
            const auto equilibrium_ratio = fluid_state_scalar.K(idx);
            state_derivative_state.setKvalue(idx, equilibrium_ratio);
        }
        const auto scalar_liquid_fraction = fluid_state_scalar.L();
        state_derivative_state.setLvalue(scalar_liquid_fraction);
        // Saturations do not enter the flash residual and remain unset here.
        using StateParamCache = typename FluidSystem::template ParameterCache<StateEval>;
        StateParamCache state_param_cache(eos_type);
        for (unsigned phase_idx = 0; phase_idx < numMisciblePhases; ++phase_idx) {
            state_param_cache.updatePhase(state_derivative_state, phase_idx);
            updateFugacityCoefficients_(state_derivative_state, state_param_cache, phase_idx);
        }

        using StateResidual = Dune::FieldVector<Scalar, num_equations>;
        using StateJacobian = Dune::FieldMatrix<Scalar, num_equations, num_state_variables>;
        StateJacobian state_jacobian;
        StateResidual state_residual;

        // This gives the state-variable Jacobian dF/ds.
        assembleNewton_<StateFluidState, StateComponentVector, num_state_variables, num_equations>(
            state_derivative_state, state_composition, state_jacobian, state_residual);

        // dF/du, from the Newton solve if it converged at exactly this state.
        using UnknownJacobian = Dune::FieldMatrix<Scalar, num_equations, numFlashUnknowns>;
        UnknownJacobian unknown_jacobian
            = (newton_jacobian && newton_jacobian->state == jacobianState_(fluid_state_scalar))
            ? newton_jacobian->jacobian
            : compositionJacobian_(fluid_state_scalar, eos_type);

        // DUNE 2.6 cannot solve directly with a matrix right-hand side, so form
        // the inverse explicitly before multiplying the state Jacobian.
        // StateJacobian xx;
        // unknown_jacobian.solve(xx, state_jacobian);
        unknown_jacobian.invert();
        state_jacobian.template leftmultiply<UnknownJacobian>(unknown_jacobian);
        // state_jacobian now stores (dF/du)^-1 dF/ds; the minus sign is
        // applied while transferring each sensitivity below.

        ComponentVector liquid_composition(numComponents);
        ComponentVector vapour_composition(numComponents);
        InputEval liquid_fraction_with_derivatives = scalar_liquid_fraction;

        // Apply du/ds = -(dF/du)^-1 dF/ds, then compose it with the
        // caller's derivatives of s.

        const auto liquid_pressure = fluid_state.pressure(FluidSystem::oilPhaseIdx);
        const auto vapour_pressure = fluid_state.pressure(FluidSystem::gasPhaseIdx);
        // at the moment, the temperature is the same for both phases
        // input_temperature is not used for non-thermal case
        [[maybe_unused]] const auto& input_temperature = fluid_state.temperature(0);

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            liquid_composition[compIdx]
                = fluid_state_scalar.moleFraction(FluidSystem::oilPhaseIdx, compIdx);
            vapour_composition[compIdx]
                = fluid_state_scalar.moleFraction(FluidSystem::gasPhaseIdx, compIdx);
        }

        // The two pressures are equal today. Keep both inputs explicit so this
        // chain rule can accommodate capillary pressure later.

        constexpr std::size_t num_derivatives = InputEval::numVars;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::vector<double> derivatives(num_derivatives, 0.);
            // derivatives from P
            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                derivatives[idx] = -state_jacobian[compIdx][0] * liquid_pressure.derivative(idx);
            }
            // derivatives from T (only for thermal case)
            if constexpr (isThermal) {
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    derivatives[idx]
                        += -state_jacobian[compIdx][1] * input_temperature.derivative(idx);
                }
            }

            for (unsigned cIdx = 0; cIdx < numComponents; ++cIdx) {
                const double composition_sensitivity
                    = -state_jacobian[compIdx][cIdx + composition_offset];
                const auto& overall_mole_fraction = z[cIdx];
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    derivatives[idx]
                        += composition_sensitivity * overall_mole_fraction.derivative(idx);
                }
            }
            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                liquid_composition[compIdx].setDerivative(idx, derivatives[idx]);
            }
            // handling y
            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                derivatives[idx]
                    = -state_jacobian[compIdx + numComponents][0] * vapour_pressure.derivative(idx);
            }
            // derivatives from T (only for thermal case)
            if constexpr (isThermal) {
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    derivatives[idx] += -state_jacobian[compIdx + numComponents][1]
                        * input_temperature.derivative(idx);
                }
            }
            for (unsigned cIdx = 0; cIdx < numComponents; ++cIdx) {
                const double composition_sensitivity
                    = -state_jacobian[compIdx + numComponents][cIdx + composition_offset];
                const auto& overall_mole_fraction = z[cIdx];
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    derivatives[idx]
                        += composition_sensitivity * overall_mole_fraction.derivative(idx);
                }
            }
            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                vapour_composition[compIdx].setDerivative(idx, derivatives[idx]);
            }

            // handling derivatives of L
            std::vector<double> liquid_fraction_derivatives(num_derivatives, 0.);
            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                liquid_fraction_derivatives[idx]
                    = -state_jacobian[2 * numComponents][0] * vapour_pressure.derivative(idx);
            }
            // derivatives from T (only for thermal case)
            if constexpr (isThermal) {
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    liquid_fraction_derivatives[idx] += -state_jacobian[2 * numComponents][1]
                        * input_temperature.derivative(idx);
                }
            }
            for (unsigned cIdx = 0; cIdx < numComponents; ++cIdx) {
                const double composition_sensitivity
                    = -state_jacobian[2 * numComponents][cIdx + composition_offset];
                const auto& overall_mole_fraction = z[cIdx];
                for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                    liquid_fraction_derivatives[idx]
                        += composition_sensitivity * overall_mole_fraction.derivative(idx);
                }
            }

            for (unsigned idx = 0; idx < num_derivatives; ++idx) {
                liquid_fraction_with_derivatives.setDerivative(idx,
                                                               liquid_fraction_derivatives[idx]);
            }
        }

        // set up the mole fractions
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluid_state.setMoleFraction(
                FluidSystem::oilPhaseIdx, compIdx, liquid_composition[compIdx]);
            fluid_state.setMoleFraction(
                FluidSystem::gasPhaseIdx, compIdx, vapour_composition[compIdx]);
        }
        fluid_state.setLvalue(liquid_fraction_with_derivatives);
    } //end updateDerivativesTwoPhase

    /*!
     * \brief A single-phase result: both phases take the mixture composition,
     *        and the liquid fraction carries no derivatives.
     */
    template <typename FlashFluidStateScalar, typename FluidState>
    static void updateDerivativesSinglePhase_(const FlashFluidStateScalar& fluid_state_scalar,
                                              FluidState& fluid_state)
    {
        using InputEval = typename FluidState::ValueType;
        // Conversion from the scalar result gives the phase label zero derivatives.
        InputEval liquid_fraction = fluid_state_scalar.L();

        // A single phase uses the overall composition under both phase labels.
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluid_state.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, fluid_state.moleFraction(compIdx) );
            fluid_state.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, fluid_state.moleFraction(compIdx) );
        }
        fluid_state.setLvalue(liquid_fraction);
    } //end updateDerivativesSinglePhase

    /*!
     * \brief The Gibbs energy of a two-phase split, per mole of feed and
     *        divided by RT, up to a constant,
     *
     * \f[ L \sum_i x_i \ln(x_i \varphi_i^{L})
     *     + (1 - L) \sum_i y_i \ln(y_i \varphi_i^{V}), \f]
     *
     * from the compositions and fugacity coefficients in the fluid state.
     * Absent components do not contribute.
     */
    template <class FlashFluidState, class ComponentVector>
    static typename ComponentVector::field_type
    splitGibbsEnergy_(const FlashFluidState& fluid_state,
                      const typename ComponentVector::field_type& L,
                      const ComponentVector& z)
    {
        typename ComponentVector::field_type gibbs = 0.;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!(Opm::getValue(z[compIdx]) > 0.)) {
                continue;
            }
            const auto x = fluid_state.moleFraction(oilPhaseIdx, compIdx);
            const auto y = fluid_state.moleFraction(gasPhaseIdx, compIdx);
            gibbs += L * x * log(x * fluid_state.fugacityCoefficient(oilPhaseIdx, compIdx))
                + (1. - L) * y * log(y * fluid_state.fugacityCoefficient(gasPhaseIdx, compIdx));
        }
        return gibbs;
    }

    /*!
     * \brief The sum of the steps that successive substitution has still to
     *        take, from its last three steps (newest first).
     *
     * Near convergence the steps follow the linear recurrence of the two
     * dominant eigenvalues of the iteration,
     * \f$ s_{n+1} + b_1 s_n + b_2 s_{n-1} = 0 \f$, with \f$b_1, b_2\f$ fitted
     * to the last three steps by least squares. Summing the recurrence gives
     *
     * \f[ \sum_{m \ge 1} s_{n+m}
     *     = -\frac{(b_1 + b_2)\, s_n + b_2\, s_{n-1}}{1 + b_1 + b_2}. \f]
     *
     * Where the two older steps are parallel, one eigenvalue
     * \f$\lambda = s_n \cdot s_n / s_{n-1} \cdot s_n\f$ is used:
     * \f$ \lambda s_n / (1 - \lambda) \f$.
     *
     * \return false where neither gives a finite sum of a contracting
     *         iteration
     */
    template <class ComponentVector>
    static bool extrapolateSubstitution_(const std::array<ComponentVector, 3>& steps,
                                         ComponentVector& extrapolation)
    {
        using FieldType = typename ComponentVector::field_type;
        const auto& [newest, middle, oldest] = steps;
        const FieldType mm = middle * middle;
        const FieldType mo = middle * oldest;
        const FieldType oo = oldest * oldest;
        const FieldType determinant = mm * oo - mo * mo;
        bool found = false;
        if (determinant > 1e-10 * mm * oo) {
            const FieldType b1 = (-(middle * newest) * oo + mo * (oldest * newest)) / determinant;
            const FieldType b2 = (-mm * (oldest * newest) + mo * (middle * newest)) / determinant;
            const FieldType denominator = 1. + b1 + b2;
            if (denominator > 1e-12) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    extrapolation[compIdx]
                        = -((b1 + b2) * newest[compIdx] + b2 * middle[compIdx]) / denominator;
                }
                found = true;
            }
        }
        if (!found) {
            const FieldType eigenvalue = (newest * newest) / (middle * newest);
            if (!(eigenvalue > 0. && eigenvalue < 1.)) {
                return false;
            }
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                extrapolation[compIdx] = newest[compIdx] * eigenvalue / (1. - eigenvalue);
            }
        }
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!std::isfinite(Opm::getValue(extrapolation[compIdx]))) {
                return false;
            }
        }
        return true;
    }

    /*!
     * \brief Successive substitution for the two-phase compositions.
     *
     * Each step recomputes the compositions from the current ratios, then
     * updates the ratios by the fugacity ratio,
     *
     * \f[ K_i \leftarrow K_i \frac{f_i^{L}}{f_i^{V}}, \f]
     *
     * until \f$\| f^{L} / f^{V} - 1 \|\f$ is below tolerance. As a
     * conditioner for Newton it takes a few steps.
     *
     * The steps in \f$\ln K\f$ shrink by a nearly constant factor, which tends
     * to one near a critical point, where the two phases become alike: there
     * plain substitution needs thousands of steps. With \p extrapolate
     * (PTFlashMethod::SsiAccelerated), every substitutionAccelerationInterval
     * steps the remaining steps are summed (general dominant eigenvalue
     * method): the last three steps give
     * the two dominant eigenvalues of the iteration, and their geometric
     * series extrapolates \f$\ln K\f$. Substitution never increases the Gibbs
     * energy of the split, see splitGibbsEnergy_(). The extrapolated ratios
     * are kept only if they lower it below its value at the last evaluated
     * ratios; otherwise the iteration continues from the plain step.
     */
    template <class FlashFluidState, class ComponentVector>
    static bool successiveSubstitutionComposition_(ComponentVector& K,
                                                   typename ComponentVector::field_type& L,
                                                   FlashFluidState& fluid_state,
                                                   const ComponentVector& z,
                                                   const bool run_newton_afterwards,
                                                   const Scalar flash_tolerance,
                                                   const EOSType& eos_type,
                                                   const int verbosity,
                                                   const bool extrapolate = false)
    {
        // Limit conditioning passes before Newton; standalone substitution runs
        // longer. With the extrapolation, the limit only stops a split that
        // does not converge.
        const int max_iterations = run_newton_afterwards ? maxHybridSubstitutionSteps : (extrapolate ? 1000 : 100);

        // Print initial guess
        if (verbosity >= 1) {
            OpmLog::debug(fmt::format("Initial guess: K = [{}] and L = {}", fmt::join(K, " "), L));
        }

        if (verbosity == 2 || verbosity == 4) {
            // Print header
            const int fugacity_width = (numComponents * 12) / 2;
            const int convergence_width = fugacity_width + 7;
            OpmLog::debug(fmt::format("{:>10}{:>{}}{:>{}}",
                                      "Iteration",
                                      "fL/fV",
                                      fugacity_width,
                                      "norm2(fL/fv-1)",
                                      convergence_width));
        }
        //
        // Successive substitution loop
        //
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!(Opm::getValue(z[compIdx]) > 0.)) {
                K[compIdx] = 1.;
            }
        }
        // The last three plain steps in ln K, newest first, and the plain
        // point to return to when an extrapolation does not lower the Gibbs
        // energy.
        using FieldType = typename ComponentVector::field_type;
        std::array<ComponentVector, 3> steps{};
        int plain_steps = 0;
        bool extrapolated = false;
        ComponentVector plain_K = K;
        FieldType plain_L = L;
        FieldType gibbs_before = 0.;
        // One parameter cache for all steps: temperature and pressure are
        // fixed, so it computes the pure-component parameters once.
        using ParamCache = typename FluidSystem::template ParameterCache<typename FlashFluidState::ValueType>;
        ParamCache paramCache(eos_type);
        Scalar previous_residual = 0.;
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            OPM_FLASH_COUNT(ssiIterations);
            // Compute (normalized) liquid and vapor mole fractions
            computeLiquidVapor_(fluid_state, L, K, z);

            // Calculate fugacity coefficient
            for (int phaseIdx=0; phaseIdx<numMisciblePhases; ++phaseIdx){
                paramCache.updatePhase(fluid_state, phaseIdx);
                updateFugacityCoefficients_(fluid_state, paramCache, phaseIdx);
            }

            // Calculate fugacity ratio
            ComponentVector fugacity_ratio;
            ComponentVector fugacity_residual;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (!(Opm::getValue(z[compIdx]) > 0.)) {
                    fugacity_ratio[compIdx] = 1.;
                    fugacity_residual[compIdx] = 0.;
                    continue;
                }
                fugacity_ratio[compIdx] = fluid_state.fugacity(oilPhaseIdx, compIdx)
                    / fluid_state.fugacity(gasPhaseIdx, compIdx);
                fugacity_residual[compIdx] = fugacity_ratio[compIdx] - 1.0;
            }

            // Print iteration info
            if (verbosity >= 2) {
                constexpr int precision = 5;
                constexpr int fugacity_width = precision + 3;
                constexpr int convergence_width = precision + 9;
                OpmLog::debug(fmt::format("{:>5}{:>{}.{}f}{:>{}.{}e}",
                                          iteration,
                                          fmt::join(fugacity_ratio, " "),
                                          fugacity_width,
                                          precision,
                                          fugacity_residual.two_norm(),
                                          convergence_width,
                                          precision));
            }

            const bool split_is_interior = L > 0. && L < 1.;
            const FieldType gibbs = extrapolate && split_is_interior ? splitGibbsEnergy_(fluid_state, L, z) : FieldType{0.};
            if (extrapolated) {
                extrapolated = false;
                if (!(split_is_interior && gibbs < gibbs_before)) {
                    // Continue from the plain step instead.
                    K = plain_K;
                    L = plain_L;
                    plain_steps = 0;
                    continue;
                }
            }

            // Check convergence
            if (fugacity_residual.two_norm() < flash_tolerance) {
                // Print info
                if (verbosity >= 1) {
                    std::vector<typename ComponentVector::field_type> liquid_composition(
                        numComponents),
                        vapour_composition(numComponents);
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                        liquid_composition[compIdx]
                            = fluid_state.moleFraction(oilPhaseIdx, compIdx);
                        vapour_composition[compIdx]
                            = fluid_state.moleFraction(gasPhaseIdx, compIdx);
                    }
                    OpmLog::debug(fmt::format("Solution converged to the following result :\n"
                                              "x = [{}]\n"
                                              "y = [{}]\n"
                                              "K = [{}]\n"
                                              "L = {}",
                                              fmt::join(liquid_composition, " "),
                                              fmt::join(vapour_composition, " "),
                                              fmt::join(K, " "),
                                              L));
                }
                return true;
            }
            //  If convergence is not met, K is updated in a successive substitution manner
            else {
                // Update K
                ComponentVector step(0.);
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    if (extrapolate) {
                        step[compIdx] = log(fugacity_ratio[compIdx]);
                    }
                    K[compIdx] *= fugacity_ratio[compIdx];
                }

                // Solve Rachford-Rice to get L from updated K
                L = RachfordRice::solve(K, z, L, 0);

                // The hand-over to Newton, see substitutionStepsBeforeNewton.
                if (run_newton_afterwards) {
                    const Scalar residual = Opm::getValue(fugacity_residual.two_norm());
                    const Scalar rate = iteration > 0 ? residual / previous_residual : Scalar{1.};
                    const Scalar predicted_steps = rate > 0. && rate < 1.
                        ? std::log(flash_tolerance / residual) / std::log(rate)
                        : std::numeric_limits<Scalar>::infinity();
                    previous_residual = residual;
                    if (!(predicted_steps <= maxPredictedSubstitutionSteps)
                        && iteration + 1 >= substitutionStepsBeforeNewton) {
                        break;
                    }
                }

                steps = {step, steps[0], steps[1]};
                ++plain_steps;
                ComponentVector extrapolation;
                if (extrapolate && split_is_interior && plain_steps >= substitutionAccelerationInterval
                    && extrapolateSubstitution_(steps, extrapolation)) {
                    plain_K = K;
                    plain_L = L;
                    gibbs_before = gibbs;
                    for (int compIdx=0; compIdx<numComponents; ++compIdx){
                        K[compIdx] *= exp(extrapolation[compIdx]);
                    }
                    try {
                        L = RachfordRice::solve(K, z, L, 0);
                        extrapolated = true;
                    }
                    catch (const NumericalProblem&) {
                        K = plain_K;
                        L = plain_L;
                    }
                    plain_steps = 0;
                }
            }
        }
        // did not get converged. check whether we will do more newton later afterward
        {
            const std::string message
                = fmt::format("Successive substitution composition update did not "
                              "converge within {} iterations.",
                              max_iterations);
            if (!run_newton_afterwards) {
                OPM_THROW_NOLOG(NumericalProblem, message);
            } else if (verbosity > 0) {
                OpmLog::debug(message);
            }
        }

        return false;
    }

};//end PTFlash

} // namespace Opm

#endif
