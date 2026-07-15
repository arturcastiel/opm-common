// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA.

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
 *
 * \brief Isenthalpic (pressure-enthalpy, "P-H") flash: given the pressure and
 *        overall composition carried by a fluid state and a specified molar
 *        enthalpy, find the temperature at which the flashed mixture attains
 *        that enthalpy, and leave the state flashed at the solution.
 *
 * The isothermal flash (PTFlash) is reused unmodified as the inner machinery:
 * at each trial temperature the mixture is flashed and its molar enthalpy
 * (MixtureEnthalpy, caloric or EoS-departure model) is compared against the
 * specified value; the temperature is bracketed and root-found.
 *
 * Contract notes:
 * - The specified enthalpy is MOLAR [J/mol] and MUST be expressed against the
 *   same reference datum as the enthalpy model (IdealGasCaloricData: H(T0) = 0) — a
 *   datum mismatch manifests as a systematically wrong temperature, not as a
 *   solver failure.
 * - Values only: after solve() the fluid state's AD derivatives (if any) are
 *   those of the final isothermal flash at the CONVERGED temperature held
 *   fixed; the sensitivity of the solution temperature itself
 *   (dT/d{H,p,z}) is NOT propagated.
 * - With the caloric model the mixture enthalpy is independent of the phase
 *   split (H = sum_i z_i h_i(T)) and strictly increasing in T, so the root is
 *   unique and smooth. With the departure model H(T) remains continuous and
 *   increasing for mixtures, with slope kinks at bubble/dew crossings and a
 *   possible jump if the single-phase label flips along the sweep — the
 *   bracketing solver converges across these, but pathological cases near
 *   phase boundaries are the province of a dedicated boundary locator, not
 *   of this solver.
 */
#ifndef OPM_PH_FLASH_HPP
#define OPM_PH_FLASH_HPP

#include <opm/material/constraintsolvers/IdealGasCaloricData.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Opm {

/*!
 * \brief Run configuration of the isenthalpic flash: the enthalpy model and
 *        the temperature search.
 *
 * The temperature bracket must contain the solution; a specified enthalpy
 * outside [H(tempMin), H(tempMax)] makes solve() return false. The bracket
 * bounds are full inner-flash evaluations, so they must lie within the
 * isothermal flash's robust operating envelope: at extreme temperatures the
 * Wilson-seeded stability/split machinery can fail to converge (and throws).
 * The defaults are chosen inside that envelope and near the cp-correlation
 * validity range; widen them deliberately, not by default.
 */
template <class Scalar, int numComponents>
struct PHFlashConfig {
    //! per-component heat-capacity polynomials. MUST be populated by the
    //! caller: the default-constructed table is all-zero and unusable —
    //! solve() rejects it (returns false) rather than "solving" H = 0.
    CpTable<Scalar, numComponents> cpTable{};
    //! enthalpy reference datum [K]. The specified enthalpy handed to
    //! solve() MUST be expressed against this same datum — a mismatch
    //! produces a systematically wrong temperature, not a solver failure.
    Scalar refTemperature = IdealGasCaloricData<Scalar>::referenceTemperature();
    //! [K] bracket bounds. Both ends are full inner-flash evaluations, so
    //! they must lie inside the isothermal flash's robust envelope AND the
    //! cp-correlation validity window (250-600 K): near-critical feeds fail
    //! to flash at extreme bounds (e.g. a 99%-methane feed at 200 K), which
    //! reports as "no solution" although the root is interior. The defaults
    //! are the field-validated window every in-tree consumer uses.
    Scalar tempMin = 270.0;
    Scalar tempMax = 460.0;
    //! convergence tolerance of the bracketing solver's TEMPERATURE
    //! interval [K]
    Scalar tolerance = 1e-6;
    //! acceptance tolerance on the ENTHALPY residual |H(T*) - hSpec|
    //! [J/mol], checked after the root-find. Where the departure-model
    //! enthalpy jumps at a phase-label flip, the bracketing loop can
    //! converge its interval onto the discontinuity although no root
    //! exists; this check is what turns that case into an honest failure.
    Scalar enthalpyTolerance = 1.0;
    int maxIterations = 100;
    EnthalpyModel model = EnthalpyModel::caloric;
};

/*!
 * \brief The isenthalpic (P-H) flash: solves H(p, T, z) = hSpec for the
 *        temperature by a bracketed root-find over the unmodified isothermal
 *        flash, using the mixture-enthalpy model selected in the config.
 *
 * The EnthalpyCalc template parameter is the enthalpy-provider seam; any
 * substitute must supply a static mixtureEnthalpy(fluidState, cpTable,
 * refTemperature, eosType, model) with MixtureEnthalpy's semantics (molar
 * enthalpy of a flashed, L-consistent state).
 */
template <class Scalar, class FluidSystem,
          class EnthalpyCalc = MixtureEnthalpy<Scalar, FluidSystem>>
struct PHFlash {
    static constexpr int numComponents = FluidSystem::numComponents;

    using EOSType = CompositionalConfig::EOSType;
    using Config = PHFlashConfig<Scalar, numComponents>;

    /*!
     * \brief Solve H(p, T, z) = hSpec for T and flash the state at the
     *        solution.
     *
     * \param fluidState carries the pressure and overall mole fractions on
     *        entry; on success it holds the flashed state at the solution
     *        temperature. The inner isothermal flash input contract (K-value
     *        and L seeding) is applied internally at every trial temperature.
     * \param hSpec specified molar mixture enthalpy [J/mol] against the
     *        config's reference datum
     * \param cfg the enthalpy model and temperature-search configuration
     * \param twoPhaseMethod inner isothermal flash iteration scheme ("ssi",
     *        "newton" or "ssi+newton" — passed through verbatim)
     * \param ptTolerance convergence tolerance of the inner isothermal flash
     *        (its fugacity-ratio residual)
     * \param eosType cubic equation-of-state variant for the inner flash and
     *        the departure enthalpy
     * \param verbosity inner-flash verbosity (passed through)
     * \return true on success — the root-find converged AND the final
     *         enthalpy residual |H(T*) - hSpec| lies within
     *         cfg.enthalpyTolerance; false when no solution is found on the
     *         bracket: hSpec outside the attainable enthalpy range, an
     *         unpopulated cp table, an inner-flash convergence failure on
     *         the sweep, or a residual exceeding the acceptance tolerance
     *         at the converged temperature. The latter is how a
     *         departure-model enthalpy JUMP at a phase-label flip reports:
     *         the bracketing loop converges its interval onto the
     *         discontinuity (it does not exhaust iterations there) and only
     *         the residual check can tell that no root exists. On false,
     *         the state is left at the last trial evaluated — treat it as
     *         invalid.
     *
     * \note std::logic_error (e.g. an unknown twoPhaseMethod string) is
     *       deliberately NOT caught: programmer errors stay loud.
     * \note The inner-flash knobs are separate arguments rather than config
     *       fields by design: they mirror the isothermal flash's own solve()
     *       signature one-to-one.
     */
    template <class FluidState>
    static bool solve(FluidState& fluidState,
                      const Scalar hSpec,
                      const Config& cfg,
                      const std::string& twoPhaseMethod,
                      const Scalar ptTolerance,
                      const EOSType& eosType,
                      const int verbosity = 0)
    {
        using Flash = Opm::PTFlash<Scalar, FluidSystem, true>;
        using ValueType = typename FluidState::ValueType;

        // reject an unpopulated cp table: all-zero polynomials make
        // H(T) identically zero, which would let hSpec = 0 "succeed" at an
        // arbitrary bracket bound
        Scalar cpMagnitude = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            const auto& c = cfg.cpTable[compIdx];
            cpMagnitude += std::abs(c.c0) + std::abs(c.c1) + std::abs(c.c2) + std::abs(c.c3);
        }
        if (!(cpMagnitude > 0.0))
            return false;

        // Warm start across NEARBY trial temperatures: a previous trial's
        // converged equilibrium ratios seed the next inner flash far better
        // than the Wilson correlation (the same answer-reuse the reference
        // implementations apply). Two hard rules keep it sound:
        //
        //   1. Only K is warm-seeded; L is ALWAYS set to the cold sentinel
        //      -1. The isothermal flash runs its phase-stability test only
        //      for a non-interior L, and the test itself starts its
        //      Michelsen trials from the incoming K — so a warm K
        //      accelerates BOTH the stability test and the split solve,
        //      while an interior warm L would silently skip the stability
        //      verdict altogether. That skip is the wrong-root trap: near a
        //      phase boundary a stale split converges the inner flash to
        //      the trivial root (x = y = z satisfies the fugacity-ratio
        //      criterion exactly), which does not throw and poisons the
        //      residual — and, on the FINAL evaluation, the state's
        //      saturations — silently. No split may ever be accepted
        //      without a stability verdict at its own (T, p, z).
        //
        //   2. The seed is reused only within a proximity window of the
        //      temperature that produced it (with the stability test always
        //      running, this is a pure performance heuristic, not a
        //      correctness guard); the bracket endpoints, evaluated back to
        //      back, therefore always run cold.
        //
        // Only a converged TWO-PHASE split with strictly positive x AND y
        // is cached (single-phase ratios are degenerate; a zero component
        // would cache K = 0, a poison seed); the cached ratios are y/x from
        // the converged state (fs.K() is an input, never an output). A
        // warm-seeded inner flash that fails retires the seed and retries
        // cold.
        bool haveWarmSeed = false;
        std::array<Scalar, numComponents> warmK{};
        Scalar warmT = 0.0;
        const Scalar maxWarmStep =
            std::min(0.1 * (cfg.tempMax - cfg.tempMin), Scalar(15.0));

        auto seedCold = [&]() {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setKvalue(compIdx, fluidState.wilsonK_(compIdx));
            fluidState.setLvalue(ValueType(-1.0));
        };

        // residual r(T) = H(p, T, z) - hSpec; evaluating it flashes the state
        // at the trial temperature (the isothermal flash's input contract)
        auto residual = [&](const Scalar T) -> Scalar {
            fluidState.setTemperature(ValueType(T));

            bool solved = false;
            if (haveWarmSeed && std::abs(T - warmT) <= maxWarmStep) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    fluidState.setKvalue(compIdx, ValueType(warmK[compIdx]));
                fluidState.setLvalue(ValueType(-1.0)); // rule 1: stability always runs
                try {
                    Flash::solve(fluidState, twoPhaseMethod, ptTolerance, eosType, verbosity);
                    solved = true;
                }
                catch (const std::runtime_error&) {
                    haveWarmSeed = false; // stale seed misled the flash — retire it
                }
            }
            if (!solved) {
                seedCold();
                Flash::solve(fluidState, twoPhaseMethod, ptTolerance, eosType, verbosity);
            }

            // cache the converged split for the next trial
            const Scalar L = Opm::getValue(fluidState.L());
            if (L > 0.0 && L < 1.0) {
                haveWarmSeed = true;
                warmT = T;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar x = Opm::getValue(
                        fluidState.moleFraction(FluidSystem::oilPhaseIdx, compIdx));
                    const Scalar y = Opm::getValue(
                        fluidState.moleFraction(FluidSystem::gasPhaseIdx, compIdx));
                    if (!(x > 0.0) || !(y > 0.0)) {
                        haveWarmSeed = false;
                        break;
                    }
                    warmK[compIdx] = y / x;
                }
            }
            else {
                haveWarmSeed = false;
            }

            const Scalar h = Opm::getValue(EnthalpyCalc::mixtureEnthalpy(
                fluidState, cfg.cpTable, cfg.refTemperature, eosType, cfg.model));
            return h - hSpec;
        };
        // The root finder's API is double-typed; bridge explicitly, and
        // memoize the two bracket-endpoint values: the root finder
        // re-evaluates f at both endpoints on entry, which would otherwise
        // duplicate the pre-check's two full (cold) inner flashes.
        double rMinCached = 0.0, rMaxCached = 0.0;
        bool endpointsCached = false;
        auto residualAsDouble = [&](const double T) -> double {
            if (endpointsCached) {
                if (T == static_cast<double>(cfg.tempMin))
                    return rMinCached;
                if (T == static_cast<double>(cfg.tempMax))
                    return rMaxCached;
            }
            return static_cast<double>(residual(static_cast<Scalar>(T)));
        };

        // Everything below reports failure as `false` per the contract:
        // an inner-flash convergence failure anywhere on the sweep, a
        // bracket without a sign change (hSpec unattainable — the bracket
        // pre-check also keeps the root finder's bracketing-failure throw
        // unreachable), or a final enthalpy residual beyond the acceptance
        // tolerance. The residual check matters because the root finder's
        // interval exit does not look at the function value: where the
        // departure-model enthalpy jumps at a phase-label flip, the loop
        // CONVERGES its interval onto the discontinuity (it does not
        // exhaust iterations there) and only the residual can tell that no
        // root exists inside the jump.
        try {
            const Scalar rMin = residual(cfg.tempMin);
            const Scalar rMax = residual(cfg.tempMax);
            if (rMin * rMax > 0.0)
                return false; // hSpec unattainable on the bracket
            rMinCached = static_cast<double>(rMin);
            rMaxCached = static_cast<double>(rMax);
            endpointsCached = true;

            int iterationsUsed = 0;
            const double temperature =
                RegulaFalsi<ThrowOnError>::solve(residualAsDouble,
                                                 static_cast<double>(cfg.tempMin),
                                                 static_cast<double>(cfg.tempMax),
                                                 cfg.maxIterations,
                                                 static_cast<double>(cfg.tolerance),
                                                 iterationsUsed);

            // final consistent state at the solution temperature — and the
            // acceptance test on the value this call produces
            const Scalar rFinal = residual(static_cast<Scalar>(temperature));
            if (!(std::abs(rFinal) <= cfg.enthalpyTolerance))
                return false; // converged onto a discontinuity, not a root
        }
        catch (const std::runtime_error&) {
            return false;
        }
        return true;
    }
};

} // namespace Opm

#endif // OPM_PH_FLASH_HPP
