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
    //! [K] lower bracket bound. Note this default extrapolates the cp
    //! correlations slightly below their nominal validity (~273 K) — accepted
    //! for the synthetic self-consistent use; narrow for quantitative work.
    Scalar tempMin = 200.0;
    Scalar tempMax = 600.0;   // [K] upper bracket bound
    //! convergence tolerance of the bracketing solver — applied to both the
    //! temperature interval [K] and the enthalpy residual [J/mol]
    Scalar tolerance = 1e-6;
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
     * \return true on success; false when no solution is found on the
     *         bracket: hSpec outside the attainable enthalpy range, an
     *         unpopulated cp table, an inner-flash convergence failure on the
     *         sweep, or root-finder iteration exhaustion (possible where the
     *         departure-model enthalpy jumps at a phase-label flip). On
     *         false, the state is left at the last trial evaluated — treat
     *         it as invalid.
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

        // residual r(T) = H(p, T, z) - hSpec; evaluating it flashes the state
        // at the trial temperature (Wilson-K + L re-seeded each trial — the
        // isothermal flash's input contract; fs.K() is an input, not output)
        auto residual = [&](const Scalar T) -> Scalar {
            fluidState.setTemperature(ValueType(T));
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setKvalue(compIdx, fluidState.wilsonK_(compIdx));
            fluidState.setLvalue(ValueType(1.0));
            Flash::solve(fluidState, twoPhaseMethod, ptTolerance, eosType, verbosity);
            const Scalar h = Opm::getValue(EnthalpyCalc::mixtureEnthalpy(
                fluidState, cfg.cpTable, cfg.refTemperature, eosType, cfg.model));
            return h - hSpec;
        };
        // the root finder's API is double-typed; bridge explicitly
        auto residualAsDouble = [&](const double T) -> double {
            return static_cast<double>(residual(static_cast<Scalar>(T)));
        };

        // Everything below reports failure as `false` per the contract:
        // an inner-flash convergence failure anywhere on the sweep, a
        // bracket without a sign change (hSpec unattainable — the bracket
        // pre-check also keeps the root finder's bracketing-failure throw
        // unreachable), or root-finder iteration exhaustion (which CAN
        // happen despite a sign change where the departure-model enthalpy
        // jumps at a single-phase label flip: a sign change then need not
        // enclose a smooth root).
        try {
            const Scalar rMin = residual(cfg.tempMin);
            const Scalar rMax = residual(cfg.tempMax);
            if (rMin * rMax > 0.0)
                return false; // hSpec unattainable on the bracket

            int iterationsUsed = 0;
            const double temperature =
                RegulaFalsi<ThrowOnError>::solve(residualAsDouble,
                                                 static_cast<double>(cfg.tempMin),
                                                 static_cast<double>(cfg.tempMax),
                                                 cfg.maxIterations,
                                                 static_cast<double>(cfg.tolerance),
                                                 iterationsUsed);

            // final consistent state at the solution temperature
            residual(static_cast<Scalar>(temperature));
        }
        catch (const std::runtime_error&) {
            return false;
        }
        return true;
    }
};

} // namespace Opm

#endif // OPM_PH_FLASH_HPP
