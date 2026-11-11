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
 * \brief The derivative pass after the scalar flash against the textbook
 *        implicit-function derivatives.
 *
 *        The reference assembles the full Jacobians of the flash residuals in
 *        the unknowns u = (x, y, L) and in the state s = (p, [T,] z) with
 *        automatic differentiation, inverts dF/du and applies
 *        du/ds = -(dF/du)^-1 dF/ds with the chain rule through the caller's
 *        derivatives. The pass of the flash must give the same derivatives to
 *        1e-8 relative, at two-phase states of a three- and a six-component
 *        mixture, isothermal and thermal, with and without the Jacobian of a
 *        Newton solve; and a single phase must keep the derivatives of the
 *        overall composition.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashDerivativePass
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/common/Exceptions.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using PTFlashMethod = Opm::PTFlashMethod;
constexpr EOSType eosType = EOSType::PR;
constexpr Scalar flashTolerance = 1e-8;
constexpr Scalar derivativeTolerance = 1e-8;

using ThreeComponentFluid = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, false>;
using SixComponentFluid = Opm::GenericOilGasWaterFluidSystem<Scalar, 6, false>;

struct Component
{
    std::string name;
    Scalar molarMass, criticalTemperature, criticalPressure, criticalVolume, acentricFactor;
};

template <class FluidSystem>
void defineFluid(const std::vector<Component>& components,
                 const std::vector<std::tuple<int, int, Scalar>>& interactions)
{
    FluidSystem::init();
    for (const auto& c : components) {
        FluidSystem::addComponent(typename FluidSystem::ComponentParam{
            c.name, c.molarMass, c.criticalTemperature, c.criticalPressure, c.criticalVolume,
            c.acentricFactor});
    }
    // The lower triangle, row by row.
    std::vector<Scalar> coefficients;
    for (int row = 1; row < FluidSystem::numComponents; ++row) {
        for (int col = 0; col < row; ++col) {
            Scalar value = 0.0;
            for (const auto& [i, j, k] : interactions) {
                if ((i == row && j == col) || (i == col && j == row)) {
                    value = k;
                }
            }
            coefficients.push_back(value);
        }
    }
    FluidSystem::setInteractionCoefficients(std::move(coefficients));
}

// CO2, methane and decane; and the six SPE5 pseudo-components.
struct Fixture
{
    Fixture()
    {
        defineFluid<ThreeComponentFluid>({{"CO2", 44.01, 304.128, 73.773e5, 0.09412, 0.22394},
                                      {"C1", 16.04, 190.564, 45.992e5, 0.09863, 0.01142},
                                      {"C10", 142.28, 617.7, 21.03e5, 0.6098, 0.4884}},
                                     {{0, 1, 0.1}, {0, 2, 0.1}, {1, 2, 0.0411}});
        defineFluid<SixComponentFluid>({{"S1", 16.04, 190.56, 4.6043e6, 0.0998, 0.013},
                                    {"S3", 44.10, 369.83, 4.2493e6, 0.2005, 0.1524},
                                    {"S6", 86.18, 507.44, 3.0123e6, 0.3698, 0.3007},
                                    {"S10", 142.29, 617.67, 2.0960e6, 0.6297, 0.4885},
                                    {"S15", 206.00, 705.56, 1.3790e6, 1.0423, 0.6500},
                                    {"S20", 282.00, 766.67, 1.1170e6, 1.3412, 0.8500}},
                                   {{0, 4, 0.05}, {0, 5, 0.05}, {1, 4, 0.005}, {1, 5, 0.005}});
    }
};

// Exposes the stages solve() chains, and the residual assembly the
// reference is built from.
template <class FluidSystem, bool isThermal>
class ExposedPtFlash : public Opm::PTFlash<Scalar, FluidSystem, isThermal>
{
    using Base = Opm::PTFlash<Scalar, FluidSystem, isThermal>;

public:
    using Base::NewtonJacobian;
    using Base::assembleNewton_;
    using Base::compositionJacobian_;
    using Base::flashScalar_;
    using Base::updateDerivatives_;
    using Base::updateFugacityCoefficients_;
};

struct Summary
{
    int twoPhase = 0;
    int singlePhase = 0;
    int reused = 0;
    Scalar worst = 0.0;
};

template <class FluidSystem, bool isThermal>
struct DerivativeCheck
{
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    // Pressure, [temperature,] and all overall mole fractions but the last.
    static constexpr int firstCompositionVar = isThermal ? 2 : 1;
    static constexpr int numVars = firstCompositionVar + numComponents - 1;
    using Evaluation = Opm::DenseAd::Evaluation<Scalar, numVars>;
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
    using ScalarFluidState = Opm::CompositionalFluidState<Scalar, FluidSystem>;
    using Flash = ExposedPtFlash<FluidSystem, isThermal>;
    using OptionalJacobian = std::optional<typename Flash::NewtonJacobian>;
    using Composition = std::array<Scalar, numComponents>;

    static FluidState makeState(const Scalar pressure, const Scalar temperature, const Composition& z)
    {
        FluidState fs;
        const auto p = Evaluation::createVariable(pressure, 0);
        fs.setPressure(oilPhaseIdx, p);
        fs.setPressure(gasPhaseIdx, p);
        if constexpr (isThermal) {
            fs.setTemperature(Evaluation::createVariable(temperature, 1));
        } else {
            fs.setTemperature(Evaluation(temperature));
        }
        Evaluation last = 1.0;
        for (int compIdx = 0; compIdx + 1 < numComponents; ++compIdx) {
            const auto zc = Evaluation::createVariable(z[compIdx], firstCompositionVar + compIdx);
            fs.setMoleFraction(compIdx, zc);
            last -= zc;
        }
        fs.setMoleFraction(numComponents - 1, last);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
        }
        fs.setLvalue(-1.0);
        return fs;
    }

    static ScalarFluidState scalarCopy(const FluidState& fs)
    {
        ScalarFluidState scalar;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            scalar.setKvalue(compIdx, Opm::getValue(fs.K(compIdx)));
            scalar.setMoleFraction(compIdx, Opm::getValue(fs.moleFraction(compIdx)));
        }
        scalar.setLvalue(Opm::getValue(fs.L()));
        scalar.setPressure(oilPhaseIdx, Opm::getValue(fs.pressure(oilPhaseIdx)));
        scalar.setPressure(gasPhaseIdx, Opm::getValue(fs.pressure(gasPhaseIdx)));
        scalar.setTemperature(Opm::getValue(fs.temperature(0)));
        return scalar;
    }

    // The derivative pass of solve().
    static FluidState flashPass(const ScalarFluidState& scalar,
                                FluidState fs,
                                const bool singlePhase,
                                const OptionalJacobian& newtonJacobian)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fs.setMoleFraction(oilPhaseIdx, compIdx, scalar.moleFraction(oilPhaseIdx, compIdx));
            fs.setMoleFraction(gasPhaseIdx, compIdx, scalar.moleFraction(gasPhaseIdx, compIdx));
        }
        Flash::updateDerivatives_(scalar, fs, eosType, singlePhase, newtonJacobian);
        return fs;
    }

    // du/ds = -(dF/du)^-1 dF/ds from the full Jacobians, with an explicit
    // inverse, composed with the caller's derivatives of s.
    static FluidState referencePass(const ScalarFluidState& scalar, FluidState fs)
    {
        constexpr int numUnknowns = 2 * numComponents + 1;
        constexpr int numStateVars = firstCompositionVar + numComponents;
        using StateEval = Opm::DenseAd::Evaluation<Scalar, numStateVars>;
        using StateFluidState = Opm::CompositionalFluidState<StateEval, FluidSystem>;
        using StateComposition = Dune::FieldVector<StateEval, numComponents>;

        StateFluidState state;
        const auto p = StateEval::createVariable(scalar.pressure(oilPhaseIdx), 0);
        state.setPressure(oilPhaseIdx, p);
        state.setPressure(gasPhaseIdx, p);
        if constexpr (isThermal) {
            state.setTemperature(StateEval::createVariable(scalar.temperature(0), 1));
        } else {
            state.setTemperature(StateEval(scalar.temperature(0)));
        }
        StateComposition z;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            z[compIdx] = StateEval::createVariable(scalar.moleFraction(compIdx),
                                                   firstCompositionVar + compIdx);
            state.setMoleFraction(oilPhaseIdx, compIdx, scalar.moleFraction(oilPhaseIdx, compIdx));
            state.setMoleFraction(gasPhaseIdx, compIdx, scalar.moleFraction(gasPhaseIdx, compIdx));
            state.setKvalue(compIdx, scalar.K(compIdx));
        }
        state.setLvalue(scalar.L());
        typename FluidSystem::template ParameterCache<StateEval> paramCache(eosType);
        for (const int phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
            paramCache.updatePhase(state, phaseIdx);
            Flash::updateFugacityCoefficients_(state, paramCache, phaseIdx);
        }
        Dune::FieldMatrix<Scalar, numUnknowns, numStateVars> stateJacobian;
        Dune::FieldVector<Scalar, numUnknowns> residual;
        Flash::template assembleNewton_<StateFluidState, StateComposition, numStateVars, numUnknowns>(
            state, z, stateJacobian, residual);

        auto inverse = Flash::compositionJacobian_(scalar, eosType);
        inverse.invert();
        stateJacobian.leftmultiply(inverse);

        std::array<Evaluation, numUnknowns> u;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            u[compIdx] = scalar.moleFraction(oilPhaseIdx, compIdx);
            u[compIdx + numComponents] = scalar.moleFraction(gasPhaseIdx, compIdx);
        }
        u[2 * numComponents] = scalar.L();
        for (int varIdx = 0; varIdx < numVars; ++varIdx) {
            std::array<Scalar, numStateVars> ds{};
            ds[0] = fs.pressure(oilPhaseIdx).derivative(varIdx);
            if constexpr (isThermal) {
                ds[1] = fs.temperature(0).derivative(varIdx);
            }
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                ds[firstCompositionVar + compIdx] = fs.moleFraction(compIdx).derivative(varIdx);
            }
            for (int unknownIdx = 0; unknownIdx < numUnknowns; ++unknownIdx) {
                Scalar derivative = 0.0;
                for (int stateIdx = 0; stateIdx < numStateVars; ++stateIdx) {
                    derivative -= stateJacobian[unknownIdx][stateIdx] * ds[stateIdx];
                }
                u[unknownIdx].setDerivative(varIdx, derivative);
            }
        }
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fs.setMoleFraction(oilPhaseIdx, compIdx, u[compIdx]);
            fs.setMoleFraction(gasPhaseIdx, compIdx, u[compIdx + numComponents]);
        }
        fs.setLvalue(u[2 * numComponents]);
        return fs;
    }

    // x, y and L of a state.
    static std::vector<Evaluation> flashUnknowns(const FluidState& fs)
    {
        std::vector<Evaluation> u;
        for (const int phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                u.push_back(fs.moleFraction(phaseIdx, compIdx));
            }
        }
        u.push_back(fs.L());
        return u;
    }

    // The largest difference of a derivative, relative to the largest
    // reference derivative of x, y and L in the same variable. Values must be
    // identical.
    static Scalar relativeDifference(const FluidState& tested, const FluidState& reference)
    {
        const auto a = flashUnknowns(tested);
        const auto b = flashUnknowns(reference);
        Scalar worst = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            BOOST_CHECK_EQUAL(a[i].value(), b[i].value());
        }
        for (int varIdx = 0; varIdx < numVars; ++varIdx) {
            Scalar scale = 0.0;
            for (const auto& v : b) {
                scale = std::max(scale, std::abs(v.derivative(varIdx)));
            }
            BOOST_CHECK(scale > 0.0);
            for (std::size_t i = 0; i < a.size(); ++i) {
                const Scalar difference = std::abs(a[i].derivative(varIdx) - b[i].derivative(varIdx));
                worst = std::max(worst, difference / scale);
            }
        }
        return worst;
    }

    static void checkSinglePhase(const FluidState& tested)
    {
        for (const int phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                const auto& x = tested.moleFraction(phaseIdx, compIdx);
                const auto& z = tested.moleFraction(compIdx);
                BOOST_CHECK_EQUAL(x.value(), z.value());
                for (int varIdx = 0; varIdx < numVars; ++varIdx) {
                    BOOST_CHECK_EQUAL(x.derivative(varIdx), z.derivative(varIdx));
                }
            }
        }
        for (int varIdx = 0; varIdx < numVars; ++varIdx) {
            BOOST_CHECK_EQUAL(tested.L().derivative(varIdx), 0.0);
        }
    }

    static Summary run(const std::vector<Scalar>& temperatures,
                       const std::vector<Composition>& compositions,
                       const std::vector<Scalar>& pressures)
    {
        using enum PTFlashMethod;
        Summary summary;
        for (const auto method : {SsiNewton, Ssi, Newton}) {
            for (const Scalar temperature : temperatures) {
                for (const auto& z : compositions) {
                    for (const Scalar pressure : pressures) {
                        BOOST_TEST_CONTEXT(numComponents << " components, "
                                           << (isThermal ? "thermal, " : "")
                                           << Opm::ptFlashMethodToString(method)
                                           << ", T = " << temperature << " K, p = "
                                           << pressure / 1e5 << " bar, z0 = " << z[0]) {
                            const auto fs = makeState(pressure, temperature, z);
                            auto scalar = scalarCopy(fs);
                            bool singlePhase = true;
                            OptionalJacobian newtonJacobian;
                            try {
                                auto result = Flash::flashScalar_(scalar, method, flashTolerance,
                                                                  eosType, 0);
                                singlePhase = result.single_phase;
                                newtonJacobian = std::move(result.newton_jacobian);
                            } catch (const Opm::NumericalProblem&) {
                                // Plain Newton has nothing to fall back on, and
                                // plain substitution can reach its iteration limit.
                                BOOST_CHECK(method != SsiNewton);
                                continue;
                            }
                            if (singlePhase) {
                                ++summary.singlePhase;
                                checkSinglePhase(flashPass(scalar, fs, true, std::nullopt));
                                continue;
                            }
                            ++summary.twoPhase;
                            const auto reference = referencePass(scalar, fs);
                            const Scalar assembled = relativeDifference(
                                flashPass(scalar, fs, false, std::nullopt), reference);
                            BOOST_CHECK_LE(assembled, derivativeTolerance);
                            summary.worst = std::max(summary.worst, assembled);
                            if (newtonJacobian) {
                                ++summary.reused;
                                const Scalar reused = relativeDifference(
                                    flashPass(scalar, fs, false, newtonJacobian), reference);
                                BOOST_CHECK_LE(reused, derivativeTolerance);
                                summary.worst = std::max(summary.worst, reused);
                            }
                        }
                    }
                }
            }
        }
        BOOST_TEST_MESSAGE(numComponents << " components" << (isThermal ? ", thermal" : "")
                           << ": " << summary.twoPhase << " two-phase states ("
                           << summary.reused << " with a Newton Jacobian), "
                           << summary.singlePhase << " single-phase; largest relative "
                           << "difference " << summary.worst);
        BOOST_CHECK(summary.twoPhase > 0);
        BOOST_CHECK(summary.singlePhase > 0);
        BOOST_CHECK(summary.reused > 0);
        return summary;
    }
};

const std::vector<std::array<Scalar, 3>> threeComponentCompositions{
    {0.1, 0.8, 0.1}, {0.3, 0.5, 0.2}, {0.6, 0.3, 0.1}, {0.2, 0.4, 0.4}, {0.4, 0.3, 0.3}};
const std::vector<std::array<Scalar, 6>> sixComponentCompositions{
    {0.635, 0.115, 0.05, 0.10, 0.075, 0.025},
    {0.5, 0.03, 0.07, 0.2, 0.15, 0.05},
    {0.70, 0.20, 0.03, 0.04, 0.02, 0.01},
    {0.3, 0.02, 0.08, 0.3, 0.2, 0.1}};
const std::vector<Scalar> pressures{10.0e5, 30.0e5, 60.0e5, 100.0e5, 150.0e5, 200.0e5, 300.0e5};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(ThreeComponents)
{
    DerivativeCheck<ThreeComponentFluid, false>::run({280.0, 350.0, 423.15},
                                                 threeComponentCompositions, pressures);
}

BOOST_AUTO_TEST_CASE(ThreeComponentsThermal)
{
    DerivativeCheck<ThreeComponentFluid, true>::run({280.0, 350.0, 423.15},
                                                threeComponentCompositions, pressures);
}

BOOST_AUTO_TEST_CASE(SixComponents)
{
    DerivativeCheck<SixComponentFluid, false>::run({300.0, 380.0, 450.0},
                                               sixComponentCompositions, pressures);
}

BOOST_AUTO_TEST_CASE(SixComponentsThermal)
{
    DerivativeCheck<SixComponentFluid, true>::run({300.0, 380.0, 450.0},
                                              sixComponentCompositions, pressures);
}
