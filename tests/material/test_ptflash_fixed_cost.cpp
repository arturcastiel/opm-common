// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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
 * \brief A flash computes what depends only on its temperature, pressure and
 *        feed once and shares it between its stages: one parameter cache,
 *        the feed on both cubic roots, one evaluation of the answer for the
 *        derivative pass.  The shortcuts must not change a bit:
 *
 *        - the values-only flash returns, with a single-phase answer, the
 *          fugacity coefficients of the feed that the flash evaluated
 *          earlier; they must be those a new evaluation of the answer gives;
 *        - solve() asks for no Newton Jacobian where the derivative pass
 *          evaluates the answer analytically anyway; its derivatives must be
 *          those of the derivative pass that reuses the Newton Jacobian, on
 *          parameter caches of its own.
 *
 *        CO2, methane and decane with volume shifts and interaction
 *        coefficients, so the volume translation of the coefficients is
 *        exercised; single- and two-phase answers, both methods.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashFixedCost
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/common/Exceptions.hpp>

#include <dune/common/fmatrix.hh>

#include <array>
#include <cstring>
#include <optional>
#include <vector>

namespace {

using Scalar = double;
constexpr int numComponents = 3;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using ScalarFluidState = Opm::CompositionalFluidState<Scalar, FluidSystem>;
using ParamCache = FluidSystem::ParameterCache<Scalar>;
using EOSType = Opm::CompositionalConfig::EOSType;
using PTFlashMethod = Opm::PTFlashMethod;
using Composition = std::array<Scalar, numComponents>;

constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr Scalar flashTolerance = 1e-8;
constexpr EOSType eosType = EOSType::PR;

class ExposedPtFlash : public Opm::PTFlash<Scalar, FluidSystem>
{
    using Base = Opm::PTFlash<Scalar, FluidSystem>;

public:
    using Base::flashScalar_;
    using Base::updateDerivatives_;
};

struct Fixture
{
    Fixture()
    {
        using CompParam = typename FluidSystem::ComponentParam;
        FluidSystem::init();
        FluidSystem::addComponent(CompParam{"CO2", 44.0, 304.128, 73.773e5, 0.09412, 0.22394, -0.05});
        FluidSystem::addComponent(CompParam{"C1", 16.04, 190.564, 45.992e5, 0.09863, 0.01142, -0.15});
        FluidSystem::addComponent(CompParam{"C10", 142.28, 617.7, 21.03e5, 0.6098, 0.4884, 0.08});
        FluidSystem::setInteractionCoefficients({0.1, 0.1, 0.0411});
    }
};

constexpr std::array<Scalar, 7> pressures{5.0e5, 20.0e5, 50.0e5, 80.0e5, 120.0e5, 200.0e5, 300.0e5};
constexpr std::array<Scalar, 8> temperatures{250.0, 280.0, 300.0, 330.0, 370.0, 420.0, 500.0, 580.0};
constexpr std::array<Composition, 3> compositions{
    Composition{0.4, 0.3, 0.3},
    Composition{0.9, 0.05, 0.05},
    Composition{0.05, 0.9, 0.05},
};
constexpr std::array<PTFlashMethod, 2> methods{PTFlashMethod::SsiNewton, PTFlashMethod::Ssi};

bool sameBits(const Scalar a, const Scalar b)
{
    return std::memcmp(&a, &b, sizeof(Scalar)) == 0;
}

bool sameBits(const Evaluation& a, const Evaluation& b)
{
    bool same = sameBits(a.value(), b.value());
    for (int i = 0; i < Evaluation::numVars; ++i) {
        same = same && sameBits(a.derivative(i), b.derivative(i));
    }
    return same;
}

ScalarFluidState coldScalarState(const Scalar pressure, const Scalar temperature, const Composition& z)
{
    ScalarFluidState fs;
    fs.setPressure(oilPhaseIdx, pressure);
    fs.setPressure(gasPhaseIdx, pressure);
    fs.setTemperature(temperature);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setMoleFraction(compIdx, z[compIdx]);
    }
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    return fs;
}

// Pressure and the first two overall mole fractions carry derivatives.
FluidState coldState(const Scalar pressure, const Scalar temperature, const Composition& z)
{
    FluidState fs;
    const auto p = Evaluation::createVariable(pressure, 0);
    fs.setPressure(oilPhaseIdx, p);
    fs.setPressure(gasPhaseIdx, p);
    fs.setTemperature(temperature);
    Evaluation last = 1.0;
    for (unsigned compIdx = 0; compIdx + 1 < numComponents; ++compIdx) {
        const auto zc = Evaluation::createVariable(z[compIdx], compIdx + 1);
        fs.setMoleFraction(compIdx, zc);
        last -= zc;
    }
    fs.setMoleFraction(numComponents - 1, last);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    return fs;
}

ScalarFluidState scalarCopy(const FluidState& fs)
{
    ScalarFluidState scalar;
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        scalar.setKvalue(compIdx, Opm::getValue(fs.K(compIdx)));
        scalar.setMoleFraction(compIdx, Opm::getValue(fs.moleFraction(compIdx)));
    }
    scalar.setLvalue(Opm::getValue(fs.L()));
    scalar.setPressure(oilPhaseIdx, Opm::getValue(fs.pressure(oilPhaseIdx)));
    scalar.setPressure(gasPhaseIdx, Opm::getValue(fs.pressure(gasPhaseIdx)));
    scalar.setTemperature(Opm::getValue(fs.temperature(0)));
    return scalar;
}

} // namespace

BOOST_FIXTURE_TEST_SUITE(PtFlashFixedCost, Fixture)

BOOST_AUTO_TEST_CASE(ValuesOnlyFlashReturnsTheCoefficientsOfItsAnswer)
{
    int single = 0;
    int split = 0;
    for (const auto method : methods) {
        for (const auto& z : compositions) {
            for (const Scalar pressure : pressures) {
                for (const Scalar temperature : temperatures) {
                    auto fs = coldScalarState(pressure, temperature, z);
                    bool isSingle = false;
                    try {
                        isSingle = ExposedPtFlash::flash_solve_scalar(fs, method, flashTolerance,
                                                                      eosType, 0);
                    } catch (const Opm::NumericalProblem&) {
                        continue;
                    } catch (const Dune::FMatrixError&) {
                        continue;
                    }
                    ++(isSingle ? single : split);
                    // A new evaluation of the answer on a cache of its own.
                    ParamCache cache(eosType);
                    for (const unsigned phaseIdx : {static_cast<unsigned>(oilPhaseIdx),
                                                    static_cast<unsigned>(gasPhaseIdx)}) {
                        cache.updatePhase(fs, phaseIdx);
                        std::array<Scalar, numComponents> expected;
                        FluidSystem::fugacityCoefficients(fs, cache, phaseIdx, expected);
                        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                            BOOST_CHECK_MESSAGE(
                                sameBits(fs.fugacityCoefficient(phaseIdx, compIdx), expected[compIdx]),
                                "phase " << phaseIdx << " component " << compIdx << " at p = "
                                         << pressure << ", T = " << temperature << ": "
                                         << fs.fugacityCoefficient(phaseIdx, compIdx) << " vs "
                                         << expected[compIdx]);
                        }
                    }
                }
            }
        }
    }
    // Both kinds of answer occur, so both closing paths are checked.
    BOOST_CHECK_GT(single, 20);
    BOOST_CHECK_GT(split, 20);
}

BOOST_AUTO_TEST_CASE(DerivativesWithoutTheNewtonJacobianAreThoseWithIt)
{
    int splits = 0;
    int reused = 0;
    for (const auto method : methods) {
        for (const auto& z : compositions) {
            for (const Scalar pressure : pressures) {
                for (const Scalar temperature : temperatures) {
                    auto fs = coldState(pressure, temperature, z);
                    auto reference = fs;
                    bool isSingle = false;
                    try {
                        isSingle = ExposedPtFlash::solve(fs, method, flashTolerance, eosType, 0);
                    } catch (const Opm::NumericalProblem&) {
                        continue;
                    } catch (const Dune::FMatrixError&) {
                        continue;
                    }
                    // The stages of solve() on their own: the scalar flash with
                    // the Newton Jacobian, the derivative pass on parameter
                    // caches of its own, reusing that Jacobian.
                    auto scalar = scalarCopy(reference);
                    auto [referenceSingle, newtonJacobian]
                        = ExposedPtFlash::flashScalar_(scalar, method, flashTolerance, eosType, 0);
                    BOOST_REQUIRE_EQUAL(isSingle, referenceSingle);
                    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                        reference.setMoleFraction(oilPhaseIdx, compIdx,
                                                  scalar.moleFraction(oilPhaseIdx, compIdx));
                        reference.setMoleFraction(gasPhaseIdx, compIdx,
                                                  scalar.moleFraction(gasPhaseIdx, compIdx));
                    }
                    ExposedPtFlash::updateDerivatives_(scalar, reference, eosType, referenceSingle,
                                                       newtonJacobian);
                    if (referenceSingle) {
                        continue;
                    }
                    ++splits;
                    reused += newtonJacobian.has_value();
                    BOOST_CHECK(sameBits(fs.L(), reference.L()));
                    for (const auto phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
                        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                            BOOST_CHECK_MESSAGE(sameBits(fs.moleFraction(phaseIdx, compIdx),
                                                         reference.moleFraction(phaseIdx, compIdx)),
                                                "phase " << phaseIdx << " component " << compIdx
                                                         << " at p = " << pressure
                                                         << ", T = " << temperature);
                        }
                    }
                }
            }
        }
    }
    BOOST_CHECK_GT(splits, 20);
    // The reference took the Newton Jacobian on most splits.
    BOOST_CHECK_GT(reused, 10);
}

BOOST_AUTO_TEST_SUITE_END()
