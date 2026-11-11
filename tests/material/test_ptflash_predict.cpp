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
 * \brief The research switch OPM_FLASH_PREDICT=predict,memo of PTFlash::solve().
 *
 *        - memo: a flash whose inputs repeat bit for bit returns the answer of
 *          the first one, values and derivatives to the last bit;
 *        - predict: a warm start a small step away from a remembered answer is
 *          accepted from its first-order prediction, and agrees with a cold
 *          flash of the same state within the flash tolerance; a large step is
 *          not accepted and still converges to the cold answer.
 *
 *        The switch is read once per process, so this test sets it before
 *        the first flash.
 */
#include "config.h"

#define OPM_FLASH_STATISTICS 1

#define BOOST_TEST_MODULE PtFlashPredict
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/FlashStatistics.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <cmath>
#include <cstdlib>

namespace {

using Scalar = double;
constexpr int numComponents = 3;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
// Pressure and the first two overall mole fractions, as in the simulator.
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using EOSType = Opm::CompositionalConfig::EOSType;
using PTFlash = Opm::PTFlash<Scalar, FluidSystem>;
using Composition = std::array<Scalar, numComponents>;

constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr Scalar flashTolerance = 1e-8;
constexpr EOSType eosType = EOSType::PR;
constexpr auto method = Opm::PTFlashMethod::SsiNewton;
constexpr Scalar temperature = 423.15;

// CO2, methane and decane with the data of the simple compositional decks,
// and the switch, before the first flash reads it.
struct Fixture
{
    Fixture()
    {
        ::setenv("OPM_FLASH_PREDICT", "predict,memo", 1);
        using CompParam = typename FluidSystem::ComponentParam;
        FluidSystem::init();
        FluidSystem::addComponent(CompParam{"CO2", 44.0, 304.128, 73.773e5, 0.09412, 0.22394});
        FluidSystem::addComponent(CompParam{"C1", 16.04, 190.564, 45.992e5, 0.09863, 0.01142});
        FluidSystem::addComponent(CompParam{"C10", 142.28, 617.7, 21.03e5, 0.6098, 0.4884});
    }
};

FluidState makeState(const Scalar pressure, const Composition& z)
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

// The state at pressure and z, warm-started from the answer of a flash, as
// the simulator's hint does.
FluidState warmState(const Scalar pressure, const Composition& z, const FluidState& from)
{
    auto fs = makeState(pressure, z);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, Opm::getValue(from.K(compIdx)));
    }
    fs.setLvalue(Opm::getValue(from.L()));
    return fs;
}

bool identical(const Evaluation& a, const Evaluation& b)
{
    if (a.value() != b.value()) {
        return false;
    }
    for (int i = 0; i < Evaluation::numVars; ++i) {
        if (a.derivative(i) != b.derivative(i)) {
            return false;
        }
    }
    return true;
}

bool identicalAnswers(const FluidState& a, const FluidState& b)
{
    bool same = identical(a.L(), b.L());
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        same = same && identical(a.moleFraction(oilPhaseIdx, compIdx), b.moleFraction(oilPhaseIdx, compIdx))
            && identical(a.moleFraction(gasPhaseIdx, compIdx), b.moleFraction(gasPhaseIdx, compIdx))
            && identical(a.K(compIdx), b.K(compIdx));
    }
    return same;
}

// The largest difference of L, x, y, and of their derivatives.
std::array<Scalar, 2> difference(const FluidState& a, const FluidState& b)
{
    std::array<Scalar, 2> d{0., 0.};
    const auto add = [&d](const Evaluation& u, const Evaluation& v) {
        d[0] = std::max(d[0], std::abs(u.value() - v.value()));
        for (int i = 0; i < Evaluation::numVars; ++i) {
            d[1] = std::max(d[1], std::abs(u.derivative(i) - v.derivative(i)));
        }
    };
    add(a.L(), b.L());
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        add(a.moleFraction(oilPhaseIdx, compIdx), b.moleFraction(oilPhaseIdx, compIdx));
        add(a.moleFraction(gasPhaseIdx, compIdx), b.moleFraction(gasPhaseIdx, compIdx));
    }
    return d;
}

long cacheHits()
{
    return Opm::FlashStatistics::get().cacheHits;
}

// Two-phase states at the decks' 150 C.
constexpr std::array<Composition, 3> compositions{
    Composition{0.1, 0.8, 0.1},
    Composition{0.3, 0.5, 0.2},
    Composition{0.6, 0.3, 0.1},
};
constexpr std::array<Scalar, 3> pressures{55.0e5, 80.0e5, 120.0e5};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(RepeatedInputsReturnTheSameAnswer)
{
    for (const auto& z : compositions) {
        for (const Scalar pressure : pressures) {
            auto first = makeState(pressure, z);
            const bool firstSingle = PTFlash::solve(first, method, flashTolerance, eosType);
            const long hits = cacheHits();
            auto second = makeState(pressure, z);
            const bool secondSingle = PTFlash::solve(second, method, flashTolerance, eosType);
            BOOST_CHECK_EQUAL(cacheHits(), hits + 1);
            BOOST_CHECK_EQUAL(firstSingle, secondSingle);
            BOOST_CHECK(identicalAnswers(first, second));
        }
    }
}

BOOST_AUTO_TEST_CASE(SmallStepIsPredicted)
{
    int twoPhase = 0;
    for (const auto& z : compositions) {
        for (const Scalar pressure : pressures) {
            auto previous = makeState(pressure, z);
            if (PTFlash::solve(previous, method, flashTolerance, eosType)) {
                continue;
            }
            ++twoPhase;
            // A late Newton iteration of a simulator: a step of 1e-6.
            const Scalar next = pressure * (1.0 + 1e-6);
            Composition moved = z;
            moved[0] += 1e-6;
            moved[1] -= 1e-6;
            const long hits = cacheHits();
            auto warm = warmState(next, moved, previous);
            BOOST_REQUIRE(!PTFlash::solve(warm, method, flashTolerance, eosType));
            BOOST_CHECK_EQUAL(cacheHits(), hits + 1);

            auto cold = makeState(next, moved);
            BOOST_REQUIRE(!PTFlash::solve(cold, method, flashTolerance, eosType));
            const auto [value, derivative] = difference(warm, cold);
            BOOST_TEST_CONTEXT("p = " << pressure / 1e5 << " bar, z0 = " << z[0]) {
                BOOST_CHECK_SMALL(value, 1e-7);
                // Derivatives in Pa^-1 and per mole fraction, evaluated at
                // states 1e-8 apart.
                BOOST_CHECK_SMALL(derivative, 1e-4);
            }
        }
    }
    BOOST_CHECK(twoPhase > 0);
}

BOOST_AUTO_TEST_CASE(LargeStepConvergesToTheColdAnswer)
{
    for (const auto& z : compositions) {
        for (const Scalar pressure : pressures) {
            auto previous = makeState(pressure, z);
            if (PTFlash::solve(previous, method, flashTolerance, eosType)) {
                continue;
            }
            const Scalar next = pressure * 1.05;
            Composition moved = z;
            moved[0] += 0.02;
            moved[1] -= 0.02;
            const long hits = cacheHits();
            auto warm = warmState(next, moved, previous);
            const bool warmSingle = PTFlash::solve(warm, method, flashTolerance, eosType);
            BOOST_CHECK_EQUAL(cacheHits(), hits);

            auto cold = makeState(next, moved);
            const bool coldSingle = PTFlash::solve(cold, method, flashTolerance, eosType);
            BOOST_CHECK_EQUAL(warmSingle, coldSingle);
            if (!warmSingle && !coldSingle) {
                BOOST_CHECK_SMALL(difference(warm, cold)[0], 1e-6);
            }
        }
    }
}
