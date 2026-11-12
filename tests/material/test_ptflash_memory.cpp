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
 * \brief PTFlash::solve() with the memory of a cell's last flash.
 *
 *        - all options off: the answer of solve() without a memory, bit for bit;
 *        - predict: the remembered sensitivities are those of the derivative
 *          pass; a warm start a small step away is accepted from its
 *          prediction and agrees with a cold flash within the flash
 *          tolerance; a large step is not accepted and reaches the cold
 *          answer; reuse gives derivatives close to the derivative pass;
 *        - phase memory: a single-phase cell that hardly moved skips its
 *          stability test with the answer of the test, bit for bit; one that
 *          moved into the two-phase region is tested and split.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashMemory
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace {

using Scalar = double;
constexpr int numComponents = 3;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
// Pressure and the first two overall mole fractions, as in the simulator.
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using EOSType = Opm::CompositionalConfig::EOSType;
using PTFlash = Opm::PTFlash<Scalar, FluidSystem>;
using Memory = PTFlash::Memory;
using Options = PTFlash::MemoryOptions;
using Composition = std::array<Scalar, numComponents>;

constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr Scalar flashTolerance = 1e-8;
constexpr EOSType eosType = EOSType::PR;
constexpr auto method = Opm::PTFlashMethod::SsiNewton;
constexpr Scalar temperature = 423.15;
// The feed of the simple compositional decks: two-phase at 75 bar, an
// undersaturated liquid at 250 bar.
constexpr Composition feed{0.1, 0.3, 0.6};

// CO2, methane and decane with the data of the simple compositional decks.
struct Fixture
{
    Fixture()
    {
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

// The largest difference of a value or derivative of x, y and L, relative to
// the largest value or derivative of the same quantity of b.
Scalar largestRelativeDifference(const FluidState& a, const FluidState& b, const bool derivatives)
{
    Scalar largest = 0.;
    const auto compare = [&](const Evaluation& u, const Evaluation& v) {
        if (!derivatives) {
            largest = std::max(largest, std::abs(u.value() - v.value()));
            return;
        }
        Scalar scale = 1e-30;
        for (int i = 0; i < Evaluation::numVars; ++i) {
            scale = std::max(scale, std::abs(v.derivative(i)));
        }
        for (int i = 0; i < Evaluation::numVars; ++i) {
            largest = std::max(largest, std::abs(u.derivative(i) - v.derivative(i)) / scale);
        }
    };
    compare(a.L(), b.L());
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        compare(a.moleFraction(oilPhaseIdx, compIdx), b.moleFraction(oilPhaseIdx, compIdx));
        compare(a.moleFraction(gasPhaseIdx, compIdx), b.moleFraction(gasPhaseIdx, compIdx));
    }
    return largest;
}

} // namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(OptionsOffGiveTheAnswerWithoutMemory)
{
    for (const Scalar pressure : {75e5, 250e5}) {
        auto plain = makeState(pressure, feed);
        auto with_memory = plain;
        const bool plain_single = PTFlash::solve(plain, method, flashTolerance, eosType);
        Memory memory;
        const bool memory_single = PTFlash::solve(with_memory, nullptr, memory, Options{}, method,
                                                  flashTolerance, eosType);
        BOOST_CHECK_EQUAL(plain_single, memory_single);
        BOOST_CHECK(identicalAnswers(plain, with_memory));
        BOOST_CHECK(!memory.has_split && !memory.has_stability);
    }
}

BOOST_AUTO_TEST_CASE(PredictionOfATwoPhaseCell)
{
    Options options;
    options.predict = true;
    auto first = makeState(75e5, feed);
    Memory memory;
    BOOST_REQUIRE(!PTFlash::solve(first, nullptr, memory, options, method, flashTolerance, eosType));
    BOOST_REQUIRE(memory.has_split);
    BOOST_CHECK(!memory.split_has_temperature);

    // The remembered sensitivities are those of the derivative pass: reused
    // at the same state, they give its derivatives to round-off.
    {
        Options reuse = options;
        reuse.reuse = true;
        auto again = warmState(75e5, feed, first);
        const long accepted = PTFlash::memoryTally_().accepted;
        Memory next;
        BOOST_REQUIRE(!PTFlash::solve(again, &memory, next, reuse, method, flashTolerance, eosType));
        BOOST_CHECK_EQUAL(PTFlash::memoryTally_().accepted, accepted + 1);
        BOOST_CHECK_LT(largestRelativeDifference(again, first, true), 1e-10);
        BOOST_CHECK_LT(largestRelativeDifference(again, first, false), 1e-8);
    }

    // A small step is accepted from its prediction and agrees with a cold
    // flash of the same state.
    const Composition moved{0.1 + 1e-7, 0.3, 0.6 - 1e-7};
    const Scalar moved_pressure = 75e5 * (1. + 1e-7);
    auto cold = makeState(moved_pressure, moved);
    BOOST_REQUIRE(!PTFlash::solve(cold, method, flashTolerance, eosType));
    for (const bool reuse : {false, true}) {
        Options step_options = options;
        step_options.reuse = reuse;
        auto warm = warmState(moved_pressure, moved, first);
        const long accepted = PTFlash::memoryTally_().accepted;
        Memory next;
        BOOST_REQUIRE(!PTFlash::solve(warm, &memory, next, step_options, method, flashTolerance, eosType));
        BOOST_CHECK_EQUAL(PTFlash::memoryTally_().accepted, accepted + 1);
        BOOST_CHECK_LT(largestRelativeDifference(warm, cold, false), 1e-7);
        BOOST_CHECK_LT(largestRelativeDifference(warm, cold, true), reuse ? 1e-4 : 1e-6);
        BOOST_CHECK(next.has_split);
    }

    // A large step is not accepted and still reaches the cold answer.
    const Composition far{0.15, 0.3, 0.55};
    auto far_cold = makeState(80e5, far);
    BOOST_REQUIRE(!PTFlash::solve(far_cold, method, flashTolerance, eosType));
    auto far_warm = warmState(80e5, far, first);
    const long accepted = PTFlash::memoryTally_().accepted;
    Memory next;
    BOOST_REQUIRE(!PTFlash::solve(far_warm, &memory, next, options, method, flashTolerance, eosType));
    BOOST_CHECK_EQUAL(PTFlash::memoryTally_().accepted, accepted);
    BOOST_CHECK_LT(largestRelativeDifference(far_warm, far_cold, false), 1e-6);
}

BOOST_AUTO_TEST_CASE(PhaseMemoryOfASinglePhaseCell)
{
    Options options;
    options.phase_memory = true;
    auto first = makeState(250e5, feed);
    Memory memory;
    BOOST_REQUIRE(PTFlash::solve(first, nullptr, memory, options, method, flashTolerance, eosType));
    BOOST_REQUIRE(memory.has_stability);

    // Hardly moved: no test, the answer of the test bit for bit.
    const Composition moved{0.1 + 1e-6, 0.3, 0.6 - 1e-6};
    const Scalar moved_pressure = 250e5 * (1. - 1e-5);
    auto tested = warmState(moved_pressure, moved, first);
    BOOST_REQUIRE(PTFlash::solve(tested, method, flashTolerance, eosType));
    auto bypassed = warmState(moved_pressure, moved, first);
    const auto& tally = PTFlash::memoryTally_();
    const long bypasses = tally.bypass_trivial + tally.bypass_margin;
    Memory next;
    BOOST_REQUIRE(PTFlash::solve(bypassed, &memory, next, options, method, flashTolerance, eosType));
    BOOST_CHECK_EQUAL(tally.bypass_trivial + tally.bypass_margin, bypasses + 1);
    BOOST_CHECK(identicalAnswers(bypassed, tested));
    // The tested state stays the reference.
    BOOST_CHECK(next.has_stability);
    BOOST_CHECK_EQUAL(next.stability_pressure, memory.stability_pressure);

    // Moved into the two-phase region: tested and split.
    auto two_phase = warmState(75e5, feed, first);
    auto two_phase_plain = two_phase;
    BOOST_REQUIRE(!PTFlash::solve(two_phase_plain, method, flashTolerance, eosType));
    Memory after_split;
    BOOST_CHECK(!PTFlash::solve(two_phase, &next, after_split, options, method, flashTolerance, eosType));
    BOOST_CHECK_EQUAL(tally.bypass_trivial + tally.bypass_margin, bypasses + 1);
    BOOST_CHECK(identicalAnswers(two_phase, two_phase_plain));
    BOOST_CHECK(!after_split.has_stability);
}
