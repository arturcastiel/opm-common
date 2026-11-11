// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
 * \brief When the P-T flash skips its stability test.
 *
 *        A cold start first tries the split; a split below the Gibbs energy
 *        of the feed proves the feed unstable, and the test is skipped. The
 *        split must be the one the stability test leads to.
 *
 *        A state that comes in single-phase with the ratios of a state whose
 *        test had both trials collapse stays single-phase without a test
 *        while it is inside the critical-distance window of that state, and
 *        is tested again outside it.
 */
#include "config.h"

#define OPM_FLASH_STATISTICS 1

#define BOOST_TEST_MODULE PtFlashSkip
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

namespace {

using Scalar = double;
constexpr int numComponents = 3;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using EOSType = Opm::CompositionalConfig::EOSType;
using Flash = Opm::PTFlash<Scalar, FluidSystem>;
using Composition = std::array<Scalar, numComponents>;

constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr Scalar flashTolerance = 1e-8;
constexpr EOSType eosType = EOSType::PR;

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

// A cell state with the given hint: ratios K and liquid fraction L.
FluidState makeState(const Scalar pressure,
                     const Scalar temperature,
                     const Composition& z,
                     const Scalar L,
                     const std::array<Scalar, numComponents>* K = nullptr)
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
        fs.setKvalue(compIdx, K ? Evaluation((*K)[compIdx]) : fs.wilsonK_(compIdx));
    }
    fs.setLvalue(L);
    return fs;
}

std::array<Scalar, numComponents> ratios(const FluidState& fs)
{
    std::array<Scalar, numComponents> K{};
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        K[compIdx] = Opm::getValue(fs.K(compIdx));
    }
    return K;
}

// Flash and return the counters of this flash alone.
Opm::FlashStatistics flash(FluidState& fs, bool& single_phase)
{
    Opm::FlashStatistics::get().reset();
    single_phase = Flash::solve(fs, Opm::PTFlashMethod::SsiNewton, flashTolerance, eosType);
    return Opm::FlashStatistics::get();
}

} // namespace

BOOST_FIXTURE_TEST_SUITE(PtFlashSkip, Fixture)

// A two-phase feed flashed cold skips the stability test, and its split is the
// one the stability test leads to (here reached through a single-phase hint,
// which always runs the test).
BOOST_AUTO_TEST_CASE(FlashFirstGivesTheSplitOfTheStabilityTest)
{
    const Composition z{0.4, 0.3, 0.3};
    const Scalar pressure = 50e5;
    const Scalar temperature = 300.0;

    bool cold_single = true;
    auto cold = makeState(pressure, temperature, z, -1.0);
    const auto cold_work = flash(cold, cold_single);
    BOOST_REQUIRE(!cold_single);
    BOOST_CHECK_EQUAL(cold_work.stabilityTests, 0);

    bool tested_single = true;
    auto tested = makeState(pressure, temperature, z, 1.0);
    const auto tested_work = flash(tested, tested_single);
    BOOST_REQUIRE(!tested_single);
    BOOST_CHECK_EQUAL(tested_work.stabilityTests, 1);

    BOOST_CHECK_SMALL(Opm::getValue(cold.L()) - Opm::getValue(tested.L()), 1e-7);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        BOOST_CHECK_SMALL(Opm::getValue(cold.moleFraction(oilPhaseIdx, compIdx))
                              - Opm::getValue(tested.moleFraction(oilPhaseIdx, compIdx)),
                          1e-7);
        BOOST_CHECK_SMALL(Opm::getValue(cold.moleFraction(gasPhaseIdx, compIdx))
                              - Opm::getValue(tested.moleFraction(gasPhaseIdx, compIdx)),
                          1e-7);
    }
}

// A single-phase feed flashed cold still gets its stability test: the split
// cannot prove anything there.
BOOST_AUTO_TEST_CASE(SinglePhaseColdStartIsTested)
{
    bool single = false;
    auto fs = makeState(20e5, 450.0, Composition{0.1, 0.85, 0.05}, -1.0);
    const auto work = flash(fs, single);
    BOOST_CHECK(single);
    BOOST_CHECK_EQUAL(work.stabilityTests, 1);
}

// The bypass: near the reference no test, far from it a test, and a hint
// with other ratios finds no reference.
BOOST_AUTO_TEST_CASE(BypassOnlyInsideTheWindowOfItsReference)
{
    const Composition z{0.1, 0.89, 0.01};
    const Scalar pressure = 20e5;
    const Scalar temperature = 500.0;

    // The reference: a deep vapour, tested, both trials collapse.
    bool single = false;
    auto reference = makeState(pressure, temperature, z, 0.0);
    const auto reference_work = flash(reference, single);
    BOOST_REQUIRE(single);
    BOOST_CHECK_EQUAL(reference_work.stabilityTests, 1);
    BOOST_CHECK_EQUAL(reference_work.cacheHits, 0);
    const auto K = ratios(reference);
    const Scalar label = Opm::getValue(reference.L());

    // A small step with the reference's hint: single-phase without a test.
    auto near = makeState(pressure * 1.001, temperature + 0.05,
                          Composition{0.1005, 0.8895, 0.01}, label, &K);
    const auto near_work = flash(near, single);
    BOOST_CHECK(single);
    BOOST_CHECK_EQUAL(near_work.stabilityTests, 0);
    BOOST_CHECK_EQUAL(near_work.cacheHits, 1);
    BOOST_CHECK_EQUAL(Opm::getValue(near.L()), label);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        BOOST_CHECK_EQUAL(Opm::getValue(near.moleFraction(oilPhaseIdx, compIdx)),
                          Opm::getValue(near.moleFraction(compIdx)));
        BOOST_CHECK_EQUAL(Opm::getValue(near.moleFraction(gasPhaseIdx, compIdx)),
                          Opm::getValue(near.moleFraction(compIdx)));
    }

    // Far in pressure with the same hint: tested.
    auto far = makeState(pressure * 3.0, temperature, z, label, &K);
    const auto far_work = flash(far, single);
    BOOST_CHECK_EQUAL(far_work.stabilityTests, 1);
    BOOST_CHECK_EQUAL(far_work.cacheHits, 0);

    // A two-phase state with the hint of a single-phase one: tested, two-phase.
    auto two_phase = makeState(50e5, 300.0, Composition{0.4, 0.3, 0.3}, 1.0);
    const auto split_work = flash(two_phase, single);
    BOOST_CHECK(!single);
    BOOST_CHECK_EQUAL(split_work.stabilityTests, 1);
}

BOOST_AUTO_TEST_SUITE_END()
