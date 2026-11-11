// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \brief The on-demand flash table of PTFlashTable.hpp in its exact mode only
 *        changes where the two-phase split starts: the answers are those of
 *        the flash without the table, with either interpolation, from a new
 *        and from a filled table. The table is cleared when the fluid
 *        system changes. In its approximate mode the answers keep the phase
 *        state and stay close to those of the flash.
 */
#include "config.h"

#define OPM_FLASH_STATISTICS 1
#define BOOST_TEST_MODULE PtFlashTable
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/FlashStatistics.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/constraintsolvers/PTFlashTable.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

using Scalar = double;
constexpr int numComponents = 3;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using Flash = Opm::PTFlash<Scalar, FluidSystem>;
using Table = Opm::PTFlashTable<Scalar, FluidSystem>;
using Settings = Opm::PTFlashTableSettings;
using EOSType = Opm::CompositionalConfig::EOSType;
using Composition = std::array<Scalar, numComponents>;

constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr Scalar flashTolerance = 1e-8;
constexpr EOSType eosType = EOSType::PR;

void makeFluid(const Scalar kij)
{
    using CompParam = typename FluidSystem::ComponentParam;
    FluidSystem::init();
    FluidSystem::addComponent(CompParam{"CO2", 44.0, 304.128, 73.773e5, 0.09412, 0.22394});
    FluidSystem::addComponent(CompParam{"C1", 16.04, 190.564, 45.992e5, 0.09863, 0.01142});
    FluidSystem::addComponent(CompParam{"C10", 142.28, 617.7, 21.03e5, 0.6098, 0.4884});
    FluidSystem::setInteractionCoefficients({0.1, 0.1, kij});
}

// A cold cell state with derivatives in pressure and the first two mole fractions.
FluidState makeState(const Scalar pressure, const Scalar temperature, const Composition& z)
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

struct Answer
{
    bool singlePhase;
    FluidState state;
};

// Two-phase and single-phase states of a CO2-methane-decane mixture.
std::vector<Answer> flashAll(const Opm::PTFlashTableMode mode)
{
    Settings::get().mode = mode;
    std::vector<Answer> answers;
    for (const Composition z : {Composition{0.4, 0.3, 0.3}, Composition{0.6, 0.3, 0.1}}) {
        for (Scalar pressure = 20e5; pressure <= 200e5; pressure += 15e5) {
            for (Scalar temperature = 280.0; temperature <= 420.0; temperature += 10.0) {
                auto fs = makeState(pressure, temperature, z);
                const bool single = Flash::solve(fs, Opm::PTFlashMethod::SsiNewton, flashTolerance, eosType);
                answers.push_back({single, fs});
            }
        }
    }
    Settings::get().mode = Opm::PTFlashTableMode::Off;
    return answers;
}

// The same answers to the flash tolerance: the splits converge from different starts.
void checkSame(const std::vector<Answer>& expected, const std::vector<Answer>& actual)
{
    BOOST_REQUIRE_EQUAL(expected.size(), actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        BOOST_TEST_CONTEXT("state " << i)
        {
            BOOST_CHECK_EQUAL(expected[i].singlePhase, actual[i].singlePhase);
            const auto& a = expected[i].state;
            const auto& b = actual[i].state;
            BOOST_CHECK_SMALL(Opm::getValue(a.L()) - Opm::getValue(b.L()), 1e-6);
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                for (const auto phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
                    const auto& xa = a.moleFraction(phaseIdx, compIdx);
                    const auto& xb = b.moleFraction(phaseIdx, compIdx);
                    BOOST_CHECK_SMALL(Opm::getValue(xa) - Opm::getValue(xb), 1e-6);
                    // The derivatives are reconstructed at the converged split.
                    for (int d = 0; d < Evaluation::numVars; ++d) {
                        BOOST_CHECK_SMALL(xa.derivative(d) - xb.derivative(d),
                                          1e-5 * (1.0 + std::abs(xa.derivative(d))));
                    }
                }
            }
        }
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(ExactTableGivesTheFlashAnswers)
{
    for (const auto interpolation : {Settings::Interpolation::Multilinear, Settings::Interpolation::Simplex}) {
        makeFluid(0.0411);
        Settings::get().interpolation = interpolation;
        Table::get().reset();
        const auto expected = flashAll(Opm::PTFlashTableMode::Off);

        Opm::FlashStatistics::get().reset();
        const auto fromNewTable = flashAll(Opm::PTFlashTableMode::Exact);
        const long newTableHits = Opm::FlashStatistics::get().cacheHits;
        checkSame(expected, fromNewTable);

        Opm::FlashStatistics::get().reset();
        const auto fromFilledTable = flashAll(Opm::PTFlashTableMode::Exact);
        checkSame(expected, fromFilledTable);
        // Every two-phase state away from the phase boundaries starts from the table.
        BOOST_CHECK_GT(Opm::FlashStatistics::get().cacheHits, 0);
        BOOST_CHECK_GE(Opm::FlashStatistics::get().cacheHits, newTableHits);
        BOOST_CHECK_GT(Table::get().vertices(), 0u);
    }
    Settings::get().interpolation = Settings::Interpolation::Automatic;
}

BOOST_AUTO_TEST_CASE(TableIsClearedWhenTheFluidChanges)
{
    makeFluid(0.0411);
    Table::get().reset();
    flashAll(Opm::PTFlashTableMode::Exact);
    const auto filled = Table::get().vertices();
    BOOST_CHECK_GT(filled, 16u);

    // One flash of another fluid: at most the 2^4 vertices of one cell.
    makeFluid(0.0);
    Settings::get().mode = Opm::PTFlashTableMode::Exact;
    auto fs = makeState(100e5, 320.0, {0.4, 0.3, 0.3});
    Flash::solve(fs, Opm::PTFlashMethod::SsiNewton, flashTolerance, eosType);
    Settings::get().mode = Opm::PTFlashTableMode::Off;
    BOOST_CHECK_LE(Table::get().vertices(), 16u);
}

BOOST_AUTO_TEST_CASE(ApproximateTableStaysClose)
{
    makeFluid(0.0411);
    Settings::get().interpolation = Settings::Interpolation::Multilinear;
    Table::get().reset();
    const auto expected = flashAll(Opm::PTFlashTableMode::Off);
    flashAll(Opm::PTFlashTableMode::Approximate);
    Opm::FlashStatistics::get().reset();
    const auto approximate = flashAll(Opm::PTFlashTableMode::Approximate);
    // From the filled table, most answers are taken without iterating.
    BOOST_CHECK_GT(2 * Opm::FlashStatistics::get().cacheHits, static_cast<long>(expected.size()));

    int changedPhaseState = 0;
    Scalar maxDifference = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (expected[i].singlePhase != approximate[i].singlePhase) {
            ++changedPhaseState;
            continue;
        }
        const auto& a = expected[i].state;
        const auto& b = approximate[i].state;
        maxDifference = std::max(maxDifference, std::abs(Opm::getValue(a.L()) - Opm::getValue(b.L())));
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (const auto phaseIdx : {oilPhaseIdx, gasPhaseIdx}) {
                const auto& x = b.moleFraction(phaseIdx, compIdx);
                maxDifference = std::max(maxDifference,
                                         std::abs(Opm::getValue(a.moleFraction(phaseIdx, compIdx))
                                                  - Opm::getValue(x)));
                for (int d = 0; d < Evaluation::numVars; ++d) {
                    BOOST_CHECK(std::isfinite(x.derivative(d)));
                }
            }
        }
    }
    BOOST_TEST_MESSAGE("approximate: phase state changed " << changedPhaseState << ", largest difference "
                       << maxDifference);
    BOOST_CHECK_EQUAL(changedPhaseState, 0);
    BOOST_CHECK_LT(maxDifference, 0.05);
    BOOST_CHECK_GT(maxDifference, 0.0);
    Settings::get().interpolation = Settings::Interpolation::Automatic;
}

BOOST_AUTO_TEST_CASE(AcceptanceTestBoundsTheApproximateAnswers)
{
    // A tight acceptance converges every two-phase answer the table cannot
    // give to that accuracy; single-phase cells still answer from the table.
    makeFluid(0.0411);
    Settings::get().interpolation = Settings::Interpolation::Multilinear;
    Settings::get().acceptance = 1e-6;
    Table::get().reset();
    const auto expected = flashAll(Opm::PTFlashTableMode::Off);
    flashAll(Opm::PTFlashTableMode::Approximate);
    const auto guarded = flashAll(Opm::PTFlashTableMode::Approximate);
    checkSame(expected, guarded);
    Settings::get().acceptance = 0.0;
    Settings::get().interpolation = Settings::Interpolation::Automatic;
}

BOOST_AUTO_TEST_CASE(NegativeFlashVerticesKeepTheAnswers)
{
    // Boundary cells: single-phase vertices next to two-phase ones store a
    // negative flash. With a tight acceptance, a state inside such a cell is
    // either split from the interpolated K, kept only with a Gibbs energy
    // below the feed's, or declared single-phase when its interpolated L is
    // beyond the margin: the answers stay those of the flash.
    makeFluid(0.0411);
    Settings::get().interpolation = Settings::Interpolation::Multilinear;
    Settings::get().acceptance = 1e-6;
    Settings::get().negativeFlash = true;
    Table::get().reset();
    const auto expected = flashAll(Opm::PTFlashTableMode::Off);
    flashAll(Opm::PTFlashTableMode::Approximate);
    const auto boundary = flashAll(Opm::PTFlashTableMode::Approximate);
    checkSame(expected, boundary);
    using Outcome = Table::Outcome;
    BOOST_TEST_MESSAGE("boundary splits " << Table::get().count(Outcome::BoundarySplit)
                       << ", boundary single-phase " << Table::get().count(Outcome::BoundarySingle));
    BOOST_CHECK_GT(Table::get().count(Outcome::BoundarySplit), 0);
    BOOST_CHECK_GT(Table::get().count(Outcome::BoundarySingle), 0);
    // The exact mode keeps only the proven splits of boundary cells.
    Table::get().reset();
    flashAll(Opm::PTFlashTableMode::Exact);
    checkSame(expected, flashAll(Opm::PTFlashTableMode::Exact));
    BOOST_CHECK_EQUAL(Table::get().count(Outcome::BoundarySingle), 0);
    Settings::get().negativeFlash = false;
    Settings::get().acceptance = 0.0;
    Settings::get().interpolation = Settings::Interpolation::Automatic;
}
