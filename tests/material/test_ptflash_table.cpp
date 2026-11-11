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

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/constraintsolvers/FlashStatistics.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PHFlash.hpp>
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
    flashAll(Opm::PTFlashTableMode::Skip);
    const auto boundary = flashAll(Opm::PTFlashTableMode::Skip);
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
    Settings::get().acceptance = 0.0;
    Settings::get().interpolation = Settings::Interpolation::Automatic;
}

BOOST_AUTO_TEST_CASE(EnthalpyOfTheVerticesStartsThePHSearch)
{
    // A P-H flash without a phase state from an earlier flash starts its
    // temperature search where the enthalpy of the table's vertices matches:
    // the same answers, from fewer isothermal flashes than from the middle
    // of the interval.
    using PhFlash = Opm::PHFlash<Scalar, FluidSystem>;
    using Enthalpy = Opm::MixtureEnthalpy<Scalar, FluidSystem>;
    makeFluid(0.0411);
    PhFlash::Config cfg;
    cfg.cpTable = {Opm::SimpleCO2<Scalar>::idealGasHeatCapacityPolynomial(),
                   Opm::C1<Scalar>::idealGasHeatCapacityPolynomial(),
                   Opm::C10<Scalar>::idealGasHeatCapacityPolynomial()};
    cfg.model = Opm::EnthalpyModel::eos_departure;
    cfg.tempMin = 250.0;
    cfg.tempMax = 600.0;
    cfg.temperatureTolerance = 1e-6;
    cfg.enthalpyTolerance = 1.0;
    cfg.maxIterations = 100;
    struct PhAnswer
    {
        bool ok;
        Scalar temperature;
        Scalar L;
    };
    const auto solveAll = [&](const Settings::EnthalpyStart start, long& evaluations) {
        Settings::get().mode = start == Settings::EnthalpyStart::Never
            ? Opm::PTFlashTableMode::Off : Opm::PTFlashTableMode::Approximate;
        Settings::get().enthalpyStart = start;
        Opm::FlashStatistics::get().reset();
        std::vector<PhAnswer> answers;
        for (Scalar pressure = 20e5; pressure <= 200e5; pressure += 30e5) {
            for (Scalar temperature = 280.0; temperature <= 420.0; temperature += 20.0) {
                const Composition z{0.4, 0.3, 0.3};
                auto reference = makeState(pressure, temperature, z);
                Flash::solve(reference, Opm::PTFlashMethod::SsiNewton, flashTolerance, eosType);
                const Scalar hSpec
                    = Enthalpy::mixtureEnthalpy(reference, cfg.cpTable, eosType, cfg.model);
                auto fs = makeState(pressure, 425.0, z);
                const bool ok = PhFlash::solve(fs, hSpec, cfg, Opm::PTFlashMethod::SsiNewton,
                                               flashTolerance, eosType);
                answers.push_back({ok, Opm::getValue(fs.temperature(0)), Opm::getValue(fs.L())});
            }
        }
        evaluations = Opm::FlashStatistics::get().phEvaluations;
        Settings::get().mode = Opm::PTFlashTableMode::Off;
        Settings::get().enthalpyStart = Settings::EnthalpyStart::Cold;
        return answers;
    };
    Settings::get().acceptance = 1e-6;
    Table::get().reset();
    long fromMiddle = 0;
    long fromTable = 0;
    const auto expected = solveAll(Settings::EnthalpyStart::Never, fromMiddle);
    const auto started = solveAll(Settings::EnthalpyStart::Cold, fromTable);
    BOOST_REQUIRE_EQUAL(expected.size(), started.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        BOOST_TEST_CONTEXT("state " << i)
        {
            BOOST_CHECK(expected[i].ok);
            BOOST_CHECK_EQUAL(expected[i].ok, started[i].ok);
            BOOST_CHECK_SMALL(expected[i].temperature - started[i].temperature, 1e-5);
            BOOST_CHECK_SMALL(expected[i].L - started[i].L, 1e-5);
        }
    }
    BOOST_TEST_MESSAGE("isothermal flashes per P-H flash: " << double(fromMiddle) / expected.size()
                       << " from the middle, " << double(fromTable) / expected.size()
                       << " from the table");
    BOOST_CHECK_LT(fromTable, fromMiddle);
    Settings::get().acceptance = 0.0;
}

BOOST_AUTO_TEST_CASE(SplitsFromAnEarlierFlashKeepTheirRatios)
{
    // A flash seeded with the split of an earlier flash (a Newton iteration
    // of a simulator: 0.1 % higher pressure) starts from its own ratios and
    // does not ask the table: the flash without the table, bit for bit.
    // Asked to, it starts from the table, with the same answers.
    makeFluid(0.0411);
    Table::get().reset();
    const auto cold = flashAll(Opm::PTFlashTableMode::Off);
    const auto warmAll = [&](const Opm::PTFlashTableMode mode, const bool fromTable) {
        Settings::get().mode = mode;
        Settings::get().warmStartFromTable = fromTable;
        std::vector<Answer> answers;
        for (const auto& answer : cold) {
            if (answer.singlePhase) {
                continue;
            }
            auto fs = answer.state;
            const auto p
                = Evaluation::createVariable(1.001 * Opm::getValue(fs.pressure(oilPhaseIdx)), 0);
            fs.setPressure(oilPhaseIdx, p);
            fs.setPressure(gasPhaseIdx, p);
            const bool single
                = Flash::solve(fs, Opm::PTFlashMethod::SsiNewton, flashTolerance, eosType);
            answers.push_back({single, fs});
        }
        Settings::get().mode = Opm::PTFlashTableMode::Off;
        Settings::get().warmStartFromTable = false;
        return answers;
    };
    const auto queries = [] {
        using Outcome = Table::Outcome;
        long sum = 0;
        for (int o = 0; o < static_cast<int>(Outcome::Count); ++o) {
            sum += Table::get().count(static_cast<Outcome>(o));
        }
        return sum;
    };
    const auto expected = warmAll(Opm::PTFlashTableMode::Off, false);
    BOOST_REQUIRE_GT(expected.size(), 50u);

    Table::get().reset();
    const long before = queries();
    const auto own = warmAll(Opm::PTFlashTableMode::Skip, false);
    BOOST_CHECK_EQUAL(queries(), before);
    BOOST_REQUIRE_EQUAL(expected.size(), own.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        BOOST_CHECK_EQUAL(expected[i].singlePhase, own[i].singlePhase);
        BOOST_CHECK_EQUAL(Opm::getValue(expected[i].state.L()), Opm::getValue(own[i].state.L()));
    }

    const auto fromTable = warmAll(Opm::PTFlashTableMode::Skip, true);
    BOOST_CHECK_GT(queries(), 0);
    checkSame(expected, fromTable);
}

BOOST_AUTO_TEST_CASE(OptionsAreOffByDefault)
{
    // Without the environment every option is off, so the flash is the
    // flash without them; the table's sub-options are on for when it is.
    const auto defaults = Settings{};
    BOOST_CHECK(defaults.mode == Opm::PTFlashTableMode::Off);
    BOOST_CHECK(!defaults.stabilityBypass);
    BOOST_CHECK(defaults.negativeFlash);
    BOOST_CHECK(defaults.enthalpyStart == Settings::EnthalpyStart::Cold);
    BOOST_CHECK(!defaults.warmStartFromTable);
}
