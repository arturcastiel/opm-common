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
 * \brief Newton must finish the split its start points to, near a critical
 *        point too.
 *
 * The SPE5 oil (six components, Peng-Robinson) close to its critical point.
 * Successive substitution hands the split to Newton after a few steps, still
 * far from convergence, and Newton's first step overshoots. With a single
 * cubic root per phase the equations do not tell the phase labels apart, so
 * Newton can end
 * - on the mirror split, the right phases with oil and gas exchanged
 *   (99.65 bar, 625.59 K: L = 0.2125 instead of 0.7875), or
 * - on coincident phases, which the flash rejects after a stability test
 *   found the feed unstable (82.88 bar, 640.72 K: the flash failed).
 * Both answers depend on the last bit of the feed: the same state with the
 * feed's last mole fraction written as one minus the others took the other
 * branch. The reference answers are those of MultiComponentFlash.jl
 * (plain substitution to 1e-8, no Newton).
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashNewtonSplitSide
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <cmath>
#include <utility>
#include <vector>

using Scalar = double;
constexpr int numComponents = 6;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
using PtFlash = Opm::PTFlash<Scalar, FluidSystem>;

namespace {

const std::array<Scalar, numComponents> feed{0.5, 0.03, 0.07, 0.20, 0.15, 0.05};

//! The SPE5 oil of MultiComponentFlash.jl's cubic_benchmark("spe5").
struct Fixture
{
    Fixture()
    {
        const Scalar psi = 6.894757293168360e+03;
        const Scalar rankine = 0.555555555555556;
        const std::array<const char*, numComponents> names{"C1", "C3", "C6", "C10", "C15", "C20"};
        const std::array<Scalar, numComponents> tcR{343, 665.7, 913.4, 1111.8, 1270.0, 1380.0};
        const std::array<Scalar, numComponents> pcPsi{667.8, 616.3, 436.9, 304.0, 200.0, 162.0};
        const std::array<Scalar, numComponents> mw{16.040, 44.100, 86.180, 142.290, 206.000, 282.000};
        const std::array<Scalar, numComponents> acentric{0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500};
        const std::array<Scalar, numComponents> zc{0.290, 0.277, 0.264, 0.257, 0.245, 0.235};
        FluidSystem::init();
        for (int i = 0; i < numComponents; ++i) {
            const Scalar tc = tcR[i] * rankine;
            const Scalar pc = pcPsi[i] * psi;
            const Scalar vc = zc[i] * 8.314 * tc / pc * 1e3;
            FluidSystem::addComponent(
                FluidSystem::ComponentParam{names[i], mw[i], tc, pc, vc, acentric[i]});
        }
        std::array<std::array<Scalar, numComponents>, numComponents> k{};
        k[0][4] = k[4][0] = 0.05;
        k[0][5] = k[5][0] = 0.05;
        k[1][4] = k[4][1] = 0.005;
        k[1][5] = k[5][1] = 0.005;
        std::vector<Scalar> bic;
        for (int row = 1; row < numComponents; ++row) {
            for (int col = 0; col < row; ++col) {
                bic.push_back(k[row][col]);
            }
        }
        FluidSystem::setInteractionCoefficients(std::move(bic));
    }
};

//! The state (ip, it) of the 1000 x 1000 map: p 1 - 400 bar, T 280 - 700 K.
std::pair<Scalar, Scalar> mapState(const int ip, const int it)
{
    const Scalar pressure = 1e5 + (4e7 - 1e5) * ip / 999;
    const Scalar temperature = 280.0 + (700.0 - 280.0) * it / 999;
    return {pressure, temperature};
}

struct Answer
{
    bool single_phase;
    Scalar L;
    std::array<Scalar, numComponents> x;
    std::array<Scalar, numComponents> y;
};

//! A cold flash (Wilson ratios, no split) of the feed, values only.
Answer flashValues(const Scalar p, const Scalar T)
{
    Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fs.setPressure(phaseIdx, p);
    }
    fs.setTemperature(T);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setMoleFraction(compIdx, feed[compIdx]);
    }
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    Answer answer{};
    answer.single_phase = PtFlash::flash_solve_scalar(fs, Opm::PTFlashMethod::SsiNewton, 1e-8,
                                                      Opm::CompositionalConfig::EOSType::PR);
    answer.L = fs.L();
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        answer.x[compIdx] = fs.moleFraction(FluidSystem::oilPhaseIdx, compIdx);
        answer.y[compIdx] = fs.moleFraction(FluidSystem::gasPhaseIdx, compIdx);
    }
    return answer;
}

//! The same flash with derivatives in p and z_0 .. z_4; z_5 = 1 - sum.
Answer flashDerivatives(const Scalar p, const Scalar T)
{
    using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents + 1>;
    Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fs.setPressure(phaseIdx, Evaluation::createVariable(p, 0));
    }
    fs.setTemperature(Evaluation(T));
    Evaluation last = 1.0;
    for (int compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
        const auto z = Evaluation::createVariable(feed[compIdx], compIdx + 1);
        fs.setMoleFraction(compIdx, z);
        last -= z;
    }
    fs.setMoleFraction(numComponents - 1, last);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    Answer answer{};
    answer.single_phase = PtFlash::solve(fs, Opm::PTFlashMethod::SsiNewton, 1e-8,
                                         Opm::CompositionalConfig::EOSType::PR);
    answer.L = fs.L().value();
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        answer.x[compIdx] = fs.moleFraction(FluidSystem::oilPhaseIdx, compIdx).value();
        answer.y[compIdx] = fs.moleFraction(FluidSystem::gasPhaseIdx, compIdx).value();
    }
    return answer;
}

//! Both paths give MultiComponentFlash.jl's split, with methane in the gas.
void checkSplit(const Answer& answer, const Scalar reference_L, const Scalar reference_x_c1)
{
    BOOST_REQUIRE(!answer.single_phase);
    BOOST_CHECK_SMALL(answer.L - reference_L, 1e-4);
    BOOST_CHECK_SMALL(answer.x[0] - reference_x_c1, 1e-4);
    BOOST_CHECK_LT(answer.x[0], answer.y[0]);
    BOOST_CHECK_GT(answer.x[5], answer.y[5]);
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(NewtonDoesNotExchangeThePhases)
{
    // 99.65 bar, 625.59 K: L = 0.78755, x_C1 = 0.47912, y_C1 = 0.57740.
    const auto [p, T] = mapState(247, 822);
    const Answer values = flashValues(p, T);
    const Answer derivatives = flashDerivatives(p, T);
    BOOST_TEST_CONTEXT("values")
    {
        checkSplit(values, 0.787547161733, 0.479120161651);
    }
    BOOST_TEST_CONTEXT("derivatives")
    {
        checkSplit(derivatives, 0.787547161733, 0.479120161651);
    }
    BOOST_CHECK_SMALL(values.L - derivatives.L, 1e-8);
}

BOOST_AUTO_TEST_CASE(NewtonDoesNotEndOnCoincidentPhases)
{
    // 82.88 bar, 640.72 K: L = 0.19836, x_C1 = 0.42971, y_C1 = 0.51739.
    const auto [p, T] = mapState(205, 858);
    Answer values{};
    Answer derivatives{};
    BOOST_REQUIRE_NO_THROW(values = flashValues(p, T));
    BOOST_REQUIRE_NO_THROW(derivatives = flashDerivatives(p, T));
    BOOST_TEST_CONTEXT("values")
    {
        checkSplit(values, 0.198362191068, 0.429713285858);
    }
    BOOST_TEST_CONTEXT("derivatives")
    {
        checkSplit(derivatives, 0.198362191068, 0.429713285858);
    }
    BOOST_CHECK_SMALL(values.L - derivatives.L, 1e-8);
}
