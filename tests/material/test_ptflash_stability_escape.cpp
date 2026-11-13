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
 * \brief Stability trials that pass close to a stationary point they do not
 *        reach must not take thousands of steps.
 *
 * The SPE5 oil (six components, Peng-Robinson) is single-phase just outside
 * its two-phase region near the critical point, at 75-110 bar and 620-650 K.
 * There a trial phase grown from Wilson's ratios slows down next to a
 * stationary point of the tangent-plane distance, or in the bottleneck where
 * one has just vanished: its substitution steps grow or shrink by a factor
 * close to one, the extrapolation of the remaining steps overshoots, and
 * plain substitution took up to 2,219 steps to leave. The trial now moves on
 * along its step while the tangent-plane distance drops. The answers are
 * those of the slow trials: single phase, labelled by the stable root.
 */
#include "config.h"

#define OPM_FLASH_STATISTICS 1

#define BOOST_TEST_MODULE PtFlashStabilityEscape
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/FlashStatistics.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <vector>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
constexpr int numComponents = 6;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, numComponents, false>;
using FluidState = Opm::CompositionalFluidState<Scalar, FluidSystem>;
using PtFlash = Opm::PTFlash<Scalar, FluidSystem>;

namespace {

// The SPE5 oil: C1, C3, C6, C10, C15, C20 (g/mol, K, Pa, m3/kmol).
struct Fixture
{
    Fixture()
    {
        const Scalar psi = 6.894757293168360e+03;
        const Scalar rankine = 0.555555555555556;
        const std::array<const char*, numComponents> names{"C1", "C3", "C6", "C10", "C15", "C20"};
        const std::array<Scalar, numComponents> tc{343, 665.7, 913.4, 1111.8, 1270.0, 1380.0};
        const std::array<Scalar, numComponents> pc{667.8, 616.3, 436.9, 304.0, 200.0, 162.0};
        const std::array<Scalar, numComponents> mw{16.040, 44.100, 86.180, 142.290, 206.000, 282.000};
        const std::array<Scalar, numComponents> omega{0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500};
        const std::array<Scalar, numComponents> zc{0.290, 0.277, 0.264, 0.257, 0.245, 0.235};
        FluidSystem::init();
        for (int i = 0; i < numComponents; ++i) {
            const Scalar t = tc[i] * rankine;
            const Scalar p = pc[i] * psi;
            FluidSystem::addComponent(FluidSystem::ComponentParam{
                names[i], mw[i], t, p, zc[i] * 8.314 * t / p * 1e3, omega[i]});
        }
        // Lower triangle, row by row: C1-C15, C1-C20 0.05; C3-C15, C3-C20 0.005.
        std::vector<Scalar> kij(numComponents * (numComponents - 1) / 2, 0.0);
        kij[6] = 0.05;  // (4, 0)
        kij[7] = 0.005; // (4, 1)
        kij[10] = 0.05;  // (5, 0)
        kij[11] = 0.005; // (5, 1)
        FluidSystem::setInteractionCoefficients(std::move(kij));
    }
};

// A state of the dense SPE5 map, 1000 x 1000 over 1-400 bar and 280-700 K.
struct MapState
{
    int ip, it;
    Scalar L; // the label of the single phase
};

// A cold flash of the SPE5 feed at a map state: single phase with the label
// state.L, after one stability test of fewer than max_steps steps.
void checkShortStableFlash(const MapState& state, const long max_steps)
{
    const std::array<Scalar, numComponents> z{0.5, 0.03, 0.07, 0.20, 0.15, 0.05};
    const Scalar p = 1e5 + (4e7 - 1e5) * state.ip / 999;
    const Scalar T = 280.0 + 420.0 * state.it / 999;
    BOOST_TEST_CONTEXT("p " << p << " Pa, T " << T << " K")
    {
        FluidState fs;
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, p);
        }
        fs.setTemperature(T);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fs.setMoleFraction(compIdx, z[compIdx]);
        }
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
        }
        fs.setLvalue(-1.0);

        auto& statistics = Opm::FlashStatistics::get();
        statistics.reset();
        const bool single
            = PtFlash::flash_solve_scalar(fs, Opm::PTFlashMethod::SsiNewton, 1e-8, EOSType::PR, 0);

        BOOST_CHECK(single);
        BOOST_CHECK_EQUAL(fs.L(), state.L);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            BOOST_CHECK_EQUAL(fs.moleFraction(FluidSystem::oilPhaseIdx, compIdx), z[compIdx]);
            BOOST_CHECK_EQUAL(fs.moleFraction(FluidSystem::gasPhaseIdx, compIdx), z[compIdx]);
        }
        // One stability test, both trials together, in fewer than max_steps steps.
        BOOST_CHECK_EQUAL(statistics.stabilityTests, 1);
        BOOST_CHECK_LT(statistics.stabilityIterations, max_steps);
    }
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(NearStationaryTrialsAreShort)
{
    // Slow vapour-like trials (2,219, 1,745, 1,325 and 1,008 steps of plain
    // substitution with extrapolation), a slow liquid-like one (1,120), and a
    // supercritical liquid (1,035).
    for (const auto& state : std::array<MapState, 6>{{{243, 831, 0.0}, {234, 839, 0.0},
                                                      {268, 807, 0.0}, {233, 840, 0.0},
                                                      {206, 861, 0.0}, {378, 673, 1.0}}}) {
        checkShortStableFlash(state, 300);
    }
}
