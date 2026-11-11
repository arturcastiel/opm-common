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
 * \brief The equilibrium ratios PTFlash::solve() returns, and a flash
 *        started from them.
 *
 * The simulator passes the K and L of a cell's last flash into its next one.
 * A two-phase result returns its converged ratios, a single phase those of
 * the stability trial that ended closest to instability. Ratios at the
 * trivial solution, as a near-critical negative flash leaves them, must not
 * stop the next stability test at once.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashWarmStart
#include <boost/test/unit_test.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <cmath>
#include <optional>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using PTFlashMethod = Opm::PTFlashMethod;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 2, false>;

namespace {

template <class Component>
void registerComponent()
{
    using CompParam = typename FluidSystem::ComponentParam;
    FluidSystem::addComponent(CompParam{Component::name(),
                                        Component::molarMass(),
                                        Component::criticalTemperature(),
                                        Component::criticalPressure(),
                                        Component::criticalVolume(),
                                        Component::acentricFactor()});
}

// Methane and decane; the component parameters are static, so they are
// registered once per process.
struct Fixture
{
    Fixture()
    {
        FluidSystem::init();
        registerComponent<Opm::C1<Scalar>>();
        registerComponent<Opm::C10<Scalar>>();
        FluidSystem::setInteractionCoefficients({0.0411});
    }
};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

namespace {

constexpr int numComponents = FluidSystem::numComponents;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents + 1>;
using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
using PtFlash = Opm::PTFlash<Scalar, FluidSystem>;
constexpr Scalar flashTolerance = 1e-8;

// A cell state with Wilson ratios and no previous split, as the simulator
// starts a cell without a hint.
FluidState makeState(const Scalar p, const Scalar T, const Scalar z0)
{
    FluidState fs;
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fs.setPressure(phaseIdx, Evaluation::createVariable(p, 0));
    }
    fs.setTemperature(Evaluation(T));
    fs.setMoleFraction(0, Evaluation::createVariable(z0, 1));
    fs.setMoleFraction(1, 1.0 - fs.moleFraction(0));
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    return fs;
}

// The K and L a flash left, as the simulator hands them to the next one.
void seed(FluidState& fs, const FluidState& from)
{
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setKvalue(compIdx, Opm::getValue(from.K(compIdx)));
    }
    fs.setLvalue(Opm::getValue(from.L()));
}

bool flash(FluidState& fs)
{
    return PtFlash::solve(fs, PTFlashMethod::SsiNewton, flashTolerance, EOSType::PR);
}

} // anonymous namespace

// A two-phase result returns its converged ratios y / x.
BOOST_AUTO_TEST_CASE(TwoPhaseReturnsConvergedRatios)
{
    FluidState fs = makeState(50e5, 400.0, 0.5);
    BOOST_REQUIRE(!flash(fs));
    const Scalar L = Opm::getValue(fs.L());
    BOOST_REQUIRE(L > 0.0 && L < 1.0);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        const Scalar x = Opm::getValue(fs.moleFraction(FluidSystem::oilPhaseIdx, compIdx));
        const Scalar y = Opm::getValue(fs.moleFraction(FluidSystem::gasPhaseIdx, compIdx));
        BOOST_CHECK_CLOSE(Opm::getValue(fs.K(compIdx)), y / x, 1e-4);
    }
}

// A single phase next to the phase boundary returns ratios that start the
// stability test of a nearby state, which then reaches the answer of a flash
// from Wilson ratios.
BOOST_AUTO_TEST_CASE(SinglePhaseRatiosStartNearbyFlash)
{
    for (const auto& [p, T] : std::array<std::array<Scalar, 2>, 3>{
             {{38.47e5, 487.3}, {199.15e5, 487.3}, {359.83e5, 392.4}}}) {
        FluidState previous = makeState(p, T, 0.9);
        BOOST_REQUIRE(flash(previous));
        Scalar log_ratio_norm = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            const Scalar K = Opm::getValue(previous.K(compIdx));
            BOOST_CHECK(std::isfinite(K) && K > 0.0);
            log_ratio_norm += std::log(K) * std::log(K);
        }
        BOOST_CHECK(log_ratio_norm > 0.0);

        const Scalar next_p = 1.001 * p;
        const Scalar next_T = T + 0.05;
        FluidState warm = makeState(next_p, next_T, 0.9);
        seed(warm, previous);
        FluidState cold = makeState(next_p, next_T, 0.9);
        const bool warm_single = flash(warm);
        const bool cold_single = flash(cold);
        BOOST_CHECK_EQUAL(warm_single, cold_single);
        BOOST_CHECK_SMALL(Opm::getValue(warm.L()) - Opm::getValue(cold.L()), 1e-6);
    }
}

// Near the critical point a negative flash returns ratios within 1e-4 of one
// (here those of a neighbouring cell). Started from them, both stability
// trials would stop at the trivial solution and miss that this state splits;
// the test must start from the Wilson ratios instead.
BOOST_AUTO_TEST_CASE(TrivialRatiosDoNotHideInstability)
{
    constexpr Scalar p = 37322033.898305088;
    constexpr Scalar T = 362.71186440677968;
    FluidState cold = makeState(p, T, 0.9);
    BOOST_REQUIRE(!flash(cold));

    FluidState warm = makeState(p, T, 0.9);
    warm.setKvalue(0, 1.00001);
    warm.setKvalue(1, 0.999935);
    warm.setLvalue(1.0);
    BOOST_REQUIRE(!flash(warm));
    BOOST_CHECK_SMALL(Opm::getValue(warm.L()) - Opm::getValue(cold.L()), 1e-6);
}
