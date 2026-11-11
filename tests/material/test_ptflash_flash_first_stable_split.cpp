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
 * \brief A cold flash that starts with the split (flash first) must end on
 *        the stable split, not on another stationary point below the feed.
 *
 * Hydrogen, propane and water at 244 bar and 327.6 K split into a water-rich
 * liquid (L = 0.0481) and a hydrogen-propane gas. Substitution from Wilson's
 * ratios instead reaches a hydrogen-propane liquid and vapour (L = 0.588).
 * That split also lies below the Gibbs energy of the feed, which proves the
 * feed unstable, but its own liquid is unstable: the water-rich split has a
 * Gibbs energy lower by 0.052 RT. The stability test of the feed finds the
 * water-rich split.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashFlashFirstStableSplit
#include <boost/test/unit_test.hpp>

#include <opm/material/components/H2.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <cmath>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, false>;

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

struct Fixture
{
    Fixture()
    {
        FluidSystem::init();
        registerComponent<Opm::H2<Scalar>>();
        // propane
        FluidSystem::addComponent(FluidSystem::ComponentParam{"C3", 44.10e-3, 369.83, 4.2493e6, 0.2005, 0.1524});
        registerComponent<Opm::H2O<Scalar>>();
        // (C3, H2), (H2O, H2), (H2O, C3)
        FluidSystem::setInteractionCoefficients({0.0, 0.13105793120443596, 0.0});
    }
};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(ColdFlashEndsOnTheStableSplit)
{
    constexpr int numComponents = FluidSystem::numComponents;
    using Evaluation = Opm::DenseAd::Evaluation<Scalar, numComponents + 1>;
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;
    using PtFlash = Opm::PTFlash<Scalar, FluidSystem>;

    const Scalar p = 24382609.35694183;
    const Scalar T = 327.6331638362285;
    const Scalar z[numComponents] = {0.48332520820349184, 0.46196719448606804, 0.054707597310439976};

    for (const auto method : {Opm::PTFlashMethod::SsiNewton, Opm::PTFlashMethod::Ssi, Opm::PTFlashMethod::Newton}) {
        BOOST_TEST_CONTEXT("method " << Opm::ptFlashMethodToString(method))
        {
            FluidState fs;
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                fs.setPressure(phaseIdx, Evaluation::createVariable(p, 0));
            }
            fs.setTemperature(T);
            Evaluation last = 1.0;
            for (int compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                const auto zc = Evaluation::createVariable(z[compIdx], compIdx + 1);
                fs.setMoleFraction(compIdx, zc);
                last -= zc;
            }
            fs.setMoleFraction(numComponents - 1, last);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
            }
            fs.setLvalue(-1.0);

            const bool single = PtFlash::solve(fs, method, 1e-8, EOSType::PR, 0);

            BOOST_CHECK(!single);
            // The water-rich liquid; the other stationary split has L = 0.588.
            BOOST_CHECK_SMALL(Opm::getValue(fs.L()) - 0.0481294890, 1e-6);
            BOOST_CHECK_GT(Opm::getValue(fs.moleFraction(FluidSystem::oilPhaseIdx, 2)), 0.999);
            BOOST_CHECK_SMALL(Opm::getValue(fs.moleFraction(FluidSystem::gasPhaseIdx, 2)) - 0.00691607, 1e-6);
        }
    }
}
