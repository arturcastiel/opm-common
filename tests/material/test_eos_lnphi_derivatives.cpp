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
 * \brief The analytic derivatives of the logarithms of the cubic
 *        equation-of-state fugacity coefficients in the phase mole fractions:
 *        - the coefficients are those of fugacityCoefficients(), bit for bit;
 *        - the derivatives are those automatic differentiation gives, to
 *          1e-10, for every equation of state, on liquid and vapour roots;
 *        - the mole-number derivatives formed from them are symmetric and
 *          satisfy the Gibbs-Duhem relation.
 */
#include "config.h"

#define BOOST_TEST_MODULE EosLnPhiDerivatives
#include <boost/test/unit_test.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, false>;
constexpr int numComponents = FluidSystem::numComponents;
using Eval = Opm::DenseAd::Evaluation<Scalar, numComponents>;
using Composition = std::array<Scalar, numComponents>;
using Derivatives = std::array<std::array<Scalar, numComponents>, numComponents>;

namespace {

template <class Component>
void registerComponent(const Scalar volumeShift)
{
    using CompParam = typename FluidSystem::ComponentParam;
    FluidSystem::addComponent(CompParam{Component::name(),
                                        Component::molarMass(),
                                        Component::criticalTemperature(),
                                        Component::criticalPressure(),
                                        Component::criticalVolume(),
                                        Component::acentricFactor(),
                                        volumeShift});
}

// Nonzero volume shifts and interaction coefficients, so that every term of
// the expressions takes part.
struct Fixture
{
    Fixture()
    {
        FluidSystem::init();
        registerComponent<Opm::SimpleCO2<Scalar>>(-0.08);
        registerComponent<Opm::C1<Scalar>>(-0.15);
        registerComponent<Opm::C10<Scalar>>(0.09);
        FluidSystem::setInteractionCoefficients({0.1, 0.1, 0.0411});
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

bool sameBits(const double a, const double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

const std::vector<EOSType> eosTypes {EOSType::PR, EOSType::PRCORR, EOSType::SRK, EOSType::RK};

// Pressure [Pa], temperature [K] and one composition, put in both phases: the
// liquid takes the smallest root and the vapour the largest.
struct State
{
    Scalar pressure;
    Scalar temperature;
    Composition composition;
};

const std::vector<State> states {
    {50e5, 300.0, {0.5, 0.3, 0.2}},
    {100e5, 350.0, {0.1, 0.8, 0.1}},
    {150e5, 423.15, {0.2, 0.4, 0.4}},
    {20e5, 280.0, {0.05, 0.05, 0.9}},
    {300e5, 450.0, {0.6, 0.3, 0.1}},
    {80e5, 250.0, {0.9, 0.05, 0.05}},
    {5e5, 400.0, {0.01, 0.98, 0.01}},
};

template <class Value>
Opm::CompositionalFluidState<Value, FluidSystem> fluidState(const State& state)
{
    Opm::CompositionalFluidState<Value, FluidSystem> fs;
    for (const auto phaseIdx : {FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx}) {
        fs.setPressure(phaseIdx, Value(state.pressure));
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Value moleFraction = state.composition[compIdx];
            if constexpr (std::is_same_v<Value, Eval>) {
                moleFraction = Eval::createVariable(state.composition[compIdx], compIdx);
            }
            fs.setMoleFraction(phaseIdx, compIdx, moleFraction);
        }
    }
    fs.setTemperature(Value(state.temperature));
    return fs;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(DerivativesAreThoseOfAutomaticDifferentiation)
{
    Scalar largestDifference = 0.0;
    for (const auto eos : eosTypes) {
        for (const auto& state : states) {
            for (const unsigned phaseIdx : {FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx}) {
                BOOST_TEST_CONTEXT("EoS " << static_cast<int>(eos) << ", p = " << state.pressure
                                   << ", T = " << state.temperature << ", phase " << phaseIdx) {
                    auto fs = fluidState<Scalar>(state);
                    FluidSystem::ParameterCache<Scalar> cache(eos);
                    cache.updatePhase(fs, phaseIdx);
                    std::array<Scalar, numComponents> reference;
                    FluidSystem::fugacityCoefficients(fs, cache, phaseIdx, reference);
                    std::array<Scalar, numComponents> coefficients;
                    Derivatives derivatives;
                    FluidSystem::fugacityCoefficientsAndDerivatives(
                        fs, cache, phaseIdx, coefficients, derivatives);

                    auto adFs = fluidState<Eval>(state);
                    FluidSystem::ParameterCache<Eval> adCache(eos);
                    adCache.updatePhase(adFs, phaseIdx);
                    std::array<Eval, numComponents> adCoefficients;
                    FluidSystem::fugacityCoefficients(adFs, adCache, phaseIdx, adCoefficients);

                    for (int i = 0; i < numComponents; ++i) {
                        BOOST_CHECK(sameBits(coefficients[i], reference[i]));
                        const Eval lnPhi = log(adCoefficients[i]);
                        for (int j = 0; j < numComponents; ++j) {
                            const Scalar ad = lnPhi.derivative(j);
                            const Scalar difference = std::abs(derivatives[i][j] - ad)
                                / std::max(1.0, std::abs(ad));
                            largestDifference = std::max(largestDifference, difference);
                            BOOST_CHECK_SMALL(difference, 1e-10);
                        }
                    }
                }
            }
        }
    }
    BOOST_TEST_MESSAGE("largest relative difference to automatic differentiation: "
                       << largestDifference);
}

BOOST_AUTO_TEST_CASE(MoleNumberDerivativesAreSymmetric)
{
    for (const auto eos : eosTypes) {
        for (const auto& state : states) {
            for (const unsigned phaseIdx : {FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx}) {
                auto fs = fluidState<Scalar>(state);
                FluidSystem::ParameterCache<Scalar> cache(eos);
                cache.updatePhase(fs, phaseIdx);
                std::array<Scalar, numComponents> coefficients;
                Derivatives derivatives;
                FluidSystem::fugacityCoefficientsAndDerivatives(
                    fs, cache, phaseIdx, coefficients, derivatives);

                // n d ln phi_i / d n_j
                Derivatives moleNumber;
                for (int i = 0; i < numComponents; ++i) {
                    Scalar weighted = 0.0;
                    for (int k = 0; k < numComponents; ++k) {
                        weighted += state.composition[k] * derivatives[i][k];
                    }
                    for (int j = 0; j < numComponents; ++j) {
                        moleNumber[i][j] = derivatives[i][j] - weighted;
                    }
                }
                for (int i = 0; i < numComponents; ++i) {
                    Scalar gibbsDuhem = 0.0;
                    for (int j = 0; j < numComponents; ++j) {
                        BOOST_CHECK_SMALL(moleNumber[i][j] - moleNumber[j][i],
                                          1e-10 * std::max(1.0, std::abs(moleNumber[i][j])));
                        gibbsDuhem += state.composition[j] * moleNumber[j][i];
                    }
                    BOOST_CHECK_SMALL(gibbsDuhem, 1e-10);
                }
            }
        }
    }
}
