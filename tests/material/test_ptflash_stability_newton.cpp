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
 * \brief The second-order steps of the stability trials.
 *
 * The composition derivatives of the fugacity coefficients behind the
 * Hessian of the tangent-plane distance are checked against central
 * differences of CubicEOS, and a stability test against the stationary-point
 * conditions of the trial phase it returns.
 */
#include "config.h"

#define BOOST_TEST_MODULE PtFlashStabilityNewton
#include <boost/test/unit_test.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/eos/CubicEOS.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <array>
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

// The component parameters are static, so they are registered once per process.
struct Fixture
{
    Fixture()
    {
        FluidSystem::init();
        registerComponent<Opm::SimpleCO2<Scalar>>();
        registerComponent<Opm::C1<Scalar>>();
        registerComponent<Opm::C10<Scalar>>();
        // CO2-C1, CO2-C10, C1-C10
        FluidSystem::setInteractionCoefficients({0.1, 0.1, 0.0411});
    }
};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

namespace {

constexpr int numComponents = FluidSystem::numComponents;
using PtFlash = Opm::PTFlash<Scalar, FluidSystem>;
using FluidState = Opm::CompositionalFluidState<Scalar, FluidSystem>;
using ParamCache = FluidSystem::ParameterCache<Scalar>;
using CubicEOS = Opm::CubicEOS<Scalar, FluidSystem>;
using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
using Matrix = Dune::FieldMatrix<Scalar, numComponents, numComponents>;

class ExposedPtFlash : public PtFlash
{
public:
    using PtFlash::checkStability_;
    using PtFlash::feedStableRoot_;
    using PtFlash::fugacityCompositionDerivatives_;
};

FluidState makeState(const Scalar pressure,
                     const Scalar temperature,
                     const unsigned phase_idx,
                     const std::array<Scalar, numComponents>& x)
{
    FluidState fs;
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fs.setPressure(phaseIdx, pressure);
    }
    fs.setTemperature(temperature);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        fs.setMoleFraction(compIdx, x[compIdx]);
        fs.setMoleFraction(phase_idx, compIdx, x[compIdx]);
    }
    return fs;
}

std::array<Scalar, numComponents> logFugacityCoefficients(const Scalar pressure,
                                                          const Scalar temperature,
                                                          const unsigned phase_idx,
                                                          std::array<Scalar, numComponents> n)
{
    Scalar total = 0.;
    for (const Scalar n_i : n) {
        total += n_i;
    }
    for (Scalar& n_i : n) {
        n_i /= total;
    }
    const FluidState fs = makeState(pressure, temperature, phase_idx, n);
    ParamCache cache(EOSType::PR);
    cache.updatePhase(fs, phase_idx);
    std::array<Scalar, numComponents> phi;
    CubicEOS::computeFugacityCoefficients(fs, cache, phase_idx, phi);
    for (Scalar& value : phi) {
        value = std::log(value);
    }
    return phi;
}

// n d ln phi_i / d n_j by central differences of the mole numbers.
Matrix finiteDifferences(const Scalar pressure,
                         const Scalar temperature,
                         const unsigned phase_idx,
                         const std::array<Scalar, numComponents>& x)
{
    Matrix derivatives;
    for (int j = 0; j < numComponents; ++j) {
        const Scalar h = 1e-6 * x[j];
        auto plus = x;
        auto minus = x;
        plus[j] += h;
        minus[j] -= h;
        const auto up = logFugacityCoefficients(pressure, temperature, phase_idx, plus);
        const auto down = logFugacityCoefficients(pressure, temperature, phase_idx, minus);
        for (int i = 0; i < numComponents; ++i) {
            derivatives[i][j] = (up[i] - down[i]) / (2. * h);
        }
    }
    return derivatives;
}

} // anonymous namespace

// Liquid and vapour roots, with and without three real roots, near-critical states.
BOOST_AUTO_TEST_CASE(FugacityDerivativesMatchFiniteDifferences)
{
    struct Point
    {
        Scalar pressure;
        Scalar temperature;
        unsigned phase_idx;
        std::array<Scalar, numComponents> x;
    };
    const std::array<Point, 6> points{{
        {50e5, 300., FluidSystem::oilPhaseIdx, {0.1, 0.2, 0.7}},
        {50e5, 300., FluidSystem::gasPhaseIdx, {0.6, 0.39, 0.01}},
        {150e5, 400., FluidSystem::oilPhaseIdx, {0.4, 0.3, 0.3}},
        {150e5, 400., FluidSystem::gasPhaseIdx, {0.4, 0.3, 0.3}},
        {250e5, 260., FluidSystem::oilPhaseIdx, {0.9, 0.05, 0.05}},
        {20e5, 500., FluidSystem::gasPhaseIdx, {0.05, 0.9, 0.05}},
    }};
    for (const auto& point : points) {
        const FluidState fs = makeState(point.pressure, point.temperature, point.phase_idx, point.x);
        ParamCache cache(EOSType::PR);
        cache.updatePhase(fs, point.phase_idx);
        Matrix analytic;
        BOOST_REQUIRE(ExposedPtFlash::fugacityCompositionDerivatives_(
            cache, point.phase_idx, point.pressure, point.temperature, point.x, analytic));
        const Matrix numeric
            = finiteDifferences(point.pressure, point.temperature, point.phase_idx, point.x);
        Scalar scale = 0.;
        for (int i = 0; i < numComponents; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                scale = std::max(scale, std::abs(numeric[i][j]));
            }
        }
        for (int i = 0; i < numComponents; ++i) {
            Scalar gibbs_duhem = 0.;
            for (int j = 0; j < numComponents; ++j) {
                BOOST_CHECK_SMALL(analytic[i][j] - numeric[i][j], 1e-6 * scale);
                // A Maxwell relation, and the Gibbs-Duhem equation.
                BOOST_CHECK_SMALL(analytic[i][j] - analytic[j][i], 1e-10 * scale);
                gibbs_duhem += point.x[j] * analytic[j][i];
            }
            BOOST_CHECK_SMALL(gibbs_duhem, 1e-10 * scale);
        }
    }
}

// Every trial that does not collapse onto the feed ends on a stationary point
// of the tangent-plane distance: ln W_i + ln phi_i(w) = ln z_i + ln phi_i(z),
// with S = sum_i W_i the trial's mole number sum and w = W / S.
BOOST_AUTO_TEST_CASE(TrialsEndOnStationaryPoints)
{
    const std::array<Scalar, numComponents> z{0.4, 0.3, 0.3};
    int non_trivial = 0;
    for (const Scalar pressure : {20e5, 60e5, 100e5, 150e5, 200e5}) {
        for (const Scalar temperature : {260., 320., 400., 480.}) {
            FluidState fs = makeState(pressure, temperature, FluidSystem::oilPhaseIdx, z);
            ComponentVector z_vector;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                z_vector[compIdx] = z[compIdx];
                fs.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, z[compIdx]);
            }
            const int reference_phase_idx
                = ExposedPtFlash::feedStableRoot_(fs, z_vector, EOSType::PR).phase_idx;
            const FluidState feed
                = makeState(pressure, temperature, reference_phase_idx, z);
            ParamCache feed_cache(EOSType::PR);
            feed_cache.updatePhase(feed, reference_phase_idx);
            std::array<Scalar, numComponents> feed_phi;
            CubicEOS::computeFugacityCoefficients(feed, feed_cache, reference_phase_idx, feed_phi);

            for (const bool is_vapour_trial : {true, false}) {
                ComponentVector K;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    K[compIdx] = fs.wilsonK_(compIdx);
                }
                bool is_trivial = false;
                ComponentVector composition;
                Scalar sum = 0.;
                ExposedPtFlash::checkStability_(fs, is_trivial, K, composition, sum, z_vector,
                                                is_vapour_trial, reference_phase_idx,
                                                EOSType::PR, 0);
                if (is_trivial) {
                    continue;
                }
                ++non_trivial;
                const unsigned trial_phase_idx
                    = is_vapour_trial ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx;
                std::array<Scalar, numComponents> w;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    w[compIdx] = composition[compIdx];
                }
                const FluidState trial = makeState(pressure, temperature, trial_phase_idx, w);
                ParamCache trial_cache(EOSType::PR);
                trial_cache.updatePhase(trial, trial_phase_idx);
                std::array<Scalar, numComponents> trial_phi;
                CubicEOS::computeFugacityCoefficients(trial, trial_cache, trial_phase_idx, trial_phi);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar gradient = std::log(sum * w[compIdx] * trial_phi[compIdx])
                        - std::log(z[compIdx] * feed_phi[compIdx]);
                    BOOST_CHECK_SMALL(gradient, 1e-4);
                }
            }
        }
    }
    // The states include stable feeds with non-trivial stationary points and
    // unstable ones; most trials of a stable feed collapse.
    BOOST_CHECK(non_trivial > 0);
}
