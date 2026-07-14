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
 * \brief Tests for the caloric mixture-enthalpy model (MvpEnthalpy/MvpCpData)
 *        — the property the isenthalpic (P-H) flash inverts for temperature.
 *
 * The cases live in the CaloricModel suite so that the EoS-consistent
 * departure extension can land as a sibling suite in this file.
 */
#include "config.h"

#define BOOST_TEST_MODULE PhMvpEnthalpy
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/MvpCpData.hpp>
#include <opm/material/constraintsolvers/MvpEnthalpy.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include "ph_mvp_fixtures.hh"

#include <array>
#include <cmath>

using Scalar = double;
using Opm::PhMvpTest::FlashCase;
using Opm::PhMvpTest::runFlash;
using Opm::PhMvpTest::f1Pressure;
using Opm::PhMvpTest::f1Z;

// F1: binary C1/nC10
using FluidSystemF1 = Opm::PhMvpTest::TwoComponentFluidSystem<Scalar>;
constexpr int numComponentsF1 = FluidSystemF1::numComponents;
using EvaluationF1 = Opm::PhMvpTest::FlashEvaluation<FluidSystemF1>;
using EnthalpyF1 = Opm::MvpEnthalpy<Scalar, FluidSystemF1>;

namespace {

const Scalar T0 = Opm::MvpCpData<Scalar>::referenceTemperature();

// unchecked probe helper: flash F1 at (P, T) and return the mixture enthalpy
// of the flashed state. Deliberately assertion-free — the calling test owns
// its expectations; do not bolt checks in here.
double mixtureEnthalpyAt(const double pressure, const double temperature)
{
    FlashCase<numComponentsF1> testCase{"enthalpy probe", pressure, temperature, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    return Opm::getValue(EnthalpyF1::mixtureEnthalpy(outcome.state, Opm::PhMvpTest::f1CpTable(), T0));
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(CaloricModel)

// Enthalpy is strictly increasing in temperature (Cp > 0). The sweep
// deliberately crosses phase-regime boundaries: with the caloric model H is
// monotone irrespective of how the split changes along the way.
BOOST_AUTO_TEST_CASE(MonotoneInTemperature)
{
    const std::array<double, 4> temperatures = {250., 300., 350., 400.};

    double previous = mixtureEnthalpyAt(f1Pressure, temperatures[0]);
    for (std::size_t i = 1; i < temperatures.size(); ++i) {
        const double current = mixtureEnthalpyAt(f1Pressure, temperatures[i]);
        BOOST_CHECK_MESSAGE(current > previous,
                            "H not monotone: H(" << temperatures[i] << ") = " << current
                                                 << " <= H(" << temperatures[i-1] << ") = " << previous);
        previous = current;
    }
}

// The reference datum: H(T0) = 0, regardless of the phase split
BOOST_AUTO_TEST_CASE(ReferenceDatum)
{
    const double H0 = mixtureEnthalpyAt(f1Pressure, T0);
    BOOST_CHECK_SMALL(H0, 1e-12);
}

// Analytic mixtureCp matches a central finite difference of the mixture
// enthalpy. Note the physics of WHY they may be compared at all: the FD probe
// re-flashes at T±h (a total derivative, split re-equilibrated) while
// mixtureCp is a frozen-split partial derivative — they agree only because
// the caloric H is flash-independent. Do NOT reuse this check unchanged for a
// departure-mode enthalpy, where the two derivatives legitimately differ.
BOOST_AUTO_TEST_CASE(CpVersusFiniteDifference)
{
    constexpr double T = 300.;
    constexpr double h = 1e-3; // [K]

    FlashCase<numComponentsF1> testCase{"cp probe", f1Pressure, T, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    const double cpAnalytic = Opm::getValue(
        EnthalpyF1::mixtureCp(outcome.state, Opm::PhMvpTest::f1CpTable()));

    const double cpFD = (mixtureEnthalpyAt(f1Pressure, T + h)
                         - mixtureEnthalpyAt(f1Pressure, T - h)) / (2.*h);

    BOOST_CHECK_CLOSE(cpAnalytic, cpFD, 1e-6); // [%]
}

// Phase-decomposed mixture enthalpy equals the feed-direct sum
// sum_i z_i*h_i(T). With the caloric model the phase split cancels exactly by
// material balance — this pins the L/x/y weighting (indexing) of the
// implementation and documents that caloric H is flash-independent.
BOOST_AUTO_TEST_CASE(PhaseDecompositionMatchesFeedSum)
{
    constexpr double T = 340.;

    FlashCase<numComponentsF1> testCase{"decomposition probe", f1Pressure, T, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    // the cancellation is only non-trivially tested on a genuine two-phase split
    BOOST_REQUIRE(!outcome.summary.single_phase);
    BOOST_REQUIRE(outcome.summary.L > 0. && outcome.summary.L < 1.);

    const auto cpTable = Opm::PhMvpTest::f1CpTable();
    const double viaPhases = Opm::getValue(EnthalpyF1::mixtureEnthalpy(outcome.state, cpTable, T0));

    double viaFeed = 0.;
    for (int compIdx = 0; compIdx < numComponentsF1; ++compIdx)
        viaFeed += f1Z[compIdx] * cpTable[compIdx].enthalpyIntegral(T, T0);

    BOOST_CHECK_CLOSE(viaPhases, viaFeed, 1e-9); // [%]
}

BOOST_AUTO_TEST_SUITE_END() // CaloricModel

// ────────────────────────────────────────────────────────────────────────────
// EoS-consistent departure (residual) enthalpy:
// H_res = -R*T^2 * sum_i w_i * dln(phi_i)/dT per phase, with the temperature
// derivative of the fugacity coefficient supplied by densead AD through the
// fluid system's own fugacityCoefficient. This is the model under which the
// enthalpy genuinely depends on the flash result (the caloric cancellation
// tested above no longer holds).
// ────────────────────────────────────────────────────────────────────────────
BOOST_AUTO_TEST_SUITE(DepartureModel)

namespace {
using EOSType = Opm::CompositionalConfig::EOSType;
}

// The AD temperature-derivative of ln(phi) matches a central finite
// difference of ln(phi(T)) evaluated at fixed pressure and frozen phase
// composition — per phase, per component.
BOOST_AUTO_TEST_CASE(AdDerivativeMatchesFiniteDifference)
{
    FlashCase<numComponentsF1> testCase{"AD probe", f1Pressure,
                                        Opm::PhMvpTest::f1Temperature, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    BOOST_REQUIRE(!outcome.summary.single_phase);

    constexpr double h = 1e-3; // FD step [K]
    for (unsigned phaseIdx : {static_cast<unsigned>(FluidSystemF1::oilPhaseIdx),
                              static_cast<unsigned>(FluidSystemF1::gasPhaseIdx)}) {
        // freeze this phase's composition from the flashed state
        std::array<double, numComponentsF1> w;
        for (int compIdx = 0; compIdx < numComponentsF1; ++compIdx)
            w[compIdx] = Opm::getValue(outcome.state.moleFraction(phaseIdx, compIdx));

        // independent, scalar evaluation path for ln(phi(T)) at fixed (P, w)
        auto lnPhi = [&](const double T, const int compIdx) {
            Opm::CompositionalFluidState<double, FluidSystemF1> fs;
            fs.setTemperature(T);
            fs.setPressure(FluidSystemF1::oilPhaseIdx, f1Pressure);
            fs.setPressure(FluidSystemF1::gasPhaseIdx, f1Pressure);
            for (int i = 0; i < numComponentsF1; ++i)
                fs.setMoleFraction(phaseIdx, i, w[i]);
            FluidSystemF1::ParameterCache<double> paramCache(EOSType::PR);
            paramCache.updatePhase(fs, phaseIdx);
            return std::log(FluidSystemF1::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx));
        };

        const double T = Opm::PhMvpTest::f1Temperature;
        for (int compIdx = 0; compIdx < numComponentsF1; ++compIdx) {
            const double ad = EnthalpyF1::phaseDLnPhiDT(outcome.state, phaseIdx, compIdx, EOSType::PR);
            const double fd = (lnPhi(T + h, compIdx) - lnPhi(T - h, compIdx)) / (2.*h);
            BOOST_CHECK_CLOSE(ad, fd, 1e-3); // [%]
        }
    }
}

// Ideal-gas limit: the gas-phase residual vanishes as P -> 0
BOOST_AUTO_TEST_CASE(IdealGasLimit)
{
    FlashCase<numComponentsF1> lowP{"near-vacuum probe", 1e3, Opm::PhMvpTest::f1Temperature, f1Z};
    FlashCase<numComponentsF1> anchor{"anchor probe", f1Pressure, Opm::PhMvpTest::f1Temperature, f1Z};

    const auto lowOutcome = runFlash<FluidSystemF1, EvaluationF1>(lowP);
    const auto anchorOutcome = runFlash<FluidSystemF1, EvaluationF1>(anchor);

    const double resLow = EnthalpyF1::phaseResidualEnthalpy(
        lowOutcome.state, FluidSystemF1::gasPhaseIdx, EOSType::PR);
    const double resAnchor = EnthalpyF1::phaseResidualEnthalpy(
        anchorOutcome.state, FluidSystemF1::gasPhaseIdx, EOSType::PR);

    BOOST_CHECK_LT(std::abs(resLow), 10.);   // [J/mol] — near-ideal at 10 mbar
    BOOST_CHECK_LT(std::abs(resLow), 0.05 * std::abs(resAnchor));
}

// With the departure term the enthalpy couples to the split: the L-weighted
// phase decomposition still holds, but the caloric feed-sum identity breaks,
// and the liquid residual is attractive (vaporization-enthalpy scale).
BOOST_AUTO_TEST_CASE(DepartureCouplesToSplit)
{
    FlashCase<numComponentsF1> testCase{"departure anchor", f1Pressure,
                                        Opm::PhMvpTest::f1Temperature, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    BOOST_REQUIRE(!outcome.summary.single_phase);
    BOOST_REQUIRE(outcome.summary.L > 0. && outcome.summary.L < 1.);

    const auto cpTable = Opm::PhMvpTest::f1CpTable();

    // (a) decomposition consistency of the model-switch overload
    const double viaMixture = Opm::getValue(EnthalpyF1::mixtureEnthalpy(
        outcome.state, cpTable, T0, EOSType::PR, Opm::EnthalpyModel::eos_departure));
    const double L = outcome.summary.L;
    const double hOil = Opm::getValue(EnthalpyF1::phaseEnthalpy(
                            outcome.state, FluidSystemF1::oilPhaseIdx, cpTable, T0))
        + EnthalpyF1::phaseResidualEnthalpy(outcome.state, FluidSystemF1::oilPhaseIdx, EOSType::PR);
    const double hGas = Opm::getValue(EnthalpyF1::phaseEnthalpy(
                            outcome.state, FluidSystemF1::gasPhaseIdx, cpTable, T0))
        + EnthalpyF1::phaseResidualEnthalpy(outcome.state, FluidSystemF1::gasPhaseIdx, EOSType::PR);
    BOOST_CHECK_CLOSE(viaMixture, L*hOil + (1. - L)*hGas, 1e-9); // [%]

    // (b) the caloric cancellation is broken: departure H differs from the
    // feed-direct ideal sum by far more than any tolerance
    double viaFeedIdeal = 0.;
    for (int compIdx = 0; compIdx < numComponentsF1; ++compIdx)
        viaFeedIdeal += f1Z[compIdx] * cpTable[compIdx].enthalpyIntegral(
            Opm::PhMvpTest::f1Temperature, T0);
    BOOST_CHECK_GT(std::abs(viaMixture - viaFeedIdeal), 100.); // [J/mol]

    // (c) the liquid's residual is negative (attractive interactions)
    BOOST_CHECK_LT(EnthalpyF1::phaseResidualEnthalpy(
        outcome.state, FluidSystemF1::oilPhaseIdx, EOSType::PR), 0.);
}

// The caloric seam is untouched: the model-switch overload with
// EnthalpyModel::caloric reproduces the original caloric overload exactly.
BOOST_AUTO_TEST_CASE(CaloricSeamUnchanged)
{
    FlashCase<numComponentsF1> testCase{"seam probe", f1Pressure,
                                        Opm::PhMvpTest::f1Temperature, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(testCase);
    const auto cpTable = Opm::PhMvpTest::f1CpTable();

    const double viaCaloric = Opm::getValue(
        EnthalpyF1::mixtureEnthalpy(outcome.state, cpTable, T0));
    const double viaSwitch = Opm::getValue(EnthalpyF1::mixtureEnthalpy(
        outcome.state, cpTable, T0, EOSType::PR, Opm::EnthalpyModel::caloric));
    BOOST_CHECK_CLOSE(viaCaloric, viaSwitch, 1e-12); // [%]
}

BOOST_AUTO_TEST_SUITE_END() // DepartureModel
