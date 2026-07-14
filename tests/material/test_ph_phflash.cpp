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
 * \brief Acceptance tests for the isenthalpic (P-H) flash: the
 *        manufactured-enthalpy round-trip. Flash at a known temperature,
 *        read the mixture enthalpy, hand that enthalpy to PHFlash — the
 *        known temperature must be recovered. Fully synthetic: the target
 *        enthalpy is manufactured from the same model that is inverted, so
 *        every test carries its own ground truth.
 *
 * Circularity limit of the trick: a self-consistent error in the enthalpy
 * model would cancel in the round-trip. The model's own correctness is
 * established independently in test_ph_enthalpy.cpp (finite-difference,
 * decomposition and ideal-gas-limit checks); this file tests the INVERSION.
 */
#include "config.h"

#define BOOST_TEST_MODULE PhMvpPhFlash
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/IdealGasCaloricData.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PHFlash.hpp>

#include <opm/material/fluidsystems/ThreeComponentFluidSystem.hh>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include "ph_mvp_fixtures.hh"

#include <cmath>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using Opm::EnthalpyModel;
using Opm::PhMvpTest::FlashCase;
using Opm::PhMvpTest::makeInitialState;
using Opm::PhMvpTest::runFlash;
using Opm::PhMvpTest::f1Pressure;
using Opm::PhMvpTest::f1Temperature;
using Opm::PhMvpTest::f1Z;
using Opm::PhMvpTest::f2Pressure;
using Opm::PhMvpTest::f2Z;

// F1: binary C1/nC10
using FluidSystemF1 = Opm::PhMvpTest::TwoComponentFluidSystem<Scalar>;
constexpr int numComponentsF1 = FluidSystemF1::numComponents;
using EvaluationF1 = Opm::PhMvpTest::FlashEvaluation<FluidSystemF1>;
using EnthalpyF1 = Opm::MixtureEnthalpy<Scalar, FluidSystemF1>;
using PhFlashF1 = Opm::PHFlash<Scalar, FluidSystemF1>;

// F2: ternary CO2/C1/nC10
using FluidSystemF2 = Opm::ThreeComponentFluidSystem<Scalar>;
constexpr int numComponentsF2 = FluidSystemF2::numComponents;
using EvaluationF2 = Opm::PhMvpTest::FlashEvaluation<FluidSystemF2>;
using EnthalpyF2 = Opm::MixtureEnthalpy<Scalar, FluidSystemF2>;
using PhFlashF2 = Opm::PHFlash<Scalar, FluidSystemF2>;

namespace {

const Scalar T0 = Opm::IdealGasCaloricData<Scalar>::referenceTemperature();

constexpr double PT_TOLERANCE = 1.e-8;        // inner isothermal flash: fugacity-ratio residual
constexpr double ROUNDTRIP_TOLERANCE = 1.e-3; // [K] outer acceptance bound on the recovered temperature

PhFlashF1::Config makeConfigF1(const EnthalpyModel model)
{
    PhFlashF1::Config cfg;
    cfg.cpTable = Opm::PhMvpTest::f1CpTable();
    cfg.model = model;
    return cfg;
}

// manufacture the target enthalpy: flash F1 at (P, Tstar) and read the
// mixture enthalpy under the given model
double manufactureHspecF1(const double pressure, const double Tstar,
                          const EnthalpyModel model)
{
    FlashCase<numComponentsF1> probe{"Hspec probe", pressure, Tstar, f1Z};
    const auto outcome = runFlash<FluidSystemF1, EvaluationF1>(probe);
    return Opm::getValue(EnthalpyF1::mixtureEnthalpy(
        outcome.state, Opm::PhMvpTest::f1CpTable(), T0, EOSType::PR, model));
}

// run the P-H flash on a fresh F1 state and return the recovered temperature
double recoverTemperatureF1(const double pressure, const double hSpec,
                            const EnthalpyModel model)
{
    // initial temperature is irrelevant to the P-H solve; use the anchor value
    auto fs = makeInitialState<FluidSystemF1, EvaluationF1>(pressure, f1Temperature, f1Z);
    const bool ok = PhFlashF1::solve(fs, hSpec, makeConfigF1(model),
                                     "ssi", PT_TOLERANCE, EOSType::PR);
    BOOST_REQUIRE_MESSAGE(ok, "PHFlash reported hSpec out of range, hSpec = " << hSpec);
    return Opm::getValue(fs.temperature(0));
}

} // anonymous namespace

// The golden round-trip, caloric model: for a spread of known temperatures
// (two-phase window plus a single-phase liquid point), manufacture the
// enthalpy at T*, then recover T* from it.
BOOST_AUTO_TEST_CASE(RoundTripCaloricF1)
{
    for (const double Tstar : {200., 260., 300., 340., 380.}) {
        const double hSpec = manufactureHspecF1(f1Pressure, Tstar, EnthalpyModel::caloric);
        const double T = recoverTemperatureF1(f1Pressure, hSpec, EnthalpyModel::caloric);
        BOOST_CHECK_MESSAGE(std::abs(T - Tstar) < ROUNDTRIP_TOLERANCE,
                            "caloric round-trip: expected T* = " << Tstar
                                                                 << " K, recovered " << T << " K");
    }
}

// The coupled-system acceptance: the same round-trip under the EoS-departure
// model, where the enthalpy genuinely depends on the phase split. This is the
// test the caloric model cannot exercise (its split-independence makes the
// inner flash irrelevant to the residual).
BOOST_AUTO_TEST_CASE(RoundTripDepartureF1)
{
    // same point spread as the caloric case, including the 200 K single-phase
    // liquid point — the departure term is largest exactly there
    for (const double Tstar : {200., 260., 300., 340., 380.}) {
        const double hSpec = manufactureHspecF1(f1Pressure, Tstar, EnthalpyModel::eos_departure);
        const double T = recoverTemperatureF1(f1Pressure, hSpec, EnthalpyModel::eos_departure);
        BOOST_CHECK_MESSAGE(std::abs(T - Tstar) < ROUNDTRIP_TOLERANCE,
                            "departure round-trip: expected T* = " << Tstar
                                                                   << " K, recovered " << T << " K");
    }
}

// Ternary round-trip (F2) at its canonical two-phase point, both models
BOOST_AUTO_TEST_CASE(RoundTripTernaryF2)
{
    PhFlashF2::Config cfg;
    cfg.cpTable = Opm::PhMvpTest::f2CpTable(); // component order compile-time asserted
    // tighter bracket around the known solution — a precaution, not a
    // measured necessity: the bounds are full inner flashes, and this
    // ternary/pressure combination is only exercised near its canonical point
    cfg.tempMin = 250.;
    cfg.tempMax = 500.;

    for (const EnthalpyModel model : {EnthalpyModel::caloric, EnthalpyModel::eos_departure}) {
        cfg.model = model;
        constexpr double Tstar = 300.;

        FlashCase<numComponentsF2> probe{"F2 Hspec probe", f2Pressure, Tstar, f2Z};
        const auto outcome = runFlash<FluidSystemF2, EvaluationF2>(probe);
        const double hSpec = Opm::getValue(EnthalpyF2::mixtureEnthalpy(
            outcome.state, cfg.cpTable, T0, EOSType::PR, model));

        auto fs = makeInitialState<FluidSystemF2, EvaluationF2>(f2Pressure, 350., f2Z);
        const bool ok = PhFlashF2::solve(fs, hSpec, cfg, "ssi", PT_TOLERANCE, EOSType::PR);
        BOOST_REQUIRE(ok);
        BOOST_CHECK_MESSAGE(std::abs(Opm::getValue(fs.temperature(0)) - Tstar) < ROUNDTRIP_TOLERANCE,
                            "F2 round-trip failed for model " << static_cast<int>(model));
    }
}

// An unattainable specified enthalpy is rejected by the bracket sign check —
// returns false, never throws.
BOOST_AUTO_TEST_CASE(OutOfRangeReturnsFalse)
{
    auto cfg = makeConfigF1(EnthalpyModel::caloric);
    const double hAboveMax = manufactureHspecF1(f1Pressure, cfg.tempMax, EnthalpyModel::caloric)
                             + 1e6; // far above H(tempMax)

    auto fs = makeInitialState<FluidSystemF1, EvaluationF1>(f1Pressure, f1Temperature, f1Z);
    const bool ok = PhFlashF1::solve(fs, hAboveMax, cfg, "ssi", PT_TOLERANCE, EOSType::PR);
    BOOST_CHECK(!ok);
}

// Isenthalpic pressure drop (Joule-Thomson): manufacture H at (50 bar, 300 K),
// then solve at 10 bar with the same H. Under the caloric model the enthalpy
// is pressure-independent, so the temperature must be unchanged; under the
// departure model the expansion cools the mixture (attractive interactions,
// well below the inversion regime for this feed).
BOOST_AUTO_TEST_CASE(JouleThomsonSign)
{
    constexpr double pLow = 10e5; // [Pa]

    const double hCaloric = manufactureHspecF1(f1Pressure, f1Temperature, EnthalpyModel::caloric);
    const double tCaloric = recoverTemperatureF1(pLow, hCaloric, EnthalpyModel::caloric);
    BOOST_CHECK_SMALL(tCaloric - f1Temperature, 1e-5);

    const double hDeparture = manufactureHspecF1(f1Pressure, f1Temperature, EnthalpyModel::eos_departure);
    const double tDeparture = recoverTemperatureF1(pLow, hDeparture, EnthalpyModel::eos_departure);
    BOOST_CHECK_MESSAGE(tDeparture < f1Temperature,
                        "expected Joule-Thomson cooling, got T = " << tDeparture
                                                                   << " K at " << pLow << " Pa");
}
