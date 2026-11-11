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
 * \brief The shortcuts of the cubic equation-of-state evaluation give the
 *        values of the full evaluation bit for bit:
 *        - CubicEOSParams::updatePure() at an unchanged temperature and
 *          pressure keeps parameters equal to a fresh evaluation;
 *        - the fugacity coefficients of a phase evaluated together equal
 *          those evaluated one component at a time, volume shift included;
 *        - cubicExtremeRoot() returns the root cubicRoots() returns at the
 *          position of the smallest or largest root, also next to double
 *          and triple roots.
 */
#include "config.h"

#define BOOST_TEST_MODULE EosEvaluationReuse
#include <boost/test/unit_test.hpp>

#include <opm/material/common/PolynomialUtils.hpp>
#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/eos/CubicEOS.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

using Scalar = double;
using EOSType = Opm::CompositionalConfig::EOSType;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, false>;
constexpr int numComponents = FluidSystem::numComponents;
using Eval = Opm::DenseAd::Evaluation<Scalar, numComponents + 1>;

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

// Nonzero volume shifts, so the translation factors take part.
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

bool sameBits(const Eval& a, const Eval& b)
{
    if (!sameBits(a.value(), b.value())) {
        return false;
    }
    for (int i = 0; i < a.size(); ++i) {
        if (!sameBits(a.derivative(i), b.derivative(i))) {
            return false;
        }
    }
    return true;
}

const std::vector<EOSType> eosTypes {EOSType::PR, EOSType::PRCORR, EOSType::SRK, EOSType::RK};

// A fluid state at (p, T) with liquid x and vapour y; pressure and the first
// mole fractions carry derivatives, as the primary variables of a simulator.
template <class Value>
Opm::CompositionalFluidState<Value, FluidSystem>
fluidState(const Scalar p, const Scalar T, const std::array<Scalar, numComponents>& x,
           const std::array<Scalar, numComponents>& y)
{
    Opm::CompositionalFluidState<Value, FluidSystem> fs;
    Value pv = p;
    if constexpr (std::is_same_v<Value, Eval>) {
        pv = Eval::createVariable(p, 0);
    }
    fs.setPressure(FluidSystem::oilPhaseIdx, pv);
    fs.setPressure(FluidSystem::gasPhaseIdx, pv);
    fs.setTemperature(Value(T));
    for (int c = 0; c < numComponents; ++c) {
        Value xc = x[c];
        Value yc = y[c];
        if constexpr (std::is_same_v<Value, Eval>) {
            if (c + 1 < numComponents) {
                xc = Eval::createVariable(x[c], c + 1);
                yc = Eval::createVariable(y[c], c + 1);
            }
        }
        fs.setMoleFraction(FluidSystem::oilPhaseIdx, c, xc);
        fs.setMoleFraction(FluidSystem::gasPhaseIdx, c, yc);
    }
    return fs;
}

template <class Value>
void checkPureParametersKept()
{
    using Params = Opm::CubicEOSParams<Value, FluidSystem, FluidSystem::oilPhaseIdx>;
    const std::array<std::pair<Scalar, Scalar>, 4> states {
        {{350.0, 8e6}, {350.0, 8e6}, {420.0, 8e6}, {420.0, 2e7}}};
    for (const auto eos : eosTypes) {
        Params reused;
        reused.setEOSType(eos);
        for (const auto& [T, p] : states) {
            Value Tv = T;
            Value pv = p;
            if constexpr (std::is_same_v<Value, Eval>) {
                pv = Eval::createVariable(p, 0);
            }
            reused.updatePure(Tv, pv);
            Params fresh;
            fresh.setEOSType(eos);
            fresh.updatePure(Tv, pv);
            for (unsigned i = 0; i < numComponents; ++i) {
                BOOST_CHECK(sameBits(reused.Ai(i), fresh.Ai(i)));
                BOOST_CHECK(sameBits(reused.Bi(i), fresh.Bi(i)));
                BOOST_CHECK(sameBits(reused.volumeTranslationFactor(i),
                                     fresh.volumeTranslationFactor(i)));
                for (unsigned j = 0; j < numComponents; ++j) {
                    BOOST_CHECK(sameBits(reused.aCache(i, j), fresh.aCache(i, j)));
                }
            }
        }
        // The other phase takes them only at the same (T, p) and equation.
        Opm::CubicEOSParams<Value, FluidSystem, FluidSystem::gasPhaseIdx> gas;
        gas.setEOSType(eos);
        Value lastPressure = 2e7;
        if constexpr (std::is_same_v<Value, Eval>) {
            lastPressure = Eval::createVariable(2e7, 0);
        }
        BOOST_CHECK(!gas.takePure(reused, Value(420.0), Value(8e6)));
        BOOST_CHECK(gas.takePure(reused, Value(420.0), lastPressure));
        for (unsigned i = 0; i < numComponents; ++i) {
            BOOST_CHECK(sameBits(gas.Ai(i), reused.Ai(i)));
            BOOST_CHECK(sameBits(gas.volumeTranslationFactor(i), reused.volumeTranslationFactor(i)));
            for (unsigned j = 0; j < numComponents; ++j) {
                BOOST_CHECK(sameBits(gas.aCache(i, j), reused.aCache(i, j)));
            }
        }
        // A change of the equation of state recomputes at the same (T, p).
        const auto other = eos == EOSType::PR ? EOSType::SRK : EOSType::PR;
        reused.setEOSType(other);
        reused.updatePure(Value(420.0), Value(2e7));
        Params fresh;
        fresh.setEOSType(other);
        fresh.updatePure(Value(420.0), Value(2e7));
        for (unsigned i = 0; i < numComponents; ++i) {
            BOOST_CHECK(sameBits(reused.Ai(i), fresh.Ai(i)));
        }
    }
}

template <class Value>
void checkFugacityCoefficientsTogether()
{
    using FluidState = Opm::CompositionalFluidState<Value, FluidSystem>;
    using ParamCache = FluidSystem::ParameterCache<Value>;
    using CubicEOS = Opm::CubicEOS<Scalar, FluidSystem>;
    const std::array<Scalar, numComponents> x {0.30, 0.20, 0.50};
    const std::array<Scalar, numComponents> y {0.45, 0.50, 0.05};
    for (const auto eos : eosTypes) {
        for (const auto& [p, T] : std::array<std::pair<Scalar, Scalar>, 3> {
                 {{5e6, 300.0}, {1.5e7, 400.0}, {3e7, 550.0}}}) {
            const FluidState fs = fluidState<Value>(p, T, x, y);
            ParamCache cache(eos);
            for (unsigned phaseIdx = 0; phaseIdx < 2; ++phaseIdx) {
                cache.updatePhase(fs, phaseIdx);
                std::array<Value, numComponents> together;
                FluidSystem::fugacityCoefficients(fs, cache, phaseIdx, together);
                std::array<Value, numComponents> eosTogether;
                CubicEOS::computeFugacityCoefficients(fs, cache, phaseIdx, eosTogether);
                for (unsigned c = 0; c < numComponents; ++c) {
                    const Value alone = FluidSystem::fugacityCoefficient(fs, cache, phaseIdx, c);
                    BOOST_CHECK(sameBits(together[c], alone));
                    const Value eosAlone
                        = CubicEOS::computeFugacityCoefficient(fs, cache, phaseIdx, c);
                    BOOST_CHECK(sameBits(eosTogether[c], eosAlone));
                    // The volume shifts are nonzero: the translation applies.
                    BOOST_CHECK(!sameBits(Opm::getValue(together[c]), Opm::getValue(eosAlone)));
                }
            }
        }
    }
}

template <class Value>
void checkExtremeRoot(const Value& a, const Value& b, const Value& c, const Value& d)
{
    Value all[3] = {0.0, 0.0, 0.0};
    unsigned numAll = 0;
    try {
        numAll = Opm::cubicRoots(all, a, b, c, d);
    }
    catch (const std::runtime_error&) {
        // A degenerate cubic cubicRoots() rejects is rejected alike.
        for (const bool largest : {false, true}) {
            Value one[3] = {0.0, 0.0, 0.0};
            BOOST_CHECK_THROW(Opm::cubicExtremeRoot(one, a, b, c, d, largest), std::runtime_error);
        }
        return;
    }
    for (const bool largest : {false, true}) {
        Value one[3] = {0.0, 0.0, 0.0};
        const unsigned numOne = Opm::cubicExtremeRoot(one, a, b, c, d, largest);
        BOOST_REQUIRE_EQUAL(numOne, numAll);
        const unsigned idx = largest ? numAll - 1 : 0;
        BOOST_CHECK(sameBits(one[idx], all[idx]));
    }
}

// Cubic coefficients from roots r0 <= r1 <= r2: (Z - r0)(Z - r1)(Z - r2).
template <class Value>
void checkFromRoots(const Value& r0, const Value& r1, const Value& r2)
{
    checkExtremeRoot<Value>(Value(1.0), -(r0 + r1 + r2), r0 * r1 + r0 * r2 + r1 * r2, -r0 * r1 * r2);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(PureParametersKeptAtSameTemperatureAndPressure)
{
    checkPureParametersKept<Scalar>();
    checkPureParametersKept<Eval>();
}

BOOST_AUTO_TEST_CASE(FugacityCoefficientsTogetherEqualOneAtATime)
{
    checkFugacityCoefficientsTogether<Scalar>();
    checkFugacityCoefficientsTogether<Eval>();
}

BOOST_AUTO_TEST_CASE(ExtremeRootEqualsSortedRoot)
{
    std::mt19937 gen(5);
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);
    // The cubic of the equation of state, m1 + m2 and m1 m2 of PR and SRK.
    for (const auto& [m_sum, m_prod] : std::array<std::pair<Scalar, Scalar>, 2> {{{2.0, -1.0}, {1.0, 0.0}}}) {
        for (int k = 0; k < 20000; ++k) {
            const Scalar A = 3.0 * unit(gen);
            const Scalar B = 0.3 * unit(gen);
            const Scalar a2 = (m_sum - 1) * B - 1;
            const Scalar a3 = A + m_prod * B * B - m_sum * B * (B + 1);
            const Scalar a4 = -A * B - m_prod * B * B * (B + 1);
            checkExtremeRoot<Scalar>(1.0, a2, a3, a4);
            const Eval Ae = Eval::createVariable(A, 0);
            const Eval Be = Eval::createVariable(B, 1);
            checkExtremeRoot<Eval>(Eval(1.0), (m_sum - 1) * Be - 1,
                                   Ae + m_prod * Be * Be - m_sum * Be * (Be + 1),
                                   -Ae * Be - m_prod * Be * Be * (Be + 1));
        }
    }
    // Next to double and triple roots, where rounding may order roots
    // differently than their exact values.
    for (const Scalar gap : {0.0, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4}) {
        for (const Scalar r : {0.003, 0.3, 1.0}) {
            checkFromRoots<Scalar>(r, r + gap, r + 0.5);
            checkFromRoots<Scalar>(r, r + 0.5, r + 0.5 + gap);
            checkFromRoots<Scalar>(r, r + gap, r + 2 * gap);
            checkFromRoots<Eval>(Eval::createVariable(r, 0), Eval(r + gap), Eval(r + 0.5));
            checkFromRoots<Eval>(Eval(r), Eval::createVariable(r + 0.5, 1), Eval(r + 0.5 + gap));
            checkFromRoots<Eval>(Eval::createVariable(r, 2), Eval(r + gap), Eval(r + 2 * gap));
        }
    }
}
