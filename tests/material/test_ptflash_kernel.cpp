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
 * \brief The P-T flash kernel (PTFlashKernel.hpp) against
 *        PTFlash::flash_solve_scalar() on the same states: the same failures,
 *        phase states, liquid fractions and phase compositions.
 *
 * Fluids: the ternary CO2/C1/nC10 (F2) with every two-phase method, cold
 * and warm starts; C1/nC20 with volume shifts (H1); and the SPE5 oil of the
 * dense-map benchmark, at states that take the Newton fallback in
 * (x, y, L) and at the one state of the 1000 x 1000 map where both fail.
 */
#include "config.h"

#define BOOST_TEST_MODULE PTFlashKernel
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/constraintsolvers/PTFlashKernel.hpp>
#include <opm/material/constraintsolvers/PTFlashKernelHost.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include "flashTestFixtures.hpp"

#include <array>
#include <cmath>
#include <exception>
#include <vector>

namespace {

using EOSType = Opm::CompositionalConfig::EOSType;
namespace K = Opm::PTFlashKernel;

// Both evaluate the same operations in the same order; the composition
// Newton's Jacobian (analytic in the kernel, automatic in PTFlash) is the
// only source of rounding differences.
constexpr double tolerance = 1e-10;

template <class FluidSystem>
struct Answer
{
    bool ok = false;
    bool single = false;
    double L = 0;
    std::array<double, FluidSystem::numComponents> x{};
    std::array<double, FluidSystem::numComponents> y{};
    std::array<double, FluidSystem::numComponents> K{};
};

//! PTFlash::flash_solve_scalar() from (K, L).
template <class FluidSystem>
Answer<FluidSystem> ptflash(double p, double T,
                            const std::array<double, FluidSystem::numComponents>& z,
                            const std::array<double, FluidSystem::numComponents>& K0, double L0,
                            Opm::PTFlashMethod method, EOSType eos)
{
    constexpr int N = FluidSystem::numComponents;
    Opm::CompositionalFluidState<double, FluidSystem> fs;
    for (unsigned ph = 0; ph < FluidSystem::numPhases; ++ph) {
        fs.setPressure(ph, p);
    }
    fs.setTemperature(T);
    for (int c = 0; c < N; ++c) {
        fs.setMoleFraction(c, z[c]);
        fs.setKvalue(c, K0[c]);
    }
    fs.setLvalue(L0);
    Answer<FluidSystem> a;
    try {
        a.single = Opm::PTFlash<double, FluidSystem>::flash_solve_scalar(fs, method, 1e-8, eos, 0);
        a.ok = std::isfinite(fs.L());
    }
    catch (const std::exception&) {
        a.ok = false;
    }
    a.L = fs.L();
    for (int c = 0; c < N; ++c) {
        a.x[c] = fs.moleFraction(FluidSystem::oilPhaseIdx, c);
        a.y[c] = fs.moleFraction(FluidSystem::gasPhaseIdx, c);
        a.K[c] = a.single ? K0[c] : a.y[c] / a.x[c];
    }
    return a;
}

//! The kernel on the same input; checks it against PTFlash and returns it.
template <class FluidSystem>
Answer<FluidSystem> checkState(double p, double T,
                               const std::array<double, FluidSystem::numComponents>& z,
                               const std::array<double, FluidSystem::numComponents>& K0, double L0,
                               Opm::PTFlashMethod method, EOSType eos)
{
    constexpr int N = FluidSystem::numComponents;
    const auto reference = ptflash<FluidSystem>(p, T, z, K0, L0, method, eos);

    const auto params = K::makeParameters<double, FluidSystem>(eos);
    K::Vec<double, N> zk, k0;
    for (int c = 0; c < N; ++c) {
        zk[c] = z[c];
        k0[c] = K0[c];
    }
    K::Settings<double> settings;
    settings.method = K::methodFrom(method);
    settings.tolerance = 1e-8;
    const auto kernel = K::flash(params, p, T, zk, k0, L0, settings);

    const bool kernel_ok = kernel.status == K::Status::ok && std::isfinite(kernel.L);
    BOOST_TEST_CONTEXT("p = " << p << " Pa, T = " << T << " K, L0 = " << L0 << ", method "
                              << static_cast<int>(method) << ", kernel status "
                              << static_cast<int>(kernel.status))
    {
        BOOST_CHECK_EQUAL(kernel_ok, reference.ok);
        if (kernel_ok && reference.ok) {
            BOOST_CHECK_EQUAL(kernel.singlePhase, reference.single);
            BOOST_CHECK_SMALL(kernel.L - reference.L, tolerance);
            for (int c = 0; c < N; ++c) {
                BOOST_CHECK_SMALL(kernel.x[c] - reference.x[c], tolerance);
                BOOST_CHECK_SMALL(kernel.y[c] - reference.y[c], tolerance);
            }
        }
    }
    Answer<FluidSystem> a;
    a.ok = kernel_ok;
    a.single = kernel.singlePhase;
    a.L = kernel.L;
    for (int c = 0; c < N; ++c) {
        a.x[c] = kernel.x[c];
        a.y[c] = kernel.y[c];
        a.K[c] = kernel.K[c];
    }
    return a;
}

template <class FluidSystem>
std::array<double, FluidSystem::numComponents> wilson(double p, double T)
{
    std::array<double, FluidSystem::numComponents> K0;
    for (int c = 0; c < FluidSystem::numComponents; ++c) {
        K0[c] = std::exp(5.37 * (1 + FluidSystem::acentricFactor(c))
                         * (1 - FluidSystem::criticalTemperature(c) / T))
            * (FluidSystem::criticalPressure(c) / p);
    }
    return K0;
}

// The SPE5 oil of bench/olav.cpp (MultiComponentFlash.jl's benchmark fluid).
using Spe5 = Opm::GenericOilGasWaterFluidSystem<double, 6, false>;

void registerSpe5()
{
    const double psi = 6.894757293168360e+03;
    const double rankine = 0.555555555555556;
    const std::array<const char*, 6> names{"C1", "C3", "C6", "C10", "C15", "C20"};
    const std::array<double, 6> tcR{343, 665.7, 913.4, 1111.8, 1270.0, 1380.0};
    const std::array<double, 6> pcPsi{667.8, 616.3, 436.9, 304.0, 200.0, 162.0};
    const std::array<double, 6> mw{16.040, 44.100, 86.180, 142.290, 206.000, 282.000};
    const std::array<double, 6> acc{0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500};
    const std::array<double, 6> zc{0.290, 0.277, 0.264, 0.257, 0.245, 0.235};
    Spe5::init();
    for (int i = 0; i < 6; ++i) {
        const double tc = tcR[i] * rankine;
        const double pc = pcPsi[i] * psi;
        const double vc = zc[i] * 8.314 * tc / pc * 1e3;
        Spe5::addComponent(typename Spe5::ComponentParam{names[i], mw[i], tc, pc, vc, acc[i]});
    }
    std::array<std::array<double, 6>, 6> k{};
    k[0][4] = k[4][0] = 0.05;
    k[0][5] = k[5][0] = 0.05;
    k[1][4] = k[4][1] = 0.005;
    k[1][5] = k[5][1] = 0.005;
    std::vector<double> bic;
    for (int row = 1; row < 6; ++row) {
        for (int col = 0; col < row; ++col) {
            bic.push_back(k[row][col]);
        }
    }
    Spe5::setInteractionCoefficients(std::move(bic));
}

struct Components
{
    Components()
    {
        Opm::FlashTest::registerF2Components<double>();
        Opm::FlashTest::registerH1Components<double>();
        registerSpe5();
    }
};

} // namespace

BOOST_GLOBAL_FIXTURE(Components);

BOOST_AUTO_TEST_CASE(TernaryEveryMethodColdAndWarm)
{
    using FS = Opm::FlashTest::F2FluidSystem<double>;
    const auto& z = Opm::FlashTest::f2Z;
    const EOSType eos = EOSType::PR;
    int two_phase = 0;
    for (const auto method : {Opm::PTFlashMethod::SsiNewton, Opm::PTFlashMethod::Ssi,
                              Opm::PTFlashMethod::SsiAccelerated, Opm::PTFlashMethod::Newton}) {
        for (const double p : {5e5, 10e5, 30e5, 60e5, 100e5}) {
            for (const double T : {250., 300., 350., 450.}) {
                const auto cold = checkState<FS>(p, T, z, wilson<FS>(p, T), -1.0, method, eos);
                two_phase += cold.ok && !cold.single;
                // A warm start from the answer at a nearby state, as the
                // simulator's next Newton iteration.
                if (cold.ok && !cold.single) {
                    checkState<FS>(p * 1.02, T + 1.0, z, cold.K, cold.L, method, eos);
                }
            }
        }
    }
    BOOST_CHECK_GT(two_phase, 0);
}

BOOST_AUTO_TEST_CASE(VolumeShiftedBinary)
{
    using FS = Opm::FlashTest::H1FluidSystem<double>;
    const auto& z = Opm::FlashTest::h1Z;
    for (const double T : {250., 270., 290., 300., 400., 600.}) {
        for (const double p : {Opm::FlashTest::h1Pressure, 150e5}) {
            checkState<FS>(p, T, z, wilson<FS>(p, T), -1.0, Opm::PTFlashMethod::SsiNewton, EOSType::PR);
        }
    }
}

BOOST_AUTO_TEST_CASE(Spe5NearCritical)
{
    // The grid of bench/olav.cpp with 1000 points per axis.
    const auto gridP = [](int ip) { return 1e5 + (4e7 - 1e5) * ip / 999; };
    const auto gridT = [](int it) { return 280.0 + 420.0 * it / 999; };
    const std::array<double, 6> z{0.5, 0.03, 0.07, 0.20, 0.15, 0.05};
    // (ip, it): ordinary states, two-phase states whose Gibbs-energy Newton
    // gives up (the Newton in (x, y, L) finishes), and the failure of both.
    const std::vector<std::array<int, 2>> states{
        {0, 0}, {100, 500}, {500, 500}, {999, 999}, {250, 700}, {205, 857}, // ordinary
        {198, 860}, {214, 854}, {243, 828}, {252, 815},                     // Newton fallback
        {205, 858}};                                                        // both fail
    for (const auto& [ip, it] : states) {
        const double p = gridP(ip);
        const double T = gridT(it);
        checkState<Spe5>(p, T, z, wilson<Spe5>(p, T), -1.0, Opm::PTFlashMethod::SsiNewton, EOSType::PR);
    }
    const double p = gridP(205);
    const double T = gridT(858);
    const auto failed = ptflash<Spe5>(p, T, z, wilson<Spe5>(p, T), -1.0, Opm::PTFlashMethod::SsiNewton,
                                      EOSType::PR);
    BOOST_CHECK(!failed.ok);
}
