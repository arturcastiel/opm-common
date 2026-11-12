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
 * \brief The host side of PTFlashKernel.hpp: the kernel's parameter struct
 *        from a fluid system, and the PTFlashMethod and EOS type in the
 *        kernel's terms.
 */
#ifndef OPM_PTFLASH_KERNEL_HOST_HPP
#define OPM_PTFLASH_KERNEL_HOST_HPP

#include <opm/material/constraintsolvers/PTFlashKernel.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <stdexcept>

namespace Opm::PTFlashKernel {

inline Eos eosFrom(CompositionalConfig::EOSType eos_type)
{
    using EOSType = CompositionalConfig::EOSType;
    switch (eos_type) {
    case EOSType::PR:
        return Eos::PR;
    case EOSType::PRCORR:
        return Eos::PRCORR;
    case EOSType::SRK:
        return Eos::SRK;
    case EOSType::RK:
        return Eos::RK;
    default:
        throw std::runtime_error("PTFlashKernel: EOS type not implemented");
    }
}

inline Method methodFrom(PTFlashMethod method)
{
    switch (method) {
    case PTFlashMethod::Ssi:
        return Method::Ssi;
    case PTFlashMethod::SsiAccelerated:
        return Method::SsiAccelerated;
    case PTFlashMethod::Newton:
        return Method::Newton;
    case PTFlashMethod::SsiNewton:
        return Method::SsiNewton;
    }
    throw std::runtime_error("PTFlashKernel: unknown method");
}

/*!
 * \brief The kernel's parameters from the static data of a fluid system
 *        (GenericOilGasWaterFluidSystem or one with its interface), with
 *        the values CubicEOSParams uses: the per-component OMEGAA as a
 *        factor on the equation's own Omega_a, OMEGAB or the equation's own,
 *        and the volume shift where the fluid system has one.
 */
template <class Scalar, class FluidSystem>
Parameters<Scalar, FluidSystem::numComponents>
makeParameters(CompositionalConfig::EOSType eos_type)
{
    constexpr int N = FluidSystem::numComponents;
    Parameters<Scalar, N> params;
    params.eos = eosFrom(eos_type);
    const Scalar constant_a = detail::constantOmegaA<Scalar>(params.eos);
    const Scalar constant_b = detail::constantOmegaB<Scalar>(params.eos);
    for (int i = 0; i < N; ++i) {
        params.criticalTemperature[i] = FluidSystem::criticalTemperature(i);
        params.criticalPressure[i] = FluidSystem::criticalPressure(i);
        params.acentricFactor[i] = FluidSystem::acentricFactor(i);
        params.criticalVolume[i] = FluidSystem::criticalVolume(i);
        params.volumeShift[i] = 0.0;
        params.omegaAScale[i] = 1.0;
        params.omegaB[i] = constant_b;
        if constexpr (requires { FluidSystem::volumeShift(0u); }) {
            params.volumeShift[i] = FluidSystem::volumeShift(i);
        }
        if constexpr (requires { FluidSystem::omegaA(0u); }) {
            if (const auto omega_a = FluidSystem::omegaA(i)) {
                params.omegaAScale[i] = *omega_a / constant_a;
            }
        }
        if constexpr (requires { FluidSystem::omegaB(0u); }) {
            if (const auto omega_b = FluidSystem::omegaB(i)) {
                params.omegaB[i] = *omega_b;
            }
        }
        for (int j = 0; j < N; ++j) {
            params.interaction[i][j] = i == j ? Scalar(0.0) : FluidSystem::interactionCoefficient(i, j);
        }
    }
    return params;
}

} // namespace Opm::PTFlashKernel

#endif // OPM_PTFLASH_KERNEL_HOST_HPP
