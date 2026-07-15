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
 * \brief Ideal-gas heat-capacity polynomials and the enthalpy reference state
 *        used by the caloric mixture-enthalpy model (MixtureEnthalpy) and the
 *        isenthalpic (P-H) flash.
 *
 * Units are SI throughout: temperature [K], molar heat capacity [J/(mol K)],
 * molar enthalpy [J/mol]. Enthalpy is zero at the reference temperature.
 */
#ifndef OPM_IDEAL_GAS_CALORIC_DATA_HPP
#define OPM_IDEAL_GAS_CALORIC_DATA_HPP

#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Opm {

/*!
 * \brief Cubic ideal-gas heat-capacity polynomial of one component:
 *        cp(T) = c0 + c1*T + c2*T^2 + c3*T^3   [J/(mol K)]
 */
template <class Scalar>
struct ComponentCp {
    Scalar c0, c1, c2, c3;

    //! cp(T) [J/(mol K)]. Generic in the evaluation type (double or AD).
    template <class Eval>
    Eval heatCapacity(const Eval& T) const
    {
        return c0 + c1*T + c2*T*T + c3*T*T*T;
    }

    /*!
     * \brief Ideal-gas enthalpy h(T) = int_{T0}^{T} cp dT' [J/mol],
     *        in closed form. h(T0) = 0 by construction.
     */
    template <class Eval>
    Eval enthalpyIntegral(const Eval& T, const Scalar T0) const
    {
        return c0*(T - T0)
             + c1/2*(T*T - T0*T0)
             + c2/3*(T*T*T - T0*T0*T0)
             + c3/4*(T*T*T*T - T0*T0*T0*T0);
    }
};

//! Per-component cp table for an N-component fluid system, indexed like the
//! fluid system's component indices.
template <class Scalar, int numComponents>
using CpTable = std::array<ComponentCp<Scalar>, numComponents>;

/*!
 * \brief The enthalpy reference state and component cp presets.
 *
 * The reference (datum) is fixed ONCE for the whole P-H stack: any specified
 * enthalpy H_spec handed to the isenthalpic flash must be expressed against
 * the same datum, H(referenceTemperature) = 0.
 */
template <class Scalar>
struct IdealGasCaloricData {
    //! reference temperature T0 [K]; enthalpy is zero here
    static constexpr Scalar referenceTemperature() { return 298.15; }

    //! reference pressure P0 [Pa] (documentation of the datum; the ideal-gas
    //! caloric enthalpy itself is pressure-independent)
    static constexpr Scalar referencePressure() { return 1e5; }

    // The coefficients below are least-squares cubic fits,
    // cp = c0 + c1*T + c2*T^2 + c3*T^3, to the ideal-gas heat capacity of
    // each fluid's reference equation of state, fitted over 250-600 K (the
    // temperature window the isenthalpic flash practically operates in):
    //   methane: Setzmann & Wagner, J. Phys. Chem. Ref. Data 20 (1991) 1061
    //            (fit RMS 0.024, max 0.063 J/(mol K))
    //   n-decane: Lemmon & Span, J. Chem. Eng. Data 51 (2006) 785
    //            (fit RMS 0.27, max 0.79 J/(mol K))
    //   CO2:     Span & Wagner, J. Phys. Chem. Ref. Data 25 (1996) 1509
    //            (fit RMS 0.003, max 0.014 J/(mol K))
    // The reference curves were sampled from the fluids' Helmholtz ideal
    // parts (CoolProp 8.0.0 as the extraction tool; CoolProp itself is not a
    // dependency). Outside 250-600 K the cubics extrapolate — refit rather
    // than trust them there. A unit test pins each preset against tabulated
    // reference values so a corrupt coefficient row cannot enter silently
    // (an earlier n-decane row of untraceable origin was ~32% low, which
    // no self-consistent round-trip test could detect).

    //! methane (C1) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> methane()
    { return {40.1503, -8.47372e-2, 2.93012e-4, -1.96125e-7}; }

    //! n-decane (nC10) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> decane()
    { return {79.4791, 3.1066e-1, 9.88317e-4, -1.00245e-6}; }

    //! carbon dioxide (CO2) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> carbonDioxide()
    { return {18.2687, 8.36359e-2, -7.75148e-5, 3.14088e-8}; }

    /*!
     * \brief Preset lookup by component name (deck-style aliases,
     *        case-insensitive).
     *
     * Throws std::runtime_error naming the component when no preset exists:
     * there is deliberately NO silent fallback — an unknown component must
     * fail loudly rather than receive somebody else's heat capacity.
     */
    static ComponentCp<Scalar> byName(const std::string_view name)
    {
        std::string n(name);
        for (auto& c : n)
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

        if (n == "C1" || n == "CH4" || n == "METHANE")
            return methane();
        if (n == "C10" || n == "NC10" || n == "DECANE" || n == "N-DECANE")
            return decane();
        if (n == "CO2" || n == "CARBONDIOXIDE" || n == "CARBON-DIOXIDE" || n == "CARBON DIOXIDE")
            return carbonDioxide();

        throw std::runtime_error(
            "IdealGasCaloricData: no ideal-gas heat-capacity preset for component '"
            + std::string(name) + "' — supply coefficients explicitly");
    }
};

} // namespace Opm

#endif // OPM_IDEAL_GAS_CALORIC_DATA_HPP
