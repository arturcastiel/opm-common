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

    // The coefficients below follow the standard ideal-gas heat-capacity
    // polynomial tabulations cp = c0 + c1*T + c2*T^2 + c3*T^3 (cf. Poling,
    // Prausnitz & O'Connell, "The Properties of Gases and Liquids"), with a
    // nominal fit validity of roughly 273-1500 K. Within the P-H stack they
    // are self-consistent by construction (the isenthalpic round-trip
    // manufactures its target enthalpy from the same table); for studies
    // where absolute enthalpy values matter, confirm the coefficients
    // against the cited tabulations.

    //! methane (C1) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> methane()
    { return {19.25, 5.213e-2, 1.197e-5, -1.132e-8}; }

    //! n-decane (nC10) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> decane()
    { return {16.35, 5.762e-1, -3.115e-4, 6.62e-8}; }

    //! carbon dioxide (CO2) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> carbonDioxide()
    { return {19.80, 7.344e-2, -5.602e-5, 1.715e-8}; }

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
