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
 * \brief Name-keyed lookup of the components' ideal-gas heat-capacity
 *        polynomials and the enthalpy reference state used by the caloric
 *        mixture-enthalpy model (MixtureEnthalpy) and the isenthalpic (P-H)
 *        flash.
 *
 * The polynomial coefficients live on the component classes (their caloric
 * identity, next to the EoS constants); the polynomial mathematics lives in
 * ComponentCp.hpp. This header owns the shared reference datum and the
 * deck-style name lookup.
 *
 * Units are SI throughout: temperature [K], molar heat capacity [J/(mol K)],
 * molar enthalpy [J/mol]. Enthalpy is zero at the reference temperature.
 */
#ifndef OPM_IDEAL_GAS_CALORIC_DATA_HPP
#define OPM_IDEAL_GAS_CALORIC_DATA_HPP

#include <opm/material/components/Air.hpp>
#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/ComponentCp.hpp>
#include <opm/material/components/H2.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/SimpleCO2.hpp>

#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Opm {

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

    // The coefficients live on the COMPONENT CLASSES (C1.hpp, C10.hpp,
    // SimpleCO2.hpp — idealGasHeatCapacityPolynomial()), each next to the
    // species' other constants and its provenance block: one identity card
    // per species. The wrappers below are the stable lookup surface; the
    // fitting recipe (least-squares cubic to the reference-EoS ideal-gas cp,
    // 250-600 K window, CoolProp 8.0.0 as extraction tool) is documented on
    // the classes. A unit test pins each preset against tabulated reference
    // values so a corrupt coefficient row cannot enter silently (an earlier
    // n-decane row of untraceable origin was ~32% low, which no
    // self-consistent round-trip test could detect).

    //! methane (C1) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> methane()
    { return C1<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! n-decane (nC10) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> decane()
    { return C10<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! carbon dioxide (CO2) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> carbonDioxide()
    { return SimpleCO2<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! nitrogen (N2) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> nitrogen()
    { return N2<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! water (H2O) ideal-gas cp polynomial [J/(mol K)] — window 275–600 K
    static constexpr ComponentCp<Scalar> water()
    { return H2O<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! hydrogen (H2) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> hydrogen()
    { return H2<Scalar>::idealGasHeatCapacityPolynomial(); }

    //! air (pseudo-component) ideal-gas cp polynomial [J/(mol K)]
    static constexpr ComponentCp<Scalar> air()
    { return Air<Scalar>::idealGasHeatCapacityPolynomial(); }

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
        if (n == "N2" || n == "NITROGEN")
            return nitrogen();
        if (n == "H2O" || n == "WATER")
            return water();
        if (n == "H2" || n == "HYDROGEN")
            return hydrogen();
        if (n == "AIR")
            return air();

        throw std::runtime_error(
            "IdealGasCaloricData: no ideal-gas heat-capacity preset for component '"
            + std::string(name) + "' — supply coefficients explicitly");
    }
};

} // namespace Opm

#endif // OPM_IDEAL_GAS_CALORIC_DATA_HPP
