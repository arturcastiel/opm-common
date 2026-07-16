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
 * \brief Cited binary-interaction coefficients (kij) for the component
 *        pairs used by this library's compositional tests and tools —
 *        written down once, consumed everywhere.
 *
 * kij is a PAIR property (it enters the cubic-EoS mixing rule
 * a_mix = sum_i sum_j x_i x_j sqrt(a_i a_j) (1 - k_ij)), so it cannot live
 * on a single component's class. This header is where the literal values
 * are authored; it changes NO machinery:
 *
 *  - fluid systems keep owning their kij at runtime
 *    (BaseFluidSystem::interactionCoefficient is the consumption interface);
 *  - deck-driven runs keep the BIC/BICS keywords as the authoritative
 *    source — this table is for the non-deck consumers (tests, opmflash)
 *    that would otherwise each carry their own copy.
 *
 * The values are the standard Peng-Robinson literature coefficients used
 * throughout this library's compositional tests. Scope note: kij values are
 * EoS-family-specific — these are PR values; a future SRK consumer extends
 * the lookup with its own cited set rather than reusing these.
 *
 * Unlisted pairs return 0.0 — the standard convention (ideal geometric-mean
 * mixing), also what every fluid system defaults to.
 */
#ifndef OPM_BINARY_INTERACTION_HPP
#define OPM_BINARY_INTERACTION_HPP

#include <string_view>

namespace Opm {

template <class Scalar>
struct BinaryInteraction {
    /*!
     * \brief Peng-Robinson kij for the (a, b) component pair, by canonical
     *        component name (the names the component classes report:
     *        "C1", "C10", "CO2"). Symmetric; unlisted pairs are 0.0.
     */
    static constexpr Scalar kij(std::string_view a, std::string_view b)
    {
        constexpr auto is = [](std::string_view x, std::string_view y,
                               std::string_view p, std::string_view q)
        { return (x == p && y == q) || (x == q && y == p); };

        if (is(a, b, "C1", "C10")) { return 0.0411; }
        if (is(a, b, "C1", "CO2")) { return 0.10; }
        if (is(a, b, "CO2", "C10")) { return 0.10; }
        return 0.0;
    }
};

} // namespace Opm

#endif // OPM_BINARY_INTERACTION_HPP
