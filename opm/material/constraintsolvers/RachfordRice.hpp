// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 NORCE.
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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
 * \copydoc Opm::RachfordRice
 */
#ifndef OPM_RACHFORD_RICE_HPP
#define OPM_RACHFORD_RICE_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <dune/common/math.hh>

#include <fmt/format.h>

#include <cmath>
#include <cstddef>
#include <optional>
#include <string_view>
#include <utility>

namespace Opm {

/*!
 * \brief The Rachford-Rice equation and the solvers for it.
 *
 * The phase split follows from the equilibrium ratios and the overall
 * composition alone, so nothing here needs a fluid system: the number of
 * components is the length of the vectors handed in.
 */
class RachfordRice
{
public:
    /*!
     * \brief The liquid fraction the equilibrium ratios and composition imply.
     *
     * The residual in the vapour fraction \f$V = 1 - L\f$,
     *
     * \f[ h(V) = \sum_i \frac{z_i (K_i - 1)}{1 + V (K_i - 1)}, \f]
     *
     * decreases monotonically between its poles, so it has exactly one root
     * there when some \f$K_i > 1\f$ and some \f$K_i < 1\f$ (negative flash:
     * the root may lie outside [0, 1]). At the root both phase compositions
     * sum to one, so no mole fraction exceeds one, which bounds the root by
     * the window of Leibovici and Neoschil (1992),
     *
     * \f[ \max_{K_i > 1} \frac{K_i z_i - 1}{K_i - 1} \le V \le
     *     \min_{K_i < 1} \frac{1 - z_i}{1 - K_i}, \f]
     *
     * tighter than the poles and finite at both ends. Newton's method is
     * kept inside the bracket, which shrinks with the sign of every
     * residual, by a bisection step whenever a Newton step would leave it, so
     * it converges for every set of ratios.
     *
     * \throws NumericalProblem when all ratios of the components present lie
     *         on one side of one: no split satisfies the equation.
     */
    template <class Vector>
    static typename Vector::field_type solve(const Vector& K, const Vector& z, int verbosity)
    {
        return solveInWindow(K, z, std::optional<typename Vector::field_type>{}, verbosity);
    }

    /*!
     * \brief As solve(K, z, verbosity), with Newton started from the liquid
     *        fraction \p initialL when it lies inside the window, e.g. the
     *        previous one of an iteration that updates K.
     */
    template <class Vector>
    static typename Vector::field_type solve(const Vector& K,
                                             const Vector& z,
                                             typename Vector::field_type initialL,
                                             int verbosity)
    {
        return solveInWindow(K, z, std::optional<typename Vector::field_type>{initialL}, verbosity);
    }

private:
    template <class Vector>
    static typename Vector::field_type
    solveInWindow(const Vector& K,
                  const Vector& z,
                  const std::optional<typename Vector::field_type>& initialL,
                  int verbosity)
    {
        using field_type = typename Vector::field_type;
        constexpr field_type tol = 1e-12;
        constexpr int max_iterations = 200;

        // The window of the root, from the components present.
        bool has_ratio_above_one = false;
        bool has_ratio_below_one = false;
        field_type lower = 0.0;
        field_type upper = 0.0;
        for (std::size_t compIdx = 0; compIdx < K.size(); ++compIdx) {
            if (!(z[compIdx] > 0.0)) {
                continue;
            }
            const field_type dK = K[compIdx] - 1.0;
            if (dK > 0.0) {
                const field_type bound = (K[compIdx] * z[compIdx] - 1.0) / dK;
                lower = has_ratio_above_one && lower > bound ? lower : bound;
                has_ratio_above_one = true;
            }
            else if (dK < 0.0) {
                const field_type bound = (1.0 - z[compIdx]) / (-dK);
                upper = has_ratio_below_one && upper < bound ? upper : bound;
                has_ratio_below_one = true;
            }
        }
        if (!has_ratio_above_one || !has_ratio_below_one) {
            OPM_THROW_NOLOG(NumericalProblem,
                            "Rachford-Rice has no root: all equilibrium ratios lie on one "
                            "side of one");
        }

        field_type V = (lower + upper) / 2;
        if (initialL) {
            const field_type initialV = 1.0 - *initialL;
            if (initialV > lower && initialV < upper) {
                V = initialV;
            }
        }
        if (verbosity == 3 || verbosity == 4) {
            OpmLog::debug(fmt::format("Initial guess {}c : V = {} and window [{}, {}]",
                                      K.size(), V, lower, upper));
            OpmLog::debug(fmt::format("{:>10}{:>16}{:>16}", "Iteration", "h(V)", "V"));
        }

        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            field_type h = 0.0;
            field_type dh = 0.0;
            for (std::size_t compIdx = 0; compIdx < K.size(); ++compIdx) {
                const field_type dK = K[compIdx] - 1.0;
                const field_type denominator = 1.0 + V * dK;
                h += z[compIdx] * dK / denominator;
                dh -= z[compIdx] * dK * dK / (denominator * denominator);
            }
            // h decreases with V: a positive residual puts the root above V.
            if (h > 0.0) {
                lower = V;
            }
            else {
                upper = V;
            }
            field_type next = V - h / dh;
            if (verbosity == 3 || verbosity == 4) {
                OpmLog::debug(fmt::format("{:>10}{:>16}{:>16}", iteration, h, next));
            }
            const field_type abs_h = h < 0.0 ? -h : h;
            if (abs_h < tol) {
                return reportSolution(1.0 - next, verbosity);
            }
            if (!(upper - lower > 1e-15 * (1.0 + (V < 0.0 ? -V : V)))) {
                return reportSolution(1.0 - V, verbosity);
            }
            if (!(next > lower && next < upper)) {
                next = (lower + upper) / 2;
            }
            V = next;
        }
        OPM_THROW_NOLOG(NumericalProblem,
                        " Rachford-Rice did not converge within maximum number of iterations");
    }

    //! \brief Log the converged liquid fraction and hand it back.
    template <class FieldType>
    static FieldType reportSolution(FieldType L, int verbosity, std::string_view method = {})
    {
        if (verbosity >= 1) {
            OpmLog::debug(
                fmt::format("Rachford-Rice{} converged to final solution L = {}", method, L));
        }
        return L;
    }
};

} // namespace Opm

#endif // OPM_RACHFORD_RICE_HPP
