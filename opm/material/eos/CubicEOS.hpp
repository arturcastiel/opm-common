// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
#ifndef CUBIC_EOS_HPP
#define CUBIC_EOS_HPP

#include <opm/material/Constants.hpp>
#include <opm/material/common/PolynomialUtils.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <array>
#include <type_traits>

namespace Opm
{
/*!
 * \brief A two-parameter cubic equation of state in the common form
 *
 * \f[ p = \frac{RT}{V - b} - \frac{a}{(V + m_1 b)(V + m_2 b)}, \f]
 *
 * with the mixing rules and the constants \f$m_1, m_2\f$ of the chosen
 * equation (Peng-Robinson, Soave-Redlich-Kwong, ...) supplied by the
 * parameter cache. The cache carries the dimensionless
 * \f$A = a p / (RT)^2\f$ and \f$B = b p / (RT)\f$, in which the cubic is
 * solved for the compressibility factor \f$Z = pV/(RT)\f$.
 */
template<class Scalar, class FluidSystem>
class CubicEOS
{
    enum { numComponents = FluidSystem::numComponents };

    static constexpr Scalar R = Constants<Scalar>::R;

    // The fugacity coefficient is clamped so a vanishing mole fraction still
    // yields a usable fugacity; the molar volume is floored against a root
    // that has collapsed.
    static constexpr Scalar maxFugacityCoefficient = 1e10;
    static constexpr Scalar minFugacityCoefficient = 1e-10;
    static constexpr Scalar minMolarVolume = 1e-7;

public:
    /*!
     * \brief The fugacity coefficient of a component in a phase.
     *
     * For the general two-parameter cubic,
     *
     * \f[ \ln \phi_i = \frac{B_i}{B}(Z - 1) - \ln(Z - B)
     *   + \frac{A}{(m_1 - m_2) B} \ln\frac{Z + m_2 B}{Z + m_1 B}
     *     \left( \frac{2 \sum_j x_j A_{ij}}{A} - \frac{B_i}{B} \right), \f]
     *
     * evaluated at the compressibility factor of the phase's cached molar
     * volume. The result is clamped to a wide range so that an unphysical
     * intermediate state during a flash still yields a finite, usable value.
     */
    template <class FluidState, class Params, class LhsEval = typename FluidState::ValueType>
    static LhsEval computeFugacityCoefficient(const FluidState& fs,
                                              const Params& params,
                                              unsigned phaseIdx,
                                              unsigned compIdx)
    {
        return componentFugacityCoefficient_(
            phaseFugacityTerms_<LhsEval>(fs, params, phaseIdx), fs, params, phaseIdx, compIdx);
    }

    /*!
     * \brief The fugacity coefficients of all components of a phase.
     *
     * The same values as computeFugacityCoefficient() for each component;
     * the terms that do not depend on the component, among them both
     * logarithms, are evaluated once for the phase.
     *
     * \param fugCoeffs Receives the coefficient of component i at index i.
     */
    template <class FluidState, class Params, class Coefficients>
    static void computeFugacityCoefficients(const FluidState& fs,
                                            const Params& params,
                                            unsigned phaseIdx,
                                            Coefficients& fugCoeffs)
    {
        using LhsEval = typename FluidState::ValueType;
        const auto terms = phaseFugacityTerms_<LhsEval>(fs, params, phaseIdx);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fugCoeffs[compIdx] = componentFugacityCoefficient_(terms, fs, params, phaseIdx, compIdx);
        }
    }

    /*!
     * \brief The fugacity coefficients of all components of a phase, as
     *        computeFugacityCoefficients() gives them, and the derivatives
     *        of their logarithms in the mole fractions of the phase.
     *
     * The mole fractions are the independent variables, as automatic
     * differentiation in them gives: with \f$S_i = \sum_j A_{ij} x_j\f$,
     * \f$\partial A / \partial x_j = 2 S_j\f$ and
     * \f$\partial B / \partial x_j = B_j\f$, the compressibility factor follows
     * the cubic \f$F(Z, A, B) = 0\f$,
     *
     * \f[ \frac{\partial Z}{\partial x_j}
     *   = -\frac{F_A\, \partial A / \partial x_j + F_B\, \partial B / \partial x_j}{F_Z}, \f]
     *
     * and \f$\partial \ln\phi_i / \partial x_j\f$ is the derivative of the
     * expression of computeFugacityCoefficient() term by term. The mole
     * fractions are taken in \f$[0, 1]\f$, where the mixing rule does not
     * limit them. A coefficient at its limits has no derivative, as for
     * automatic differentiation.
     *
     * Mole-number derivatives follow from
     * \f$n \partial \ln\phi_i / \partial n_j
     *   = \partial \ln\phi_i / \partial x_j - \sum_k x_k \partial \ln\phi_i / \partial x_k\f$.
     *
     * \param fugCoeffs Receives the coefficient of component i at index i.
     * \param lnPhiDerivatives Receives \f$\partial \ln\phi_i / \partial x_j\f$ at [i][j].
     */
    template <class FluidState, class Params, class Coefficients, class Derivatives>
    static void computeFugacityCoefficientsAndDerivatives(const FluidState& fs,
                                                          const Params& params,
                                                          unsigned phaseIdx,
                                                          Coefficients& fugCoeffs,
                                                          Derivatives& lnPhiDerivatives)
    {
        Coefficients* noStateDerivatives = nullptr;
        coefficientsAndDerivatives_(fs, params, phaseIdx, fugCoeffs, lnPhiDerivatives,
                                    noStateDerivatives, noStateDerivatives);
    }

    /*!
     * \brief computeFugacityCoefficientsAndDerivatives(), and the derivatives
     *        of the logarithms in pressure at fixed temperature and phase
     *        composition.
     *
     * \f$A_{ij}\f$, \f$B_i\f$ and so \f$A\f$ and \f$B\f$ are proportional to
     * pressure, so \f$\partial A / \partial p = A / p\f$,
     * \f$\partial B / \partial p = B / p\f$, and only \f$Z\f$ and the
     * logarithms move:
     *
     * \f[ \frac{\partial \ln\phi_i}{\partial p}
     *   = \frac{B_i}{B} Z_p - \frac{Z_p - B_p}{Z - B}
     *   + \frac{A}{(m_1 - m_2) B}
     *     \left(\frac{Z_p + m_2 B_p}{Z + m_2 B} - \frac{Z_p + m_1 B_p}{Z + m_1 B}\right)
     *     \left(\frac{2 \sum_j x_j A_{ij}}{A} - \frac{B_i}{B}\right), \f]
     *
     * with \f$Z_p\f$ from the cubic, or \f$Z / p\f$ where the molar volume is
     * at its floor. A coefficient at its limits has no derivative.
     *
     * \param lnPhiPressureDerivatives Receives \f$\partial \ln\phi_i / \partial p\f$ at i.
     */
    template <class FluidState, class Params, class Coefficients, class Derivatives>
    static void computeFugacityCoefficientsAndDerivatives(const FluidState& fs,
                                                          const Params& params,
                                                          unsigned phaseIdx,
                                                          Coefficients& fugCoeffs,
                                                          Derivatives& lnPhiDerivatives,
                                                          Coefficients& lnPhiPressureDerivatives)
    {
        Coefficients* noTemperatureDerivatives = nullptr;
        coefficientsAndDerivatives_(fs, params, phaseIdx, fugCoeffs, lnPhiDerivatives,
                                    &lnPhiPressureDerivatives, noTemperatureDerivatives);
    }

    /*!
     * \brief computeFugacityCoefficientsAndDerivatives() with the pressure
     *        derivatives, and the derivatives of the logarithms in temperature
     *        at fixed pressure and phase composition.
     *
     * With \f$a_i = \partial \ln A_i / \partial T\f$ of the parameter cache,
     * \f$\partial A_{ij} / \partial T = A_{ij} (a_i + a_j) / 2\f$, so
     * \f$A_T = \sum_i x_i a_i S_i\f$ and
     * \f$\partial S_i / \partial T = (a_i S_i + \sum_j A_{ij} x_j a_j) / 2\f$
     * with \f$S_i = \sum_j A_{ij} x_j\f$; \f$B\f$ and every \f$B_i\f$ are
     * proportional to \f$1 / T\f$. \f$Z_T\f$ follows from the cubic, or is
     * \f$-Z / T\f$ where the molar volume is at its floor. A coefficient at
     * its limits has no derivative.
     *
     * \param lnPhiTemperatureDerivatives Receives \f$\partial \ln\phi_i / \partial T\f$ at i.
     */
    template <class FluidState, class Params, class Coefficients, class Derivatives>
    static void computeFugacityCoefficientsAndDerivatives(const FluidState& fs,
                                                          const Params& params,
                                                          unsigned phaseIdx,
                                                          Coefficients& fugCoeffs,
                                                          Derivatives& lnPhiDerivatives,
                                                          Coefficients& lnPhiPressureDerivatives,
                                                          Coefficients& lnPhiTemperatureDerivatives)
    {
        coefficientsAndDerivatives_(fs, params, phaseIdx, fugCoeffs, lnPhiDerivatives,
                                    &lnPhiPressureDerivatives, &lnPhiTemperatureDerivatives);
    }

    /*!
     * \brief The temperature derivative at fixed pressure and composition of
     *        \f$\sum_i x_i \ln\phi_i = Z - 1 - \ln(Z - B)
     *        + \frac{A}{(m_1 - m_2) B} \ln\frac{Z + m_2 B}{Z + m_1 B}\f$,
     *        the sum without the limits of the coefficients and without a
     *        volume translation, on the phase's cached molar volume.
     *
     * \f$A_T\f$, \f$B_T\f$ and \f$Z_T\f$ as in
     * computeFugacityCoefficientsAndDerivatives() with the temperature
     * derivatives.
     */
    template <class FluidState, class Params>
    static typename FluidState::ValueType
    sumLnPhiTemperatureDerivative(const FluidState& fs, const Params& params, unsigned phaseIdx)
    {
        using LhsEval = typename FluidState::ValueType;
        static_assert(std::is_floating_point_v<LhsEval>,
                      "The derivative is that of a value-only evaluation");
        const LhsEval T = fs.temperature(phaseIdx);
        const LhsEval p = fs.pressure(phaseIdx);
        const LhsEval A = params.A(phaseIdx);
        const LhsEval B = params.B(phaseIdx);
        const LhsEval m1 = params.m1(phaseIdx);
        const LhsEval m2 = params.m2(phaseIdx);
        const LhsEval Z = p * params.molarVolume(phaseIdx) / (R * T);

        std::array<LhsEval, numComponents> xa; // x_j a_j
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            xa[compIdx] = fs.moleFraction(phaseIdx, compIdx)
                * params.lnAiTemperatureDerivative(phaseIdx, compIdx);
        }
        LhsEval AT = 0.0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            LhsEval S = 0.0;
            for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                S += params.aCache(phaseIdx, compIdx, compJIdx) * fs.moleFraction(phaseIdx, compJIdx);
            }
            AT += xa[compIdx] * S;
        }
        const LhsEval BT = -B / T;
        const LhsEval c2 = (m1 + m2 - 1) * B - 1;
        const LhsEval c1 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1);
        const LhsEval FZ = (3 * Z + 2 * c2) * Z + c1;
        const LhsEval FA = Z - B;
        const LhsEval FB = (m1 + m2 - 1) * Z * Z + (2 * m1 * m2 * B - (m1 + m2) * (2 * B + 1)) * Z
            - (A + m1 * m2 * (3 * B + 2) * B);
        const bool volumeAtFloor = !(params.molarVolume(phaseIdx) > minMolarVolume);
        const LhsEval ZT = volumeAtFloor ? -Z / T : -(FA * AT + FB * BT) / FZ;

        const LhsEval P1 = Z + m1 * B;
        const LhsEval P2 = Z + m2 * B;
        const LhsEval Lambda = log(P2 / P1);
        const LhsEval LambdaT = (ZT + m2 * BT) / P2 - (ZT + m1 * BT) / P1;
        const LhsEval invB = 1 / B;
        return ZT - (ZT - BT) / (Z - B)
            + ((AT * invB - A * BT * invB * invB) * Lambda + A * invB * LambdaT) / (m1 - m2);
    }

private:
    template <class FluidState, class Params, class Coefficients, class Derivatives>
    static void coefficientsAndDerivatives_(const FluidState& fs,
                                            const Params& params,
                                            unsigned phaseIdx,
                                            Coefficients& fugCoeffs,
                                            Derivatives& lnPhiDerivatives,
                                            Coefficients* lnPhiPressureDerivatives,
                                            Coefficients* lnPhiTemperatureDerivatives)
    {
        using LhsEval = typename FluidState::ValueType;
        static_assert(std::is_floating_point_v<LhsEval>,
                      "The derivatives are those of a value-only evaluation");
        const auto terms = phaseFugacityTerms_<LhsEval>(fs, params, phaseIdx);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fugCoeffs[compIdx] = componentFugacityCoefficient_(terms, fs, params, phaseIdx, compIdx);
        }

        const LhsEval A = params.A(phaseIdx);
        const LhsEval B = params.B(phaseIdx);
        const LhsEval m1 = params.m1(phaseIdx);
        const LhsEval m2 = params.m2(phaseIdx);
        const LhsEval Z = terms.Z;
        std::array<LhsEval, numComponents> S;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            S[compIdx] = 0.0;
            for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                S[compIdx] += params.aCache(phaseIdx, compIdx, compJIdx)
                    * fs.moleFraction(phaseIdx, compJIdx);
            }
        }

        // The cubic of computeMolarVolume() and its partial derivatives; a
        // molar volume at its floor does not move.
        const LhsEval c2 = (m1 + m2 - 1) * B - 1;
        const LhsEval c1 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1);
        const LhsEval FZ = (3 * Z + 2 * c2) * Z + c1;
        const LhsEval FA = Z - B;
        const LhsEval FB = (m1 + m2 - 1) * Z * Z + (2 * m1 * m2 * B - (m1 + m2) * (2 * B + 1)) * Z
            - (A + m1 * m2 * (3 * B + 2) * B);
        const bool volumeAtFloor = !(params.molarVolume(phaseIdx) > minMolarVolume);

        const LhsEval ZB = Z - B;
        const LhsEval P1 = Z + m1 * B;
        const LhsEval P2 = Z + m2 * B;
        const LhsEval Lambda = log(P2 / P1);
        const LhsEval Delta = m1 - m2;
        const LhsEval beta = terms.beta;
        const LhsEval Zm1 = Z - 1;
        const LhsEval invA = 1 / A;
        const LhsEval invB = 1 / B;
        std::array<LhsEval, numComponents> gamma;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            gamma[compIdx] = 2 * invA * S[compIdx] - params.Bi(phaseIdx, compIdx) * invB;
        }

        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
            const LhsEval Bj = params.Bi(phaseIdx, compJIdx);
            const LhsEval Aj = 2 * S[compJIdx];
            const LhsEval Zj = volumeAtFloor ? LhsEval{0.0} : -(FA * Aj + FB * Bj) / FZ;
            const LhsEval LambdaJ = (Zj + m2 * Bj) / P2 - (Zj + m1 * Bj) / P1;
            const LhsEval betaJ = ((Aj * invB - A * Bj * invB * invB) * Lambda + A * invB * LambdaJ)
                / Delta;
            // The parts common to every component i.
            const LhsEval common = -(Zj - Bj) / ZB;
            const LhsEval perBi = Zj * invB - Bj * Zm1 * invB * invB + beta * Bj * invB * invB;
            const LhsEval perS = -2 * beta * Aj * invA * invA;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                lnPhiDerivatives[compIdx][compJIdx] = common
                    + params.Bi(phaseIdx, compIdx) * perBi + betaJ * gamma[compIdx]
                    + 2 * beta * invA * params.aCache(phaseIdx, compIdx, compJIdx)
                    + perS * S[compIdx];
            }
        }
        if (lnPhiPressureDerivatives) {
            const LhsEval p = fs.pressure(phaseIdx);
            const LhsEval Ap = A / p;
            const LhsEval Bp = B / p;
            const LhsEval Zp = volumeAtFloor ? Z / p : -(FA * Ap + FB * Bp) / FZ;
            const LhsEval LambdaP = (Zp + m2 * Bp) / P2 - (Zp + m1 * Bp) / P1;
            const LhsEval common = -(Zp - Bp) / ZB;
            const LhsEval betaP = A * invB * LambdaP / Delta;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                (*lnPhiPressureDerivatives)[compIdx] = common
                    + params.Bi(phaseIdx, compIdx) * invB * Zp + betaP * gamma[compIdx];
            }
        }
        if (lnPhiTemperatureDerivatives) {
            const LhsEval T = fs.temperature(phaseIdx);
            std::array<LhsEval, numComponents> a;
            std::array<LhsEval, numComponents> xa; // x_j a_j
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                a[compIdx] = params.lnAiTemperatureDerivative(phaseIdx, compIdx);
                xa[compIdx] = fs.moleFraction(phaseIdx, compIdx) * a[compIdx];
            }
            LhsEval AT = 0.0;
            std::array<LhsEval, numComponents> ST;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                LhsEval weighted = 0.0;
                for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                    weighted += params.aCache(phaseIdx, compIdx, compJIdx) * xa[compJIdx];
                }
                ST[compIdx] = 0.5 * (a[compIdx] * S[compIdx] + weighted);
                AT += xa[compIdx] * S[compIdx];
            }
            const LhsEval BT = -B / T;
            const LhsEval ZT = volumeAtFloor ? -Z / T : -(FA * AT + FB * BT) / FZ;
            const LhsEval LambdaT = (ZT + m2 * BT) / P2 - (ZT + m1 * BT) / P1;
            const LhsEval common = -(ZT - BT) / ZB;
            const LhsEval betaT
                = ((AT * invB - A * BT * invB * invB) * Lambda + A * invB * LambdaT) / Delta;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const LhsEval gammaT = 2 * (ST[compIdx] * invA - S[compIdx] * AT * invA * invA);
                (*lnPhiTemperatureDerivatives)[compIdx] = common
                    + params.Bi(phaseIdx, compIdx) * invB * ZT + betaT * gamma[compIdx]
                    + beta * gammaT;
            }
        }
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (!(fugCoeffs[compIdx] < maxFugacityCoefficient)
                || !(fugCoeffs[compIdx] > minFugacityCoefficient)) {
                for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                    lnPhiDerivatives[compIdx][compJIdx] = 0.0;
                }
                if (lnPhiPressureDerivatives) {
                    (*lnPhiPressureDerivatives)[compIdx] = 0.0;
                }
                if (lnPhiTemperatureDerivatives) {
                    (*lnPhiTemperatureDerivatives)[compIdx] = 0.0;
                }
            }
        }
    }

public:

    /*!
     * \brief The molar volume of a phase from the cubic in \f$Z\f$.
     *
     * \f[ Z^3 + \big((m_1 + m_2 - 1) B - 1\big) Z^2
     *   + \big(A + m_1 m_2 B^2 - (m_1 + m_2) B (B + 1)\big) Z
     *   - \big(A B + m_1 m_2 B^2 (B + 1)\big) = 0. \f]
     *
     * With three real roots the largest belongs to the vapour and the
     * smallest to the liquid; with one root both phases share it. The
     * volume \f$V = Z R T / p\f$ is floored so a collapsed root cannot
     * propagate.
     *
     * \param fs The fluid state holding the phase's pressure and temperature.
     * \param params The parameter cache holding \f$A\f$, \f$B\f$, \f$m_1\f$ and \f$m_2\f$.
     * \param phaseIdx The phase.
     * \param isGasPhase Whether to take the vapour root when there are three.
     */
    template <class FluidState, class Params>
    static typename FluidState::ValueType computeMolarVolume(const FluidState& fs,
                                                          Params& params,
                                                          unsigned phaseIdx,
                                                          bool isGasPhase)
    {
        Valgrind::CheckDefined(fs.temperature(phaseIdx));
        Valgrind::CheckDefined(fs.pressure(phaseIdx));

        using Evaluation = typename FluidState::ValueType;

        // extract variables
        const Evaluation& T = fs.temperature(phaseIdx);
        const Evaluation& p = fs.pressure(phaseIdx);
        const Evaluation& A = params.A(phaseIdx);
        const Evaluation& B = params.B(phaseIdx);
        const Evaluation& m1 = params.m1(phaseIdx);
        const Evaluation& m2 = params.m2(phaseIdx);

        // coefficients in cubic solver
        const Evaluation a1 = 1.0;  // cubic term
        const Evaluation a2 = (m1 + m2 - 1) * B - 1;  // quadratic term
        const Evaluation a3 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1);  // linear term
        const Evaluation a4 = -A * B - m1 * m2 * B * B * (B + 1);  // constant term
        Valgrind::CheckDefined(a1);
        Valgrind::CheckDefined(a2);
        Valgrind::CheckDefined(a3);
        Valgrind::CheckDefined(a4);

        // root of cubic equation
        Evaluation Vm = 0;
        Valgrind::SetUndefined(Vm);
        Evaluation Z[3] = {0.0, 0.0, 0.0};
        // Only the root of the phase is used: the largest for the vapour,
        // the smallest for the liquid.
        int numSol = cubicExtremeRoot(Z, a1, a2, a3, a4, isGasPhase);

        // pick correct root
        const Evaluation RT_p = R * T / p;
        if (numSol == 3) {
            // the EOS has three intersections with the pressure,
            // i.e. the molar volume of gas is the largest one and the
            // molar volume of liquid is the smallest one
            if (isGasPhase) {
                Vm = max(minMolarVolume, Z[2] * RT_p);
            } else {
                Vm = max(minMolarVolume, Z[0] * RT_p);
            }
        }
        else if (numSol == 1) {
            // Only one EOS root exists, so both phase labels use it.
            Vm = max(minMolarVolume, Z[0] * RT_p);
        }

        Valgrind::CheckDefined(Vm);
        assert(std::isfinite(scalarValue(Vm)));
        assert(Vm > 0);
        return Vm;

    }

private:
    //! The terms of \f$\ln \phi_i\f$ shared by all components of a phase.
    template <class LhsEval>
    struct PhaseFugacityTerms
    {
        LhsEval B; //!< \f$B\f$
        LhsEval Z; //!< compressibility factor of the cached molar volume
        LhsEval negLogZB; //!< \f$-\ln(Z - B)\f$
        LhsEval Zm1; //!< \f$Z - 1\f$
        LhsEval beta; //!< \f$\frac{A}{(m_1 - m_2) B} \ln\frac{Z + m_2 B}{Z + m_1 B}\f$
        LhsEval twoOverA; //!< \f$2 / A\f$
    };

    template <class LhsEval, class FluidState, class Params>
    static PhaseFugacityTerms<LhsEval> phaseFugacityTerms_(const FluidState& fs,
                                                           const Params& params,
                                                           unsigned phaseIdx)
    {
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.

        // extract variables
        LhsEval Vm = params.molarVolume(phaseIdx);
        LhsEval T = fs.temperature(phaseIdx);
        LhsEval p = fs.pressure(phaseIdx); // molar volume in [bar]
        LhsEval A = params.A(phaseIdx);
        LhsEval B = params.B(phaseIdx);
        LhsEval m1 = params.m1(phaseIdx);
        LhsEval m2 = params.m2(phaseIdx);

        // Calculate the compressibility factor
        LhsEval RT = R * T;
        LhsEval Z = p * Vm / RT; // compressibility factor

        PhaseFugacityTerms<LhsEval> terms;
        terms.negLogZB = -log(Z - B);
        terms.Zm1 = Z - 1;
        terms.beta = log((Z + m2 * B) / (Z + m1 * B)) * A / ((m1 - m2) * B);
        terms.twoOverA = 2 / A;
        terms.B = B;
        terms.Z = Z;
        return terms;
    }

    template <class LhsEval, class FluidState, class Params>
    static LhsEval componentFugacityCoefficient_(const PhaseFugacityTerms<LhsEval>& terms,
                                                 const FluidState& fs,
                                                 const Params& params,
                                                 unsigned phaseIdx,
                                                 unsigned compIdx)
    {
        LhsEval Bi = params.Bi(phaseIdx, compIdx);

        // Calculate Bi / B
        LhsEval Bi_B = Bi / terms.B;

        // calculate sum(A_ij * x_j)
        LhsEval A_s = 0.0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
            A_s += params.aCache(phaseIdx, compIdx, compJIdx) * fs.moleFraction(phaseIdx, compJIdx);
        }

        // calculate fugacity coefficient
        LhsEval alpha;
        LhsEval gamma;
        LhsEval ln_phi;
        LhsEval fugCoeff;

        alpha = terms.negLogZB + Bi_B * terms.Zm1;
        gamma = terms.twoOverA * A_s - Bi_B;
        ln_phi = alpha + (terms.beta * gamma);

        fugCoeff = exp(ln_phi);

        ////////
        // limit the fugacity coefficient to a reasonable range:
        //
        // on one side, we want the mole fraction to be at
        // least 10^-3 if the fugacity is at the current pressure
        //
        fugCoeff = min(maxFugacityCoefficient, fugCoeff);
        //
        // on the other hand, if the mole fraction of the component is 100%, we want the
        // fugacity to be at least 10^-3 Pa
        //
        fugCoeff = max(minFugacityCoefficient, fugCoeff);
        ///////////

        return fugCoeff;
    }

};  // class CubicEOS

}  // namespace Opm

#endif
