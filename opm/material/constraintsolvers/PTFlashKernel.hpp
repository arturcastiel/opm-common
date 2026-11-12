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
 * \brief The isothermal P-T flash of PTFlash::flash_solve_scalar() as a
 *        kernel: one free function per state on fixed-size arrays.
 *
 * The algorithm is that of PTFlash (flash first, the two-trial stability
 * test with its extrapolation, successive substitution, the Gibbs-energy
 * Newton and its polish step, the composition Newton as fallback, the
 * negative-flash, coincidence and lower-split checks, the single-phase
 * label), with the operations in the same order, so that both give the
 * same answers. The form is the one a GPU or a SIMD loop needs:
 *
 * - the fluid is a plain parameter struct (Parameters), the number of
 *   components a template parameter;
 * - the parameters that depend on pressure and temperature only are
 *   computed once per flash (Pure);
 * - every state lives in fixed-size arrays (Vec): no allocation, no
 *   virtual call, no fluid-system or fluid-state class, no exception;
 * - errors are status codes (Status), propagated where PTFlash throws and
 *   absorbed where PTFlash catches;
 * - every function is OPM_PTFLASH_KERNEL_FUNCTION, OPM_HOST_DEVICE by
 *   default.
 *
 * Values only: the derivative pass of PTFlash::solve() is not part of it.
 * PTFlashKernelHost.hpp builds the parameters from a fluid system.
 */
#ifndef OPM_PTFLASH_KERNEL_HPP
#define OPM_PTFLASH_KERNEL_HPP

#include <opm/common/utility/gpuDecorators.hpp>

#include <cmath>

#ifndef OPM_PTFLASH_KERNEL_FUNCTION
#define OPM_PTFLASH_KERNEL_FUNCTION OPM_HOST_DEVICE
#endif

namespace Opm::PTFlashKernel {

#define OPM_PTFK OPM_PTFLASH_KERNEL_FUNCTION inline

//! The math functions of the kernel; the device versions where compiled for one.
namespace math {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
// CUDA and HIP give these in the global namespace for device code.
#define OPM_PTFK_MATH(name) template <class T> OPM_PTFK T name(T x) { return ::name(x); }
#else
#define OPM_PTFK_MATH(name) template <class T> OPM_PTFK T name(T x) { using std::name; return name(x); }
#endif
OPM_PTFK_MATH(log)
OPM_PTFK_MATH(exp)
OPM_PTFK_MATH(sqrt)
OPM_PTFK_MATH(acos)
OPM_PTFK_MATH(cos)
OPM_PTFK_MATH(acosh)
OPM_PTFK_MATH(cosh)
OPM_PTFK_MATH(asinh)
OPM_PTFK_MATH(sinh)
#undef OPM_PTFK_MATH
template <class T> OPM_PTFK T abs(T x) { return x < T(0) ? -x : x; }
//! std::isfinite without the library: false for NaN and both infinities.
template <class T> OPM_PTFK bool isfinite(T x) { return x - x == T(0); }
//! std::min and std::max, which the cubic EOS uses, with their NaN behaviour.
template <class T> OPM_PTFK T min(T a, T b) { return (b < a) ? b : a; }
template <class T> OPM_PTFK T max(T a, T b) { return (a < b) ? b : a; }
template <class T> OPM_PTFK T infinity() { return static_cast<T>(HUGE_VAL); }
constexpr double pi = 3.141592653589793238462643383279502884;
} // namespace math

//! A fixed-size array that device code can copy and index.
template <class T, int n>
struct Vec
{
    T v[n];
    OPM_PTFK T& operator[](int i) { return v[i]; }
    OPM_PTFK const T& operator[](int i) const { return v[i]; }
};

template <class T, int n>
OPM_PTFK Vec<T, n> filled(T value)
{
    Vec<T, n> result;
    for (int i = 0; i < n; ++i) {
        result[i] = value;
    }
    return result;
}

//! The cubic equation of state.
enum class Eos : int { PR = 0, PRCORR = 1, SRK = 2, RK = 3 };

//! The two-phase method, as PTFlashMethod.
enum class Method : int { Ssi = 0, SsiAccelerated = 1, Newton = 2, SsiNewton = 3 };

/*!
 * \brief How a flash ended. Where PTFlash throws a NumericalProblem or a
 *        Dune::FMatrixError the kernel returns the matching code, and it
 *        continues where PTFlash catches one; cubicDegenerate is the
 *        std::runtime_error of the cubic solver, which no catch of PTFlash
 *        absorbs.
 */
enum class Status : int {
    ok = 0,
    nonFiniteParameters,      //!< pure or mixture EOS parameter not finite
    stabilityNotConverged,    //!< a stability trial at its iteration limit
    rachfordRiceNoRoot,       //!< all ratios on one side of one
    rachfordRiceNotConverged, //!< Rachford-Rice at its iteration limit
    substitutionNotConverged, //!< substitution at its iteration limit
    newtonTrivialStart,       //!< composition Newton started on coincident phases
    newtonNoProgress,         //!< composition Newton did not halve its residual
    newtonNotConverged,       //!< composition Newton at its iteration limit
    newtonNonPhysical,        //!< composition Newton converged outside the bounds
    singularJacobian,         //!< the composition Newton's Jacobian is singular
    splitNotConverged,        //!< the two-phase method did not converge
    splitNonPhysical,         //!< the two-phase answer is not physical
    warmRatiosNoSplit,        //!< warm-start ratios imply no split
    coincidentPhases,         //!< a cold split converged to coincident phases
    retryFailed,              //!< the split after a new stability test failed
    cubicDegenerate           //!< p = 0 in the cubic solver
};

//! Whether PTFlash's catch sites absorb an error of this kind.
OPM_PTFK bool recoverable(Status s)
{
    return s != Status::ok && s != Status::cubicDegenerate;
}

//! The work counters of one flash, as FlashStatistics counts them.
struct Counters
{
    int stabilityTests = 0;
    int stabilityIterations = 0;
    int ssiIterations = 0;
    int newtonIterations = 0;
    int rachfordRiceSolves = 0;
    int splitRetries = 0;
};

/*!
 * \brief The fluid: what GenericOilGasWaterFluidSystem and the EOS
 *        parameter classes give the flash, as plain numbers.
 */
template <class Scalar, int N>
struct Parameters
{
    Eos eos = Eos::PR;
    Vec<Scalar, N> criticalTemperature;  //!< [K]
    Vec<Scalar, N> criticalPressure;     //!< [Pa]
    Vec<Scalar, N> acentricFactor;       //!< [-]
    Vec<Scalar, N> criticalVolume;       //!< [m^3/kmol], as the fluid systems give it
    Vec<Scalar, N> volumeShift;          //!< SSHIFT s_i [-]
    //! OMEGAA over the equation's own Omega_a; one without OMEGAA.
    Vec<Scalar, N> omegaAScale;
    //! OMEGAB, or the equation's own Omega_b.
    Vec<Scalar, N> omegaB;
    Vec<Vec<Scalar, N>, N> interaction;  //!< k_ij, symmetric, zero diagonal
};

//! What the flash is asked to do.
template <class Scalar>
struct Settings
{
    Method method = Method::SsiNewton;
    Scalar tolerance = 1e-8;
};

//! The answer of one flash.
template <class Scalar, int N>
struct Answer
{
    Status status = Status::ok;
    bool singlePhase = false;
    Scalar L = 0;
    Vec<Scalar, N> x;  //!< liquid (oil) composition
    Vec<Scalar, N> y;  //!< vapour (gas) composition
    Vec<Scalar, N> K;  //!< the equilibrium ratios the flash ended with
    Counters counters;
};

//! Wilson's estimate of the equilibrium ratio of component i (wilsonK_() of the fluid states).
template <class Scalar, int N>
OPM_PTFK Scalar wilsonK(const Parameters<Scalar, N>& params, Scalar p, Scalar T, int i)
{
    const Scalar acf = params.acentricFactor[i];
    const Scalar Tc = params.criticalTemperature[i];
    const Scalar pc = params.criticalPressure[i];
    return math::exp(5.37 * (1 + acf) * (1 - Tc / T)) * (pc / p);
}

namespace detail {

template <class Scalar>
struct Constants
{
    static constexpr Scalar R = static_cast<Scalar>(8.314472);
    // PTFlash
    static constexpr Scalar stabilityTolerance = 1e-8;
    static constexpr Scalar substitutionTolerance = 1e-10;
    static constexpr Scalar trivialSolutionTolerance = 1e-5;
    static constexpr int substitutionAccelerationInterval = 5;
    static constexpr int substitutionStepsBeforeNewton = 5;
    static constexpr int maxHybridSubstitutionSteps = 20;
    static constexpr Scalar maxPredictedSubstitutionSteps = 8;
    static constexpr int stabilityAccelerationInterval = 5;
    static constexpr Scalar maxStabilityExtrapolation = 1.0;
    static constexpr int flashFirstIterations = 2;
    static constexpr int flashFirstMaxIterations = 3;
    static constexpr Scalar maxFugacityCoefficient = 1e10;
    static constexpr Scalar minFugacityCoefficient = 1e-10;
    static constexpr int newtonMaxIterations = 1000;
    static constexpr int newtonMaxStepsWithoutProgress = 10;
    static constexpr Scalar newtonProgressFactor = 0.5;
    static constexpr Scalar firstHessianShift = 1e-4;
    static constexpr int maxHessianShifts = 8;
    static constexpr int maxStabilityIterations = 20000;
    // CubicEOS
    static constexpr Scalar minMolarVolume = 1e-7;
};

constexpr int liquid = 0; // FluidSystem::oilPhaseIdx
constexpr int vapour = 1; // FluidSystem::gasPhaseIdx

template <class Scalar, int N>
using V = Vec<Scalar, N>;

template <class Scalar, int N>
using M = Vec<Vec<Scalar, N>, N>;

#define OPM_PTFK_TRY(expr)                                                                         \
    do {                                                                                           \
        const ::Opm::PTFlashKernel::Status opm_ptfk_status_ = (expr);                              \
        if (opm_ptfk_status_ != ::Opm::PTFlashKernel::Status::ok) {                                \
            return opm_ptfk_status_;                                                               \
        }                                                                                          \
    } while (false)

//! The EOS parameters of one (p, T): CubicEOSParams::updatePure, once per flash.
template <class Scalar, int N>
struct Pure
{
    Scalar p;
    Scalar T;
    Scalar m1;
    Scalar m2;
    V<Scalar, N> A;
    V<Scalar, N> B;
    V<Scalar, N> translation; //!< exp(-s_i B_i)
    M<Scalar, N> a;           //!< sqrt(A_i A_j) (1 - k_ij)
};

//! What the flash carries between its stages, as the scalar fluid state.
template <class Scalar, int N>
struct Phases
{
    V<Scalar, N> x[2];   //!< mole fractions, liquid and vapour
    V<Scalar, N> phi[2]; //!< fugacity coefficients of the last update
};

//! The flash state of PTFlash::flashScalar_(): fixed inputs and counters.
template <class Scalar, int N>
struct Context
{
    const Parameters<Scalar, N>& params;
    Pure<Scalar, N> pure;
    V<Scalar, N> incomingK; //!< the fluid state's K on entry
    Method method;
    Scalar tolerance;
    Counters& counters;
};

template <class Scalar>
OPM_PTFK Scalar constantOmegaA(Eos eos)
{
    return (eos == Eos::PR || eos == Eos::PRCORR) ? Scalar(0.457235529) : Scalar(0.4274802);
}

template <class Scalar>
OPM_PTFK Scalar constantOmegaB(Eos eos)
{
    return (eos == Eos::PR || eos == Eos::PRCORR) ? Scalar(0.077796074) : Scalar(0.08664035);
}

template <class Scalar, int N>
OPM_PTFK Scalar omegaA(const Parameters<Scalar, N>& params, Scalar T, int i)
{
    const Scalar Tr = T / params.criticalTemperature[i];
    const Scalar omega = params.acentricFactor[i];
    switch (params.eos) {
    case Eos::RK:
        return Scalar(0.4274802) / math::sqrt(Tr);
    case Eos::SRK: {
        const Scalar f_omega = 0.48 + omega * (1.574 + omega * (-0.176));
        const Scalar tmp = 1 + f_omega * (1 - math::sqrt(Tr));
        return Scalar(0.4274802) * tmp * tmp;
    }
    default: {
        Scalar f_omega;
        if (params.eos != Eos::PRCORR || omega <= 0.49) {
            f_omega = 0.37464 + omega * (1.54226 + omega * (-0.26992));
        }
        else {
            f_omega = 0.379642 + omega * (1.48503 + omega * (-0.164423 + omega * 0.016666));
        }
        const Scalar tmp = 1 + f_omega * (1 - math::sqrt(Tr));
        return Scalar(0.457235529) * tmp * tmp;
    }
    }
}

//! CubicEOSParams::updatePure() and updateACache_().
template <class Scalar, int N>
OPM_PTFK Status makePure(const Parameters<Scalar, N>& params, Scalar p, Scalar T, Pure<Scalar, N>& pure)
{
    pure.p = p;
    pure.T = T;
    if (params.eos == Eos::PR || params.eos == Eos::PRCORR) {
        pure.m1 = Scalar(1 + math::sqrt(2.0));
        pure.m2 = Scalar(1 - math::sqrt(2.0));
    }
    else {
        pure.m1 = Scalar(0.0);
        pure.m2 = Scalar(1.0);
    }
    for (int i = 0; i < N; ++i) {
        const Scalar pr = p / params.criticalPressure[i];
        const Scalar Tr = T / params.criticalTemperature[i];
        Scalar OmegaA = omegaA(params, T, i);
        OmegaA *= params.omegaAScale[i];
        const Scalar OmegaB = params.omegaB[i];
        const Scalar newA = OmegaA * pr / (Tr * Tr);
        const Scalar newB = OmegaB * pr / Tr;
        if (!math::isfinite(newA) || !math::isfinite(newB)) {
            return Status::nonFiniteParameters;
        }
        pure.A[i] = newA;
        pure.B[i] = newB;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            const Scalar Psi = i == j ? Scalar(0.0) : params.interaction[i][j];
            pure.a[i][j] = math::sqrt(pure.A[i] * pure.A[j]) * (1 - Psi);
            pure.a[j][i] = pure.a[i][j];
        }
    }
    for (int i = 0; i < N; ++i) {
        pure.translation[i] = math::exp(-params.volumeShift[i] * pure.B[i]);
    }
    return Status::ok;
}

//! Sorts three numbers ascending.
template <class Scalar>
OPM_PTFK void sort3(Scalar* s)
{
    if (s[1] < s[0]) { const Scalar t = s[0]; s[0] = s[1]; s[1] = t; }
    if (s[2] < s[1]) { const Scalar t = s[1]; s[1] = s[2]; s[2] = t; }
    if (s[1] < s[0]) { const Scalar t = s[0]; s[0] = s[1]; s[1] = t; }
}

/*!
 * \brief cubicExtremeRoot() for a = 1: the smallest (liquid) or largest
 *        (vapour) real root, at sol[0] or sol[numSol - 1].
 */
template <class Scalar>
OPM_PTFK Status cubicExtremeRoot(Scalar* sol, int& numSol, Scalar b, Scalar c, Scalar d, bool largest)
{
    const Scalar a = 1.0;
    const Scalar p = (3.0 * a * c - b * b) / (3.0 * a * a);
    const Scalar q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * d * a * a) / (27.0 * a * a * a);
    const Scalar discr = 4.0 * p * p * p + 27.0 * q * q;
    if (discr < 0.0) {
        const Scalar cos_3theta = ((3.0 * q) / (2.0 * p)) * math::sqrt(-3.0 / p);
        if (math::abs(cos_3theta) < 1.0) {
            // cubicExtremeOfThreeRoots_()
            const Scalar theta = (1.0 / 3.0) * math::acos(((3.0 * q) / (2.0 * p)) * math::sqrt(-3.0 / p));
            const Scalar radius = 2.0 * math::sqrt(-p / 3.0);
            const Scalar shift = b / (3.0 * a);
            constexpr double angleMargin = 1e-6;
            const double angle = theta;
            const bool separated = radius > 1e-6 * math::abs(shift);
            numSol = 3;
            if (separated && largest && angle < math::pi / 3.0 - angleMargin) {
                sol[2] = radius * math::cos(theta) - shift;
                return Status::ok;
            }
            if (separated && !largest && angle > angleMargin) {
                sol[0] = radius * math::cos(theta - ((4.0 * math::pi) / 3.0)) - shift;
                return Status::ok;
            }
            sol[0] = radius * math::cos(theta) - shift;
            sol[1] = radius * math::cos(theta - ((2.0 * math::pi) / 3.0)) - shift;
            sol[2] = radius * math::cos(theta - ((4.0 * math::pi) / 3.0)) - shift;
            sort3(sol);
            return Status::ok;
        }
    }
    else if (discr > 0.0) {
        if (p < 0) {
            const Scalar cosh_3theta = ((-3.0 * math::abs(q)) / (2.0 * p)) * math::sqrt(-3.0 / p);
            if (cosh_3theta > 1.0) {
                const Scalar theta = (1.0 / 3.0) * math::acosh(cosh_3theta);
                const Scalar t = ((-2.0 * math::abs(q)) / q) * math::sqrt(-p / 3.0) * math::cosh(theta);
                sol[0] = t - b / (3.0 * a);
                numSol = 1;
                return Status::ok;
            }
        }
        else if (p > 0) {
            const Scalar theta = (1.0 / 3.0) * math::asinh(((3.0 * q) / (2.0 * p)) * math::sqrt(3.0 / p));
            const Scalar t = -2.0 * math::sqrt(p / 3.0) * math::sinh(theta);
            sol[0] = t - b / (3.0 * a);
            numSol = 1;
            return Status::ok;
        }
        else {
            return Status::cubicDegenerate;
        }
    }
    numSol = 3;
    if (p == 0) {
        sol[0] = sol[1] = sol[2] = 0.0 - b / (3.0 * a);
    }
    else {
        sol[0] = (3.0 * q / p) - b / (3.0 * a);
        sol[1] = sol[2] = (-3.0 * q) / (2.0 * p) - b / (3.0 * a);
        sort3(sol);
    }
    return Status::ok;
}

//! The terms of ln phi a phase's components share (CubicEOS::PhaseFugacityTerms).
template <class Scalar>
struct PhaseTerms
{
    Scalar A;
    Scalar B;
    Scalar Vm;
    Scalar Z;
    Scalar negLogZB;
    Scalar Zm1;
    Scalar beta;
    Scalar twoOverA;
};

/*!
 * \brief The mixing rule, the root of the phase label and the shared terms:
 *        PTFlashParameterCache::updatePhase() and phaseFugacityTerms_().
 */
template <class Scalar, int N>
OPM_PTFK Status phaseTerms(const Pure<Scalar, N>& pure, const V<Scalar, N>& w, bool gas, PhaseTerms<Scalar>& t)
{
    // CubicEOSParams::updateMix(): mole fractions limited to [0, 1].
    V<Scalar, N> xc;
    for (int i = 0; i < N; ++i) {
        xc[i] = math::max(Scalar(0.0), math::min(Scalar(1.0), w[i]));
    }
    Scalar newA = 0;
    Scalar newB = 0;
    for (int i = 0; i < N; ++i) {
        const Scalar xi = xc[i];
        for (int j = 0; j < N; ++j) {
            newA += xi * xc[j] * pure.a[i][j];
        }
        newB += math::max(Scalar(0.0), xi) * pure.B[i];
    }
    if (!math::isfinite(newA) || !math::isfinite(newB)) {
        return Status::nonFiniteParameters;
    }
    const Scalar A = newA;
    const Scalar B = newB;
    const Scalar m1 = pure.m1;
    const Scalar m2 = pure.m2;

    // CubicEOS::computeMolarVolume()
    const Scalar a2 = (m1 + m2 - 1) * B - 1;
    const Scalar a3 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1);
    const Scalar a4 = -A * B - m1 * m2 * B * B * (B + 1);
    Scalar Zr[3] = {0.0, 0.0, 0.0};
    int numSol = 0;
    OPM_PTFK_TRY(cubicExtremeRoot(Zr, numSol, a2, a3, a4, gas));
    const Scalar RT_p = Constants<Scalar>::R * pure.T / pure.p;
    Scalar Vm = 0;
    if (numSol == 3) {
        Vm = math::max(Constants<Scalar>::minMolarVolume, (gas ? Zr[2] : Zr[0]) * RT_p);
    }
    else if (numSol == 1) {
        Vm = math::max(Constants<Scalar>::minMolarVolume, Zr[0] * RT_p);
    }

    // CubicEOS::phaseFugacityTerms_()
    const Scalar RT = Constants<Scalar>::R * pure.T;
    const Scalar Z = pure.p * Vm / RT;
    t.A = A;
    t.B = B;
    t.Vm = Vm;
    t.Z = Z;
    t.negLogZB = -math::log(Z - B);
    t.Zm1 = Z - 1;
    t.beta = math::log((Z + m2 * B) / (Z + m1 * B)) * A / ((m1 - m2) * B);
    t.twoOverA = 2 / A;
    return Status::ok;
}

//! CubicEOS::componentFugacityCoefficient_(), times the translation if asked.
template <class Scalar, int N>
OPM_PTFK void fugacityCoefficients(const Pure<Scalar, N>& pure,
                                   const PhaseTerms<Scalar>& t,
                                   const V<Scalar, N>& w,
                                   bool translate,
                                   V<Scalar, N>& phi)
{
    for (int i = 0; i < N; ++i) {
        const Scalar Bi_B = pure.B[i] / t.B;
        Scalar A_s = 0.0;
        for (int j = 0; j < N; ++j) {
            A_s += pure.a[i][j] * w[j];
        }
        const Scalar alpha = t.negLogZB + Bi_B * t.Zm1;
        const Scalar gamma = t.twoOverA * A_s - Bi_B;
        const Scalar ln_phi = alpha + (t.beta * gamma);
        Scalar fugCoeff = math::exp(ln_phi);
        fugCoeff = math::min(Constants<Scalar>::maxFugacityCoefficient, fugCoeff);
        fugCoeff = math::max(Constants<Scalar>::minFugacityCoefficient, fugCoeff);
        phi[i] = translate ? fugCoeff * pure.translation[i] : fugCoeff;
    }
}

/*!
 * \brief The fugacity coefficients of a phase (translated) and the
 *        derivatives of their logarithms in its mole fractions, as
 *        GenericOilGasWaterFluidSystem::fugacityCoefficientsAndDerivatives().
 */
template <class Scalar, int N>
OPM_PTFK Status fugacityCoefficientsAndDerivatives(const Pure<Scalar, N>& pure,
                                                   const V<Scalar, N>& w,
                                                   bool gas,
                                                   V<Scalar, N>& phi,
                                                   M<Scalar, N>& d)
{
    PhaseTerms<Scalar> t;
    OPM_PTFK_TRY(phaseTerms(pure, w, gas, t));
    V<Scalar, N> raw;
    fugacityCoefficients(pure, t, w, /*translate=*/false, raw);

    const Scalar A = t.A;
    const Scalar B = t.B;
    const Scalar m1 = pure.m1;
    const Scalar m2 = pure.m2;
    const Scalar Z = t.Z;
    V<Scalar, N> S;
    for (int i = 0; i < N; ++i) {
        S[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            S[i] += pure.a[i][j] * w[j];
        }
    }
    const Scalar c2 = (m1 + m2 - 1) * B - 1;
    const Scalar c1 = A + m1 * m2 * B * B - (m1 + m2) * B * (B + 1);
    const Scalar FZ = (3 * Z + 2 * c2) * Z + c1;
    const Scalar FA = Z - B;
    const Scalar FB = (m1 + m2 - 1) * Z * Z + (2 * m1 * m2 * B - (m1 + m2) * (2 * B + 1)) * Z
        - (A + m1 * m2 * (3 * B + 2) * B);
    const bool volumeAtFloor = !(t.Vm > Constants<Scalar>::minMolarVolume);

    const Scalar ZB = Z - B;
    const Scalar P1 = Z + m1 * B;
    const Scalar P2 = Z + m2 * B;
    const Scalar Lambda = math::log(P2 / P1);
    const Scalar Delta = m1 - m2;
    const Scalar beta = t.beta;
    const Scalar Zm1 = Z - 1;
    const Scalar invA = 1 / A;
    const Scalar invB = 1 / B;
    V<Scalar, N> gamma;
    for (int i = 0; i < N; ++i) {
        gamma[i] = 2 * invA * S[i] - pure.B[i] * invB;
    }
    for (int j = 0; j < N; ++j) {
        const Scalar Bj = pure.B[j];
        const Scalar Aj = 2 * S[j];
        const Scalar Zj = volumeAtFloor ? Scalar{0.0} : -(FA * Aj + FB * Bj) / FZ;
        const Scalar LambdaJ = (Zj + m2 * Bj) / P2 - (Zj + m1 * Bj) / P1;
        const Scalar betaJ = ((Aj * invB - A * Bj * invB * invB) * Lambda + A * invB * LambdaJ) / Delta;
        const Scalar common = -(Zj - Bj) / ZB;
        const Scalar perBi = Zj * invB - Bj * Zm1 * invB * invB + beta * Bj * invB * invB;
        const Scalar perS = -2 * beta * Aj * invA * invA;
        for (int i = 0; i < N; ++i) {
            d[i][j] = common + pure.B[i] * perBi + betaJ * gamma[i] + 2 * beta * invA * pure.a[i][j]
                + perS * S[i];
        }
    }
    for (int i = 0; i < N; ++i) {
        if (!(raw[i] < Constants<Scalar>::maxFugacityCoefficient)
            || !(raw[i] > Constants<Scalar>::minFugacityCoefficient)) {
            for (int j = 0; j < N; ++j) {
                d[i][j] = 0.0;
            }
        }
        phi[i] = raw[i] * pure.translation[i];
    }
    return Status::ok;
}

//! updatePhase() and updateFugacityCoefficients_() of one phase of a state.
template <class Scalar, int N>
OPM_PTFK Status updatePhase(const Pure<Scalar, N>& pure, Phases<Scalar, N>& s, int phase)
{
    PhaseTerms<Scalar> t;
    OPM_PTFK_TRY(phaseTerms(pure, s.x[phase], phase == vapour, t));
    fugacityCoefficients(pure, t, s.x[phase], /*translate=*/true, s.phi[phase]);
    return Status::ok;
}

//! The fugacity of a component in a phase of a state, p phi x.
template <class Scalar, int N>
OPM_PTFK Scalar fugacity(const Pure<Scalar, N>& pure, const Phases<Scalar, N>& s, int phase, int i)
{
    return pure.p * s.phi[phase][i] * s.x[phase][i];
}

template <class Scalar, int N>
OPM_PTFK Scalar dot(const V<Scalar, N>& a, const V<Scalar, N>& b)
{
    Scalar result(0);
    for (int i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template <class Scalar, int n>
OPM_PTFK Scalar twoNorm(const Vec<Scalar, n>& a)
{
    Scalar result(0);
    for (int i = 0; i < n; ++i) {
        result += a[i] * a[i];
    }
    return math::sqrt(result);
}

//! RachfordRice::solveInWindow()
template <class Scalar, int N>
OPM_PTFK Status rachfordRice(const V<Scalar, N>& K, const V<Scalar, N>& z, bool hasInitial, Scalar initialL, Scalar& L)
{
    constexpr Scalar tol = 1e-12;
    constexpr int max_iterations = 200;
    bool has_ratio_above_one = false;
    bool has_ratio_below_one = false;
    Scalar lower = 0.0;
    Scalar upper = 0.0;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.0)) {
            continue;
        }
        const Scalar dK = K[i] - 1.0;
        if (dK > 0.0) {
            const Scalar bound = (K[i] * z[i] - 1.0) / dK;
            lower = has_ratio_above_one && lower > bound ? lower : bound;
            has_ratio_above_one = true;
        }
        else if (dK < 0.0) {
            const Scalar bound = (1.0 - z[i]) / (-dK);
            upper = has_ratio_below_one && upper < bound ? upper : bound;
            has_ratio_below_one = true;
        }
    }
    if (!has_ratio_above_one || !has_ratio_below_one) {
        return Status::rachfordRiceNoRoot;
    }
    Scalar Vf = (lower + upper) / 2;
    if (hasInitial) {
        const Scalar initialV = 1.0 - initialL;
        if (initialV > lower && initialV < upper) {
            Vf = initialV;
        }
    }
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Scalar h = 0.0;
        Scalar dh = 0.0;
        for (int i = 0; i < N; ++i) {
            const Scalar dK = K[i] - 1.0;
            const Scalar denominator = 1.0 + Vf * dK;
            h += z[i] * dK / denominator;
            dh -= z[i] * dK * dK / (denominator * denominator);
        }
        if (h > 0.0) {
            lower = Vf;
        }
        else {
            upper = Vf;
        }
        Scalar next = Vf - h / dh;
        const Scalar abs_h = h < 0.0 ? -h : h;
        if (abs_h < tol) {
            L = 1.0 - next;
            return Status::ok;
        }
        if (!(upper - lower > 1e-15 * (1.0 + (Vf < 0.0 ? -Vf : Vf)))) {
            L = 1.0 - Vf;
            return Status::ok;
        }
        if (!(next > lower && next < upper)) {
            next = (lower + upper) / 2;
        }
        Vf = next;
    }
    return Status::rachfordRiceNotConverged;
}

//! PTFlash::trivialSolutionMeasure_() and isTrivialSolution_()
template <class Scalar, int N>
OPM_PTFK bool isTrivialSolution(const V<Scalar, N>& K)
{
    for (int i = 0; i < N; ++i) {
        if (!(K[i] >= 0.99 && K[i] <= 1.01)) {
            return false;
        }
    }
    Scalar measure = 0.;
    for (int i = 0; i < N; ++i) {
        const Scalar k = K[i];
        if (!(k > 0.) || !math::isfinite(k)) {
            return false; // an infinite measure
        }
        const Scalar log_ratio = math::log(k);
        measure += log_ratio * log_ratio;
    }
    return measure < Constants<Scalar>::trivialSolutionTolerance;
}

//! PTFlash::phasesCoincide_()
template <class Scalar, int N>
OPM_PTFK bool phasesCoincide(const Phases<Scalar, N>& s)
{
    V<Scalar, N> ratios = filled<Scalar, N>(1.);
    int active = 0;
    for (int i = 0; i < N; ++i) {
        const Scalar xl = s.x[liquid][i];
        const Scalar xv = s.x[vapour][i];
        if (!(xl > 0.) || !(xv > 0.)) {
            continue;
        }
        ++active;
        ratios[i] = xv / xl;
    }
    return active >= 2 && isTrivialSolution(ratios);
}

//! PTFlash::phaseCompositionsArePhysical_()
template <class Scalar, int N>
OPM_PTFK bool phaseCompositionsArePhysical(const V<Scalar, N>& xl, const V<Scalar, N>& xv)
{
    const Scalar slack = 1.e-8;
    const auto out_of_range = [slack](Scalar value) {
        return !math::isfinite(value) || value < -slack || value > Scalar{1} + slack;
    };
    for (int i = 0; i < N; ++i) {
        if (out_of_range(xl[i]) || out_of_range(xv[i])) {
            return false;
        }
    }
    return true;
}

//! PTFlash::computeLiquidVapor_()
template <class Scalar, int N>
OPM_PTFK void computeLiquidVapor(Phases<Scalar, N>& s, Scalar L, const V<Scalar, N>& K, const V<Scalar, N>& z)
{
    V<Scalar, N> x;
    V<Scalar, N> y;
    Scalar sumx = 0;
    Scalar sumy = 0;
    for (int i = 0; i < N; ++i) {
        x[i] = z[i] / (L + (1 - L) * K[i]);
        sumx += x[i];
        y[i] = (K[i] * z[i]) / (L + (1 - L) * K[i]);
        sumy += y[i];
    }
    for (int i = 0; i < N; ++i) {
        x[i] /= sumx;
        y[i] /= sumy;
    }
    s.x[liquid] = x;
    s.x[vapour] = y;
}

//! PTFlash::extrapolateSubstitution_()
template <class Scalar, int N>
OPM_PTFK bool extrapolateSubstitution(const V<Scalar, N> (&steps)[3], V<Scalar, N>& extrapolation)
{
    const V<Scalar, N>& newest = steps[0];
    const V<Scalar, N>& middle = steps[1];
    const V<Scalar, N>& oldest = steps[2];
    const Scalar mm = dot(middle, middle);
    const Scalar mo = dot(middle, oldest);
    const Scalar oo = dot(oldest, oldest);
    const Scalar determinant = mm * oo - mo * mo;
    bool found = false;
    if (determinant > 1e-10 * mm * oo) {
        const Scalar b1 = (-dot(middle, newest) * oo + mo * dot(oldest, newest)) / determinant;
        const Scalar b2 = (-mm * dot(oldest, newest) + mo * dot(middle, newest)) / determinant;
        const Scalar denominator = 1. + b1 + b2;
        if (denominator > 1e-12) {
            for (int i = 0; i < N; ++i) {
                extrapolation[i] = -((b1 + b2) * newest[i] + b2 * middle[i]) / denominator;
            }
            found = true;
        }
    }
    if (!found) {
        const Scalar eigenvalue = dot(newest, newest) / dot(middle, newest);
        if (!(eigenvalue > 0. && eigenvalue < 1.)) {
            return false;
        }
        for (int i = 0; i < N; ++i) {
            extrapolation[i] = newest[i] * eigenvalue / (1. - eigenvalue);
        }
    }
    for (int i = 0; i < N; ++i) {
        if (!math::isfinite(extrapolation[i])) {
            return false;
        }
    }
    return true;
}

//! The stable cubic root of a feed (PTFlash::FeedRoot).
template <class Scalar>
struct FeedRoot
{
    int phase;
    bool hasTwoRoots;
    Scalar molarVolume;
};

/*!
 * \brief PTFlash::feedStableRoot_(): the root of the feed with the lower
 *        Gibbs energy and the untranslated fugacity coefficients on it.
 */
template <class Scalar, int N>
OPM_PTFK Status feedStableRoot(const Pure<Scalar, N>& pure,
                               const V<Scalar, N>& z,
                               FeedRoot<Scalar>& root,
                               V<Scalar, N>& stable)
{
    PhaseTerms<Scalar> tl;
    PhaseTerms<Scalar> tv;
    OPM_PTFK_TRY(phaseTerms(pure, z, /*gas=*/false, tl));
    OPM_PTFK_TRY(phaseTerms(pure, z, /*gas=*/true, tv));
    const Scalar liquid_volume = tl.Vm;
    const Scalar vapour_volume = tv.Vm;
    if (!(liquid_volume < (1.0 - 1e-6) * vapour_volume)) {
        fugacityCoefficients(pure, tl, z, /*translate=*/false, stable);
        root = {liquid, false, vapour_volume};
        return Status::ok;
    }
    V<Scalar, N> lc;
    fugacityCoefficients(pure, tv, z, /*translate=*/false, stable);
    fugacityCoefficients(pure, tl, z, /*translate=*/false, lc);
    Scalar sum_v = 0.0;
    Scalar sum_l = 0.0;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            continue;
        }
        sum_v += z[i] * math::log(stable[i]);
    }
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            continue;
        }
        sum_l += z[i] * math::log(lc[i]);
    }
    const bool vapour_is_stable = sum_v < sum_l;
    if (!vapour_is_stable) {
        stable = lc;
    }
    root = {vapour_is_stable ? vapour : liquid, true, vapour_is_stable ? vapour_volume : liquid_volume};
    return Status::ok;
}

/*!
 * \brief PTFlash::checkStability_(): one trial phase of Michelsen's test,
 *        with the extrapolation every stabilityAccelerationInterval steps.
 */
template <class Scalar, int N>
OPM_PTFK Status checkStability(const Pure<Scalar, N>& pure,
                               Counters& counters,
                               bool& isTrivial,
                               V<Scalar, N>& K,
                               V<Scalar, N>& trial_composition,
                               Scalar& trial_sum,
                               const V<Scalar, N>& z,
                               bool is_vapour_trial,
                               const V<Scalar, N>& reference_coefficients)
{
    using C = Constants<Scalar>;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            K[i] = 1.;
        }
    }
    V<Scalar, N> steps[3] = {filled<Scalar, N>(0.), filled<Scalar, N>(0.), filled<Scalar, N>(0.)};
    int plain_steps = 0;
    bool extrapolated = false;
    V<Scalar, N> plain_K = K;
    Scalar tm_before = 0.;
    V<Scalar, N> reference_fugacity;
    for (int i = 0; i < N; ++i) {
        reference_fugacity[i] = pure.p * reference_coefficients[i] * z[i];
    }

    for (int iteration = 0; iteration < C::maxStabilityIterations; ++iteration) {
        ++counters.stabilityIterations;
        for (int i = 0; i < N; ++i) {
            trial_composition[i] = is_vapour_trial ? K[i] * z[i] : z[i] / K[i];
        }
        trial_sum = Scalar(0.0);
        for (int i = 0; i < N; ++i) {
            trial_sum = trial_sum + trial_composition[i];
        }
        for (int i = 0; i < N; ++i) {
            trial_composition[i] /= trial_sum;
        }
        PhaseTerms<Scalar> t;
        OPM_PTFK_TRY(phaseTerms(pure, trial_composition, is_vapour_trial, t));
        V<Scalar, N> trial_phi;
        fugacityCoefficients(pure, t, trial_composition, /*translate=*/false, trial_phi);

        V<Scalar, N> R;
        for (int i = 0; i < N; ++i) {
            if (!(z[i] > 0.)) {
                R[i] = 1.;
                continue;
            }
            const Scalar trial_fugacity = pure.p * trial_phi[i] * trial_composition[i];
            R[i] = is_vapour_trial ? reference_fugacity[i] / trial_fugacity / trial_sum
                                   : trial_fugacity / reference_fugacity[i] * trial_sum;
        }

        V<Scalar, N> step = filled<Scalar, N>(0.);
        Scalar tm = 1. - trial_sum;
        if (extrapolated || plain_steps + 3 >= C::stabilityAccelerationInterval) {
            for (int i = 0; i < N; ++i) {
                if (!(z[i] > 0.)) {
                    continue;
                }
                step[i] = math::log(R[i]);
                const Scalar mole_number = trial_sum * trial_composition[i];
                tm += is_vapour_trial ? -mole_number * step[i] : mole_number * step[i];
            }
        }
        if (extrapolated) {
            extrapolated = false;
            if (!(tm < tm_before) || !math::isfinite(tm)) {
                K = plain_K;
                plain_steps = 0;
                continue;
            }
        }

        for (int i = 0; i < N; ++i) {
            K[i] *= R[i];
        }
        Scalar substitution_residual = 0.0;
        for (int i = 0; i < N; ++i) {
            if (!(z[i] > 0.)) {
                continue;
            }
            const Scalar e = R[i] - 1.0;
            substitution_residual += e * e;
        }

        isTrivial = isTrivialSolution(K);
        if (isTrivial || substitution_residual < C::substitutionTolerance) {
            return Status::ok;
        }

        steps[2] = steps[1];
        steps[1] = steps[0];
        steps[0] = step;
        ++plain_steps;
        V<Scalar, N> extrapolation;
        if (plain_steps >= C::stabilityAccelerationInterval && extrapolateSubstitution(steps, extrapolation)) {
            plain_K = K;
            tm_before = tm;
            // Dune's infinity_norm(); a zero or NaN norm leaves the step unshortened, as there.
            Scalar largest = math::abs(extrapolation[0]);
            for (int i = 1; i < N; ++i) {
                largest = math::max(math::abs(extrapolation[i]), largest);
            }
            const Scalar shortening = largest > C::maxStabilityExtrapolation
                ? C::maxStabilityExtrapolation / largest : Scalar{1.};
            for (int i = 0; i < N; ++i) {
                K[i] *= math::exp(shortening * extrapolation[i]);
            }
            extrapolated = true;
            plain_steps = 0;
        }
    }
    return Status::stabilityNotConverged;
}

/*!
 * \brief PTFlash::phaseStabilityTest_(): both trials against the feed on
 *        its stable root; for an unstable feed, K from the trials.
 */
template <class Scalar, int N>
OPM_PTFK Status phaseStabilityTest(const Pure<Scalar, N>& pure,
                                   Counters& counters,
                                   bool& isStable,
                                   V<Scalar, N>& K,
                                   Phases<Scalar, N>& s,
                                   const V<Scalar, N>& z,
                                   FeedRoot<Scalar>& feed_root)
{
    using C = Constants<Scalar>;
    ++counters.stabilityTests;
    bool liquid_trivial = false;
    bool vapour_trivial = false;
    V<Scalar, N> liquid_composition;
    V<Scalar, N> vapour_composition;
    Scalar liquid_sum = 0;
    Scalar vapour_sum = 0;
    V<Scalar, N> vapour_k = K;
    V<Scalar, N> liquid_k = K;
    V<Scalar, N> reference;
    OPM_PTFK_TRY(feedStableRoot(pure, z, feed_root, reference));
    OPM_PTFK_TRY(checkStability(pure, counters, vapour_trivial, vapour_k, vapour_composition, vapour_sum, z,
                                /*is_vapour_trial=*/true, reference));
    const bool vapour_no_instability = (vapour_sum < (1.0 + C::stabilityTolerance)) || vapour_trivial;
    OPM_PTFK_TRY(checkStability(pure, counters, liquid_trivial, liquid_k, liquid_composition, liquid_sum, z,
                                /*is_vapour_trial=*/false, reference));
    const bool liquid_no_instability = (liquid_sum < (1.0 + C::stabilityTolerance)) || liquid_trivial;

    isStable = liquid_no_instability && vapour_no_instability;
    if (isStable) {
        s.x[vapour] = z;
        s.x[liquid] = z;
    }
    else {
        const V<Scalar, N> starting_k = K;
        for (int i = 0; i < N; ++i) {
            K[i] = z[i] > 0. ? vapour_composition[i] / liquid_composition[i] : Scalar(1.);
        }
        if (isTrivialSolution(K)) {
            const V<Scalar, N>& trial = vapour_sum > liquid_sum ? vapour_composition : liquid_composition;
            Scalar lightness = 0.;
            for (int i = 0; i < N; ++i) {
                const Scalar ratio = starting_k[i];
                if (z[i] > 0. && ratio > 0.) {
                    lightness += (trial[i] - z[i]) * math::log(ratio);
                }
            }
            for (int i = 0; i < N; ++i) {
                if (z[i] > 0.) {
                    K[i] = lightness > 0. ? trial[i] / z[i] : z[i] / trial[i];
                }
            }
        }
    }
    return Status::ok;
}

//! PTFlash::splitGibbsEnergy_(fluid_state, L, z): from the stored coefficients.
template <class Scalar, int N>
OPM_PTFK Scalar splitGibbsEnergy(const Phases<Scalar, N>& s, Scalar L, const V<Scalar, N>& z)
{
    Scalar gibbs = 0.;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            continue;
        }
        const Scalar x = s.x[liquid][i];
        const Scalar y = s.x[vapour][i];
        gibbs += L * x * math::log(x * s.phi[liquid][i]) + (1. - L) * y * math::log(y * s.phi[vapour][i]);
    }
    return gibbs;
}

//! PTFlash::splitGibbsEnergy_(fluid_state, L, eos_type): updates both phases first.
template <class Scalar, int N>
OPM_PTFK Status splitGibbsEnergyUpdated(const Pure<Scalar, N>& pure, Phases<Scalar, N>& s, Scalar L, Scalar& energy)
{
    energy = 0.0;
    for (int phase = 0; phase < 2; ++phase) {
        OPM_PTFK_TRY(updatePhase(pure, s, phase));
        const Scalar amount = phase == liquid ? L : 1.0 - L;
        for (int i = 0; i < N; ++i) {
            const Scalar w = s.x[phase][i];
            if (w > 0.) {
                energy += amount * w * math::log(w * s.phi[phase][i]);
            }
        }
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::successiveSubstitutionComposition_(); \p converged false
 *        is the hand-over to Newton of the hybrid method.
 */
template <class Scalar, int N>
OPM_PTFK Status successiveSubstitution(Context<Scalar, N>& ctx,
                                       V<Scalar, N>& K,
                                       Scalar& L,
                                       Phases<Scalar, N>& s,
                                       const V<Scalar, N>& z,
                                       bool run_newton_afterwards,
                                       bool extrapolate,
                                       bool& converged)
{
    using C = Constants<Scalar>;
    const Scalar flash_tolerance = ctx.tolerance;
    const int max_iterations
        = run_newton_afterwards ? C::maxHybridSubstitutionSteps : (extrapolate ? 1000 : 100);
    converged = false;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            K[i] = 1.;
        }
    }
    V<Scalar, N> steps[3] = {filled<Scalar, N>(0.), filled<Scalar, N>(0.), filled<Scalar, N>(0.)};
    int plain_steps = 0;
    bool extrapolated = false;
    V<Scalar, N> plain_K = K;
    Scalar plain_L = L;
    Scalar gibbs_before = 0.;
    Scalar previous_residual = 0.;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        ++ctx.counters.ssiIterations;
        computeLiquidVapor(s, L, K, z);
        for (int phase = 0; phase < 2; ++phase) {
            OPM_PTFK_TRY(updatePhase(ctx.pure, s, phase));
        }
        V<Scalar, N> fugacity_ratio;
        V<Scalar, N> fugacity_residual;
        for (int i = 0; i < N; ++i) {
            if (!(z[i] > 0.)) {
                fugacity_ratio[i] = 1.;
                fugacity_residual[i] = 0.;
                continue;
            }
            fugacity_ratio[i] = fugacity(ctx.pure, s, liquid, i) / fugacity(ctx.pure, s, vapour, i);
            fugacity_residual[i] = fugacity_ratio[i] - 1.0;
        }
        const bool split_is_interior = L > 0. && L < 1.;
        const Scalar gibbs = extrapolate && split_is_interior ? splitGibbsEnergy(s, L, z) : Scalar{0.};
        if (extrapolated) {
            extrapolated = false;
            if (!(split_is_interior && gibbs < gibbs_before)) {
                K = plain_K;
                L = plain_L;
                plain_steps = 0;
                continue;
            }
        }
        const Scalar residual = twoNorm(fugacity_residual);
        if (residual < flash_tolerance) {
            converged = true;
            return Status::ok;
        }
        V<Scalar, N> step = filled<Scalar, N>(0.);
        for (int i = 0; i < N; ++i) {
            if (extrapolate) {
                step[i] = math::log(fugacity_ratio[i]);
            }
            K[i] *= fugacity_ratio[i];
        }
        OPM_PTFK_TRY(rachfordRice(K, z, true, L, L));

        if (run_newton_afterwards) {
            if (iteration + 1 >= C::substitutionStepsBeforeNewton) {
                const Scalar rate = iteration > 0 ? residual / previous_residual : Scalar{1.};
                const Scalar predicted_steps = rate > 0. && rate < 1.
                    ? math::log(flash_tolerance / residual) / math::log(rate)
                    : math::infinity<Scalar>();
                if (!(predicted_steps <= C::maxPredictedSubstitutionSteps)) {
                    break;
                }
            }
            previous_residual = residual;
        }

        steps[2] = steps[1];
        steps[1] = steps[0];
        steps[0] = step;
        ++plain_steps;
        V<Scalar, N> extrapolation;
        if (extrapolate && split_is_interior && plain_steps >= C::substitutionAccelerationInterval
            && extrapolateSubstitution(steps, extrapolation)) {
            plain_K = K;
            plain_L = L;
            gibbs_before = gibbs;
            for (int i = 0; i < N; ++i) {
                K[i] *= math::exp(extrapolation[i]);
            }
            Scalar extrapolated_L = L;
            const Status rr = rachfordRice(K, z, true, L, extrapolated_L);
            if (rr == Status::ok) {
                L = extrapolated_L;
                extrapolated = true;
            }
            else {
                K = plain_K;
                L = plain_L;
            }
            plain_steps = 0;
        }
    }
    if (!run_newton_afterwards) {
        return Status::substitutionNotConverged;
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::solveModifiedCholesky_(): (H + s diag(d)) step = rhs on
 *        the leading n x n block with the smallest shift that factorises.
 */
template <class Scalar, int N>
OPM_PTFK bool solveModifiedCholesky(const M<Scalar, N>& hessian,
                                    const V<Scalar, N>& diagonal,
                                    const V<Scalar, N>& rhs,
                                    int n,
                                    V<Scalar, N>& solution)
{
    using C = Constants<Scalar>;
    M<Scalar, N> factor;
    for (int i = 0; i < N; ++i) {
        factor[i] = filled<Scalar, N>(0.);
    }
    Scalar shift = 0.;
    for (int attempt = 0; attempt < C::maxHessianShifts; ++attempt) {
        bool positive_definite = true;
        for (int j = 0; j < n && positive_definite; ++j) {
            Scalar pivot = hessian[j][j] + shift * diagonal[j];
            for (int k = 0; k < j; ++k) {
                pivot -= factor[j][k] * factor[j][k];
            }
            if (!(pivot > 1e-12 * diagonal[j])) {
                positive_definite = false;
                break;
            }
            factor[j][j] = math::sqrt(pivot);
            for (int i = j + 1; i < n; ++i) {
                Scalar entry = hessian[i][j];
                for (int k = 0; k < j; ++k) {
                    entry -= factor[i][k] * factor[j][k];
                }
                factor[i][j] = entry / factor[j][j];
            }
        }
        if (positive_definite) {
            for (int i = 0; i < n; ++i) {
                Scalar entry = rhs[i];
                for (int k = 0; k < i; ++k) {
                    entry -= factor[i][k] * solution[k];
                }
                solution[i] = entry / factor[i][i];
            }
            for (int i = n; i-- > 0;) {
                Scalar entry = solution[i];
                for (int k = i + 1; k < n; ++k) {
                    entry -= factor[k][i] * solution[k];
                }
                solution[i] = entry / factor[i][i];
            }
            return true;
        }
        shift = attempt == 0 ? C::firstHessianShift : 10. * shift;
    }
    return false;
}

/*!
 * \brief PTFlash::gibbsNewtonComposition_(): Newton on the Gibbs energy in
 *        the vapour mole numbers, with the final polish step.
 *
 * \param gaveUp set where Newton gives up (PTFlash's std::nullopt)
 */
template <class Scalar, int N>
OPM_PTFK Status gibbsNewton(Context<Scalar, N>& ctx,
                            V<Scalar, N>& K,
                            Scalar& L,
                            Phases<Scalar, N>& s,
                            const V<Scalar, N>& z,
                            bool& gaveUp)
{
    using C = Constants<Scalar>;
    gaveUp = true;
    int present[N];
    int num_present = 0;
    V<Scalar, N> vap = filled<Scalar, N>(0.);
    V<Scalar, N> feed = filled<Scalar, N>(0.);
    if (!(L > 0. && L < 1.)) {
        return Status::ok;
    }
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            continue;
        }
        const Scalar moles = z[i] * (1. - L) * K[i] / (L + (1. - L) * K[i]);
        if (!(moles > 0. && moles < z[i])) {
            return Status::ok;
        }
        present[num_present] = i;
        vap[num_present] = moles;
        feed[num_present] = z[i];
        ++num_present;
    }
    if (num_present < 2) {
        return Status::ok;
    }

    constexpr int num_equations = 2 * N + 1;
    Vec<Scalar, num_equations> residual;
    V<Scalar, N> liquid_composition = filled<Scalar, N>(0.);
    V<Scalar, N> vapour_composition = filled<Scalar, N>(0.);
    Scalar liquid_moles = 0.;
    Scalar vapour_moles = 0.;
    Scalar liquid_fraction = 0.;
    const auto set_compositions = [&]() {
        liquid_moles = 0.;
        vapour_moles = 0.;
        for (int a = 0; a < num_present; ++a) {
            vapour_moles += vap[a];
            liquid_moles += feed[a] - vap[a];
        }
        for (int a = 0; a < num_present; ++a) {
            liquid_composition[present[a]] = (feed[a] - vap[a]) / liquid_moles;
            vapour_composition[present[a]] = vap[a] / vapour_moles;
        }
        s.x[liquid] = liquid_composition;
        s.x[vapour] = vapour_composition;
        liquid_fraction = liquid_moles / (liquid_moles + vapour_moles);
    };
    set_compositions();
    if (phasesCoincide(s)) {
        return Status::ok;
    }

    const Scalar pressure = ctx.pure.p;
    M<Scalar, N> derivatives[2];
    const auto evaluate = [&](Scalar& norm) -> Status {
        ++ctx.counters.newtonIterations;
        set_compositions();
        for (int phase = 0; phase < 2; ++phase) {
            OPM_PTFK_TRY(fugacityCoefficientsAndDerivatives(ctx.pure, s.x[phase], phase == vapour, s.phi[phase],
                                                            derivatives[phase]));
        }
        Scalar liquid_sum = 0.;
        Scalar vapour_sum = 0.;
        for (int i = 0; i < N; ++i) {
            residual[i] = -z[i] + liquid_fraction * liquid_composition[i]
                + (1 - liquid_fraction) * vapour_composition[i];
            residual[i + N] = (fugacity(ctx.pure, s, liquid, i) - fugacity(ctx.pure, s, vapour, i)) / pressure;
            liquid_sum += liquid_composition[i];
            vapour_sum += vapour_composition[i];
        }
        residual[num_equations - 1] = liquid_sum - vapour_sum;
        const Scalar n2 = twoNorm(residual);
        norm = math::isfinite(n2) ? n2 : math::infinity<Scalar>();
        return Status::ok;
    };
    const auto newton_step = [&](V<Scalar, N>& step) {
        V<Scalar, N> gradient = filled<Scalar, N>(0.);
        V<Scalar, N> ideal_diagonal = filled<Scalar, N>(0.);
        M<Scalar, N> hessian;
        for (int i = 0; i < N; ++i) {
            hessian[i] = filled<Scalar, N>(0.);
        }
        for (int phase = 0; phase < 2; ++phase) {
            const V<Scalar, N>& composition = phase == liquid ? liquid_composition : vapour_composition;
            const M<Scalar, N>& d = derivatives[phase];
            const Scalar moles = phase == liquid ? liquid_moles : vapour_moles;
            for (int a = 0; a < num_present; ++a) {
                const int i = present[a];
                Scalar weighted = 0.;
                for (int b = 0; b < num_present; ++b) {
                    weighted += composition[present[b]] * d[i][present[b]];
                }
                for (int b = 0; b < num_present; ++b) {
                    hessian[a][b] += (d[i][present[b]] - weighted - 1.) / moles;
                }
                const Scalar ideal = 1. / (moles * composition[i]);
                hessian[a][a] += ideal;
                ideal_diagonal[a] += ideal;
            }
        }
        for (int a = 0; a < num_present; ++a) {
            const int i = present[a];
            gradient[a] = math::log(fugacity(ctx.pure, s, liquid, i) / fugacity(ctx.pure, s, vapour, i));
        }
        return solveModifiedCholesky(hessian, ideal_diagonal, gradient, num_present, step);
    };

    bool converged = false;
    int iteration = 0;
    int steps_without_progress = 0;
    Scalar residual_norm = 0;
    OPM_PTFK_TRY(evaluate(residual_norm));
    Scalar progress_residual = residual_norm;
    while (iteration < C::newtonMaxIterations) {
        converged = residual_norm < ctx.tolerance;
        if (converged) {
            break;
        }
        if (residual_norm < C::newtonProgressFactor * progress_residual) {
            progress_residual = residual_norm;
            steps_without_progress = 0;
        }
        else if (++steps_without_progress >= C::newtonMaxStepsWithoutProgress
                 || !math::isfinite(residual_norm)) {
            return Status::ok;
        }
        V<Scalar, N> step = filled<Scalar, N>(0.);
        if (!newton_step(step)) {
            return Status::ok;
        }
        for (int a = 0; a < num_present; ++a) {
            vap[a] += step[a];
            if (!(vap[a] > 0. && vap[a] < feed[a])) {
                return Status::ok;
            }
        }
        OPM_PTFK_TRY(evaluate(residual_norm));
        ++iteration;
    }
    if (!converged) {
        return Status::ok;
    }
    V<Scalar, N> last_step = filled<Scalar, N>(0.);
    if (newton_step(last_step)) {
        bool inside = true;
        for (int a = 0; a < num_present; ++a) {
            inside = inside && vap[a] + last_step[a] > 0. && vap[a] + last_step[a] < feed[a];
        }
        if (inside) {
            for (int a = 0; a < num_present; ++a) {
                vap[a] += last_step[a];
            }
            set_compositions();
        }
    }
    for (int i = 0; i < N; ++i) {
        K[i] = vapour_composition[i] / liquid_composition[i];
    }
    L = liquid_fraction;
    gaveUp = false;
    return Status::ok;
}

/*!
 * \brief Solves J dx = r in place of r by LU with partial pivoting, as
 *        Dune::DenseMatrix::solve() for n > 3. The factors overwrite A.
 */
template <class Scalar, int n>
OPM_PTFK Status luSolve(Vec<Vec<Scalar, n>, n>& A, Vec<Scalar, n>& rhs)
{
    for (int i = 0; i < n; ++i) {
        Scalar pivmax = math::abs(A[i][i]);
        int imax = i;
        for (int k = i + 1; k < n; ++k) {
            const Scalar a = math::abs(A[k][i]);
            if (a > pivmax) {
                pivmax = a;
                imax = k;
            }
        }
        if (imax != i) {
            const Vec<Scalar, n> row = A[i];
            A[i] = A[imax];
            A[imax] = row;
        }
        const Scalar t = rhs[i];
        rhs[i] = rhs[imax];
        rhs[imax] = t;
        if (!(pivmax != Scalar(0))) {
            return Status::singularJacobian;
        }
        for (int k = i + 1; k < n; ++k) {
            const Scalar factor = A[k][i] / A[i][i];
            A[k][i] = factor;
            for (int j = i + 1; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            rhs[k] -= factor * rhs[i];
        }
    }
    Vec<Scalar, n>& x = rhs;
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; ++j) {
            rhs[i] -= A[i][j] * x[j];
        }
        x[i] = rhs[i] / A[i][i];
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::newtonComposition_(): the Gibbs-energy Newton, and where
 *        it gives up Newton in (x, y, L). The Jacobian of the latter is the
 *        analytic one of PTFlash::analyticCompositionJacobian_() where
 *        PTFlash differentiates automatically: the same matrix to rounding.
 */
template <class Scalar, int N>
OPM_PTFK Status newtonComposition(Context<Scalar, N>& ctx,
                                  V<Scalar, N>& K,
                                  Scalar& L,
                                  Phases<Scalar, N>& s,
                                  const V<Scalar, N>& z)
{
    using C = Constants<Scalar>;
    bool gaveUp = true;
    OPM_PTFK_TRY(gibbsNewton(ctx, K, L, s, z, gaveUp));
    if (!gaveUp) {
        return Status::ok;
    }

    constexpr int n = 2 * N + 1;
    computeLiquidVapor(s, L, K, z);
    V<Scalar, N> lc = s.x[liquid];
    V<Scalar, N> vc = s.x[vapour];
    Scalar lf = L;
    if (phasesCoincide(s)) {
        return Status::newtonTrivialStart;
    }
    const Scalar pref = ctx.pure.p;
    Vec<Scalar, n> residual;
    Vec<Vec<Scalar, n>, n> jacobian;
    const auto evaluate = [&](Scalar& norm) -> Status {
        ++ctx.counters.newtonIterations;
        s.x[liquid] = lc;
        s.x[vapour] = vc;
        M<Scalar, N> d[2];
        for (int phase = 0; phase < 2; ++phase) {
            OPM_PTFK_TRY(fugacityCoefficientsAndDerivatives(ctx.pure, s.x[phase], phase == vapour, s.phi[phase],
                                                            d[phase]));
        }
        for (int r = 0; r < n; ++r) {
            jacobian[r] = filled<Scalar, n>(0.);
        }
        for (int i = 0; i < N; ++i) {
            residual[i] = -z[i] + lf * lc[i] + (1 - lf) * vc[i];
            residual[i + N] = (fugacity(ctx.pure, s, liquid, i) - fugacity(ctx.pure, s, vapour, i)) / pref;
            jacobian[i][i] = lf;
            jacobian[i][N + i] = 1. - lf;
            jacobian[i][n - 1] = lc[i] - vc[i];
            jacobian[n - 1][i] = 1.;
            jacobian[n - 1][N + i] = -1.;
            for (int phase = 0; phase < 2; ++phase) {
                const Scalar scaled = s.phi[phase][i];
                const Scalar w = s.x[phase][i];
                const Scalar sign = phase == liquid ? Scalar(1) : Scalar(-1);
                for (int j = 0; j < N; ++j) {
                    jacobian[N + i][phase * N + j] = sign * (scaled * w * d[phase][i][j]);
                }
                jacobian[N + i][phase * N + i] += sign * scaled;
            }
        }
        Scalar liquid_sum = 0.;
        Scalar vapour_sum = 0.;
        for (int i = 0; i < N; ++i) {
            liquid_sum += lc[i];
            vapour_sum += vc[i];
        }
        residual[n - 1] = liquid_sum - vapour_sum;
        const Scalar n2 = twoNorm(residual);
        norm = math::isfinite(n2) ? n2 : math::infinity<Scalar>();
        return Status::ok;
    };

    bool converged = false;
    int iteration = 0;
    int steps_without_progress = 0;
    Scalar residual_norm = 0;
    OPM_PTFK_TRY(evaluate(residual_norm));
    Scalar progress_residual = residual_norm;
    while (iteration < C::newtonMaxIterations) {
        converged = residual_norm < ctx.tolerance;
        if (converged) {
            break;
        }
        if (residual_norm < C::newtonProgressFactor * progress_residual) {
            progress_residual = residual_norm;
            steps_without_progress = 0;
        }
        else if (++steps_without_progress >= C::newtonMaxStepsWithoutProgress
                 || !math::isfinite(residual_norm)) {
            return Status::newtonNoProgress;
        }
        Vec<Scalar, n> step = residual;
        OPM_PTFK_TRY(luSolve(jacobian, step));
        for (int i = 0; i < N; ++i) {
            lc[i] -= step[i];
            vc[i] -= step[i + N];
        }
        lf -= step[n - 1];
        OPM_PTFK_TRY(evaluate(residual_norm));
        ++iteration;
    }
    if (!converged) {
        return Status::newtonNotConverged;
    }
    if (!math::isfinite(lf) || !phaseCompositionsArePhysical(lc, vc)) {
        return Status::newtonNonPhysical;
    }
    s.x[liquid] = lc;
    s.x[vapour] = vc;
    for (int i = 0; i < N; ++i) {
        K[i] = vc[i] / lc[i];
    }
    L = lf;
    return Status::ok;
}

//! PTFlash::flash_2ph(): the two-phase compositions by the requested method.
template <class Scalar, int N>
OPM_PTFK Status flash2ph(Context<Scalar, N>& ctx,
                         const V<Scalar, N>& z,
                         V<Scalar, N>& K,
                         Scalar& L,
                         Phases<Scalar, N>& s)
{
    bool converged = false;
    switch (ctx.method) {
    case Method::Newton:
        OPM_PTFK_TRY(newtonComposition(ctx, K, L, s, z));
        converged = true;
        break;
    case Method::Ssi:
        OPM_PTFK_TRY(successiveSubstitution(ctx, K, L, s, z, false, false, converged));
        break;
    case Method::SsiAccelerated:
        OPM_PTFK_TRY(successiveSubstitution(ctx, K, L, s, z, false, true, converged));
        break;
    case Method::SsiNewton:
        OPM_PTFK_TRY(successiveSubstitution(ctx, K, L, s, z, true, false, converged));
        if (converged) {
            // The Newton polish; a failed one returns the substitution answer.
            const V<Scalar, N> K_substitution = K;
            const Scalar L_substitution = L;
            const Phases<Scalar, N> s_substitution = s;
            const Status polish = newtonComposition(ctx, K, L, s, z);
            if (polish != Status::ok) {
                if (!recoverable(polish)) {
                    return polish;
                }
                K = K_substitution;
                L = L_substitution;
                s = s_substitution;
            }
        }
        if (!converged) {
            const Status newton = newtonComposition(ctx, K, L, s, z);
            if (newton == Status::ok) {
                converged = true;
            }
            else if (recoverable(newton)) {
                OPM_PTFK_TRY(successiveSubstitution(ctx, K, L, s, z, false, true, converged));
            }
            else {
                return newton;
            }
        }
        break;
    }
    if (!converged) {
        return Status::splitNotConverged;
    }
    if (!math::isfinite(L) || !phaseCompositionsArePhysical(s.x[liquid], s.x[vapour])) {
        return Status::splitNonPhysical;
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::feedTangentPlane_(): d_i = ln z_i + ln phi_i(z) on the
 *        root of the lower Gibbs energy, with translated coefficients.
 *
 * \param exists false where a coefficient is at its clamp
 */
template <class Scalar, int N>
OPM_PTFK Status feedTangentPlane(const Pure<Scalar, N>& pure, const V<Scalar, N>& z, bool& exists, V<Scalar, N>& plane)
{
    using C = Constants<Scalar>;
    Phases<Scalar, N> feed;
    feed.x[liquid] = z;
    feed.x[vapour] = z;
    exists = false;
    Scalar plane_gibbs = 0.0;
    for (int phase = 0; phase < 2; ++phase) {
        OPM_PTFK_TRY(updatePhase(pure, feed, phase));
        V<Scalar, N> phase_plane = filled<Scalar, N>(0.0);
        Scalar phase_gibbs = 0.0;
        for (int i = 0; i < N; ++i) {
            const Scalar z_i = z[i];
            if (!(z_i > 0.)) {
                continue;
            }
            const Scalar phi = feed.phi[phase][i];
            if (!(phi > C::minFugacityCoefficient && phi < C::maxFugacityCoefficient)) {
                exists = false;
                return Status::ok;
            }
            phase_plane[i] = math::log(z_i * phi);
            phase_gibbs += z_i * phase_plane[i];
        }
        if (!exists || phase_gibbs < plane_gibbs) {
            plane = phase_plane;
            plane_gibbs = phase_gibbs;
            exists = true;
        }
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::liquidPhaseIsStable_(): a liquid-like trial grown from the
 *        liquid of a split, from the incoming ratios.
 */
template <class Scalar, int N>
OPM_PTFK Status liquidPhaseIsStable(Context<Scalar, N>& ctx, const Phases<Scalar, N>& s, bool& stable)
{
    const V<Scalar, N> liquid_x = s.x[liquid];
    V<Scalar, N> trial_K = ctx.incomingK;
    FeedRoot<Scalar> root;
    V<Scalar, N> reference;
    OPM_PTFK_TRY(feedStableRoot(ctx.pure, liquid_x, root, reference));
    bool is_trivial = false;
    V<Scalar, N> trial_composition;
    Scalar trial_sum = 0.0;
    OPM_PTFK_TRY(checkStability(ctx.pure, ctx.counters, is_trivial, trial_K, trial_composition, trial_sum,
                                liquid_x, /*is_vapour_trial=*/false, reference));
    stable = is_trivial || trial_sum < 1.0 + Constants<Scalar>::stabilityTolerance;
    return Status::ok;
}

/*!
 * \brief PTFlash::replaceByLowerSplit_(): the split from the stability test
 *        of the feed in place of the given one when its Gibbs energy is lower.
 */
template <class Scalar, int N>
OPM_PTFK Status replaceByLowerSplit(Context<Scalar, N>& ctx,
                                    V<Scalar, N>& K,
                                    Scalar& L,
                                    Phases<Scalar, N>& s,
                                    const Phases<Scalar, N>& incoming,
                                    const V<Scalar, N>& z)
{
    Phases<Scalar, N> other = incoming;
    V<Scalar, N> other_K = ctx.incomingK;
    bool feed_is_stable = false;
    FeedRoot<Scalar> root;
    OPM_PTFK_TRY(phaseStabilityTest(ctx.pure, ctx.counters, feed_is_stable, other_K, other, z, root));
    if (feed_is_stable) {
        return Status::ok;
    }
    // The rest is inside PTFlash's try: its errors keep the given split.
    const auto attempt = [&]() -> Status {
        ++ctx.counters.rachfordRiceSolves;
        Scalar other_L = 0;
        OPM_PTFK_TRY(rachfordRice(other_K, z, false, Scalar(0), other_L));
        if (!(other_L > 0. && other_L < 1.)) {
            return Status::ok;
        }
        OPM_PTFK_TRY(flash2ph(ctx, z, other_K, other_L, other));
        if (!(other_L > 0. && other_L < 1.) || phasesCoincide(other)) {
            return Status::ok;
        }
        Phases<Scalar, N> given = s;
        Scalar other_energy = 0;
        Scalar given_energy = 0;
        OPM_PTFK_TRY(splitGibbsEnergyUpdated(ctx.pure, other, other_L, other_energy));
        OPM_PTFK_TRY(splitGibbsEnergyUpdated(ctx.pure, given, L, given_energy));
        if (other_energy < given_energy - Constants<Scalar>::stabilityTolerance) {
            K = other_K;
            L = other_L;
            s = other;
        }
        return Status::ok;
    };
    const Status status = attempt();
    if (status != Status::ok && !recoverable(status)) {
        return status;
    }
    return Status::ok;
}

/*!
 * \brief PTFlash::flashFirst_(): substitution steps of the split from the
 *        incoming ratios; a split below the feed's tangent plane makes the
 *        stability test unnecessary.
 */
template <class Scalar, int N>
OPM_PTFK Status flashFirst(Context<Scalar, N>& ctx,
                           V<Scalar, N>& K,
                           Scalar& L,
                           Phases<Scalar, N>& s,
                           const V<Scalar, N>& z,
                           bool& accepted)
{
    using C = Constants<Scalar>;
    accepted = false;
    if (L == 0. || L == 1.) {
        return Status::ok;
    }
    const Phases<Scalar, N> incoming = s;
    V<Scalar, N> trial_K = K;
    for (int i = 0; i < N; ++i) {
        if (!(z[i] > 0.)) {
            trial_K[i] = 1.;
        }
    }
    Scalar trial_L = 0.;
    bool has_plane = false;
    V<Scalar, N> feed_plane;
    bool feed_is_unstable = false;
    for (int iteration = 0; iteration < C::flashFirstMaxIterations; ++iteration) {
        ++ctx.counters.rachfordRiceSolves;
        Scalar rr_L = trial_L;
        if (rachfordRice(trial_K, z, false, Scalar(0), rr_L) != Status::ok) {
            break;
        }
        trial_L = rr_L;
        if (!(trial_L > 0. && trial_L < 1.)) {
            break;
        }
        ++ctx.counters.ssiIterations;
        computeLiquidVapor(s, trial_L, trial_K, z);
        bool is_clamped = false;
        for (int phase = 0; phase < 2; ++phase) {
            OPM_PTFK_TRY(updatePhase(ctx.pure, s, phase));
            for (int i = 0; i < N; ++i) {
                const Scalar phi = s.phi[phase][i];
                is_clamped = is_clamped || !(phi > C::minFugacityCoefficient && phi < C::maxFugacityCoefficient);
            }
        }
        if (is_clamped) {
            break;
        }
        if (!has_plane) {
            OPM_PTFK_TRY(feedTangentPlane(ctx.pure, z, has_plane, feed_plane));
            if (!has_plane) {
                break;
            }
        }
        Scalar distance[2] = {0.0, 0.0};
        for (int phase = 0; phase < 2; ++phase) {
            for (int i = 0; i < N; ++i) {
                const Scalar w = s.x[phase][i];
                if (z[i] > 0. && w > 0.) {
                    distance[phase] += w * (math::log(w * s.phi[phase][i]) - feed_plane[i]);
                }
            }
        }
        const Scalar liquid_fraction = trial_L;
        if (liquid_fraction * distance[liquid] + (1. - liquid_fraction) * distance[vapour]
            < -C::stabilityTolerance) {
            feed_is_unstable = true;
            break;
        }
        const bool phase_below_plane
            = distance[liquid] < -C::stabilityTolerance || distance[vapour] < -C::stabilityTolerance;
        if (!phase_below_plane && iteration + 1 >= C::flashFirstIterations) {
            break;
        }
        for (int i = 0; i < N; ++i) {
            if (z[i] > 0.) {
                trial_K[i] *= fugacity(ctx.pure, s, liquid, i) / fugacity(ctx.pure, s, vapour, i);
            }
        }
    }

    if (feed_is_unstable) {
        // PTFlash's try: an error restores the incoming state below, with
        // whatever K and L were taken by then.
        const auto attempt = [&]() -> Status {
            OPM_PTFK_TRY(flash2ph(ctx, z, trial_K, trial_L, s));
            const Scalar liquid_fraction = trial_L;
            if (liquid_fraction > 0. && liquid_fraction < 1. && !phasesCoincide(s)) {
                K = trial_K;
                L = trial_L;
                bool liquid_stable = true;
                OPM_PTFK_TRY(liquidPhaseIsStable(ctx, s, liquid_stable));
                if (!liquid_stable) {
                    OPM_PTFK_TRY(replaceByLowerSplit(ctx, K, L, s, incoming, z));
                }
                accepted = true;
            }
            return Status::ok;
        };
        const Status status = attempt();
        if (status != Status::ok && !recoverable(status)) {
            return status;
        }
        if (accepted) {
            return Status::ok;
        }
    }
    s = incoming;
    return Status::ok;
}

//! Li's single-phase label (PTFlash::li_single_phase_label_()).
template <class Scalar, int N>
OPM_PTFK Scalar liLabel(const Parameters<Scalar, N>& params, Scalar T, const V<Scalar, N>& z)
{
    Scalar sumVz = 0.0;
    for (int i = 0; i < N; ++i) {
        sumVz += (params.criticalVolume[i] * z[i]);
    }
    Scalar Tc_est = 0.0;
    for (int i = 0; i < N; ++i) {
        Tc_est += (params.criticalVolume[i] * params.criticalTemperature[i] * z[i] / sumVz);
    }
    return T < Tc_est ? Scalar(1.0) : Scalar(0.0);
}

/*!
 * \brief PTFlash::flashScalar_() without a memory.
 */
template <class Scalar, int N>
OPM_PTFK Status flashScalar(Context<Scalar, N>& ctx,
                            const V<Scalar, N>& z,
                            V<Scalar, N>& K,
                            Scalar& L,
                            Phases<Scalar, N>& s,
                            bool& is_single_phase)
{
    bool has_feed_root = false;
    FeedRoot<Scalar> stable_feed_root{liquid, false, Scalar(0)};
    is_single_phase = false;
    bool is_negative_flash = false;

    const bool is_cold_start = L <= 0 || L >= 1;
    bool split_from_flash_first = false;
    if (is_cold_start) {
        OPM_PTFK_TRY(flashFirst(ctx, K, L, s, z, split_from_flash_first));
    }
    if (is_cold_start && !split_from_flash_first) {
        OPM_PTFK_TRY(phaseStabilityTest(ctx.pure, ctx.counters, is_single_phase, K, s, z, stable_feed_root));
        has_feed_root = true;
    }
    if (!is_single_phase && !split_from_flash_first) {
        bool is_warm_split = !is_cold_start;
        const auto try_solve_split = [&]() -> Status {
            ++ctx.counters.rachfordRiceSolves;
            Scalar rr_L = L;
            OPM_PTFK_TRY(rachfordRice(K, z, true, L, rr_L));
            L = rr_L;
            if (is_warm_split && !(L > 0. && L < 1.)) {
                return Status::warmRatiosNoSplit;
            }
            return flash2ph(ctx, z, K, L, s);
        };
        const auto classify_negative_flash = [&]() {
            const Scalar liquid_fraction = L;
            is_negative_flash = liquid_fraction <= 0. || liquid_fraction >= 1.;
            is_single_phase = is_negative_flash;
            if (is_negative_flash) {
                L = liquid_fraction <= 0. ? 0. : 1.;
            }
        };

        const Status split = try_solve_split();
        if (split != Status::ok && !recoverable(split)) {
            return split;
        }
        const bool split_succeeded = split == Status::ok;
        if (split_succeeded) {
            classify_negative_flash();
        }
        const bool phases_coincide
            = split_succeeded && (!is_single_phase || !is_cold_start) && phasesCoincide(s);
        if (phases_coincide && is_cold_start) {
            return Status::coincidentPhases;
        }
        if (!split_succeeded || (phases_coincide && !is_cold_start)) {
            is_single_phase = false;
            is_negative_flash = false;
            for (int i = 0; i < N; ++i) {
                K[i] = wilsonK(ctx.params, ctx.pure.p, ctx.pure.T, i);
            }
            OPM_PTFK_TRY(phaseStabilityTest(ctx.pure, ctx.counters, is_single_phase, K, s, z, stable_feed_root));
            has_feed_root = true;
            if (!is_single_phase) {
                ++ctx.counters.splitRetries;
                is_warm_split = false;
                const Status retry = try_solve_split();
                if (retry != Status::ok) {
                    return recoverable(retry) ? Status::retryFailed : retry;
                }
                classify_negative_flash();
                if (!is_single_phase && phasesCoincide(s)) {
                    return Status::coincidentPhases;
                }
            }
        }
    }
    if (is_single_phase) {
        s.x[vapour] = z;
        s.x[liquid] = z;
        if (!is_negative_flash) {
            FeedRoot<Scalar> root = stable_feed_root;
            if (!has_feed_root) {
                V<Scalar, N> unused;
                OPM_PTFK_TRY(feedStableRoot(ctx.pure, z, root, unused));
            }
            if (root.hasTwoRoots) {
                L = root.phase == liquid ? 1.0 : 0.0;
            }
            else {
                Scalar critical_volume = 0.0;
                for (int i = 0; i < N; ++i) {
                    critical_volume += z[i] * ctx.params.criticalVolume[i] / 1.0e3;
                }
                if (root.molarVolume > critical_volume) {
                    L = 0.0;
                }
                else {
                    L = liLabel(ctx.params, ctx.pure.T, z);
                }
            }
        }
    }
    return Status::ok;
}

} // namespace detail

/*!
 * \brief The isothermal flash of one state: PTFlash::flash_solve_scalar()
 *        on (p, T, z) with the incoming ratios \p K and liquid fraction
 *        \p L (L outside (0, 1): a cold start, with the stability test).
 *
 * \return The answer: status, phase state, L, x, y, the final ratios and
 *         the work counters. For a single phase x = y = z and L is the
 *         label (1 liquid, 0 vapour).
 */
template <class Scalar, int N>
OPM_PTFK Answer<Scalar, N> flash(const Parameters<Scalar, N>& params,
                                 Scalar p,
                                 Scalar T,
                                 const Vec<Scalar, N>& z,
                                 const Vec<Scalar, N>& K,
                                 Scalar L,
                                 const Settings<Scalar>& settings = {})
{
    Answer<Scalar, N> answer;
    answer.K = K;
    answer.L = L;
    detail::Context<Scalar, N> ctx{params, {}, K, settings.method, settings.tolerance, answer.counters};
    answer.status = detail::makePure(params, p, T, ctx.pure);
    if (answer.status != Status::ok) {
        return answer;
    }
    detail::Phases<Scalar, N> s;
    for (int phase = 0; phase < 2; ++phase) {
        s.x[phase] = filled<Scalar, N>(0.);
        s.phi[phase] = filled<Scalar, N>(0.);
    }
    answer.status = detail::flashScalar(ctx, z, answer.K, answer.L, s, answer.singlePhase);
    answer.x = s.x[detail::liquid];
    answer.y = s.x[detail::vapour];
    return answer;
}

#undef OPM_PTFK_TRY
#undef OPM_PTFK

} // namespace Opm::PTFlashKernel

#endif // OPM_PTFLASH_KERNEL_HPP
