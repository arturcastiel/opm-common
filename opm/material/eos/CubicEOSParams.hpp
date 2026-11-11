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
#ifndef CUBIC_EOS_PARAMS_HPP
#define CUBIC_EOS_PARAMS_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/material/Constants.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/material/eos/PRParams.hpp>
#include <opm/material/eos/RKParams.hpp>
#include <opm/material/eos/SRKParams.hpp>

#include <array>
#include <cmath>

#include <fmt/format.h>

namespace Opm
{

namespace detail {
//! The pure-component parameters of the last updatePure() on this thread
//! that computed them, for one scalar type and fluid system.
template <class Scalar, class FluidSystem>
struct CubicEOSPureMemo
{
    static constexpr unsigned numComponents = FluidSystem::numComponents;

    bool valid = false;
    unsigned generation = 0;
    CompositionalConfig::EOSType eosType{};
    Scalar temperature{};
    Scalar pressure{};
    std::array<Scalar, numComponents> Ai{};
    std::array<Scalar, numComponents> Bi{};
    std::array<std::array<Scalar, numComponents>, numComponents> aCache{};
    std::array<Scalar, numComponents> translation{};

    static CubicEOSPureMemo& get()
    {
        static thread_local CubicEOSPureMemo memo;
        return memo;
    }
};
} // namespace detail

template <class Scalar, class FluidSystem, unsigned phaseIdx>
class CubicEOSParams
{
    enum { numComponents = FluidSystem::numComponents };

    static constexpr Scalar R = Constants<Scalar>::R;

    using EOSType = CompositionalConfig::EOSType;
    using PR = Opm::PRParams<Scalar, FluidSystem>;
    using RK = Opm::RKParams<Scalar, FluidSystem>;
    using SRK = Opm::SRKParams<Scalar, FluidSystem>;

public:

    void setEOSType(const EOSType eos_type)
    {
        EosType_= eos_type;
        pureUpToDate_ = false;
    }

    /*!
     * \brief The pure-component parameters \f$A_i\f$, \f$B_i\f$, the mixing
     *        matrix \f$A_{ij}\f$ and the volume translation factors at a
     *        temperature and pressure.
     *
     * They depend on temperature and pressure only. A repeated call at the
     * temperature and pressure of the last completed one keeps them, so the
     * composition updates of a flash, which all share its temperature and
     * pressure, compute them once.
     */
    void updatePure(Scalar temperature, Scalar pressure)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(pressure);

        if (pureIsCurrent(temperature, pressure)) {
            return;
        }
        pureUpToDate_ = false;

        // The parameter caches of one flash, and of the phases of a cache,
        // share the temperature and pressure: the last computation on this
        // thread holds them if the fluid system has not changed since.
        using Memo = detail::CubicEOSPureMemo<Scalar, FluidSystem>;
        constexpr bool hasGeneration = requires { FluidSystem::parameterGeneration(); };
        if constexpr (hasGeneration) {
            const Memo& memo = Memo::get();
            if (memo.valid && memo.generation == FluidSystem::parameterGeneration()
                && memo.eosType == EosType_ && memo.temperature == temperature
                && memo.pressure == pressure) {
                Ai_ = memo.Ai;
                Bi_ = memo.Bi;
                aCache_ = memo.aCache;
                translation_ = memo.translation;
                pureTemperature_ = temperature;
                purePressure_ = pressure;
                pureUpToDate_ = true;
                return;
            }
        }

        // calculate Ai and Bi
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar pr = pressure / FluidSystem::criticalPressure(compIdx);
            Scalar Tr = temperature / FluidSystem::criticalTemperature(compIdx);
            Scalar OmegaA = OmegaA_(temperature, compIdx);
            Scalar OmegaB = OmegaB_();
            // OMEGAA and OMEGAB replace the constants of the equation per component.
            // Fluid systems without support for them keep the constants.
            if constexpr (requires { FluidSystem::omegaA(0u); }) {
                if (const auto omegaA = FluidSystem::omegaA(compIdx)) {
                    OmegaA *= *omegaA / constantOmegaA_();
                }
            }
            if constexpr (requires { FluidSystem::omegaB(0u); }) {
                if (const auto omegaB = FluidSystem::omegaB(compIdx)) {
                    OmegaB = *omegaB;
                }
            }

            Scalar newA = OmegaA * pr / (Tr * Tr);
            Scalar newB = OmegaB * pr / Tr;
            if (!std::isfinite(scalarValue(newA)) || !std::isfinite(scalarValue(newB))) {
                OPM_THROW(NumericalProblem,
                          "CubicEOSParams::updatePure: non-finite pure-component EoS "
                          "parameter (check the temperature/pressure input)");
            }

            setAi(newA, compIdx);
            setBi(newB, compIdx);
            Valgrind::CheckDefined(Ai(compIdx));
            Valgrind::CheckDefined(Bi(compIdx));
        }

        // Update Aij
        updateACache_();

        // The SSHIFT factor of the fugacity coefficients, exp(-s_i B_i).
        if constexpr (requires { FluidSystem::volumeShift(0u); }) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                translation_[compIdx] = exp(-FluidSystem::volumeShift(compIdx) * Bi(compIdx));
            }
        }

        pureTemperature_ = temperature;
        purePressure_ = pressure;
        pureUpToDate_ = true;

        if constexpr (hasGeneration) {
            Memo& memo = Memo::get();
            memo.valid = true;
            memo.generation = FluidSystem::parameterGeneration();
            memo.eosType = EosType_;
            memo.temperature = temperature;
            memo.pressure = pressure;
            memo.Ai = Ai_;
            memo.Bi = Bi_;
            memo.aCache = aCache_;
            memo.translation = translation_;
        }
    }

    /*!
     * \brief Whether the pure-component parameters are those of this
     *        temperature and pressure.
     */
    bool pureIsCurrent(const Scalar& temperature, const Scalar& pressure) const
    {
        return pureUpToDate_ && temperature == pureTemperature_ && pressure == purePressure_;
    }

    /*!
     * \brief Takes the pure-component parameters of another phase when that
     *        phase holds them for this temperature, pressure and equation of
     *        state; they do not depend on the phase.
     *
     * \return Whether they were taken.
     */
    template <unsigned otherPhaseIdx>
    bool takePure(const CubicEOSParams<Scalar, FluidSystem, otherPhaseIdx>& other,
                  const Scalar& temperature,
                  const Scalar& pressure)
    {
        if (!other.pureUpToDate_ || other.EosType_ != EosType_
            || !(temperature == other.pureTemperature_) || !(pressure == other.purePressure_)) {
            return false;
        }
        Ai_ = other.Ai_;
        Bi_ = other.Bi_;
        aCache_ = other.aCache_;
        translation_ = other.translation_;
        pureTemperature_ = other.pureTemperature_;
        purePressure_ = other.purePressure_;
        pureUpToDate_ = true;
        return true;
    }

    template <class FluidState>
    void updateMix(const FluidState& fs)
    {
        using FlashEval = typename FluidState::ValueType;

        // The mole fractions limited to [0, 1], once per component.
        std::array<FlashEval, numComponents> x;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const FlashEval moleFrac = fs.moleFraction(phaseIdx, compIdx);
            x[compIdx] = max(0.0, min(1.0, moleFrac));
            Valgrind::CheckDefined(x[compIdx]);
        }

        FlashEval newA = 0;
        FlashEval newB = 0;
        for (unsigned compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            const FlashEval& xi = x[compIIdx];

            for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                // Calculate A
                newA +=  xi * x[compJIdx] * aCache_[compIIdx][compJIdx];
            }

            // Calculate B
            newB += max(0.0, xi) * Bi(compIIdx);
        }

        // check once after the double loop rather than on every accumulation step
        if (!std::isfinite(scalarValue(newA)) || !std::isfinite(scalarValue(newB))) {
            OPM_THROW(NumericalProblem,
                      fmt::format("CubicEOSParams::updateMix: non-finite mixture parameters "
                                  "for phase {} (A = {}, B = {}); typically a diverged "
                                  "composition update feeding the EoS",
                                  phaseIdx, scalarValue(newA), scalarValue(newB)));
        }

        // assign A and B
        setA(decay<Scalar>(newA));
        setB(decay<Scalar>(newB));
        Valgrind::CheckDefined(A());
        Valgrind::CheckDefined(B());
    }

    template <class FluidState>
    void updateSingleMoleFraction(const FluidState& fs,
                                  unsigned /*compIdx*/)
    {
        updateMix(fs);
    }

    Scalar aCache(unsigned compIIdx, unsigned compJIdx ) const
    {
        return aCache_[compIIdx][compJIdx];
    }

    void setAi(Scalar value, unsigned compIdx)
    {
        Ai_[compIdx] = value;
        pureUpToDate_ = false;
    }

    void setBi(Scalar value, unsigned compIdx)
    {
        Bi_[compIdx] = value;
        pureUpToDate_ = false;
    }

    /*!
     * \brief The volume translation factor \f$\exp(-s_i B_i)\f$ of a
     *        component's fugacity coefficient, for fluid systems with a
     *        volume shift \f$s_i\f$ (SSHIFT).
     */
    Scalar volumeTranslationFactor(unsigned compIdx) const
    {
        return translation_[compIdx];
    }

    Scalar Ai(unsigned compIdx) const
    {
        return Ai_[compIdx];
    }

    Scalar Bi(unsigned compIdx) const
    {
        return Bi_[compIdx];
    }

    /*!
     * \brief \f$\partial \ln A_i / \partial T\f$ at the temperature and
     *        pressure of the pure-component parameters:
     *        \f$A_i \propto \Omega_{a,i}(T)\, p / T^2\f$, and OMEGAA scales
     *        \f$\Omega_{a,i}\f$ by a constant. \f$B_i \propto p / T\f$.
     */
    Scalar lnAiTemperatureDerivative(unsigned compIdx) const
    {
        const Scalar& T = pureTemperature_;
        Scalar lnOmegaA;
        switch (EosType_) {
            case EOSType::PRCORR:
                lnOmegaA = PR::calcLnOmegaATemperatureDerivative(T, compIdx, /*modified=*/true);
                break;
            case EOSType::PR:
                lnOmegaA = PR::calcLnOmegaATemperatureDerivative(T, compIdx, /*modified=*/false);
                break;
            case EOSType::RK:
                lnOmegaA = RK::calcLnOmegaATemperatureDerivative(T, compIdx);
                break;
            case EOSType::SRK:
                lnOmegaA = SRK::calcLnOmegaATemperatureDerivative(T, compIdx);
                break;
            default:
                throw std::runtime_error("EOS type not implemented!");
        }
        return lnOmegaA - 2 / T;
    }

    void setA(Scalar value)
    {
        A_ = value;
    }

    void setB(Scalar value)
    {
        B_ = value;
    }

    Scalar A() const
    {
        return A_;
    }

    Scalar B() const
    {
        return B_;
    }

    Scalar m1() const
    {
        switch (EosType_) {
            case EOSType::PRCORR:
            case EOSType::PR:
                return PR::calcm1();
            case EOSType::RK:
                return RK::calcm1();
            case EOSType::SRK:
                return SRK::calcm1();
            default:
                throw std::runtime_error("EOS type not implemented!");
        }
    }

    Scalar m2() const
    {
        switch (EosType_) {
            case EOSType::PRCORR:
            case EOSType::PR:
                return PR::calcm2();
            case EOSType::RK:
                return RK::calcm2();
            case EOSType::SRK:
                return SRK::calcm2();
            default:
                throw std::runtime_error("EOS type not implemented!");
        }
    }

protected:
    std::array<Scalar, numComponents> Ai_;
    std::array<Scalar, numComponents> Bi_;
    Scalar A_;
    Scalar B_;
    std::array<std::array<Scalar, numComponents>, numComponents> aCache_;
    std::array<Scalar, numComponents> translation_;

    EOSType EosType_;

    // The temperature and pressure of the last completed updatePure().
    Scalar pureTemperature_;
    Scalar purePressure_;
    bool pureUpToDate_ = false;

private:
    void updateACache_()
    {
        // A_ij is symmetric, since A_i A_j = A_j A_i in floating point too.
        for (unsigned compIIdx = 0; compIIdx < numComponents; ++ compIIdx) {
            for (unsigned compJIdx = compIIdx; compJIdx < numComponents; ++ compJIdx) {
                // interaction coefficient as given in SPE5
                Scalar Psi = FluidSystem::interactionCoefficient(compIIdx, compJIdx);

                aCache_[compIIdx][compJIdx] = sqrt(Ai(compIIdx) * Ai(compJIdx)) * (1 - Psi);
                aCache_[compJIdx][compIIdx] = aCache_[compIIdx][compJIdx];
            }
        }
    }

    Scalar OmegaA_(Scalar temperature, unsigned compIdx)
    {
        switch (EosType_) {
            case EOSType::PRCORR:
                return PR::calcOmegaA(temperature, compIdx, /*modified=*/true);
            case EOSType::PR:
                return PR::calcOmegaA(temperature, compIdx, /*modified=*/false);
            case EOSType::RK:
                return RK::calcOmegaA(temperature, compIdx);
            case EOSType::SRK:
                return SRK::calcOmegaA(temperature, compIdx);
            default:
                throw std::runtime_error("EOS type not implemented!");
        }
    }

    Scalar OmegaB_()
    {
        switch (EosType_) {
            case EOSType::PRCORR:
            case EOSType::PR:
                return PR::calcOmegaB();
            case EOSType::RK:
                return RK::calcOmegaB();
            case EOSType::SRK:
                return SRK::calcOmegaB();
            default:
                throw std::runtime_error("EOS type not implemented!");
        }
    }

    //! \brief The \f$\Omega_a\f$ that OmegaA_() includes, before temperature scaling.
    Scalar constantOmegaA_() const
    {
        switch (EosType_) {
        case EOSType::PRCORR:
        case EOSType::PR:
            return PR::omegaA();
        case EOSType::RK:
            return RK::omegaA();
        case EOSType::SRK:
            return SRK::omegaA();
        default:
            throw std::runtime_error("EOS type not implemented!");
        }
    }

    template <class, class, unsigned>
    friend class CubicEOSParams;
};  // class CubicEOSParams

}  // namespace Opm

#endif
