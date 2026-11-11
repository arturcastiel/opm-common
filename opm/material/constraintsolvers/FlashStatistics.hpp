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
 * \brief Work counters of the P-T and P-H flashes, for benchmarks.
 *
 * The counters cost nothing unless OPM_FLASH_STATISTICS is defined before the
 * flash headers are included.
 */
#ifndef OPM_FLASH_STATISTICS_HPP
#define OPM_FLASH_STATISTICS_HPP

namespace Opm {

//! Work done by the flashes of one thread since the last reset().
struct FlashStatistics
{
    long stabilityTests = 0;      //!< calls of the two-trial stability test
    long stabilityIterations = 0; //!< substitution steps of the stability trials
    long ssiIterations = 0;       //!< substitution steps of the two-phase split
    long newtonIterations = 0;    //!< Newton steps of the two-phase split
    long rachfordRiceSolves = 0;  //!< Rachford-Rice solves
    long splitRetries = 0;        //!< splits redone after a new stability test
    long eosUpdates = 0;          //!< phase updates of the EOS parameter cache
    long phEvaluations = 0;       //!< isothermal flashes of the P-H temperature search
    long cacheHits = 0;           //!< answers or starts taken from a cache

    static FlashStatistics& get()
    {
        static thread_local FlashStatistics statistics;
        return statistics;
    }

    void reset()
    {
        *this = FlashStatistics{};
    }
};

} // namespace Opm

#if defined(OPM_FLASH_STATISTICS)
#define OPM_FLASH_COUNT(field) (++::Opm::FlashStatistics::get().field)
#else
#define OPM_FLASH_COUNT(field) ((void)0)
#endif

#endif // OPM_FLASH_STATISTICS_HPP
