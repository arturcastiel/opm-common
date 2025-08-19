/*
  Copyright 2014 Statoil IT
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
#include "config.h"
#include <cstddef>

#define BOOST_TEST_MODULE Restart_File_IO

#include <boost/test/unit_test.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/AggregateAquiferData.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/EclIOdata.hpp>
#include <opm/io/eclipse/OutputStream.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FIPRegionStatistics.hpp>
#include <opm/input/eclipse/EclipseState/Grid/RegionSetMatcher.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Eqldims.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/MSW/SegmentMatcher.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/input/eclipse/Utility/Functional.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <tests/WorkArea.hpp>

using namespace Opm;

namespace {

data::GroupAndNetworkValues mkGroups()
{
    return {};
}

data::Wells mkWells()
{
    data::Rates r1, r2, rc1, rc2, rc3;
    r1.set( data::Rates::opt::wat, 5.67 );
    r1.set( data::Rates::opt::oil, 6.78 );
    r1.set( data::Rates::opt::gas, 7.89 );

    r2.set( data::Rates::opt::wat, 8.90 );
    r2.set( data::Rates::opt::oil, 9.01 );
    r2.set( data::Rates::opt::gas, 10.12 );

    rc1.set( data::Rates::opt::wat, 20.41 );
    rc1.set( data::Rates::opt::oil, 21.19 );
    rc1.set( data::Rates::opt::gas, 22.41 );

    rc2.set( data::Rates::opt::wat, 23.19 );
    rc2.set( data::Rates::opt::oil, 24.41 );
    rc2.set( data::Rates::opt::gas, 25.19 );

    rc3.set( data::Rates::opt::wat, 26.41 );
    rc3.set( data::Rates::opt::oil, 27.19 );
    rc3.set( data::Rates::opt::gas, 28.41 );

    data::Well w1, w2;
    w1.rates = r1;
    w1.thp = 1.0;
    w1.bhp = 1.23;
    w1.temperature = 3.45;
    w1.control = 1;

    /*
     *  the completion keys (active indices) and well names correspond to the
     *  input deck. All other entries in the well structures are arbitrary.
     */
    Opm::data::ConnectionFiltrate con_filtrate {0.1, 1, 3, 0.4, 1.e-9, 0.2, 0.05, 10.}; // values are not used in this test
    w1.connections.push_back( { 88, rc1, 30.45, 123.4, 543.21, 0.62, 0.15, 1.0e3, 1.234, 0.0, 1.23, con_filtrate, } );
    w1.connections.push_back( { 288, rc2, 33.19, 123.4, 432.1, 0.26, 0.45, 2.56, 2.345, 0.0, 0.98, con_filtrate } );

    w2.rates = r2;
    w2.thp = 2.0;
    w2.bhp = 2.34;
    w2.temperature = 4.56;
    w2.control = 2;
    w2.connections.push_back( { 188, rc3, 36.22, 123.4, 256.1, 0.55, 0.0125, 314.15, 3.456, 0.0, 2.46, con_filtrate } );

    {
        data::Wells wellRates;

        wellRates["OP_1"] = w1;
        wellRates["OP_2"] = w2;

        return wellRates;
    }
}

data::Wells mkWellsLGR_Global()
{
    // This function creates a Wells object with two wells, each having two connections matching the one in the LGR_BASESIM2WELLS.DATA
    data::Rates r1, r2, rc1, rc2, rc3;
    r1.set( data::Rates::opt::wat, 5.67 );
    r1.set( data::Rates::opt::oil, 6.78 );
    r1.set( data::Rates::opt::gas, 7.89 );

    r2.set( data::Rates::opt::wat, 8.90 );
    r2.set( data::Rates::opt::oil, 9.01 );
    r2.set( data::Rates::opt::gas, 10.12 );

    rc1.set( data::Rates::opt::wat, 20.41 );
    rc1.set( data::Rates::opt::oil, 21.19 );
    rc1.set( data::Rates::opt::gas, 22.41 );

    rc2.set( data::Rates::opt::wat, 23.19 );
    rc2.set( data::Rates::opt::oil, 24.41 );
    rc2.set( data::Rates::opt::gas, 25.19 );



    data::Well w1, w2;
    w1.rates = r1;
    w1.thp = 1.0;
    w1.bhp = 1.23;
    w1.temperature = 3.45;
    w1.control = 1;

    /*
     *  the completion keys (active indices) and well names correspond to the
     *  input deck. All other entries in the well structures are arbitrary.
     */
    Opm::data::ConnectionFiltrate con_filtrate {0.1, 1, 3, 0.4, 1.e-9, 0.2, 0.05, 10.}; // values are not used in this test
    w1.connections.push_back( { 2, rc1, 30.45, 123.4, 543.21, 0.62, 0.15, 1.0e3, 1.234, 0.0, 1.23, con_filtrate, 1 } );

    w2.rates = r2;
    w2.thp = 2.0;
    w2.bhp = 2.34;
    w2.temperature = 4.56;
    w2.control = 2;
    w2.connections.push_back( { 0, rc2, 36.22, 123.4, 256.1, 0.55, 0.0125, 314.15, 3.456, 0.0, 2.46, con_filtrate,2 } );

    {
        data::Wells wellRates;

        wellRates["PROD"] = w1;
        wellRates["INJ"] = w2;

        return wellRates;
    }
}



data::Wells mkWellsLGR_LGR1()
{
    // This function creates a Wells object with two wells, each having two connections matching the one in the LGR_BASESIM2WELLS.DATA
    data::Rates r1, r2, rc1, rc2, rc3;
    r1.set( data::Rates::opt::wat, 5.67 );
    r1.set( data::Rates::opt::oil, 6.78 );
    r1.set( data::Rates::opt::gas, 7.89 );

    rc1.set( data::Rates::opt::wat, 20.41 );
    rc1.set( data::Rates::opt::oil, 21.19 );
    rc1.set( data::Rates::opt::gas, 22.41 );

    data::Well w1;
    w1.rates = r1;
    w1.thp = 1.0;
    w1.bhp = 1.23;
    w1.temperature = 3.45;
    w1.control = 1;

    /*
     *  the completion keys (active indices) and well names correspond to the
     *  input deck. All other entries in the well structures are arbitrary.
     */
    Opm::data::ConnectionFiltrate con_filtrate {0.1, 1, 3, 0.4, 1.e-9, 0.2, 0.05, 10.}; // values are not used in this test
    w1.connections.push_back( { 8, rc1, 30.45, 123.4, 543.21, 0.62, 0.15, 1.0e3, 1.234, 0.0, 1.23, con_filtrate, 1 } );

    {
        data::Wells wellRates;

        wellRates["PROD"] = w1;

        return wellRates;
    }
}


data::Wells mkWellsLGR_LGR2()
{
    // This function creates a Wells object with two wells, each having two connections matching the one in the LGR_BASESIM2WELLS.DATA
    data::Rates r2,  rc2;

    r2.set( data::Rates::opt::wat, 8.90 );
    r2.set( data::Rates::opt::oil, 9.01 );
    r2.set( data::Rates::opt::gas, 10.12 );


    rc2.set( data::Rates::opt::wat, 23.19 );
    rc2.set( data::Rates::opt::oil, 24.41 );
    rc2.set( data::Rates::opt::gas, 25.19 );



    data::Well  w2;


    /*
     *  the completion keys (active indices) and well names correspond to the
     *  input deck. All other entries in the well structures are arbitrary.
     */
    Opm::data::ConnectionFiltrate con_filtrate {0.1, 1, 3, 0.4, 1.e-9, 0.2, 0.05, 10.}; // values are not used in this test

    w2.rates = r2;
    w2.thp = 2.0;
    w2.bhp = 2.34;
    w2.temperature = 4.56;
    w2.control = 2;
    w2.connections.push_back( { 1, rc2, 36.22, 123.4, 256.1, 0.55, 0.0125, 314.15, 3.456, 0.0, 2.46, con_filtrate,2 } );

    {
        data::Wells wellRates;

        wellRates["INJ"] = w2;

        return wellRates;
    }
}





data::Solution mkSolution(int numCells)
{
    using measure = UnitSystem::measure;

    auto sol = data::Solution {
        { "PRESSURE", data::CellData { measure::pressure,    {}, data::TargetType::RESTART_SOLUTION } },
        { "TEMP",     data::CellData { measure::temperature, {}, data::TargetType::RESTART_SOLUTION } },
        { "SWAT",     data::CellData { measure::identity,    {}, data::TargetType::RESTART_SOLUTION } },
        { "SGAS",     data::CellData { measure::identity,    {}, data::TargetType::RESTART_SOLUTION } },
    };

    sol.data<double>("PRESSURE").assign( numCells, 6.0 );
    sol.data<double>("TEMP").assign( numCells, 7.0 );
    sol.data<double>("SWAT").assign( numCells, 8.0 );
    sol.data<double>("SGAS").assign( numCells, 9.0 );

    fun::iota rsi( 300.0, 300.0 + numCells );
    fun::iota rvi( 400.0, 400.0 + numCells );

    sol.insert("RS", measure::identity,
               std::vector<double>{ rsi.begin(), rsi.end() },
               data::TargetType::RESTART_SOLUTION);
    sol.insert("RV", measure::identity,
               std::vector<double>{ rvi.begin(), rvi.end() },
               data::TargetType::RESTART_SOLUTION);

    return sol;
}

data::Solution mkSolutionFIP(const int numCells)
{
    using measure = UnitSystem::measure;

    auto sol = data::Solution {
        { "PRESSURE", data::CellData { measure::pressure, {}, data::TargetType::RESTART_SOLUTION } },
        { "SWAT",     data::CellData { measure::identity, {}, data::TargetType::RESTART_SOLUTION } },
        { "SGAS",     data::CellData { measure::identity, {}, data::TargetType::RESTART_SOLUTION } },
        { "FIPOIL",   data::CellData { measure::identity, {}, data::TargetType::RESTART_SOLUTION } },
        { "FIPWAT",   data::CellData { measure::identity, {}, data::TargetType::RESTART_SOLUTION } },
        { "FIPGAS",   data::CellData { measure::identity, {}, data::TargetType::RESTART_SOLUTION } },
    };

    sol.data<double>("PRESSURE").assign(numCells, 6.0);
    sol.data<double>("SWAT").assign(numCells, 8.0);
    sol.data<double>("SGAS").assign(numCells, 9.0);
    sol.data<double>("FIPOIL").assign(numCells, 10.0);
    sol.data<double>("FIPWAT").assign(numCells, 11.0);
    sol.data<double>("FIPGAS").assign(numCells, 12.0);

    return sol;
}

Opm::SummaryState sim_stateLGR(const Opm::Schedule& sched)
{
    auto state = Opm::SummaryState {
        TimeService::now(),
        sched.back().udq().params().undefinedValue()
    };

    for (const auto& well : sched.getWellsatEnd()) {
        for (const auto& connection : well.getConnections()) {
            state.update_conn_var(well.name(), "CPR", connection.global_index() + 1, 111);
            if (well.isInjector()) {
                state.update_conn_var(well.name(), "COIR", connection.global_index() + 1, 222);
                state.update_conn_var(well.name(), "CGIR", connection.global_index() + 1, 333);
                state.update_conn_var(well.name(), "CWIR", connection.global_index() + 1, 444);
                state.update_conn_var(well.name(), "CVIR", connection.global_index() + 1, 555);

                state.update_conn_var(well.name(), "COIT", connection.global_index() + 1, 222 * 2.0);
                state.update_conn_var(well.name(), "CGIT", connection.global_index() + 1, 333 * 2.0);
                state.update_conn_var(well.name(), "CWIT", connection.global_index() + 1, 444 * 2.0);
                state.update_conn_var(well.name(), "CWIT", connection.global_index() + 1, 555 * 2.0);
            } else {
                state.update_conn_var(well.name(), "COPR", connection.global_index() + 1, 666);
                state.update_conn_var(well.name(), "CGPR", connection.global_index() + 1, 777);
                state.update_conn_var(well.name(), "CWPR", connection.global_index() + 1, 888);
                state.update_conn_var(well.name(), "CVPR", connection.global_index() + 1, 999);

                state.update_conn_var(well.name(), "CGOR", connection.global_index() + 1, 777.0 / 666.0);

                state.update_conn_var(well.name(), "COPT", connection.global_index() + 1, 555 * 2.0);
                state.update_conn_var(well.name(), "CGPT", connection.global_index() + 1, 666 * 2.0);
                state.update_conn_var(well.name(), "CWPT", connection.global_index() + 1, 777 * 2.0);
                state.update_conn_var(well.name(), "CVPT", connection.global_index() + 1, 999 * 2.0);
            }
        }
    }

    state.update_well_var("INJ", "WOPR", 1.0);
    state.update_well_var("INJ", "WWPR", 2.0);
    state.update_well_var("INJ", "WGPR", 3.0);
    state.update_well_var("INJ", "WVPR", 4.0);
    state.update_well_var("INJ", "WOPT", 10.0);
    state.update_well_var("INJ", "WWPT", 20.0);
    state.update_well_var("INJ", "WGPT", 30.0);
    state.update_well_var("INJ", "WVPT", 40.0);
    state.update_well_var("INJ", "WWIR", 0.0);
    state.update_well_var("INJ", "WGIR", 0.0);
    state.update_well_var("INJ", "WWIT", 0.0);
    state.update_well_var("INJ", "WGIT", 0.0);
    state.update_well_var("INJ", "WVIT", 0.0);
    state.update_well_var("INJ", "WWCT", 0.625);
    state.update_well_var("INJ", "WGOR", 234.5);
    state.update_well_var("INJ", "WBHP", 314.15);
    state.update_well_var("INJ", "WTHP", 123.45);
    state.update_well_var("INJ", "WOPTH", 345.6);
    state.update_well_var("INJ", "WWPTH", 456.7);
    state.update_well_var("INJ", "WGPTH", 567.8);
    state.update_well_var("INJ", "WWITH", 0.0);
    state.update_well_var("INJ", "WGITH", 0.0);
    state.update_well_var("INJ", "WGVIR", 0.0);
    state.update_well_var("INJ", "WWVIR", 0.0);

    state.update_well_var("PROD", "WOPR", 0.0);
    state.update_well_var("PROD", "WWPR", 0.0);
    state.update_well_var("PROD", "WGPR", 0.0);
    state.update_well_var("PROD", "WVPR", 0.0);
    state.update_well_var("PROD", "WOPT", 0.0);
    state.update_well_var("PROD", "WWPT", 0.0);
    state.update_well_var("PROD", "WGPT", 0.0);
    state.update_well_var("PROD", "WVPT", 0.0);
    state.update_well_var("PROD", "WWIR", 100.0);
    state.update_well_var("PROD", "WGIR", 200.0);
    state.update_well_var("PROD", "WWIT", 1000.0);
    state.update_well_var("PROD", "WGIT", 2000.0);
    state.update_well_var("PROD", "WVIT", 1234.5);
    state.update_well_var("PROD", "WWCT", 0.0);
    state.update_well_var("PROD", "WGOR", 0.0);
    state.update_well_var("PROD", "WBHP", 400.6);
    state.update_well_var("PROD", "WTHP", 234.5);
    state.update_well_var("PROD", "WOPTH", 0.0);
    state.update_well_var("PROD", "WWPTH", 0.0);
    state.update_well_var("PROD", "WGPTH", 0.0);
    state.update_well_var("PROD", "WWITH", 1515.0);
    state.update_well_var("PROD", "WGITH", 3030.0);
    state.update_well_var("PROD", "WGVIR", 1234.0);
    state.update_well_var("PROD", "WWVIR", 4321.0);



    state.update_group_var("G1", "GOPR" ,     110.0);
    state.update_group_var("G1", "GWPR" ,     120.0);
    state.update_group_var("G1", "GGPR" ,     130.0);
    state.update_group_var("G1", "GVPR" ,     140.0);
    state.update_group_var("G1", "GOPT" ,    1100.0);
    state.update_group_var("G1", "GWPT" ,    1200.0);
    state.update_group_var("G1", "GGPT" ,    1300.0);
    state.update_group_var("G1", "GVPT" ,    1400.0);
    state.update_group_var("G1", "GWIR" , -   256.0);
    state.update_group_var("G1", "GGIR" , - 65536.0);
    state.update_group_var("G1", "GWIT" ,   31415.9);
    state.update_group_var("G1", "GGIT" ,   27182.8);
    state.update_group_var("G1", "GVIT" ,   44556.6);
    state.update_group_var("G1", "GWCT" ,       0.625);
    state.update_group_var("G1", "GGOR" ,    1234.5);
    state.update_group_var("G1", "GGVIR",     123.45);
    state.update_group_var("G1", "GWVIR",    1234.56);
    state.update_group_var("G1", "GOPTH",    5678.90);
    state.update_group_var("G1", "GWPTH",    6789.01);
    state.update_group_var("G1", "GGPTH",    7890.12);
    state.update_group_var("G1", "GWITH",    8901.23);
    state.update_group_var("G1", "GGITH",    9012.34);

    state.update("FOPR" ,     1100.0);
    state.update("FWPR" ,     1200.0);
    state.update("FGPR" ,     1300.0);
    state.update("FVPR" ,     1400.0);
    state.update("FOPT" ,    11000.0);
    state.update("FWPT" ,    12000.0);
    state.update("FGPT" ,    13000.0);
    state.update("FVPT" ,    14000.0);
    state.update("FWIR" , -   2560.0);
    state.update("FGIR" , - 655360.0);
    state.update("FWIT" ,   314159.2);
    state.update("FGIT" ,   271828.1);
    state.update("FVIT" ,   445566.77);
    state.update("FWCT" ,        0.625);
    state.update("FGOR" ,     1234.5);
    state.update("FOPTH",    56789.01);
    state.update("FWPTH",    67890.12);
    state.update("FGPTH",    78901.23);
    state.update("FWITH",    89012.34);
    state.update("FGITH",    90123.45);
    state.update("FGVIR",     1234.56);
    state.update("FWVIR",    12345.67);

    return state;
}


Opm::SummaryState sim_state(const Opm::Schedule& sched)
{
    auto state = Opm::SummaryState {
        TimeService::now(),
        sched.back().udq().params().undefinedValue()
    };

    for (const auto& well : sched.getWellsatEnd()) {
        for (const auto& connection : well.getConnections()) {
            state.update_conn_var(well.name(), "CPR", connection.global_index() + 1, 111);
            if (well.isInjector()) {
                state.update_conn_var(well.name(), "COIR", connection.global_index() + 1, 222);
                state.update_conn_var(well.name(), "CGIR", connection.global_index() + 1, 333);
                state.update_conn_var(well.name(), "CWIR", connection.global_index() + 1, 444);
                state.update_conn_var(well.name(), "CVIR", connection.global_index() + 1, 555);

                state.update_conn_var(well.name(), "COIT", connection.global_index() + 1, 222 * 2.0);
                state.update_conn_var(well.name(), "CGIT", connection.global_index() + 1, 333 * 2.0);
                state.update_conn_var(well.name(), "CWIT", connection.global_index() + 1, 444 * 2.0);
                state.update_conn_var(well.name(), "CWIT", connection.global_index() + 1, 555 * 2.0);
            } else {
                state.update_conn_var(well.name(), "COPR", connection.global_index() + 1, 666);
                state.update_conn_var(well.name(), "CGPR", connection.global_index() + 1, 777);
                state.update_conn_var(well.name(), "CWPR", connection.global_index() + 1, 888);
                state.update_conn_var(well.name(), "CVPR", connection.global_index() + 1, 999);

                state.update_conn_var(well.name(), "CGOR", connection.global_index() + 1, 777.0 / 666.0);

                state.update_conn_var(well.name(), "COPT", connection.global_index() + 1, 555 * 2.0);
                state.update_conn_var(well.name(), "CGPT", connection.global_index() + 1, 666 * 2.0);
                state.update_conn_var(well.name(), "CWPT", connection.global_index() + 1, 777 * 2.0);
                state.update_conn_var(well.name(), "CVPT", connection.global_index() + 1, 999 * 2.0);
            }
        }
    }

    state.update_well_var("OP_1", "WOPR", 1.0);
    state.update_well_var("OP_1", "WWPR", 2.0);
    state.update_well_var("OP_1", "WGPR", 3.0);
    state.update_well_var("OP_1", "WVPR", 4.0);
    state.update_well_var("OP_1", "WOPT", 10.0);
    state.update_well_var("OP_1", "WWPT", 20.0);
    state.update_well_var("OP_1", "WGPT", 30.0);
    state.update_well_var("OP_1", "WVPT", 40.0);
    state.update_well_var("OP_1", "WWIR", 0.0);
    state.update_well_var("OP_1", "WGIR", 0.0);
    state.update_well_var("OP_1", "WWIT", 0.0);
    state.update_well_var("OP_1", "WGIT", 0.0);
    state.update_well_var("OP_1", "WVIT", 0.0);
    state.update_well_var("OP_1", "WWCT", 0.625);
    state.update_well_var("OP_1", "WGOR", 234.5);
    state.update_well_var("OP_1", "WBHP", 314.15);
    state.update_well_var("OP_1", "WTHP", 123.45);
    state.update_well_var("OP_1", "WOPTH", 345.6);
    state.update_well_var("OP_1", "WWPTH", 456.7);
    state.update_well_var("OP_1", "WGPTH", 567.8);
    state.update_well_var("OP_1", "WWITH", 0.0);
    state.update_well_var("OP_1", "WGITH", 0.0);
    state.update_well_var("OP_1", "WGVIR", 0.0);
    state.update_well_var("OP_1", "WWVIR", 0.0);

    state.update_well_var("OP_2", "WOPR", 0.0);
    state.update_well_var("OP_2", "WWPR", 0.0);
    state.update_well_var("OP_2", "WGPR", 0.0);
    state.update_well_var("OP_2", "WVPR", 0.0);
    state.update_well_var("OP_2", "WOPT", 0.0);
    state.update_well_var("OP_2", "WWPT", 0.0);
    state.update_well_var("OP_2", "WGPT", 0.0);
    state.update_well_var("OP_2", "WVPT", 0.0);
    state.update_well_var("OP_2", "WWIR", 100.0);
    state.update_well_var("OP_2", "WGIR", 200.0);
    state.update_well_var("OP_2", "WWIT", 1000.0);
    state.update_well_var("OP_2", "WGIT", 2000.0);
    state.update_well_var("OP_2", "WVIT", 1234.5);
    state.update_well_var("OP_2", "WWCT", 0.0);
    state.update_well_var("OP_2", "WGOR", 0.0);
    state.update_well_var("OP_2", "WBHP", 400.6);
    state.update_well_var("OP_2", "WTHP", 234.5);
    state.update_well_var("OP_2", "WOPTH", 0.0);
    state.update_well_var("OP_2", "WWPTH", 0.0);
    state.update_well_var("OP_2", "WGPTH", 0.0);
    state.update_well_var("OP_2", "WWITH", 1515.0);
    state.update_well_var("OP_2", "WGITH", 3030.0);
    state.update_well_var("OP_2", "WGVIR", 1234.0);
    state.update_well_var("OP_2", "WWVIR", 4321.0);

    state.update_well_var("OP_3", "WOPR", 11.0);
    state.update_well_var("OP_3", "WWPR", 12.0);
    state.update_well_var("OP_3", "WGPR", 13.0);
    state.update_well_var("OP_3", "WVPR", 14.0);
    state.update_well_var("OP_3", "WOPT", 110.0);
    state.update_well_var("OP_3", "WWPT", 120.0);
    state.update_well_var("OP_3", "WGPT", 130.0);
    state.update_well_var("OP_3", "WVPT", 140.0);
    state.update_well_var("OP_3", "WWIR", 0.0);
    state.update_well_var("OP_3", "WGIR", 0.0);
    state.update_well_var("OP_3", "WWIT", 0.0);
    state.update_well_var("OP_3", "WGIT", 0.0);
    state.update_well_var("OP_3", "WVIT", 0.0);
    state.update_well_var("OP_3", "WWCT", 0.0625);
    state.update_well_var("OP_3", "WGOR", 1234.5);
    state.update_well_var("OP_3", "WBHP", 314.15);
    state.update_well_var("OP_3", "WTHP", 246.9);
    state.update_well_var("OP_3", "WOPTH", 2345.6);
    state.update_well_var("OP_3", "WWPTH", 3456.7);
    state.update_well_var("OP_3", "WGPTH", 4567.8);
    state.update_well_var("OP_3", "WWITH", 0.0);
    state.update_well_var("OP_3", "WGITH", 0.0);
    state.update_well_var("OP_3", "WGVIR", 0.0);
    state.update_well_var("OP_3", "WWVIR", 43.21);

    state.update_group_var("OP", "GOPR" ,     110.0);
    state.update_group_var("OP", "GWPR" ,     120.0);
    state.update_group_var("OP", "GGPR" ,     130.0);
    state.update_group_var("OP", "GVPR" ,     140.0);
    state.update_group_var("OP", "GOPT" ,    1100.0);
    state.update_group_var("OP", "GWPT" ,    1200.0);
    state.update_group_var("OP", "GGPT" ,    1300.0);
    state.update_group_var("OP", "GVPT" ,    1400.0);
    state.update_group_var("OP", "GWIR" , -   256.0);
    state.update_group_var("OP", "GGIR" , - 65536.0);
    state.update_group_var("OP", "GWIT" ,   31415.9);
    state.update_group_var("OP", "GGIT" ,   27182.8);
    state.update_group_var("OP", "GVIT" ,   44556.6);
    state.update_group_var("OP", "GWCT" ,       0.625);
    state.update_group_var("OP", "GGOR" ,    1234.5);
    state.update_group_var("OP", "GGVIR",     123.45);
    state.update_group_var("OP", "GWVIR",    1234.56);
    state.update_group_var("OP", "GOPTH",    5678.90);
    state.update_group_var("OP", "GWPTH",    6789.01);
    state.update_group_var("OP", "GGPTH",    7890.12);
    state.update_group_var("OP", "GWITH",    8901.23);
    state.update_group_var("OP", "GGITH",    9012.34);

    state.update("FOPR" ,     1100.0);
    state.update("FWPR" ,     1200.0);
    state.update("FGPR" ,     1300.0);
    state.update("FVPR" ,     1400.0);
    state.update("FOPT" ,    11000.0);
    state.update("FWPT" ,    12000.0);
    state.update("FGPT" ,    13000.0);
    state.update("FVPT" ,    14000.0);
    state.update("FWIR" , -   2560.0);
    state.update("FGIR" , - 655360.0);
    state.update("FWIT" ,   314159.2);
    state.update("FGIT" ,   271828.1);
    state.update("FVIT" ,   445566.77);
    state.update("FWCT" ,        0.625);
    state.update("FGOR" ,     1234.5);
    state.update("FOPTH",    56789.01);
    state.update("FWPTH",    67890.12);
    state.update("FGPTH",    78901.23);
    state.update("FWITH",    89012.34);
    state.update("FGITH",    90123.45);
    state.update("FGVIR",     1234.56);
    state.update("FWVIR",    12345.67);

    return state;
}

struct Setup
{
    EclipseState es;
    const EclipseGrid& grid;
    Schedule schedule;
    SummaryConfig summary_config;

    explicit Setup(const char* path)
        : Setup { Parser{}.parseFile(path) }
    {}

    explicit Setup(const Deck& deck)
        : es             { deck }
        , grid           { es.getInputGrid() }
        , schedule       { deck, es, std::make_shared<Python>() }
        , summary_config { deck, schedule, es.fieldProps(), es.aquifer() }
    {
        es.getIOConfig().setEclCompatibleRST(false);
    }
};

RestartValue
first_sim(const Setup&         setup,
          const Action::State& action_state,
          SummaryState&        st,
          UDQState&            udq_state,
          bool                 write_double)
{
    WellTestState wtest_state;
    EclipseIO eclWriter(setup.es, setup.grid, setup.schedule, setup.summary_config);

    const auto num_cells = setup.grid.getNumActive( );
    const int report_step = 1;
    const auto start_time = setup.schedule.getStartTime();
    const auto first_step = setup.schedule.simTime(report_step);

    const auto sol = mkSolution(num_cells);
    const auto wells = mkWells();
    const auto groups = mkGroups();
    const auto& udq = setup.schedule.getUDQConfig(report_step);
    auto segmentMatcherFactory = []() { return std::make_unique<SegmentMatcher>(ScheduleState{}); };
    auto regionSetMatcherFactory = []() { return std::make_unique<RegionSetMatcher>(FIPRegionStatistics {}); };

    udq.eval(report_step,
             setup.schedule.wellMatcher(report_step),
             setup.schedule[report_step].group_order(),
             segmentMatcherFactory,
             regionSetMatcherFactory,
             st, udq_state);

    RestartValue restart_value(sol, wells, groups, {});
    eclWriter.writeTimeStep(action_state,
                            wtest_state,
                            st,
                            udq_state,
                            report_step,
                            false,
                            std::difftime(first_step, start_time),
                            restart_value,
                            write_double);

    return restart_value;
}

RestartValue
second_sim(const Setup&                   setup,
           Action::State&                 action_state,
           SummaryState&                  summary_state,
           const std::vector<RestartKey>& solution_keys)
{
    EclipseIO writer(setup.es, setup.grid, setup.schedule, setup.summary_config);
    return writer.loadRestart(action_state, summary_state, solution_keys);
}

void compare(const RestartValue&            fst,
             const RestartValue&            snd,
             const std::vector<RestartKey>& solution_keys)
{
    for (const auto& value : solution_keys) {
        auto tol = 0.00001;
        const auto& key = value.key;

        if (key == "TEMP") {
            tol *= 10.0;
        }

        auto first = fst.solution.data<double>(key).begin();
        auto second = snd.solution.data<double>(key).begin();

        for (; first != fst.solution.data<double>(key).end(); ++first, ++second) {
            BOOST_CHECK_CLOSE(*first, *second, tol);
        }
    }
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(EclipseReadWriteWellStateData)
{
    const std::vector<RestartKey> keys {
        {"PRESSURE" , UnitSystem::measure::pressure},
        {"SWAT" , UnitSystem::measure::identity},
        {"SGAS" , UnitSystem::measure::identity},
        {"TEMP" , UnitSystem::measure::temperature},
    };

    WorkArea test_area("test_restart");
    test_area.copyIn("BASE_SIM.DATA");
    test_area.copyIn("RESTART_SIM.DATA");

    Setup base_setup("BASE_SIM.DATA");
    auto st = sim_state(base_setup.schedule);
    Action::State action_state;
    UDQState udq_state(19);
    const auto state1 = first_sim( base_setup , action_state, st, udq_state, false );

    Setup restart_setup("RESTART_SIM.DATA");
    const auto state2 = second_sim( restart_setup , action_state, st , keys );
    compare(state1, state2 , keys);

    BOOST_CHECK_THROW( second_sim( restart_setup, action_state, st, {{"SOIL", UnitSystem::measure::pressure}} ) , std::runtime_error );
    BOOST_CHECK_THROW( second_sim( restart_setup, action_state, st, {{"SOIL", UnitSystem::measure::pressure, true}}) , std::runtime_error );
}


BOOST_AUTO_TEST_CASE(ECL_FORMATTED)
{
    namespace OS = ::Opm::EclIO::OutputStream;

    WorkArea test_area("test_Restart");
    test_area.copyIn("BASE_SIM.DATA");

    Setup base_setup("BASE_SIM.DATA");
    auto& io_config = base_setup.es.getIOConfig();
    {
        auto num_cells = base_setup.grid.getNumActive( );
        auto cells = mkSolution( num_cells );
        auto wells = mkWells();
        auto groups = mkGroups();
        auto sumState = sim_state(base_setup.schedule);
        auto udqState = UDQState{1};
        auto aquiferData = std::optional<Opm::RestartIO::Helpers::AggregateAquiferData>{std::nullopt};
        Action::State action_state;
        WellTestState wtest_state;
        {
            RestartValue restart_value(cells, wells, groups, {});

            io_config.setEclCompatibleRST( false );
            restart_value.addExtra("EXTRA", UnitSystem::measure::pressure, std::vector<double>{10.0,1.0,2.0,3.0});

            const auto outputDir = test_area.currentWorkingDirectory();

            {
                const auto seqnum = 1;
                auto rstFile = OS::Restart {
                    OS::ResultSet{ outputDir, "OPM_FILE" }, seqnum,
                    OS::Formatted{ false }, OS::Unified{ true }
                };

                RestartIO::save(rstFile, seqnum,
                                100,
                                restart_value,
                                base_setup.es,
                                base_setup.grid,
                                base_setup.schedule,
                                action_state,
                                wtest_state,
                                sumState,
                                udqState,
                                aquiferData,
                                true);
            }

            {
                const auto rstFile = ::Opm::EclIO::OutputStream::
                    outputFileName({outputDir, "OPM_FILE"}, "UNRST");

                EclIO::ERst rst{ rstFile };

                BOOST_CHECK_MESSAGE(rst.hasKey("SWAT"), "Restart file must have SWAT vector");
                BOOST_CHECK_MESSAGE(rst.hasKey("EXTRA"), "Restart file must have EXTRA vector");
            }

            io_config.setEclCompatibleRST( true );
            {
                const auto seqnum = 1;
                auto rstFile = OS::Restart {
                    OS::ResultSet{ outputDir, "ECL_FILE" }, seqnum,
                    OS::Formatted{ false }, OS::Unified{ true }
                };

                RestartIO::save(rstFile, seqnum,
                                100,
                                restart_value,
                                base_setup.es,
                                base_setup.grid,
                                base_setup.schedule,
                                action_state,
                                wtest_state,
                                sumState,
                                udqState,
                                aquiferData,
                                true);
            }

            {
                const auto rstFile = ::Opm::EclIO::OutputStream::
                    outputFileName({outputDir, "ECL_FILE"}, "UNRST");

                EclIO::ERst rst{ rstFile };

                BOOST_CHECK_MESSAGE(rst.hasKey("SWAT"), "Restart file must have SWAT vector");
                BOOST_CHECK_MESSAGE(!rst.hasKey("EXTRA"), "Restart file must NOT have EXTRA vector");
            }
        }
    }
}


BOOST_AUTO_TEST_CASE(ECL_LGRFORMATTED)
{
    namespace OS = ::Opm::EclIO::OutputStream;

    WorkArea test_area("test_Restart");
    test_area.copyIn("LGR_BASESIM2WELLS.DATA");

    Setup base_setup("LGR_BASESIM2WELLS.DATA");
    auto& io_config = base_setup.es.getIOConfig();
    {
        const auto& lgr_labels = base_setup.grid.get_all_lgr_labels();
        auto num_lgr_cells = lgr_labels.size();
        std::vector <std::size_t> num_cells(num_lgr_cells+1);
        num_cells[0] = base_setup.grid.getNumActive();
        std::transform(lgr_labels.begin(), lgr_labels.end(),
                        num_cells.begin() + 1 , [&base_setup](const std::string& lgr_tag) {
                                                return base_setup.grid.getLGRCell(lgr_tag).getNumActive();})
                                                ;

        std::vector<data::Solution> cells(num_lgr_cells+1);
        std::transform(num_cells.begin(), num_cells.end(), cells.begin(),
                        [](int n) { return mkSolution(n); });

        std::vector<data::Wells> wells(num_lgr_cells+1);
        wells.push_back(mkWellsLGR_Global());
        wells.push_back(mkWellsLGR_LGR1());
        wells.push_back(mkWellsLGR_LGR2());

        auto groups = mkGroups();
        auto sumState = sim_stateLGR(base_setup.schedule);
        auto udqState = UDQState{1};
        auto aquiferData = std::optional<Opm::RestartIO::Helpers::AggregateAquiferData>{std::nullopt};
        Action::State action_state;
        WellTestState wtest_state;
        {

            std::vector<RestartValue> restart_value;

            for (std::size_t i = 0; i < num_lgr_cells + 1; ++i) {
                restart_value.emplace_back(cells[i], wells[i], groups, data::Aquifers{});
            }


            io_config.setEclCompatibleRST( false );
            std::transform(restart_value.begin(), restart_value.end(),
                            restart_value.begin(),
                            [](RestartValue& rv) {
                                rv.addExtra("EXTRA", UnitSystem::measure::pressure, std::vector<double>{10.0,1.0,2.0,3.0});
                                return rv;
                            });

            const auto outputDir = test_area.currentWorkingDirectory();

            {
                const auto seqnum = 1;
                auto rstFile = OS::Restart {
                    OS::ResultSet{ outputDir, "OPM_FILE" }, seqnum,
                    OS::Formatted{ false }, OS::Unified{ true }
                };

                RestartIO::save(rstFile, seqnum,
                                100,
                                restart_value,
                                base_setup.es,
                                base_setup.grid,
                                base_setup.schedule,
                                action_state,
                                wtest_state,
                                sumState,
                                udqState,
                                aquiferData,
                                true);
            }

            {
                const auto rstFile = ::Opm::EclIO::OutputStream::
                    outputFileName({outputDir, "OPM_FILE"}, "UNRST");

                EclIO::ERst rst{ rstFile };

                BOOST_CHECK_MESSAGE(rst.hasKey("SWAT"), "Restart file must have SWAT vector");
                BOOST_CHECK_MESSAGE(rst.hasKey("EXTRA"), "Restart file must have EXTRA vector");
            }

            io_config.setEclCompatibleRST( true );
            {
                const auto seqnum = 1;
                auto rstFile = OS::Restart {
                    OS::ResultSet{ outputDir, "ECL_FILE" }, seqnum,
                    OS::Formatted{ false }, OS::Unified{ true }
                };

                RestartIO::save(rstFile, seqnum,
                                100,
                                restart_value,
                                base_setup.es,
                                base_setup.grid,
                                base_setup.schedule,
                                action_state,
                                wtest_state,
                                sumState,
                                udqState,
                                aquiferData,
                                true);
            }

            {
                const auto rstFile = ::Opm::EclIO::OutputStream::
                    outputFileName({outputDir, "ECL_FILE"}, "UNRST");

                EclIO::ERst rst{ rstFile };

                BOOST_CHECK_MESSAGE(rst.hasKey("SWAT"), "Restart file must have SWAT vector");
                BOOST_CHECK_MESSAGE(!rst.hasKey("EXTRA"), "Restart file must NOT have EXTRA vector");
            }
        }
    }
}
