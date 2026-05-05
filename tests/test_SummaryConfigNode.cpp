/*
  Copyright 2026 Equinor

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

#define BOOST_TEST_MODULE SummaryConfigNode

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/io/eclipse/SummaryNode.hpp>

#include <initializer_list>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace Opm::EclIO {
    template <typename CharT, typename Traits>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
               const SummaryNode::Type            t)
    {
        switch (t) {
        case SummaryNode::Type::Rate:
            return os << "Rate";

        case SummaryNode::Type::Total:
            return os << "Total";

        case SummaryNode::Type::Ratio:
            return os << "Ratio";

        case SummaryNode::Type::Pressure:
            return os << "Pressure";

        case SummaryNode::Type::Count:
            return os << "Count";

        case SummaryNode::Type::Mode:
            return os << "Mode";

        case SummaryNode::Type::ProdIndex:
            return os << "ProdIndex";

        case SummaryNode::Type::Undefined:
            return os << "Undefined";
        }

        return os << "<Unknown ("
                  << static_cast<std::underlying_type_t<EclIO::SummaryNode::Type>>(t)
                  << ")>";
    }
} // namespace Opm::EclIO

BOOST_AUTO_TEST_SUITE(ParseKeywords)

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Type)

BOOST_AUTO_TEST_SUITE(Total)

BOOST_AUTO_TEST_CASE(BOPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("BOPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(COPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("COPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(FOPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("FOPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(GOPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("GOPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(ROPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("ROPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(WOPT)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("WOPT"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_CASE(CGMITL)
{
    BOOST_CHECK_EQUAL(Opm::parseKeywordType("CGMITL"),
                      Opm::SummaryConfigNode::Type::Total);
}

BOOST_AUTO_TEST_SUITE_END()     // Total

BOOST_AUTO_TEST_SUITE_END()     // Type

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Category)

BOOST_AUTO_TEST_CASE(LWOPR_is_Well)
{
    using Cat = Opm::SummaryConfigNode::Category;
    BOOST_CHECK(Opm::parseKeywordCategory("LWOPR") == Cat::Well);
}

BOOST_AUTO_TEST_CASE(LCOFR_is_Connection)
{
    using Cat = Opm::SummaryConfigNode::Category;
    BOOST_CHECK(Opm::parseKeywordCategory("LCOFR") == Cat::Connection);
}

BOOST_AUTO_TEST_CASE(LBPR_is_Block)
{
    using Cat = Opm::SummaryConfigNode::Category;
    BOOST_CHECK(Opm::parseKeywordCategory("LBPR") == Cat::Block);
}

BOOST_AUTO_TEST_CASE(WOPR_is_Well)
{
    using Cat = Opm::SummaryConfigNode::Category;
    BOOST_CHECK(Opm::parseKeywordCategory("WOPR") == Cat::Well);
}

BOOST_AUTO_TEST_CASE(unknown_L_is_Misc)
{
    using Cat = Opm::SummaryConfigNode::Category;
    BOOST_CHECK(Opm::parseKeywordCategory("LXYZ") == Cat::Miscellaneous);
}

BOOST_AUTO_TEST_SUITE_END()     // Category

BOOST_AUTO_TEST_SUITE_END()     // ParseKeywords

// =====================================================================

BOOST_AUTO_TEST_SUITE(LGR)

BOOST_AUTO_TEST_CASE(lgr_field_roundtrip)
{
    using Node = Opm::SummaryConfigNode;
    using Cat  = Opm::SummaryConfigNode::Category;

    // Build a node and attach LGR info
    auto node = Node("LWOPR", Cat::Well, Opm::KeywordLocation{});
    node.lgr(Opm::EclIO::lgr_info{ "LGR1", {0, 0, 0} });

    BOOST_CHECK(node.lgr().has_value());
    BOOST_CHECK_EQUAL(node.lgr()->name,   "LGR1");
    BOOST_CHECK_EQUAL(node.lgr()->ijk[0], 0);
    BOOST_CHECK_EQUAL(node.lgr()->ijk[1], 0);
    BOOST_CHECK_EQUAL(node.lgr()->ijk[2], 0);

    // Conversion to EclIO::SummaryNode preserves lgr
    Opm::EclIO::SummaryNode sn = node;
    BOOST_CHECK(sn.lgr.has_value());
    BOOST_CHECK_EQUAL(sn.lgr->name,   "LGR1");
    BOOST_CHECK_EQUAL(sn.lgr->ijk[0], 0);

    // Non-LGR node has empty optional
    auto global = Node("WOPR", Cat::Well, Opm::KeywordLocation{});
    BOOST_CHECK(!global.lgr().has_value());
    Opm::EclIO::SummaryNode sn2 = global;
    BOOST_CHECK(!sn2.lgr.has_value());

    // LB* node with specific ijk
    auto blk = Node("LBPR", Cat::Block, Opm::KeywordLocation{});
    blk.lgr(Opm::EclIO::lgr_info{ "LGR1", {2, 3, 4} });
    Opm::EclIO::SummaryNode sn3 = blk;
    BOOST_CHECK_EQUAL(sn3.lgr->ijk[0], 2);
    BOOST_CHECK_EQUAL(sn3.lgr->ijk[1], 3);
    BOOST_CHECK_EQUAL(sn3.lgr->ijk[2], 4);
}

BOOST_AUTO_TEST_CASE(lgr_dedup_distinct_lgrs)
{
    // Regression: two nodes with same keyword+well but different LGRs
    // must NOT compare equal (silent dedup drop prevention)
    using Node = Opm::SummaryConfigNode;
    using Cat  = Opm::SummaryConfigNode::Category;

    auto n1 = Node("LWOPR", Cat::Well, Opm::KeywordLocation{});
    n1.namedEntity("PROD1");
    n1.lgr(Opm::EclIO::lgr_info{ "LGR1", {0, 0, 0} });

    auto n2 = Node("LWOPR", Cat::Well, Opm::KeywordLocation{});
    n2.namedEntity("PROD1");
    n2.lgr(Opm::EclIO::lgr_info{ "LGR2", {0, 0, 0} });

    BOOST_CHECK(n1 != n2);
    BOOST_CHECK(n1 < n2 || n2 < n1); // strict weak ordering holds

    // Same LGR — must be equal
    auto n3 = Node("LWOPR", Cat::Well, Opm::KeywordLocation{});
    n3.namedEntity("PROD1");
    n3.lgr(Opm::EclIO::lgr_info{ "LGR1", {0, 0, 0} });

    BOOST_CHECK(n1 == n3);
    BOOST_CHECK(!(n1 < n3) && !(n3 < n1));
}

BOOST_AUTO_TEST_CASE(lgr_no_lgr_sorts_before_lgr)
{
    // non-LGR node sorts before LGR node with same keyword+entity
    using Node = Opm::SummaryConfigNode;
    using Cat  = Opm::SummaryConfigNode::Category;

    auto global = Node("WOPR", Cat::Well, Opm::KeywordLocation{});
    global.namedEntity("PROD1");

    auto lgr = Node("WOPR", Cat::Well, Opm::KeywordLocation{});
    lgr.namedEntity("PROD1");
    lgr.lgr(Opm::EclIO::lgr_info{ "LGR1", {0, 0, 0} });

    BOOST_CHECK(global < lgr);
    BOOST_CHECK(!(lgr < global));
}

BOOST_AUTO_TEST_SUITE_END() // LGR

// =====================================================================

BOOST_AUTO_TEST_SUITE(NoSumLgr)

BOOST_AUTO_TEST_CASE(default_is_false)
{
    // A default-constructed SummaryConfig has noSumLgr_ == false
    Opm::SummaryConfig sc;
    BOOST_CHECK(!sc.noSumLgr());
}

BOOST_AUTO_TEST_CASE(serialization_test_object_sets_flag)
{
    // serializationTestObject() sets noSumLgr_ = true to exercise serialization
    const auto sc = Opm::SummaryConfig::serializationTestObject();
    BOOST_CHECK(sc.noSumLgr());
}

BOOST_AUTO_TEST_SUITE_END() // NoSumLgr

// =====================================================================
// Deck-parsing integration tests using the two LGR .DATA files.
// Exercises keywordLW / keywordLC / keywordLB + handleKW L-prefix dispatch
// added in PR01. Tests that SummaryConfig correctly populates lgr_ on nodes
// produced from real deck keywords.
// =====================================================================

BOOST_AUTO_TEST_SUITE(LGR_SummaryConfig_Deck)

namespace {

Opm::SummaryConfig makeSummaryConfig(const std::string& path)
{
    const auto deck  = Opm::Parser{}.parseFile(path);
    const auto es    = Opm::EclipseState { deck };
    const auto sched = Opm::Schedule { deck, es, std::make_shared<Opm::Python>() };
    return Opm::SummaryConfig { deck, sched, es.fieldProps(), es.aquifer() };
}

} // anonymous namespace

// --- 1LGR deck: INJ in LGR1, PROD in global grid ---

BOOST_AUTO_TEST_CASE(lgr_1lgr_lw_node_populated)
{
    // LWBHP 'LGR1' 'INJ' / — should produce one Well node with lgr_.name=="LGR1"
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-1LGR.DATA");
    const auto hits = cfg.keywords("LWBHP");

    BOOST_REQUIRE_EQUAL(hits.size(), 1U);
    BOOST_REQUIRE(hits[0].lgr().has_value());
    BOOST_CHECK_EQUAL(hits[0].lgr()->name,   "LGR1");
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[0], 0);   // LW* sentinel: whole well
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[1], 0);
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[2], 0);
    BOOST_CHECK_EQUAL(hits[0].namedEntity(), "INJ");
    BOOST_CHECK(hits[0].category() == Opm::SummaryConfigNode::Category::Well);
}

BOOST_AUTO_TEST_CASE(lgr_1lgr_lc_node_populated)
{
    // LCGFR 'LGR1' 'INJ' 2 2 1 / — should produce one Connection node
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-1LGR.DATA");
    const auto hits = cfg.keywords("LCGFR");

    BOOST_REQUIRE_EQUAL(hits.size(), 1U);
    BOOST_REQUIRE(hits[0].lgr().has_value());
    BOOST_CHECK_EQUAL(hits[0].lgr()->name,   "LGR1");
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[0], 2);   // explicit coords from deck
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[1], 2);
    BOOST_CHECK_EQUAL(hits[0].lgr()->ijk[2], 1);
    BOOST_CHECK_EQUAL(hits[0].namedEntity(), "INJ");
    BOOST_CHECK(hits[0].category() == Opm::SummaryConfigNode::Category::Connection);
}

BOOST_AUTO_TEST_CASE(lgr_1lgr_lb_node_populated)
{
    // LBPR 'LGR1' 2 2 1 / — should produce one Block node
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-1LGR.DATA");
    const auto hits = cfg.keywords("LBPR");

    // 1LGR deck requests LBPR for multiple cells; at minimum one at (2,2,1)
    BOOST_CHECK(!hits.empty());
    const auto it = std::find_if(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.lgr().has_value()
                && n.lgr()->name   == "LGR1"
                && n.lgr()->ijk[0] == 2
                && n.lgr()->ijk[1] == 2
                && n.lgr()->ijk[2] == 1;
        });
    BOOST_CHECK(it != hits.end());
    BOOST_CHECK(it->category() == Opm::SummaryConfigNode::Category::Block);
}

BOOST_AUTO_TEST_CASE(lgr_1lgr_global_well_has_no_lgr)
{
    // PROD is on the global grid — its WOPR node must have no lgr (regression)
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-1LGR.DATA");
    const auto hits = cfg.keywords("WOPR");

    const auto it = std::find_if(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.namedEntity() == "PROD";
        });
    BOOST_REQUIRE(it != hits.end());
    BOOST_CHECK(!it->lgr().has_value());
}

// --- 2LGR deck: INJ in LGR1, PROD in LGR2 — Fix4 dedup stress ---

BOOST_AUTO_TEST_CASE(lgr_2lgr_lw_two_nodes_not_deduped)
{
    // LWBHP has 'LGR1'/'INJ' AND 'LGR2'/'PROD' — must produce 2 distinct nodes
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-2LGR.DATA");
    const auto hits = cfg.keywords("LWBHP");

    BOOST_REQUIRE_EQUAL(hits.size(), 2U);
    BOOST_CHECK(hits[0] != hits[1]);

    // One for each LGR
    const bool has_lgr1 = std::any_of(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.lgr().has_value() && n.lgr()->name == "LGR1"
                && n.namedEntity() == "INJ";
        });
    const bool has_lgr2 = std::any_of(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.lgr().has_value() && n.lgr()->name == "LGR2"
                && n.namedEntity() == "PROD";
        });
    BOOST_CHECK(has_lgr1);
    BOOST_CHECK(has_lgr2);
}

BOOST_AUTO_TEST_CASE(lgr_2lgr_lb_dedup_by_lgr_name)
{
    // LBPR 'LGR1' 2 2 1 / AND 'LGR2' 2 2 1 /
    // Same keyword + same ijk — ONLY lgr_.name differs.
    // Without lgr_ in Block operator==, second node is silently dropped.
    const auto cfg  = makeSummaryConfig("LGR-WELL-3X3-2LGR.DATA");
    const auto hits = cfg.keywords("LBPR");

    // Expect at least the two (2,2,1) nodes — one per LGR
    const auto lgr1_221 = std::count_if(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.lgr().has_value() && n.lgr()->name == "LGR1"
                && n.lgr()->ijk[0] == 2 && n.lgr()->ijk[1] == 2 && n.lgr()->ijk[2] == 1;
        });
    const auto lgr2_221 = std::count_if(hits.begin(), hits.end(),
        [](const Opm::SummaryConfigNode& n) {
            return n.lgr().has_value() && n.lgr()->name == "LGR2"
                && n.lgr()->ijk[0] == 2 && n.lgr()->ijk[1] == 2 && n.lgr()->ijk[2] == 1;
        });
    BOOST_CHECK_EQUAL(lgr1_221, 1);
    BOOST_CHECK_EQUAL(lgr2_221, 1);
}

BOOST_AUTO_TEST_SUITE_END() // LGR_SummaryConfig_Deck

// =====================================================================
// Inline schema tests — parse LW*/LC*/LB* from a minimal string deck.
// These catch schema bugs (wrong "size", missing items, bad regex) without
// needing .DATA files on disk. Any keyword-recognition or record-format
// regression shows up here before Jenkins.
// =====================================================================

BOOST_AUTO_TEST_SUITE(LGR_Schema_Inline)

namespace {

// Parse a minimal complete deck from an inline string and return SummaryConfig.
// The RUNSPEC/GRID/PROPS/SCHEDULE sections are minimal stubs so EclipseState
// can be constructed without needing any external files.
Opm::SummaryConfig parseSummarySection(const std::string& summary_body)
{
    const std::string deck_str =
        "RUNSPEC\n"
        "TITLE\n"
        " LGR_SCHEMA_TEST /\n"
        "DIMENS\n"
        " 3 3 1 /\n"
        "WELLDIMS\n"
        " 4 1 1 4 /\n"
        "START\n"
        " 1 JAN 2020 /\n"
        "GRID\n"
        "DXV\n"
        " 3*100.0 /\n"
        "DYV\n"
        " 3*100.0 /\n"
        "DZV\n"
        " 1*10.0 /\n"
        "DEPTHZ\n"
        " 16*2000.0 /\n"
        "PORO\n"
        " 9*0.3 /\n"
        "PROPS\n"
        "SOLUTION\n"
        "SUMMARY\n"
        + summary_body +
        "\nSCHEDULE\n"
        // WELSPECL: WELL GROUP LGR I J DEPTH PHASE ...
        "WELSPECL\n"
        " 'WELL1'  'G'  'LGR1'  2 2  2000.0  'OIL' /\n"
        " 'WELL2'  'G'  'LGR2'  3 3  2000.0  'OIL' /\n"
        " 'INJ'    'G'  'LGR1'  1 1  2000.0  'WAT' /\n"
        " 'PROD'   'G'  'LGR2'  3 3  2000.0  'OIL' /\n"
        "/\n"
        "DATES\n"
        " 1 FEB 2020 /\n"
        "/\n"
        "END\n";

    const auto deck  = Opm::Parser{}.parseString(deck_str);
    const auto es    = Opm::EclipseState { deck };
    const auto sched = Opm::Schedule { deck, es, std::make_shared<Opm::Python>() };
    return Opm::SummaryConfig { deck, sched, es.fieldProps(), es.aquifer() };
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(lw_single_record_produces_node)
{
    // Verify that LWBHP with one LGR/well record produces exactly one node.
    // Catches: "size":1 regression (would throw before reaching keywordLW),
    // AND the lgr_well_tag matching (would silently produce zero nodes).
    const auto cfg = parseSummarySection(
        "LWBHP\n"
        " 'LGR1'  'WELL1' /\n"
        "/\n"
    );

    const auto nodes = cfg.keywords("LWBHP");
    BOOST_REQUIRE_EQUAL(nodes.size(), 1u);
    BOOST_CHECK_EQUAL(nodes[0].namedEntity(), "WELL1");
    BOOST_REQUIRE(nodes[0].lgr().has_value());
    BOOST_CHECK_EQUAL(nodes[0].lgr()->name, "LGR1");
}

BOOST_AUTO_TEST_CASE(lw_multi_record_produces_two_nodes)
{
    // Verify that LWBHP with two records (LGR1/WELL1, LGR2/WELL2) produces two nodes.
    // Catches: "size":1 bug — schema must be table format (no "size" field).
    const auto cfg = parseSummarySection(
        "LWBHP\n"
        " 'LGR1'  'WELL1' /\n"
        " 'LGR2'  'WELL2' /\n"
        "/\n"
    );

    const auto nodes = cfg.keywords("LWBHP");
    BOOST_REQUIRE_EQUAL(nodes.size(), 2u);
    BOOST_CHECK_EQUAL(nodes[0].lgr()->name, "LGR1");
    BOOST_CHECK_EQUAL(nodes[1].lgr()->name, "LGR2");
}

BOOST_AUTO_TEST_CASE(lc_multi_record_produces_nodes)
{
    // Verify that LCOFR with multiple records (LGR + well + I J K) produces nodes.
    const auto cfg = parseSummarySection(
        "LCOFR\n"
        " 'LGR1'  'WELL1'  1 1 1 /\n"
        " 'LGR1'  'WELL1'  2 1 1 /\n"
        "/\n"
    );

    const auto nodes = cfg.keywords("LCOFR");
    BOOST_CHECK_EQUAL(nodes.size(), 2u);
    for (const auto& node : nodes) {
        BOOST_CHECK(node.lgr().has_value());
        BOOST_CHECK_EQUAL(node.lgr()->name, "LGR1");
    }
}

BOOST_AUTO_TEST_CASE(lb_multi_record_produces_nodes)
{
    // Verify that LBPR with multiple records (LGR + I J K) produces nodes.
    // keywordLB has no well lookup — tests pure schema parsing.
    const auto cfg = parseSummarySection(
        "LBPR\n"
        " 'LGR1'  1 1 1 /\n"
        " 'LGR1'  2 2 1 /\n"
        "/\n"
    );

    const auto nodes = cfg.keywords("LBPR");
    BOOST_REQUIRE_EQUAL(nodes.size(), 2u);
    BOOST_CHECK(nodes[0].lgr().has_value());
    BOOST_CHECK_EQUAL(nodes[0].lgr()->name, "LGR1");
}

BOOST_AUTO_TEST_CASE(lw_regex_matches_lwstat)
{
    // Verify the LW.+ regex covers LWSTAT and produces a node.
    const auto cfg = parseSummarySection(
        "LWSTAT\n"
        " 'LGR1'  'INJ' /\n"
        "/\n"
    );

    const auto nodes = cfg.keywords("LWSTAT");
    BOOST_REQUIRE_EQUAL(nodes.size(), 1u);
    BOOST_CHECK_EQUAL(nodes[0].lgr()->name, "LGR1");
}

BOOST_AUTO_TEST_SUITE_END() // LGR_Schema_Inline
