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

#include <opm/io/eclipse/SummaryNode.hpp>

#include <initializer_list>
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
