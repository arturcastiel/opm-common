/*
   Copyright 2020 Equinor ASA.

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM.  If not, see <http://www.gnu.org/licenses/>.
   */

#include <opm/io/eclipse/ESmry.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/output/eclipse/WStat.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <list>
#include <optional>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

    constexpr std::size_t column_width { 8 } ;
    constexpr std::size_t column_space { 5 } ;

    constexpr std::size_t column_count { 10 } ;

    constexpr std::size_t total_column { column_width + column_space } ;
    constexpr std::size_t total_width  { total_column * column_count } ;

    const std::string block_separator_line { } ;

    // the fact that the dashed header line has 127 rather than 130 dashes has no provenance
    const std::string divider_line { std::string(total_width - 3, '-') } ;

    const std::vector<std::string> month_names = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

    std::string format_date(const Opm::time_point& tp) {
        auto ts = Opm::TimeStampUTC( Opm::TimeService::to_time_t(tp) );
        char buffer[12];
        std::snprintf(buffer, 12, "%2d-%3s-%4d", ts.day(), month_names[ts.month() - 1].c_str(), ts.year());
        return std::string(buffer, 11);
    }

    std::string block_header_line(const std::string& run_name) {
        std::string date_string = format_date(Opm::TimeService::from_time_t(std::time(nullptr)));
        return "SUMMARY OF RUN " + run_name + " at: " + date_string;
    }

    void write_line(std::ostream& os, const std::string& line, char prefix = ' ') {
        os << prefix << std::setw(total_width) << std::left << line << '\n';
    }

    void write_padding(std::ostream& os, std::size_t columns) {
        for (std::size_t icol = columns + 1; icol < column_count; icol++)
            os << std::setw(column_width + column_space) << " ";
    }

    void print_text_element(std::ostream& os, const std::string& element) {
        os << std::setw(8) << std::left << element << std::setw(5) << "";
    }

    void print_time_element(std::ostream& os, const std::string& element) {
        os << std::setw(11) << element << std::setw(2) << "";
    }

    std::string format_float_element(float element) {
        std::stringstream ss;
        static const std::regex integer_regex { "\\.0*$" };
        auto element_string = std::to_string(element);
        if (element_string.size() > 8) {
            element_string = element_string.substr(0, 8);
        }
        element_string = std::regex_replace(element_string, integer_regex, "");
        ss << std::setw(8) << std::right << element_string;
        return ss.str();
    }

    void print_float_element(std::ostream& os, float element) {
        os << format_float_element(element) << std::setw(5) << "";
    }

    void write_header_columns(std::ostream& os, const std::string& time_column, const std::vector<Opm::EclIO::SummaryNode>& vectors, std::function<void(std::ostream&, const Opm::EclIO::SummaryNode&)> print_element, char prefix = ' ') {
        os << prefix;

        print_text_element(os, time_column);
        for (const auto& vector : vectors) {
            print_element(os, vector);
        }
        write_padding(os, vectors.size());

        os << '\n';
    }

    std::string convert_wstat(double numeric_wstat) {
        static const std::unordered_map<int, std::string> wstat_map = {
            {Opm::WStat::numeric::UNKNOWN, Opm::WStat::symbolic::UNKNOWN},
            {Opm::WStat::numeric::PROD,    Opm::WStat::symbolic::PROD},
            {Opm::WStat::numeric::INJ,     Opm::WStat::symbolic::INJ},
            {Opm::WStat::numeric::SHUT,    Opm::WStat::symbolic::SHUT},
            {Opm::WStat::numeric::STOP,    Opm::WStat::symbolic::STOP},
            {Opm::WStat::numeric::PSHUT,   Opm::WStat::symbolic::PSHUT},
            {Opm::WStat::numeric::PSTOP,   Opm::WStat::symbolic::PSTOP},
        };
        return wstat_map.at(static_cast<int>(numeric_wstat));
    }

    void write_data_row(std::ostream& os, const std::vector<std::string>& time_column, const std::vector<Opm::EclIO::SummaryNode>& summary_nodes, const std::vector<std::pair<std::vector<float>, int>>& data, std::size_t time_index, char prefix = ' ') {
        os << prefix;

        print_time_element( os, time_column[time_index] );
        for (std::size_t row_index = 0; row_index < data.size(); row_index++) {
            const auto& [time_series, scale_factor] = data[row_index];
            const auto& summary_node = summary_nodes[row_index];

            if (summary_node.keyword == "WSTAT")
                print_text_element(os, convert_wstat(time_series[time_index]));
            else
                print_float_element(os, time_series[time_index] * std::pow(10.0, -scale_factor));
        }
        write_padding(os, data.size());

        os << '\n';
    }

    void write_scale_columns(std::ostream& os,
                             const std::vector<std::pair<std::vector<float>, int>>& data,
                             char prefix = ' ')
    {
        os << prefix;

        print_text_element(os, "");
        for (const auto& [_ , scale_factor] : data) {
            (void)_;
            if (scale_factor) {
                print_text_element(os, "*10**" + std::to_string(scale_factor));
            } else {
                print_text_element(os, "");
            }
        }

        os << '\n';
    }

}

namespace Opm { namespace EclIO {

void ESmry::write_block(std::ostream& os,
                        bool write_dates,
                        const std::vector<std::string>& time_column,
                        const std::vector<SummaryNode>& vectors) const
{
    write_line(os, block_separator_line, '1');
    write_line(os, divider_line);
    write_line(os, block_header_line(inputFileName.stem().generic_string()));
    write_line(os, divider_line);

    std::vector<std::pair<std::vector<float>, int>> data;

    bool has_scale_factors { false } ;
    for (const auto& vector : vectors) {
        const auto& vector_data { this->get(vector) } ;

        auto max = *std::max_element(vector_data.begin(), vector_data.end());
        // log10 for 0 is undefined and log10 for negative values yields nan.
        // We skip the scale factor in these cases to prevent undefined behavior
        int scale_factor {
            max <= 0 ? 0 :
            std::max(0, 3 * static_cast<int>(std::floor(( std::log10(max) - 4 ) / 3 ))) } ;

        // Make sure that 10**scale_factor is less than 13 character
        scale_factor = std::min(99999999, scale_factor);

        if (scale_factor) {
            has_scale_factors = true;
        }

        data.emplace_back(vector_data, scale_factor);
    }

    {
        std::size_t rows { data[0].first.size() };
        std::string time_header = "TIME";
        std::string time_unit = "DAYS";

        if (write_dates) {
            time_header = "DATE";
            time_unit = "";
        }
        write_header_columns(os, time_header, vectors, [](std::ostream& oss, const SummaryNode& node) { print_text_element(oss, node.keyword); });
        write_header_columns(os, time_unit, vectors, [this](std::ostream& oss, const SummaryNode& node) { print_text_element(oss, this->get_unit(node)); });
        if (has_scale_factors) {
            write_scale_columns(os, data);
        }
        write_header_columns(os, "", vectors, [](std::ostream& oss, const SummaryNode& node) { print_text_element(oss, node.display_name().value_or("")); });
        write_header_columns(os, "", vectors, [](std::ostream& oss, const SummaryNode& node) { print_text_element(oss, node.display_number().value_or("")); });

        write_line(os, divider_line);

        for (std::size_t i { 0 } ; i < rows; i++) {
            write_data_row(os, time_column, vectors, data, i);
        }
    }

    os << std::flush;
}

void ESmry::write_rsm(std::ostream& os) const
{
    bool write_dates = false;
    std::vector<SummaryNode> data_vectors;
    std::remove_copy_if(summaryNodes.begin(), summaryNodes.end(), std::back_inserter(data_vectors), [](const SummaryNode& node){
        if (node.keyword == "TIME" || node.keyword == "DAY" || node.keyword == "MONTH" || node.keyword == "YEAR") {
            return true;
        } else {
            return false;
        }
    });

    /*
      Ensure that the YEARS vector is the first in the data_vectors; could in
      principle embark on a more general sorting here.
    */
    if (data_vectors[0].keyword != "YEARS") {
        auto years_iter = std::find_if(data_vectors.begin(), data_vectors.end(), [](const SummaryNode& node) { return (node.keyword == "YEARS"); });
        if (years_iter != data_vectors.end())
            std::swap(data_vectors[0], *years_iter);
    }

    std::vector<std::list<SummaryNode>> data_vector_blocks;
    constexpr std::size_t data_column_count { column_count - 1 } ;
    for (std::size_t i { 0 } ; i < data_vectors.size(); i += data_column_count) {
        auto last = std::min(data_vectors.size(), i + data_column_count);
        data_vector_blocks.emplace_back(data_vectors.begin() + i, data_vectors.begin() + last);
    }

    this->loadData();
    std::vector<std::string> time_column;
    if (this->hasKey("DAY") && this->hasKey("MONTH") && this->hasKey("YEAR")) {
        write_dates = true;
        const auto dates = this->dates();
        std::transform(dates.begin(), dates.end(),
                       std::back_inserter(time_column),
                       [](const auto& t)
                       {
                           return format_date(t);
                       });
    } else {
        const auto& time_data = this->get("TIME");
        std::transform(time_data.begin(), time_data.end(),
                       std::back_inserter(time_column),
                       [](const auto& t)
                       {
                           return format_float_element(t);
                       });
    }

    for (const auto& data_vector_block : data_vector_blocks) {
        write_block(os, write_dates, time_column, {data_vector_block.begin(), data_vector_block.end() });
    }
}

void ESmry::write_rsm_file(std::optional<std::filesystem::path> filename) const
{
    std::filesystem::path summary_file_name { filename.value_or(inputFileName) } ;
    summary_file_name.replace_extension("RSM");

    std::ofstream rsm_file { summary_file_name } ;

    if (!rsm_file.is_open()) {
        OPM_THROW(std::runtime_error, "Could not open file " + summary_file_name.generic_string() + " for writing");
    }

    write_rsm(rsm_file);

    rsm_file.close();
}

}} // namespace Opm::EclIO
