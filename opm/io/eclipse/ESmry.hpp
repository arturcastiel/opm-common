/*
   Copyright 2019 Equinor ASA.

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

#ifndef OPM_IO_ESMRY_HPP
#define OPM_IO_ESMRY_HPP

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>
#include <map>
#include <stdint.h>

#include <opm/common/utility/TimeService.hpp>
#include <opm/io/eclipse/SummaryNode.hpp>

namespace Opm { namespace EclIO {

using ArrSourceEntry = std::tuple<std::string, std::string, int, uint64_t>;
using TimeStepEntry = std::tuple<int, int, uint64_t>;
using RstEntry = std::tuple<std::string, int>;

class ESmry
{
public:

    // input is smspec (or fsmspec file)
    explicit ESmry(const std::string& filename, bool loadBaseRunData=false);

    int numberOfVectors() const { return nVect; }

    bool hasKey(const std::string& key) const;

    const std::vector<float>& get(const std::string& name) const;
    const std::vector<float>& get(const SummaryNode& node) const;
    std::vector<time_point> dates() const;

    std::vector<float> get_at_rstep(const std::string& name) const;
    std::vector<float> get_at_rstep(const SummaryNode& node) const;
    std::vector<time_point> dates_at_rstep() const;

    void loadData(const std::vector<std::string>& vectList) const;
    void loadData() const;

    bool make_esmry_file();

    time_point startdate() const { return tp_startdat; }
    const std::vector<int>& start_v() const { return start_vect; }

    const std::vector<std::string>& keywordList() const;
    std::vector<std::string> keywordList(const std::string& pattern) const;
    const std::vector<SummaryNode>& summaryNodeList() const;

    int timestepIdxAtReportstepStart(const int reportStep) const;

    size_t numberOfTimeSteps() const { return nTstep; }

    const std::string& get_unit(const std::string& name) const;
    const std::string& get_unit(const SummaryNode& node) const;

    void write_rsm(std::ostream&) const;
    void write_rsm_file(std::optional<std::filesystem::path> = std::nullopt) const;

    bool all_steps_available();
    std::string rootname() { return inputFileName.stem().generic_string(); }
    std::tuple<double, double> get_io_elapsed() const;

private:
    std::filesystem::path inputFileName;
    RstEntry restart_info;

    int nI, nJ, nK, nSpecFiles;
    bool fromSingleRun;
    size_t nVect, nTstep;

    std::vector<bool> formattedFiles;
    std::vector<std::string> dataFileList;
    mutable std::vector<std::vector<float>> vectorData;
    mutable std::vector<bool> vectorLoaded;
    std::vector<TimeStepEntry> timeStepList;
    std::vector<TimeStepEntry> miniStepList;
    std::vector<std::map<int, int>> arrayPos;
    std::vector<std::string> keyword;
    std::map<std::string, int> keyword_index;
    std::vector<int> nParamsSpecFile;

    std::vector<std::vector<std::string>> keywordListSpecFile;

    std::vector<int> seqIndex;
    std::vector<int> mini_steps;

    std::vector<std::string> ignore_keyword_list = {"TNAVHEAD", "TNAVTIME"};

    void ijk_from_global_index(int glob, int &i, int &j, int &k) const;

    std::vector<SummaryNode> summaryNodes;
    std::unordered_map<std::string, std::string> kwunits;

    time_point tp_startdat;
    std::vector<int> start_vect;

    mutable double m_io_opening;
    mutable double m_io_loading;

    std::vector<std::string> checkForMultipleResultFiles(const std::filesystem::path& rootN, bool formatted) const;

    void getRstString(const std::vector<std::string>& restartArray,
                      std::filesystem::path& pathRst,
                      std::filesystem::path& rootN) const;

    void updatePathAndRootName(std::filesystem::path& dir, std::filesystem::path& rootN) const;


    std::string makeKeyString(const std::string& keyword, const std::string& wgname, int num,
                              const std::optional<Opm::EclIO::lgr_info> lgr_info) const;

    std::string unpackNumber(const SummaryNode&) const;
    std::string lookupKey(const SummaryNode&) const;


    void write_block(std::ostream &, bool write_dates, const std::vector<std::string>& time_column, const std::vector<SummaryNode>&) const;

    template <typename T>
    std::vector<T> rstep_vector(const std::vector<T>& full_vector) const {
        std::vector<T> result;
        result.reserve(seqIndex.size());

        std::transform(seqIndex.begin(), seqIndex.end(),
                       std::back_inserter(result),
                       [&full_vector](const auto& ind)
                       {
                           return full_vector[ind];
                       });

        return result;
    }

    std::vector<std::tuple <std::string, uint64_t>>
    getListOfArrays(const std::string& filename, bool formatted);

    std::vector<int> makeKeywPosVector(int speInd) const;
    std::string read_string_from_disk(std::fstream& fileH, uint64_t size) const;

    void read_ministeps_from_disk();
    int read_ministep_formatted(std::fstream& fileH);
};

}} // namespace Opm::EclIO

inline std::ostream& operator<<(std::ostream& os, const Opm::EclIO::ESmry& smry) {
    smry.write_rsm(os);

    return os;
}

#endif // OPM_IO_ESMRY_HPP
