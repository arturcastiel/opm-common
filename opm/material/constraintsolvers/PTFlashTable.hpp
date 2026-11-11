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
 * \brief An on-demand table of P-T flash results on a uniform mesh of the
 *        state space, in the manner of operator-based linearization.
 *
 * The state (p, T, z) is mapped to the axes (ln p, 1/T, z without one
 * dependent component). Each axis is cut into intervals of equal width, so
 * the state space is a lattice of hypercubes ("cells") with 2^d vertices,
 * d = numComponents + 1. A vertex is flashed exactly the first time a cell
 * that has it is visited, and both vertices and cells are kept in hash maps
 * for the rest of the run. Inside a cell whose vertices are all two-phase,
 * ln K is the multilinear interpolant of the vertex values. From four
 * components on, the d + 1 vertices of the Kuhn simplex that holds the state
 * are used instead (no cell store), as open-DARTS does for many dimensions.
 * Over (ln p, 1/T) the Wilson estimate of ln K is linear, so the
 * interpolation error there comes only from the non-ideal part of ln K.
 *
 * Reference: D.V. Voskov, Operator-based linearization approach for modeling
 * of multiphase multi-component flow in porous media, J. Comput. Phys. 337
 * (2017); M. Khait and D.V. Voskov, Adaptive parameterization for solving of
 * thermal/compositional nonlinear flow and transport with buoyancy, SPE J. 23
 * (2018).
 *
 * PTFlash consults the table when PTFlashTableSettings::get().mode is not Off,
 * for a flash without a split from an earlier flash (L outside (0, 1); with
 * PTFlashTableSettings::warmStartFromTable, for every flash):
 * - Exact: in an all-two-phase cell the interpolated K start the split and
 *   the stability test is skipped; a split that does not converge to a
 *   distinct interior split is redone by the flash without the table.
 * - Skip: as Exact, and in a cell whose vertices are all single-phase with
 *   one label the answer is that phase, without a stability test.
 * - Approximate: in an all-two-phase cell the answer is Rachford-Rice on the
 *   interpolated K; in a cell whose vertices are all single-phase with one
 *   label the answer is that phase; elsewhere the flash without the table.
 *
 * With negative-flash vertices (PTFlashTableSettings::negativeFlash), a cell
 * whose vertices are partly two-phase and partly single-phase is kept as a
 * boundary cell when every single-phase vertex has a negative flash: a split
 * with L outside [0, 1] that continues the tie-lines of the two-phase region
 * across its boundary (Whitson and Michelsen 1989). ln K then interpolates
 * over the whole cell, and Rachford-Rice at the state's own composition
 * places the boundary inside it.
 *
 * The table is thread-local and is cleared when the fluid system parameters
 * or the equation of state change.
 */
#ifndef OPM_PTFLASH_TABLE_HPP
#define OPM_PTFLASH_TABLE_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {

//! What the on-demand flash table does for a P-T flash.
enum class PTFlashTableMode {
    Off,         //!< every flash is solved without the table
    Exact,       //!< all-two-phase cells give the start of the split
    Skip,        //!< Exact, and all-single-phase cells answer single-phase
    Approximate, //!< cells with one phase state give the answer
};

/*!
 * \brief The optional speed-ups of PTFlash and PHFlash, read once per
 *        process.
 *
 * Everything is off by default: the flash is then the flash without these
 * options. get() reads the environment at its first call; a program can
 * also set the fields of get() before its first flash. Nothing here is read
 * again in the flash itself.
 *
 * | variable | values | default | what it does |
 * |---|---|---|---|
 * | OPM_FLASH_OBL | off, exact, skip, approx | off | the table of PTFlashTable, see PTFlashTableMode |
 * | OPM_FLASH_OBL_N | intervals per axis | 32 | the table's resolution |
 * | OPM_FLASH_OBL_NEGATIVE | 0, 1 | 1 | boundary cells from negative-flash vertices |
 * | OPM_FLASH_OBL_MARGIN | ln K error | 0.05 | margin of a single-phase verdict in a boundary cell |
 * | OPM_FLASH_OBL_PH | never, cold, always | cold | P-H search starts from the table's enthalpy |
 * | OPM_FLASH_OBL_WARM | 0, 1 | 0 | a split from an earlier flash also starts from the table |
 * | OPM_FLASH_OBL_ACCEPT | tolerance | 0 | approx: largest |ln(f_L/f_V)| of an interpolated split |
 * | OPM_FLASH_OBL_INTERP | multilinear, simplex | by size | interpolation; simplex from four components |
 * | OPM_FLASH_OBL_STATS | 0, 1 | 0 | table size and use on stderr |
 * | OPM_FLASH_BYPASS | 0, 1 | 0 | stability bypass of PTFlash::stabilityBypass_() |
 *
 * Recommended: OPM_FLASH_OBL=skip OPM_FLASH_BYPASS=1. The sub-options of
 * the table act only when OPM_FLASH_OBL is not off. Skip and approx are
 * approximate in principle (a single-phase answer without a stability
 * test); exact changes only where splits start.
 */
struct PTFlashTableSettings
{
    PTFlashTableMode mode = PTFlashTableMode::Off;
    int intervals = 32; //!< intervals per axis over the nominal ranges below
    double pressureMin = 1e5; //!< [Pa], the ln p axis
    double pressureMax = 500e5;
    double temperatureMin = 250.0; //!< [K], the 1/T axis
    double temperatureMax = 700.0;
    bool printStatistics = false; //!< table sizes on stderr when a table is cleared
    //! Boundary cells: single-phase vertices next to two-phase ones also
    //! store a negative flash, see PTFlashTable.
    bool negativeFlash = true;
    //! Boundary cells: a state whose interpolated split stays outside
    //! [0, 1] when every ln K is off by this much is single-phase without a
    //! stability test (not in the exact mode).
    double margin = 0.05;
    //! The P-H flash starts its temperature search where the enthalpy of
    //! the table's vertices matches: Never, Cold (a state without a phase
    //! state from an earlier flash, L outside [0, 1]) or Always.
    enum class EnthalpyStart { Never, Cold, Always };
    EnthalpyStart enthalpyStart = EnthalpyStart::Cold;
    //! Approximate mode: the largest |ln(f_L / f_V)| a two-phase answer may
    //! have; above it the split is converged from there. Zero: no test.
    double acceptance = 0.0;
    //! Interpolate on the d + 1 vertices of the Kuhn simplex of a state
    //! instead of the 2^d vertices of its cell; by default from d = 5
    //! (four components) on, where a cell has 32 or more vertices.
    enum class Interpolation { Automatic, Multilinear, Simplex };
    Interpolation interpolation = Interpolation::Automatic;
    //! A P-T flash with a split from an earlier flash (0 < L < 1) starts
    //! from its own ratios, which a simulator keeps between Newton
    //! iterations and are closer than the table's; with this, from the
    //! table too.
    bool warmStartFromTable = false;
    //! PTFlash: a single-phase state near a reference state whose stability
    //! test was decisive stays single-phase without a test (the critical
    //! distance of MultiComponentFlash.jl). Not a table option.
    bool stabilityBypass = false;

    static PTFlashTableSettings& get()
    {
        static PTFlashTableSettings settings = fromEnvironment();
        return settings;
    }

    static PTFlashTableSettings fromEnvironment()
    {
        PTFlashTableSettings settings;
        const auto read = [](const char* name) -> std::string {
            const char* value = std::getenv(name);
            return value ? std::string(value) : std::string();
        };
        const std::string mode = read("OPM_FLASH_OBL");
        if (mode == "exact") {
            settings.mode = PTFlashTableMode::Exact;
        }
        else if (mode == "skip") {
            settings.mode = PTFlashTableMode::Skip;
        }
        else if (mode == "approx") {
            settings.mode = PTFlashTableMode::Approximate;
        }
        if (const std::string n = read("OPM_FLASH_OBL_N"); !n.empty()) {
            settings.intervals = std::max(1, std::atoi(n.c_str()));
        }
        const std::string interpolation = read("OPM_FLASH_OBL_INTERP");
        if (interpolation == "simplex") {
            settings.interpolation = Interpolation::Simplex;
        }
        else if (interpolation == "multilinear") {
            settings.interpolation = Interpolation::Multilinear;
        }
        if (const std::string acceptance = read("OPM_FLASH_OBL_ACCEPT"); !acceptance.empty()) {
            settings.acceptance = std::atof(acceptance.c_str());
        }
        settings.printStatistics = read("OPM_FLASH_OBL_STATS") == "1";
        if (const std::string negative = read("OPM_FLASH_OBL_NEGATIVE"); !negative.empty()) {
            settings.negativeFlash = negative == "1";
        }
        if (const std::string margin = read("OPM_FLASH_OBL_MARGIN"); !margin.empty()) {
            settings.margin = std::atof(margin.c_str());
        }
        const std::string start = read("OPM_FLASH_OBL_PH");
        if (start == "never") {
            settings.enthalpyStart = EnthalpyStart::Never;
        }
        else if (start == "always") {
            settings.enthalpyStart = EnthalpyStart::Always;
        }
        settings.warmStartFromTable = read("OPM_FLASH_OBL_WARM") == "1";
        settings.stabilityBypass = read("OPM_FLASH_BYPASS") == "1";
        return settings;
    }
};

/*!
 * \brief The reports of OPM_FLASH_OBL_STATS = 1, one per thread and store.
 *
 * A thread's table and the stability memory of PTFlash register their
 * report here, and report when their thread ends. The threads of a pool
 * (OpenMP) do not end before the process, so what is still registered
 * reports when the process exits. Nothing registers without the statistics.
 */
class PTFlashReports
{
public:
    static PTFlashReports& get()
    {
        // Never destroyed: a thread may end after the static destructors.
        static PTFlashReports* reports = [] {
            auto* r = new PTFlashReports;
            std::atexit([] { PTFlashReports::get().reportAll(); });
            return r;
        }();
        return *reports;
    }

    //! The number of the calling thread in the reports, in order of first use.
    static int threadNumber()
    {
        static std::atomic<int> next{0};
        thread_local const int number = next++;
        return number;
    }

    void add(const void* owner, std::function<void()> report)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        entries_.emplace_back(owner, std::move(report));
    }

    //! Reports for \p owner and forgets it.
    void remove(const void* owner)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto it = entries_.begin(); it != entries_.end(); ++it) {
            if (it->first == owner) {
                it->second();
                entries_.erase(it);
                return;
            }
        }
    }

private:
    void reportAll()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        for (const auto& entry : entries_) {
            entry.second();
        }
        entries_.clear();
    }

    std::mutex mutex_;
    std::vector<std::pair<const void*, std::function<void()>>> entries_;
};

/*!
 * \brief The vertex and cell stores of one fluid system and thread.
 */
template <class Scalar, class FluidSystem>
class PTFlashTable
{
public:
    static constexpr int numComponents = FluidSystem::numComponents;
    //! ln p, 1/T and every mole fraction but the dependent one
    static constexpr int numAxes = numComponents + 1;
    static constexpr int numCorners = 1 << numAxes;
    static_assert(numAxes <= 16, "the lattice key packs at most 16 axes");

    //! Unused: a cell that stopped flashing its corners once it could not be used.
    //! Boundary: two-phase and single-phase vertices, all with a split.
    enum class Phase : std::uint8_t { Failed, TwoPhase, Liquid, Vapour, Mixed, Unused, Boundary };

    //! Whether a single-phase vertex has a negative flash; not tried until a
    //! boundary cell needs it.
    enum class Split : std::uint8_t { Unknown, Yes, No };

    //! An exact flash at a lattice point. L and ln K hold the split of a
    //! two-phase vertex, or the negative flash of a single-phase one.
    struct Vertex
    {
        Phase phase = Phase::Failed;
        Split negative = Split::Unknown;
        Scalar L = 0.0;
        std::array<Scalar, numComponents> lnK{};

        bool hasSplit() const
        {
            return phase == Phase::TwoPhase || negative == Split::Yes;
        }
    };

    //! What a query led to, for the statistics.
    enum class Outcome : std::uint8_t {
        TwoPhaseStart,  //!< all-two-phase vertices: split from the table's K, kept
        TwoPhaseRetry,  //!< all-two-phase vertices: that split not kept
        SinglePhase,    //!< all-single-phase vertices: answered single-phase
        BoundarySplit,  //!< boundary cell, interior L: split kept
        BoundaryRetry,  //!< boundary cell, interior L: split not kept
        BoundarySingle, //!< boundary cell, L beyond the margin: single-phase, no test
        BoundaryBand,   //!< boundary cell, L inside the margin: the flash
        Flash,          //!< any other: the flash without the table
        Count
    };

    //! A hypercube: the common phase state of its vertices, Mixed when they
    //! differ, Failed when one cannot be flashed; ln K of every corner
    //! (corner-major) when all are two-phase.
    struct Cell
    {
        Phase phase = Phase::Failed;
        std::vector<Scalar> lnK;
    };

    //! A state's cell and its position inside it, each fraction in [0, 1).
    struct Location
    {
        std::array<int, numAxes> index{};
        std::array<Scalar, numAxes> fraction{};
    };

    static PTFlashTable& get()
    {
        static thread_local PTFlashTable table;
        return table;
    }

    PTFlashTable()
    {
        if (PTFlashTableSettings::get().printStatistics) {
            thread_ = PTFlashReports::threadNumber();
            PTFlashReports::get().add(this, [this] { report(); });
        }
    }

    ~PTFlashTable()
    {
        if (PTFlashTableSettings::get().printStatistics) {
            PTFlashReports::get().remove(this);
        }
    }

    PTFlashTable(const PTFlashTable&) = delete;
    PTFlashTable& operator=(const PTFlashTable&) = delete;

    //! Empties the table at the next flash, which also reads the settings again.
    void reset()
    {
        initialized_ = false;
    }

    //! Whether vertices are being flashed; their flashes must not use the table.
    bool building() const
    {
        return building_;
    }

    /*!
     * \brief A cheap hash of the fluid a vertex flash depends on: the
     *        equation of state and the component parameters.
     *
     * Not the method or the tolerance of the call: the table gives starts
     * and phase states, and every answer is converged by the caller's
     * method. A simulator that flashes its cells and its wells with
     * different methods keeps one table.
     */
    static std::uint64_t fingerprint(int eos)
    {
        std::uint64_t h = 0x9e3779b97f4a7c15ull;
        const auto mix = [&h](Scalar v) {
            std::uint64_t bits = 0;
            const double d = static_cast<double>(v);
            std::memcpy(&bits, &d, sizeof(bits));
            h = mix64(h ^ bits);
        };
        mix(eos);
        for (int i = 0; i < numComponents; ++i) {
            mix(FluidSystem::criticalTemperature(i));
            mix(FluidSystem::criticalPressure(i));
            mix(FluidSystem::acentricFactor(i));
            for (int j = 0; j < i; ++j) {
                mix(FluidSystem::interactionCoefficient(i, j));
            }
        }
        return h;
    }

    /*!
     * \brief Clears the table when the fingerprint changed; a new table takes
     *        the largest component of the first state as the dependent one.
     */
    template <class Composition>
    void prepare(std::uint64_t fingerprint, const Composition& z)
    {
        if (initialized_ && fingerprint == fingerprint_) {
            return;
        }
        report();
        vertices_.clear();
        cells_.clear();
        enthalpies_.clear();
        enthalpyQueries_ = enthalpyStarts_ = 0;
        queries_ = hits_ = vertexFlashes_ = negativeFlashes_ = 0;
        outcomes_.fill(0);
        fingerprint_ = fingerprint;
        initialized_ = true;
        const auto& s = PTFlashTableSettings::get();
        intervals_ = s.intervals;
        negative_ = s.negativeFlash;
        using Interpolation = PTFlashTableSettings::Interpolation;
        simplex_ = s.interpolation == Interpolation::Simplex
            || (s.interpolation == Interpolation::Automatic && numAxes >= 5);
        origin_[0] = std::log(s.pressureMin);
        step_[0] = (std::log(s.pressureMax) - origin_[0]) / s.intervals;
        origin_[1] = 1.0 / s.temperatureMax;
        step_[1] = (1.0 / s.temperatureMin - origin_[1]) / s.intervals;
        dependent_ = 0;
        for (int c = 1; c < numComponents; ++c) {
            if (z[c] > z[dependent_]) {
                dependent_ = c;
            }
        }
        for (int a = 2; a < numAxes; ++a) {
            origin_[a] = 0.0;
            step_[a] = 1.0 / s.intervals;
        }
    }

    //! The cell of a state; false outside the range of the lattice key.
    template <class Composition>
    bool locate(Scalar p, Scalar T, const Composition& z, Location& location)
    {
        ++queries_;
        std::array<Scalar, numAxes> u;
        u[0] = std::log(p);
        u[1] = 1.0 / T;
        for (int c = 0, a = 2; c < numComponents; ++c) {
            if (c != dependent_) {
                u[a++] = z[c];
            }
        }
        for (int a = 0; a < numAxes; ++a) {
            const Scalar s = (u[a] - origin_[a]) / step_[a];
            if (!std::isfinite(s)) {
                return false;
            }
            const Scalar f = std::floor(s);
            if (f < -keyOffset || f >= keyOffset - 1) {
                return false;
            }
            location.index[a] = static_cast<int>(f);
            location.fraction[a] = s - f;
        }
        return true;
    }

    /*!
     * \brief The phase state of the vertices that interpolate a state and,
     *        when they are all two-phase, the interpolated ln K.
     *
     * The vertices are the 2^d corners of the cell (multilinear), or the
     * d + 1 corners of the Kuhn simplex that holds the state: with the
     * fractions sorted, f_(1) >= ... >= f_(d), the simplex walks from the
     * lower corner one step along each axis in that order, and the weights
     * are 1 - f_(1), f_(1) - f_(2), ..., f_(d).
     */
    template <class VertexFlash, class NegativeFlash, class Useful>
    Phase evaluate(const Location& location, const VertexFlash& flash,
                   const NegativeFlash& negativeFlash, const Useful& useful,
                   std::array<Scalar, numComponents>& lnK)
    {
        if (!simplex_) {
            const Cell& c = cell(location, flash, negativeFlash, useful);
            if (c.phase == Phase::TwoPhase || c.phase == Phase::Boundary) {
                interpolate(c, location, lnK);
            }
            return c.phase;
        }
        if (negative_) {
            return evaluateSimplexWithBoundary(location, flash, negativeFlash, useful, lnK);
        }
        std::array<int, numAxes> order;
        for (int a = 0; a < numAxes; ++a) {
            order[a] = a;
        }
        std::sort(order.begin(), order.end(), [&location](int a, int b) {
            return location.fraction[a] > location.fraction[b];
        });
        std::array<int, numAxes> index = location.index;
        Phase phase = Phase::Failed;
        bool present = true;
        lnK.fill(0.0);
        building_ = true;
        for (int k = 0; k <= numAxes; ++k) {
            if (k > 0) {
                ++index[order[k - 1]];
            }
            const Scalar weight = k == 0 ? 1.0 - location.fraction[order[0]]
                : (k == numAxes ? location.fraction[order[k - 1]]
                                : location.fraction[order[k - 1]] - location.fraction[order[k]]);
            const Vertex& v = vertex(index, flash, &present);
            if (k == 0) {
                phase = v.phase;
            }
            else if (v.phase != phase) {
                phase = (v.phase == Phase::Failed || phase == Phase::Failed) ? Phase::Failed : Phase::Mixed;
            }
            if (phase == Phase::Failed || phase == Phase::Mixed) {
                break;
            }
            if (!useful(phase)) {
                phase = Phase::Unused;
                break;
            }
            for (int c = 0; c < numComponents; ++c) {
                lnK[c] += weight * v.lnK[c];
            }
        }
        building_ = false;
        hits_ += present;
        return phase;
    }

    /*!
     * \brief The cell at a location, built on first use.
     *
     * \param flash returns the Vertex of an exact flash at (p, T, z).
     * \param useful whether a cell of the phase state would be used; a cell
     *        stops flashing its corners once it cannot be.
     */
    template <class VertexFlash, class NegativeFlash, class Useful>
    const Cell& cell(const Location& location, const VertexFlash& flash,
                     const NegativeFlash& negativeFlash, const Useful& useful)
    {
        const std::uint64_t key = pack(location.index);
        const auto found = cells_.find(key);
        if (found != cells_.end()) {
            ++hits_;
            return found->second;
        }
        if (negative_) {
            return cellWithBoundary(key, location, flash, negativeFlash, useful);
        }
        building_ = true;
        Cell cell;
        std::vector<Scalar> lnK(numCorners * numComponents);
        for (int corner = 0; corner < numCorners; ++corner) {
            std::array<int, numAxes> index = location.index;
            for (int a = 0; a < numAxes; ++a) {
                index[a] += (corner >> a) & 1;
            }
            const Vertex& v = vertex(index, flash);
            if (corner == 0) {
                cell.phase = v.phase;
            }
            else if (v.phase != cell.phase) {
                cell.phase = (v.phase == Phase::Failed || cell.phase == Phase::Failed)
                    ? Phase::Failed : Phase::Mixed;
            }
            if (cell.phase == Phase::Failed || cell.phase == Phase::Mixed) {
                break;
            }
            if (!useful(cell.phase)) {
                cell.phase = Phase::Unused;
                break;
            }
            std::copy(v.lnK.begin(), v.lnK.end(), lnK.begin() + corner * numComponents);
        }
        building_ = false;
        if (cell.phase == Phase::TwoPhase) {
            cell.lnK = std::move(lnK);
        }
        return cells_.emplace(key, std::move(cell)).first->second;
    }

    //! The multilinear interpolant of ln K in an all-two-phase cell.
    void interpolate(const Cell& cell, const Location& location,
                     std::array<Scalar, numComponents>& lnK) const
    {
        std::array<Scalar, numCorners> weight;
        weight[0] = 1.0;
        for (int a = 0, size = 1; a < numAxes; ++a, size *= 2) {
            const Scalar t = location.fraction[a];
            for (int corner = 0; corner < size; ++corner) {
                weight[corner + size] = weight[corner] * t;
                weight[corner] *= 1.0 - t;
            }
        }
        lnK.fill(0.0);
        for (int corner = 0; corner < numCorners; ++corner) {
            const Scalar* v = cell.lnK.data() + corner * numComponents;
            for (int c = 0; c < numComponents; ++c) {
                lnK[c] += weight[corner] * v[c];
            }
        }
    }

    /*!
     * \brief As cell(), when boundary cells are kept: every corner is flashed,
     *        and a cell with two-phase and single-phase corners is a boundary
     *        cell when each single-phase corner has a negative flash, started
     *        from ln K of the nearest two-phase corner.
     */
    template <class VertexFlash, class NegativeFlash, class Useful>
    const Cell& cellWithBoundary(std::uint64_t key, const Location& location,
                                 const VertexFlash& flash, const NegativeFlash& negativeFlash,
                                 const Useful& useful)
    {
        building_ = true;
        std::array<Vertex*, numCorners> corners{};
        std::array<std::array<int, numAxes>, numCorners> indices{};
        for (int corner = 0; corner < numCorners; ++corner) {
            indices[corner] = location.index;
            for (int a = 0; a < numAxes; ++a) {
                indices[corner][a] += (corner >> a) & 1;
            }
            corners[corner] = &vertex(indices[corner], flash);
        }
        const auto distance = [](int a, int b) { return std::popcount(static_cast<unsigned>(a ^ b)); };
        Cell cell;
        cell.phase = classify(corners, indices, distance, negativeFlash);
        building_ = false;
        if (cell.phase != Phase::Boundary && cell.phase != Phase::Mixed
            && cell.phase != Phase::Failed && !useful(cell.phase)) {
            cell.phase = Phase::Unused;
        }
        if (cell.phase == Phase::TwoPhase || cell.phase == Phase::Boundary) {
            cell.lnK.resize(numCorners * numComponents);
            for (int corner = 0; corner < numCorners; ++corner) {
                std::copy(corners[corner]->lnK.begin(), corners[corner]->lnK.end(),
                          cell.lnK.begin() + corner * numComponents);
            }
        }
        return cells_.emplace(key, std::move(cell)).first->second;
    }

    /*!
     * \brief As the simplex part of evaluate(), when boundary cells are kept:
     *        the d + 1 vertices of the simplex, a boundary when two-phase and
     *        single-phase vertices all have a split.
     */
    template <class VertexFlash, class NegativeFlash, class Useful>
    Phase evaluateSimplexWithBoundary(const Location& location, const VertexFlash& flash,
                                      const NegativeFlash& negativeFlash, const Useful& useful,
                                      std::array<Scalar, numComponents>& lnK)
    {
        std::array<int, numAxes> order;
        for (int a = 0; a < numAxes; ++a) {
            order[a] = a;
        }
        std::sort(order.begin(), order.end(), [&location](int a, int b) {
            return location.fraction[a] > location.fraction[b];
        });
        std::array<Vertex*, numAxes + 1> vertices{};
        std::array<std::array<int, numAxes>, numAxes + 1> indices{};
        std::array<Scalar, numAxes + 1> weights{};
        bool present = true;
        building_ = true;
        std::array<int, numAxes> index = location.index;
        for (int k = 0; k <= numAxes; ++k) {
            if (k > 0) {
                ++index[order[k - 1]];
            }
            weights[k] = k == 0 ? 1.0 - location.fraction[order[0]]
                : (k == numAxes ? location.fraction[order[k - 1]]
                                : location.fraction[order[k - 1]] - location.fraction[order[k]]);
            indices[k] = index;
            vertices[k] = &vertex(index, flash, &present);
        }
        // Along the walk, vertex k and j differ in |k - j| coordinates.
        const auto distance = [](int a, int b) { return a > b ? a - b : b - a; };
        Phase phase = classify(vertices, indices, distance, negativeFlash);
        building_ = false;
        hits_ += present;
        if (phase != Phase::Boundary && phase != Phase::Mixed && phase != Phase::Failed
            && !useful(phase)) {
            return Phase::Unused;
        }
        if (phase == Phase::TwoPhase || phase == Phase::Boundary) {
            lnK.fill(0.0);
            for (int k = 0; k <= numAxes; ++k) {
                for (int c = 0; c < numComponents; ++c) {
                    lnK[c] += weights[k] * vertices[k]->lnK[c];
                }
            }
        }
        return phase;
    }

    /*!
     * \brief A start for the temperature search of a P-H flash: where the
     *        enthalpy of the vertices equals \p hSpec.
     *
     * At each temperature level of the lattice the enthalpy is interpolated
     * multilinearly over the other axes (ln p and z) at the state, and
     * between two levels linearly in T. The search begins at the levels
     * around \p T0 and moves one level towards \p hSpec at a time, so only
     * the vertices near the answer are flashed. A vertex enthalpy comes from
     * \p enthalpy(p, T, z, vertex) and is kept for the rest of the run;
     * \p enthalpyKey identifies the enthalpy model, a change empties them.
     *
     * \return nothing outside [tMin, tMax], where a vertex has no enthalpy,
     *         or where the enthalpy does not increase between two levels
     */
    template <class Composition, class VertexFlash, class Enthalpy>
    std::optional<Scalar> temperature(Scalar p, Scalar T0, const Composition& z, Scalar hSpec,
                                      Scalar tMin, Scalar tMax, std::uint64_t enthalpyKey,
                                      const VertexFlash& flash, const Enthalpy& enthalpy)
    {
        if (enthalpyKey != enthalpyKey_) {
            enthalpies_.clear();
            enthalpyKey_ = enthalpyKey;
        }
        Location location;
        if (!locate(p, T0, z, location)) {
            return std::nullopt;
        }
        ++enthalpyQueries_;
        // The enthalpy at level j of the 1/T axis, over the other axes.
        const auto level = [&](int j) -> std::optional<Scalar> {
            Scalar h = 0.0;
            for (int corner = 0; corner < numCorners / 2; ++corner) {
                std::array<int, numAxes> index = location.index;
                index[1] = j;
                Scalar weight = 1.0;
                for (int a = 0, bit = 0; a < numAxes; ++a) {
                    if (a == 1) {
                        continue;
                    }
                    const int upper = (corner >> bit++) & 1;
                    index[a] += upper;
                    weight *= upper ? location.fraction[a] : 1.0 - location.fraction[a];
                }
                if (weight == 0.0) {
                    continue;
                }
                const Scalar hv = vertexEnthalpy(index, flash, enthalpy);
                if (!std::isfinite(hv)) {
                    return std::nullopt;
                }
                h += weight * hv;
            }
            return h;
        };
        const auto levelTemperature = [this](int j) { return 1.0 / (origin_[1] + j * step_[1]); };
        // Level j is the warmer end of the interval, j + 1 the colder.
        int j = location.index[1];
        for (int move = 0; move < 2 * intervals_; ++move) {
            const Scalar tWarm = levelTemperature(j);
            const Scalar tCold = levelTemperature(j + 1);
            if (!(tWarm > tMin && tCold < tMax)) {
                return std::nullopt;
            }
            const auto hWarm = level(j);
            const auto hCold = level(j + 1);
            if (!hWarm || !hCold || !(*hWarm > *hCold)) {
                return std::nullopt;
            }
            if (hSpec > *hWarm) {
                --j;
            }
            else if (hSpec < *hCold) {
                ++j;
            }
            else {
                ++enthalpyStarts_;
                const Scalar T = tCold + (hSpec - *hCold) / (*hWarm - *hCold) * (tWarm - tCold);
                return std::clamp(T, tMin, tMax);
            }
        }
        return std::nullopt;
    }

    //! Counts what a query led to.
    void note(Outcome outcome)
    {
        ++outcomes_[static_cast<int>(outcome)];
    }

    //! How many queries led to an outcome since the table was last emptied.
    long count(Outcome outcome) const
    {
        return outcomes_[static_cast<int>(outcome)];
    }

    std::size_t vertices() const
    {
        return vertices_.size();
    }

    //! Starts for the P-H flash given, and asked for.
    long enthalpyStarts() const
    {
        return enthalpyStarts_;
    }

    std::size_t cells() const
    {
        return cells_.size();
    }

    //! Table sizes and use on stderr, when OPM_FLASH_OBL_STATS = 1.
    void report() const
    {
        if (!PTFlashTableSettings::get().printStatistics || queries_ == 0) {
            return;
        }
        std::size_t phases[7] = {0, 0, 0, 0, 0, 0, 0};
        std::size_t cellBytes = 0;
        for (const auto& [key, c] : cells_) {
            ++phases[static_cast<int>(c.phase)];
            cellBytes += sizeof(Cell) + sizeof(key) + 2 * sizeof(void*) + c.lnK.size() * sizeof(Scalar);
        }
        std::size_t failedVertices = 0;
        for (const auto& [key, v] : vertices_) {
            failedVertices += v.phase == Phase::Failed;
        }
        const std::size_t vertexBytes
            = vertices_.size() * (sizeof(Vertex) + sizeof(std::uint64_t) + 2 * sizeof(void*));
        std::cerr << "OBL-TABLE thread " << thread_ << " components " << numComponents
                  << " intervals " << intervals_
                  << (simplex_ ? " simplex" : " multilinear")
                  << " queries " << queries_ << " present " << hits_ << " vertex-flashes "
                  << vertexFlashes_ << " vertices " << vertices_.size() << " (failed " << failedVertices
                  << ") cells " << cells_.size() << " (two-phase " << phases[1] << ", liquid " << phases[2]
                  << ", vapour " << phases[3] << ", mixed " << phases[4] << ", unused " << phases[5]
                  << ", failed " << phases[0] << ", boundary " << phases[6]
                  << ") negative-flashes " << negativeFlashes_
                  << " bytes " << vertexBytes + cellBytes << std::endl;
        if (std::accumulate(outcomes_.begin(), outcomes_.end(), 0L) > 0) {
            static constexpr const char* names[] = {"two-phase-start", "two-phase-retry",
                "single-phase", "boundary-split", "boundary-retry", "boundary-single",
                "boundary-band", "flash"};
            std::cerr << "OBL-OUTCOMES thread " << thread_;
            for (int o = 0; o < static_cast<int>(Outcome::Count); ++o) {
                std::cerr << ' ' << names[o] << ' ' << outcomes_[o];
            }
            std::cerr << std::endl;
        }
        if (enthalpyQueries_ > 0) {
            std::cerr << "OBL-ENTHALPY thread " << thread_ << " queries " << enthalpyQueries_
                      << " starts " << enthalpyStarts_ << " vertex-enthalpies "
                      << enthalpies_.size() << std::endl;
        }
    }

private:
    static constexpr int keyBits = std::min(16, 64 / numAxes);
    static constexpr Scalar keyOffset = Scalar(1 << (keyBits - 1));

    static std::uint64_t mix64(std::uint64_t x)
    {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ull;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebull;
        x ^= x >> 31;
        return x;
    }

    struct KeyHash
    {
        std::size_t operator()(std::uint64_t key) const
        {
            return static_cast<std::size_t>(mix64(key));
        }
    };

    static std::uint64_t pack(const std::array<int, numAxes>& index)
    {
        std::uint64_t key = 0;
        for (int a = 0; a < numAxes; ++a) {
            const auto shifted = static_cast<std::uint64_t>(index[a] + (1 << (keyBits - 1)));
            key = (key << keyBits) | (shifted & ((std::uint64_t(1) << keyBits) - 1));
        }
        return key;
    }

    /*!
     * \brief The state of a lattice point; false for a point outside the
     *        simplex. A mole fraction on the lower face is a trace, so that
     *        its K is defined.
     */
    bool point(const std::array<int, numAxes>& index, Scalar& p, Scalar& T,
               std::array<Scalar, numComponents>& z) const
    {
        constexpr Scalar trace = 1e-6;
        Scalar sum = 0.0;
        for (int c = 0, a = 2; c < numComponents; ++c) {
            if (c != dependent_) {
                z[c] = std::max(origin_[a] + index[a] * step_[a], trace);
                sum += z[c];
                ++a;
            }
        }
        const Scalar zDependent = 1.0 - sum;
        if (!(zDependent > -1e-12)) {
            return false;
        }
        z[dependent_] = std::max(zDependent, trace);
        sum += z[dependent_];
        for (auto& zc : z) {
            zc /= sum;
        }
        p = std::exp(origin_[0] + index[0] * step_[0]);
        T = 1.0 / (origin_[1] + index[1] * step_[1]);
        return true;
    }

    template <class VertexFlash>
    Vertex& vertex(const std::array<int, numAxes>& index, const VertexFlash& flash,
                   bool* present = nullptr)
    {
        const std::uint64_t key = pack(index);
        const auto found = vertices_.find(key);
        if (found != vertices_.end()) {
            return found->second;
        }
        if (present) {
            *present = false;
        }
        Vertex v;
        Scalar p = 0.0;
        Scalar T = 0.0;
        std::array<Scalar, numComponents> z{};
        if (point(index, p, T, z)) {
            ++vertexFlashes_;
            v = flash(p, T, z);
            if (v.phase == Phase::TwoPhase) {
                v.negative = Split::Yes;
            }
        }
        return vertices_.emplace(key, v).first->second;
    }

    //! The enthalpy of a vertex, NaN where it has none.
    template <class VertexFlash, class Enthalpy>
    Scalar vertexEnthalpy(const std::array<int, numAxes>& index, const VertexFlash& flash,
                          const Enthalpy& enthalpy)
    {
        const std::uint64_t key = pack(index);
        const auto found = enthalpies_.find(key);
        if (found != enthalpies_.end()) {
            return found->second;
        }
        building_ = true;
        const Vertex& v = vertex(index, flash);
        building_ = false;
        Scalar h = std::numeric_limits<Scalar>::quiet_NaN();
        Scalar p = 0.0;
        Scalar T = 0.0;
        std::array<Scalar, numComponents> z{};
        if (v.phase != Phase::Failed && point(index, p, T, z)) {
            if (const auto value = enthalpy(p, T, z, v)) {
                h = *value;
            }
        }
        return enthalpies_.emplace(key, h).first->second;
    }

    /*!
     * \brief The phase state of the vertices of a cell or simplex, with the
     *        negative flashes a boundary needs.
     *
     * \param distance the number of coordinates in which two vertices differ;
     *        a negative flash starts from ln K of the nearest two-phase vertex.
     */
    template <std::size_t size, class Distance, class NegativeFlash>
    Phase classify(const std::array<Vertex*, size>& vertices,
                   const std::array<std::array<int, numAxes>, size>& indices,
                   const Distance& distance, const NegativeFlash& negativeFlash)
    {
        bool same = true;
        bool twoPhase = false;
        for (std::size_t k = 0; k < size; ++k) {
            if (vertices[k]->phase == Phase::Failed) {
                return Phase::Failed;
            }
            same = same && vertices[k]->phase == vertices[0]->phase;
            twoPhase = twoPhase || vertices[k]->phase == Phase::TwoPhase;
        }
        if (same) {
            return vertices[0]->phase;
        }
        if (!twoPhase) {
            return Phase::Mixed;
        }
        for (std::size_t k = 0; k < size; ++k) {
            Vertex& v = *vertices[k];
            if (v.negative == Split::Unknown) {
                std::size_t nearest = size;
                for (std::size_t j = 0; j < size; ++j) {
                    if (vertices[j]->phase == Phase::TwoPhase
                        && (nearest == size
                            || distance(static_cast<int>(j), static_cast<int>(k))
                                < distance(static_cast<int>(nearest), static_cast<int>(k)))) {
                        nearest = j;
                    }
                }
                Scalar p = 0.0;
                Scalar T = 0.0;
                std::array<Scalar, numComponents> z{};
                point(indices[k], p, T, z);
                ++negativeFlashes_;
                v.negative = negativeFlash(p, T, z, vertices[nearest]->lnK, v.L, v.lnK)
                    ? Split::Yes : Split::No;
            }
            if (v.negative != Split::Yes) {
                return Phase::Mixed;
            }
        }
        return Phase::Boundary;
    }

    bool initialized_ = false;
    int thread_ = 0;
    bool simplex_ = false;
    bool negative_ = false;
    std::uint64_t fingerprint_ = 0;
    int dependent_ = 0;
    int intervals_ = 0;
    std::array<Scalar, numAxes> origin_{};
    std::array<Scalar, numAxes> step_{};
    bool building_ = false;
    long queries_ = 0;
    long hits_ = 0;
    long vertexFlashes_ = 0;
    long negativeFlashes_ = 0;
    long enthalpyQueries_ = 0;
    long enthalpyStarts_ = 0;
    std::uint64_t enthalpyKey_ = 0;
    std::unordered_map<std::uint64_t, Scalar, KeyHash> enthalpies_;
    std::array<long, static_cast<int>(Outcome::Count)> outcomes_{};
    std::unordered_map<std::uint64_t, Vertex, KeyHash> vertices_;
    std::unordered_map<std::uint64_t, Cell, KeyHash> cells_;
};

} // namespace Opm

#endif // OPM_PTFLASH_TABLE_HPP
