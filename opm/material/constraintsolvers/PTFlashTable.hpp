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
 * PTFlash consults the table when PTFlashTableSettings::get().mode is not Off:
 * - Exact: in an all-two-phase cell the interpolated K start the split and
 *   the stability test is skipped; a split that does not converge to a
 *   distinct interior split is redone by the flash without the table.
 * - Approximate: in an all-two-phase cell the answer is Rachford-Rice on the
 *   interpolated K; in a cell whose vertices are all single-phase with one
 *   label the answer is that phase; elsewhere the flash without the table.
 *
 * The table is thread-local and is cleared when the fluid system parameters,
 * the equation of state, the method or the tolerance change.
 */
#ifndef OPM_PTFLASH_TABLE_HPP
#define OPM_PTFLASH_TABLE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {

//! What the on-demand flash table does for a P-T flash.
enum class PTFlashTableMode {
    Off,         //!< every flash is solved without the table
    Exact,       //!< all-two-phase cells give the start of the split
    Approximate, //!< cells with one phase state give the answer
};

//! The table settings, read once from the environment.
struct PTFlashTableSettings
{
    //! Off unless asked for: OPM_FLASH_OBL, or set here before the first flash.
    PTFlashTableMode mode = PTFlashTableMode::Off;
    int intervals = 32; //!< intervals per axis over the nominal ranges below
    double pressureMin = 1e5; //!< [Pa], the ln p axis
    double pressureMax = 500e5;
    double temperatureMin = 250.0; //!< [K], the 1/T axis
    double temperatureMax = 700.0;
    bool printStatistics = false; //!< table sizes on stderr when a table is cleared
    //! Approximate mode: the largest |ln(f_L / f_V)| a two-phase answer may
    //! have; above it the split is converged from there. Zero: no test.
    double acceptance = 0.0;
    //! Interpolate on the d + 1 vertices of the Kuhn simplex of a state
    //! instead of the 2^d vertices of its cell; by default from d = 5
    //! (four components) on, where a cell has 32 or more vertices.
    enum class Interpolation { Automatic, Multilinear, Simplex };
    Interpolation interpolation = Interpolation::Automatic;

    //! OPM_FLASH_OBL = off | exact | approx, OPM_FLASH_OBL_N = intervals,
    //! OPM_FLASH_OBL_INTERP = multilinear | simplex, OPM_FLASH_OBL_ACCEPT =
    //! acceptance, OPM_FLASH_OBL_STATS = 1.
    static PTFlashTableSettings& get()
    {
        static PTFlashTableSettings settings = fromEnvironment();
        return settings;
    }

    static PTFlashTableSettings fromEnvironment()
    {
        PTFlashTableSettings settings;
        if (const char* mode = std::getenv("OPM_FLASH_OBL")) {
            const std::string value(mode);
            if (value == "off" || value == "0") {
                settings.mode = PTFlashTableMode::Off;
            }
            else if (value == "exact") {
                settings.mode = PTFlashTableMode::Exact;
            }
            else if (value == "approx") {
                settings.mode = PTFlashTableMode::Approximate;
            }
        }
        if (const char* n = std::getenv("OPM_FLASH_OBL_N")) {
            settings.intervals = std::max(1, std::atoi(n));
        }
        if (const char* interpolation = std::getenv("OPM_FLASH_OBL_INTERP")) {
            const std::string value(interpolation);
            if (value == "simplex") {
                settings.interpolation = Interpolation::Simplex;
            }
            else if (value == "multilinear") {
                settings.interpolation = Interpolation::Multilinear;
            }
        }
        if (const char* acceptance = std::getenv("OPM_FLASH_OBL_ACCEPT")) {
            settings.acceptance = std::atof(acceptance);
        }
        if (const char* stats = std::getenv("OPM_FLASH_OBL_STATS")) {
            settings.printStatistics = std::string(stats) == "1";
        }
        return settings;
    }
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
    enum class Phase : std::uint8_t { Failed, TwoPhase, Liquid, Vapour, Mixed, Unused };

    //! An exact flash at a lattice point.
    struct Vertex
    {
        Phase phase = Phase::Failed;
        Scalar L = 0.0;
        std::array<Scalar, numComponents> lnK{};
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

    ~PTFlashTable()
    {
        report();
    }

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

    //! A cheap hash of everything a vertex flash depends on but the state.
    static std::uint64_t fingerprint(int eos, int method, Scalar tolerance)
    {
        std::uint64_t h = 0x9e3779b97f4a7c15ull;
        const auto mix = [&h](Scalar v) {
            std::uint64_t bits = 0;
            const double d = static_cast<double>(v);
            std::memcpy(&bits, &d, sizeof(bits));
            h = mix64(h ^ bits);
        };
        mix(eos);
        mix(method);
        mix(tolerance);
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
        queries_ = hits_ = vertexFlashes_ = 0;
        fingerprint_ = fingerprint;
        initialized_ = true;
        const auto& s = PTFlashTableSettings::get();
        intervals_ = s.intervals;
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
    template <class VertexFlash, class Useful>
    Phase evaluate(const Location& location, const VertexFlash& flash, const Useful& useful,
                   std::array<Scalar, numComponents>& lnK)
    {
        if (!simplex_) {
            const Cell& c = cell(location, flash, useful);
            if (c.phase == Phase::TwoPhase) {
                interpolate(c, location, lnK);
            }
            return c.phase;
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
    template <class VertexFlash, class Useful>
    const Cell& cell(const Location& location, const VertexFlash& flash, const Useful& useful)
    {
        const std::uint64_t key = pack(location.index);
        const auto found = cells_.find(key);
        if (found != cells_.end()) {
            ++hits_;
            return found->second;
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

    std::size_t vertices() const
    {
        return vertices_.size();
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
        std::size_t phases[6] = {0, 0, 0, 0, 0, 0};
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
        std::cerr << "OBL-TABLE components " << numComponents << " intervals " << intervals_
                  << (simplex_ ? " simplex" : " multilinear")
                  << " queries " << queries_ << " present " << hits_ << " vertex-flashes "
                  << vertexFlashes_ << " vertices " << vertices_.size() << " (failed " << failedVertices
                  << ") cells " << cells_.size() << " (two-phase " << phases[1] << ", liquid " << phases[2]
                  << ", vapour " << phases[3] << ", mixed " << phases[4] << ", unused " << phases[5]
                  << ", failed " << phases[0]
                  << ") bytes " << vertexBytes + cellBytes << std::endl;
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

    template <class VertexFlash>
    const Vertex& vertex(const std::array<int, numAxes>& index, const VertexFlash& flash,
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
        // The lattice point; a mole fraction on the lower face is flashed as a
        // trace, so its K is defined, and a point outside the simplex is Failed.
        constexpr Scalar trace = 1e-6;
        std::array<Scalar, numComponents> z{};
        Scalar sum = 0.0;
        for (int c = 0, a = 2; c < numComponents; ++c) {
            if (c != dependent_) {
                z[c] = std::max(origin_[a] + index[a] * step_[a], trace);
                sum += z[c];
                ++a;
            }
        }
        Vertex v;
        const Scalar zDependent = 1.0 - sum;
        if (zDependent > -1e-12) {
            z[dependent_] = std::max(zDependent, trace);
            sum += z[dependent_];
            for (auto& zc : z) {
                zc /= sum;
            }
            const Scalar p = std::exp(origin_[0] + index[0] * step_[0]);
            const Scalar T = 1.0 / (origin_[1] + index[1] * step_[1]);
            ++vertexFlashes_;
            v = flash(p, T, z);
        }
        return vertices_.emplace(key, v).first->second;
    }

    bool initialized_ = false;
    bool simplex_ = false;
    std::uint64_t fingerprint_ = 0;
    int dependent_ = 0;
    int intervals_ = 0;
    std::array<Scalar, numAxes> origin_{};
    std::array<Scalar, numAxes> step_{};
    bool building_ = false;
    long queries_ = 0;
    long hits_ = 0;
    long vertexFlashes_ = 0;
    std::unordered_map<std::uint64_t, Vertex, KeyHash> vertices_;
    std::unordered_map<std::uint64_t, Cell, KeyHash> cells_;
};

} // namespace Opm

#endif // OPM_PTFLASH_TABLE_HPP
