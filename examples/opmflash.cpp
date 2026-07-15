/*
  Copyright 2026, Equinor ASA

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
/*
  opmflash — flash a compositional state from the command line.

  Runs the isothermal (P-T) or isenthalpic (P-H) flash on a synthetic
  mixture, through the same fluid system and solvers the compositional
  simulator uses. Intended for manual validation, quick checks and
  benchmark-curve generation (--sweep + --format csv), without writing
  C++ or building a simulator deck.

    opmflash --spec pt --p 50 --t 300 --z 0.5,0.5 --components C1,C10
    opmflash --spec ph --p 50 --h -25000 --z 0.5,0.5 --components C1,C10
    opmflash --spec ph --p 50 --z 0.5,0.5 --components C1,C10 \
             --sweep h=-40000:5000:1000 --format csv
    opmflash --list-components
    opmflash --selftest

  Units: --p bar · --t K · --h J/mol (molar enthalpy against the caloric
  datum H(298.15 K) = 0). Output enthalpies are J/mol.
*/

#include "config.h"

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/constraintsolvers/IdealGasCaloricData.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PHFlash.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using EOSType = Opm::CompositionalConfig::EOSType;

// ── component database ─────────────────────────────────────────────────────
// EoS parameters from the component classes of this library (C1.hpp,
// C10.hpp, SimpleCO2.hpp), converted to the fluid-system conventions
// (molar mass g/mol, critical volume m3/kmol). Heat-capacity data comes
// from the IdealGasCaloricData presets of the same names. Binary
// interaction coefficients are the standard PR literature values used
// throughout this library's compositional tests: C1/nC10 0.0411;
// C1/CO2 and CO2/nC10 0.10.
struct ComponentEntry {
    std::string_view name;    // canonical, understood by IdealGasCaloricData
    double molarMassGramPerMol;
    double criticalTemperature;
    double criticalPressure;
    double criticalVolume;    // [m3/kmol]
    double acentricFactor;
};

template <class Comp>
ComponentEntry makeEntry()
{
    return {Comp::name(), Comp::molarMass() * 1e3, Comp::criticalTemperature(),
            Comp::criticalPressure(), Comp::criticalVolume(), Comp::acentricFactor()};
}

const std::vector<ComponentEntry>& componentDatabase()
{
    static const std::vector<ComponentEntry> db{
        makeEntry<Opm::C1<double>>(),
        makeEntry<Opm::C10<double>>(),
        makeEntry<Opm::SimpleCO2<double>>(),
    };
    return db;
}

const ComponentEntry& lookupComponent(const std::string& name)
{
    std::string upper;
    for (char c : name) { upper.push_back(std::toupper(static_cast<unsigned char>(c))); }
    // accept the deck-style aliases IdealGasCaloricData understands
    const auto canonical = [&]() -> std::string {
        if (upper == "C1" || upper == "CH4" || upper == "METHANE") { return "C1"; }
        if (upper == "C10" || upper == "NC10" || upper == "DECANE" || upper == "N-DECANE") { return "C10"; }
        if (upper == "CO2" || upper == "CARBONDIOXIDE") { return "CO2"; }
        return upper;
    }();
    for (const auto& entry : componentDatabase()) {
        if (entry.name == canonical) {
            return entry;
        }
    }
    throw std::runtime_error("unknown component '" + name +
                             "' — run: opmflash --list-components");
}

double pairKij(std::string_view a, std::string_view b)
{
    const auto is = [](std::string_view x, std::string_view y,
                       std::string_view p, std::string_view q)
    { return (x == p && y == q) || (x == q && y == p); };
    if (is(a, b, "C1", "C10")) { return 0.0411; }
    if (is(a, b, "C1", "CO2")) { return 0.10; }
    if (is(a, b, "CO2", "C10")) { return 0.10; }
    return 0.0;
}

// ── command line ───────────────────────────────────────────────────────────
struct Sweep {
    char var = '\0';          // 'p' | 't' | 'h'
    double start = 0.0, stop = 0.0, step = 0.0;
};

struct Cli {
    std::string spec = "pt";              // pt | ph
    double pBar = 50.0;
    std::optional<double> tKelvin;        // pt: required
    std::optional<double> hMolar;         // ph: required (unless sweeping h)
    std::vector<std::string> components{"C1", "C10"};
    std::vector<double> z{0.5, 0.5};
    std::string eos = "PR";
    std::string enthalpyModel = "eos_departure";
    double tempMin = 270.0;
    double tempMax = 460.0;
    double tolerance = 1e-6;
    int maxIterations = 100;
    double ptTolerance = 1e-8;
    std::string twoPhaseMethod = "ssi";
    std::optional<Sweep> sweep;
    std::string format = "table";         // table | csv
    int verbosity = 0;
    bool selftest = false;
};

std::vector<std::string> splitCommas(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) { out.push_back(item); }
    return out;
}

Sweep parseSweep(const std::string& s)
{
    // var=start:stop:step
    const auto eq = s.find('=');
    const auto c1 = s.find(':', eq);
    const auto c2 = s.find(':', c1 + 1);
    if (eq != 1 || c1 == std::string::npos || c2 == std::string::npos ||
        (s[0] != 'p' && s[0] != 't' && s[0] != 'h')) {
        throw std::runtime_error("--sweep expects p|t|h=start:stop:step, got '" + s + "'");
    }
    return {s[0],
            std::stod(s.substr(eq + 1, c1 - eq - 1)),
            std::stod(s.substr(c1 + 1, c2 - c1 - 1)),
            std::stod(s.substr(c2 + 1))};
}

void listComponents()
{
    std::cout <<
        "Available components (parameters from this library's component\n"
        "classes; heat capacity from the IdealGasCaloricData presets, cubic\n"
        "fits to the reference-EoS ideal-gas cp, valid 250-600 K):\n\n"
        "  name  aliases          MW[g/mol]   Tc[K]    Pc[bar]  Vc[m3/kmol]  omega\n";
    std::cout << std::fixed;
    for (const auto& e : componentDatabase()) {
        std::string aliases;
        if (e.name == "C1") { aliases = "CH4, METHANE"; }
        else if (e.name == "C10") { aliases = "nC10, DECANE"; }
        else if (e.name == "CO2") { aliases = "CARBONDIOXIDE"; }
        std::cout << "  " << std::setw(5) << std::left << e.name
                  << std::setw(17) << aliases << std::right
                  << std::setw(9) << std::setprecision(3) << e.molarMassGramPerMol
                  << std::setw(9) << std::setprecision(2) << e.criticalTemperature
                  << std::setw(9) << e.criticalPressure / 1e5
                  << std::setw(11) << std::setprecision(4) << e.criticalVolume
                  << std::setw(8) << std::setprecision(3) << e.acentricFactor
                  << '\n';
    }
    std::cout <<
        "\nBinary interaction coefficients (applied automatically per pair;\n"
        "standard PR literature values, the same ones this library's\n"
        "compositional tests use):\n"
        "  C1/C10   0.0411\n"
        "  C1/CO2   0.10\n"
        "  CO2/C10  0.10\n"
        "  (any pair not listed: 0)\n\n"
        "Select with --components, e.g.  --components C1,C10  or\n"
        "--components CO2,C1,C10  (2 or 3 components; --z must match).\n";
}

void printUsage()
{
    std::cout <<
        "opmflash — flash a compositional state from the command line\n"
        "\n"
        "Runs the isothermal (P-T) or isenthalpic (P-H) flash on a synthetic\n"
        "mixture through the same fluid system and solvers the compositional\n"
        "simulator uses. No deck, no C++: state in, split out.\n"
        "\n"
        "usage:\n"
        "  opmflash --spec pt --p <bar> --t <K>     --z <z,z[,z]> --components <c,c[,c]> [options]\n"
        "  opmflash --spec ph --p <bar> --h <J/mol> --z <z,z[,z]> --components <c,c[,c]> [options]\n"
        "  opmflash --list-components\n"
        "  opmflash --selftest\n"
        "\n"
        "specifications:\n"
        "  pt   isothermal: flash at (p, T, z); reports the phase split and the\n"
        "       molar enthalpy of the flashed state (both enthalpy models) —\n"
        "       also the way to LOOK UP the enthalpy for a later ph run\n"
        "  ph   isenthalpic: finds the temperature at which the flashed mixture\n"
        "       attains --h, then reports the same quantities\n"
        "\n"
        "options:\n"
        "  --components      comma list; see --list-components       [C1,C10]\n"
        "  --z               feed mole fractions, must sum to 1      [0.5,0.5]\n"
        "  --eos             PR | SRK | ...                          [PR]\n"
        "  --enthalpy-model  caloric | eos_departure                 [eos_departure]\n"
        "                    caloric = ideal-gas cp integral only;\n"
        "                    eos_departure adds the EoS residual (real-fluid)\n"
        "  --tmin, --tmax    P-H temperature search bracket [K]      [270, 460]\n"
        "  --tol             P-H solver tolerance                    [1e-6]\n"
        "  --maxiter         P-H solver iteration cap                [100]\n"
        "  --pt-tol          inner isothermal flash tolerance        [1e-8]\n"
        "  --method          inner two-phase method                  [ssi]\n"
        "  --sweep           p|t|h=start:stop:step — one variable;\n"
        "                    t-sweeps need --spec pt, h-sweeps --spec ph\n"
        "  --format          table | csv (csv: machine-readable,\n"
        "                    provenance in '#' header lines)         [table]\n"
        "  --verbosity       0..5 (passed to the solvers)            [0]\n"
        "\n"
        "units: --p bar · --t K · --h J/mol. All enthalpies are MOLAR [J/mol]\n"
        "against the caloric datum H(298.15 K) = 0.\n"
        "\n"
        "examples:\n"
        "  # what phase is this state in, and what is its enthalpy?\n"
        "  opmflash --spec pt --p 50 --t 300 --z 0.5,0.5 --components C1,C10\n"
        "\n"
        "  # invert: which temperature gives that enthalpy back?\n"
        "  opmflash --spec ph --p 50 --h -25048 --z 0.5,0.5 --components C1,C10\n"
        "\n"
        "  # benchmark curve H -> T as CSV (for plots or reference comparison)\n"
        "  opmflash --spec ph --p 50 --z 0.5,0.5 --components C1,C10 \\\n"
        "           --sweep h=-40000:5000:1000 --format csv > curve.csv\n"
        "\n"
        "  # ternary with CO2\n"
        "  opmflash --spec pt --p 50 --t 320 --z 0.4,0.2,0.4 --components CO2,C1,C10\n";
}

Cli parseCli(int argc, char** argv)
{
    Cli cli;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const auto value = [&]() -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value after " + arg);
            }
            return argv[++i];
        };
        if (arg == "--help" || arg == "-h") { printUsage(); std::exit(EXIT_SUCCESS); }
        else if (arg == "--list-components") { listComponents(); std::exit(EXIT_SUCCESS); }
        else if (arg == "--selftest") { cli.selftest = true; }
        else if (arg == "--spec") { cli.spec = value(); }
        else if (arg == "--p") { cli.pBar = std::stod(value()); }
        else if (arg == "--t") { cli.tKelvin = std::stod(value()); }
        else if (arg == "--h") { cli.hMolar = std::stod(value()); }
        else if (arg == "--components") { cli.components = splitCommas(value()); }
        else if (arg == "--z") {
            cli.z.clear();
            for (const auto& v : splitCommas(value())) { cli.z.push_back(std::stod(v)); }
        }
        else if (arg == "--eos") { cli.eos = value(); }
        else if (arg == "--enthalpy-model") { cli.enthalpyModel = value(); }
        else if (arg == "--tmin") { cli.tempMin = std::stod(value()); }
        else if (arg == "--tmax") { cli.tempMax = std::stod(value()); }
        else if (arg == "--tol") { cli.tolerance = std::stod(value()); }
        else if (arg == "--maxiter") { cli.maxIterations = std::stoi(value()); }
        else if (arg == "--pt-tol") { cli.ptTolerance = std::stod(value()); }
        else if (arg == "--method") { cli.twoPhaseMethod = value(); }
        else if (arg == "--sweep") { cli.sweep = parseSweep(value()); }
        else if (arg == "--format") { cli.format = value(); }
        else if (arg == "--verbosity") { cli.verbosity = std::stoi(value()); }
        else {
            throw std::runtime_error("unknown option '" + arg + "' (see --help)");
        }
    }
    return cli;
}

// ── the flash, for a compile-time component count ──────────────────────────
template <int N>
struct FlashRunner {
    using FluidSystem = Opm::GenericOilGasWaterFluidSystem<double, N, false>;
    using FlashEval = Opm::DenseAd::Evaluation<double, N + 1>;
    using FluidState = Opm::CompositionalFluidState<FlashEval, FluidSystem>;
    using PtFlash = Opm::PTFlash<double, FluidSystem>;
    using PhFlash = Opm::PHFlash<double, FluidSystem>;
    using Enthalpy = Opm::MixtureEnthalpy<double, FluidSystem>;

    struct Result {
        bool ok;
        bool singlePhase;
        double temperature;
        double L;
        std::array<double, N> x, y, K;
        double hCaloric, hDeparture;
    };

    static void setup(const Cli& cli)
    {
        FluidSystem::init();
        std::vector<std::string_view> names;
        for (const auto& c : cli.components) {
            const auto& e = lookupComponent(c);
            names.push_back(e.name);
            FluidSystem::addComponent(typename FluidSystem::ComponentParam{
                std::string(e.name), e.molarMassGramPerMol, e.criticalTemperature,
                e.criticalPressure, e.criticalVolume, e.acentricFactor});
        }
        // packed lower-triangle order: (1,0), (2,0), (2,1), ...
        std::vector<double> bic;
        for (int row = 1; row < N; ++row) {
            for (int col = 0; col < row; ++col) {
                bic.push_back(pairKij(names[row], names[col]));
            }
        }
        FluidSystem::setInteractionCoefficients(std::move(bic));
        FluidSystem::initEnthalpyFromComponentNames();
    }

    static FluidState makeState(double p, double t, const std::vector<double>& z)
    {
        FluidState fs;
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, p);
        }
        fs.setTemperature(t);
        for (int compIdx = 0; compIdx < N; ++compIdx) {
            fs.setMoleFraction(compIdx, z[compIdx]);
        }
        for (int compIdx = 0; compIdx < N; ++compIdx) {
            fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
        }
        fs.setLvalue(-1.0);
        return fs;
    }

    static typename PhFlash::Config phConfig(const Cli& cli)
    {
        typename PhFlash::Config cfg;
        cfg.cpTable = cpTable();
        cfg.model = Opm::enthalpyModelFromString(cli.enthalpyModel);
        cfg.tempMin = cli.tempMin;
        cfg.tempMax = cli.tempMax;
        cfg.tolerance = cli.tolerance;
        cfg.maxIterations = cli.maxIterations;
        return cfg;
    }

    static Opm::CpTable<double, N> cpTable()
    {
        Opm::CpTable<double, N> table;
        for (int compIdx = 0; compIdx < N; ++compIdx) {
            table[compIdx] = Opm::IdealGasCaloricData<double>::byName(
                std::string(FluidSystem::componentName(compIdx)));
        }
        return table;
    }

    static Result collect(FluidState& fs, bool singlePhase, bool ok,
                          const EOSType& eos)
    {
        Result r{};
        r.ok = ok;
        r.singlePhase = singlePhase;
        r.temperature = Opm::getValue(fs.temperature(0));
        r.L = Opm::getValue(fs.L());
        for (int c = 0; c < N; ++c) {
            r.x[c] = Opm::getValue(fs.moleFraction(FluidSystem::oilPhaseIdx, c));
            r.y[c] = Opm::getValue(fs.moleFraction(FluidSystem::gasPhaseIdx, c));
            r.K[c] = r.x[c] > 0.0 ? r.y[c] / r.x[c] : 0.0;
        }
        if (ok) {
            const auto table = cpTable();
            constexpr double T0 = 298.15;
            r.hCaloric = Opm::getValue(Enthalpy::mixtureEnthalpy(
                fs, table, T0, eos, Opm::EnthalpyModel::caloric));
            r.hDeparture = Opm::getValue(Enthalpy::mixtureEnthalpy(
                fs, table, T0, eos, Opm::EnthalpyModel::eos_departure));
        }
        return r;
    }

    static Result runPt(const Cli& cli, double pPa, double t)
    {
        const auto eos = Opm::CompositionalConfig::eosTypeFromString(cli.eos);
        FluidState fs = makeState(pPa, t, cli.z);
        const bool singlePhase =
            PtFlash::solve(fs, cli.twoPhaseMethod, cli.ptTolerance, eos, cli.verbosity);
        return collect(fs, singlePhase, true, eos);
    }

    static Result runPh(const Cli& cli, double pPa, double h)
    {
        const auto eos = Opm::CompositionalConfig::eosTypeFromString(cli.eos);
        // the initial temperature only seeds Wilson K; PHFlash brackets over T
        FluidState fs = makeState(pPa, 0.5 * (cli.tempMin + cli.tempMax), cli.z);
        const bool ok = PhFlash::solve(fs, h, phConfig(cli), cli.twoPhaseMethod,
                                       cli.ptTolerance, eos, cli.verbosity);
        const bool singlePhase = !(Opm::getValue(fs.L()) > 0.0 &&
                                   Opm::getValue(fs.L()) < 1.0);
        return collect(fs, singlePhase, ok, eos);
    }
};

// ── output ─────────────────────────────────────────────────────────────────
template <int N>
void printHeader(const Cli& cli)
{
    const char* prefix = (cli.format == "csv") ? "# " : "";
    std::cout << prefix << "opmflash · spec " << cli.spec << " · eos " << cli.eos
              << " · enthalpy-model " << cli.enthalpyModel
              << " · datum H(298.15 K) = 0 · enthalpies J/mol\n";
    if (cli.format == "csv") {
        std::cout << "p_bar,T_K,h_spec";
        for (int c = 0; c < N; ++c) { std::cout << ",z" << c; }
        std::cout << ",single_phase,L";
        for (int c = 0; c < N; ++c) { std::cout << ",x" << c; }
        for (int c = 0; c < N; ++c) { std::cout << ",y" << c; }
        for (int c = 0; c < N; ++c) { std::cout << ",K" << c; }
        std::cout << ",H_caloric,H_departure,ok\n";
    }
}

template <int N>
void printRow(const Cli& cli, double pPa, std::optional<double> hSpec,
              const typename FlashRunner<N>::Result& r)
{
    if (cli.format == "csv") {
        std::cout << std::setprecision(10)
                  << pPa / 1e5 << ',' << r.temperature << ','
                  << (hSpec ? std::to_string(*hSpec) : "");
        for (int c = 0; c < N; ++c) { std::cout << ',' << cli.z[c]; }
        std::cout << ',' << (r.singlePhase ? 1 : 0) << ',' << r.L;
        for (int c = 0; c < N; ++c) { std::cout << ',' << r.x[c]; }
        for (int c = 0; c < N; ++c) { std::cout << ',' << r.y[c]; }
        for (int c = 0; c < N; ++c) { std::cout << ',' << r.K[c]; }
        std::cout << ',' << r.hCaloric << ',' << r.hDeparture
                  << ',' << (r.ok ? 1 : 0) << '\n';
        return;
    }

    std::cout << std::fixed << std::setprecision(4);
    if (!r.ok) {
        std::cout << "  NO SOLUTION on the temperature bracket [" << cli.tempMin
                  << ", " << cli.tempMax << "] K for h = "
                  << (hSpec ? *hSpec : 0.0) << " J/mol at p = " << pPa / 1e5
                  << " bar — widen --tmin/--tmax or check the input\n";
        return;
    }
    std::cout << "  p = " << pPa / 1e5 << " bar";
    if (hSpec) { std::cout << "   h_spec = " << *hSpec << " J/mol"; }
    std::cout << "\n  T = " << r.temperature << " K   "
              << (r.singlePhase ? "single-phase" : "two-phase")
              << "   L = " << r.L << "\n  ";
    for (int c = 0; c < N; ++c) {
        std::cout << "x" << c << " = " << r.x[c] << "  ";
    }
    for (int c = 0; c < N; ++c) {
        std::cout << "y" << c << " = " << r.y[c] << "  ";
    }
    std::cout << "\n  H caloric = " << r.hCaloric
              << "   H eos_departure = " << r.hDeparture << "  [J/mol]\n\n";
}

// ── drivers ────────────────────────────────────────────────────────────────
template <int N>
int runCases(const Cli& cli)
{
    FlashRunner<N>::setup(cli);
    printHeader<N>(cli);

    const double pPa = cli.pBar * 1e5;
    bool allOk = true;

    const auto runOne = [&](double p, std::optional<double> t, std::optional<double> h) {
        if (cli.spec == "pt") {
            if (!t) { throw std::runtime_error("--spec pt requires --t (or a t-sweep)"); }
            const auto r = FlashRunner<N>::runPt(cli, p, *t);
            printRow<N>(cli, p, std::nullopt, r);
            allOk = allOk && r.ok;
        }
        else {
            if (!h) { throw std::runtime_error("--spec ph requires --h (or an h-sweep)"); }
            const auto r = FlashRunner<N>::runPh(cli, p, *h);
            printRow<N>(cli, p, h, r);
            allOk = allOk && r.ok;
        }
    };

    if (!cli.sweep) {
        runOne(pPa, cli.tKelvin, cli.hMolar);
        return allOk ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    const auto& sw = *cli.sweep;
    if ((sw.stop - sw.start) * sw.step <= 0.0) {
        throw std::runtime_error("--sweep step does not advance from start to stop");
    }
    if ((cli.spec == "pt" && sw.var == 'h') || (cli.spec == "ph" && sw.var == 't')) {
        throw std::runtime_error(std::string("--sweep over '") + sw.var +
                                 "' does not apply to --spec " + cli.spec);
    }
    for (double v = sw.start; (sw.step > 0.0) ? v <= sw.stop : v >= sw.stop; v += sw.step) {
        switch (sw.var) {
        case 'p': runOne(v * 1e5, cli.tKelvin, cli.hMolar); break;
        case 't': runOne(pPa, v, cli.hMolar); break;
        case 'h': runOne(pPa, cli.tKelvin, v); break;
        }
    }
    return allOk ? EXIT_SUCCESS : EXIT_FAILURE;
}

int dispatchOnComponentCount(const Cli& cli)
{
    if (cli.components.size() != cli.z.size()) {
        throw std::runtime_error("--components and --z must have the same length");
    }
    double sum = 0.0;
    for (double v : cli.z) { sum += v; }
    if (std::abs(sum - 1.0) > 1e-9) {
        throw std::runtime_error("--z must sum to 1");
    }

    switch (cli.components.size()) {
    case 2: return runCases<2>(cli);
    case 3: return runCases<3>(cli);
    default:
        throw std::runtime_error("supported component counts: 2, 3");
    }
}

// The canonical F1 round-trip as an executable check of the WHOLE CLI
// surface (parsing aside): flash at the anchor temperature, hand the
// resulting enthalpy to the isenthalpic solver, require the temperature
// back within 1e-4 K; then exercise the F2 ternary the same way.
int selftest()
{
    int failures = 0;
    const auto check = [&](const char* name, const Cli& cli, double tAnchor) {
        const auto verdict = [&](bool pass, double tBack) {
            std::cout << (pass ? "  PASS  " : "  FAIL  ") << name
                      << "  T = " << tAnchor << " -> " << tBack << " K\n";
            if (!pass) { ++failures; }
        };
        if (cli.components.size() == 2) {
            FlashRunner<2>::setup(cli);
            const auto pt = FlashRunner<2>::runPt(cli, cli.pBar * 1e5, tAnchor);
            const auto ph = FlashRunner<2>::runPh(cli, cli.pBar * 1e5, pt.hDeparture);
            verdict(ph.ok && std::abs(ph.temperature - tAnchor) < 1e-4, ph.temperature);
        }
        else {
            FlashRunner<3>::setup(cli);
            const auto pt = FlashRunner<3>::runPt(cli, cli.pBar * 1e5, tAnchor);
            const auto ph = FlashRunner<3>::runPh(cli, cli.pBar * 1e5, pt.hDeparture);
            verdict(ph.ok && std::abs(ph.temperature - tAnchor) < 1e-4, ph.temperature);
        }
    };

    Cli f1;
    f1.components = {"C1", "C10"};
    f1.z = {0.5, 0.5};
    check("F1 binary round-trip (C1/C10, 50 bar, 300 K)", f1, 300.0);

    Cli f2;
    f2.components = {"CO2", "C1", "C10"};
    f2.z = {0.4, 0.2, 0.4};
    check("F2 ternary round-trip (CO2/C1/C10, 50 bar, 320 K)", f2, 320.0);

    std::cout << (failures == 0 ? "opmflash selftest: ALL GREEN\n"
                                : "opmflash selftest: FAILURES\n");
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // anonymous namespace

int main(int argc, char** argv)
{
    try {
        if (argc <= 1) {
            printUsage();
            return EXIT_FAILURE;
        }
        const Cli cli = parseCli(argc, argv);
        return cli.selftest ? selftest() : dispatchOnComponentCount(cli);
    }
    catch (const std::exception& e) {
        std::cerr << "opmflash: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
