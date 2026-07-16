// Integration tests: run mesh_generator on each model listed in
// reference_hashes.txt and check the SHA-256 of the generated file against a
// stored reference. Because the exact predicates make the pipeline fully
// deterministic, these hashes must be identical on every OS/compiler -- so a
// mismatch means a portability/nondeterminism regression.
//
// The reference hashes are produced by tests/compute_hashes.py (which uses the
// same SHA-256), so adding a new model is: drop the .off in models/, run the
// python helper, paste the line into reference_hashes.txt.
//
// Models listed in the manifest but not present on disk are skipped (with a
// warning) rather than failing -- this lets CI run whatever subset of models is
// committed to the repo.

#include <catch2/catch_test_macros.hpp>

#include "sha256.h"
#include "tet_orientation.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Entry {
    std::string file;   // e.g. "sphere.off"
    std::string mode;   // "default" | "blackfaces" | "skin"
    std::string sha256; // expected 64-hex digest of the generated file
};

std::vector<Entry> load_manifest() {
    std::vector<Entry> entries;
    std::ifstream in(VRTEST_HASHES_FILE);
    std::string line;
    while (std::getline(in, line)) {
        const auto p = line.find_first_not_of(" \t\r\n");
        if (p == std::string::npos || line[p] == '#')
            continue;
        std::istringstream ss(line);
        Entry e;
        ss >> e.file >> e.mode >> e.sha256;
        if (!e.file.empty() && !e.mode.empty() && !e.sha256.empty())
            entries.push_back(e);
    }
    return entries;
}

std::string flags_for(const std::string& mode) {
    if (mode == "blackfaces") return "-b";
    if (mode == "skin") return "-s";
    if (mode == "tet") return "-t";
    return ""; // default -> volume.msh
}

std::string output_for(const std::string& mode) {
    if (mode == "blackfaces") return "black_faces.off";
    if (mode == "skin") return "skin.off";
    if (mode == "tet") return "volume.tet";
    return "volume.msh";
}

// Run mesh_generator on 'model' in a private temp dir (the tool writes its
// output to the current directory with a fixed name) and return the SHA-256 of
// the generated file. Empty string if no output was produced.
std::string run_and_hash(const std::string& model, const std::string& mode) {
    static int counter = 0;
    fs::path work = fs::temp_directory_path() / ("vrtest_" + std::to_string(counter++));
    std::error_code ec;
    fs::remove_all(work, ec);
    fs::create_directories(work, ec);

    const std::string bin = VRTEST_MESH_GENERATOR;
    const std::string flags = flags_for(mode);
    std::string cmd;
#ifdef _WIN32
    cmd = "cd /d \"" + work.string() + "\" && \"" + bin + "\" \"" + model + "\" " + flags + " > nul 2>&1";
#else
    cmd = "cd '" + work.string() + "' && '" + bin + "' '" + model + "' " + flags + " > /dev/null 2>&1";
#endif
    std::system(cmd.c_str());

    const fs::path out = work / output_for(mode);
    std::string digest = fs::exists(out) ? vrtest::sha256_file(out.string()) : std::string();
    // Tetrahedralization output must be a valid simplicial complex: no overlapping /
    // inconsistently-oriented tets. Combinatorial, so cheap even on the big models.
    if (mode == "tet" && fs::exists(out)) {
        const vrtest::OrientationResult orient = vrtest::check_tet_orientation_file(out.string());
        INFO("orientation " << model << ": same_winding=" << orient.same_winding
                            << " over_shared=" << orient.over_shared << " (" << orient.error << ")");
        CHECK(orient.ok);
    }
    fs::remove_all(work, ec);
    return digest;
}

// Run `mesh_generator -t -r` then the independent verify_tracking on 'model' and
// return the last line of verify_tracking's output ("TRACKING OK" on success). The
// exact-rational verifier is the source of truth for the face-provenance tracking.
std::string run_tracking_check(const std::string& model, bool strict) {
    static int counter = 1000000;
    fs::path work = fs::temp_directory_path() / ("vrtrack_" + std::to_string(counter++));
    std::error_code ec;
    fs::remove_all(work, ec);
    fs::create_directories(work, ec);

    const std::string bin = VRTEST_MESH_GENERATOR;
    const std::string vt = VRTEST_VERIFY_TRACKING;
    const std::string strictflag = strict ? " --strict" : "";
    std::string gen, chk;
    // -a (keep_all_cells): skip the in/out min-cut and tetrahedralize the whole domain,
    // so every input triangle is realized as a face between two kept cells and the
    // surface tracking is exact (no thin feature shaved off, no near-coplanar loss).
#ifdef _WIN32
    gen = "cd /d \"" + work.string() + "\" && \"" + bin + "\" \"" + model + "\" -t -r -a > nul 2>&1";
    chk = "cd /d \"" + work.string() + "\" && \"" + vt + "\" \"" + model + "\" volume.tet" +
          strictflag + " > vt.out 2>&1";
#else
    gen = "cd '" + work.string() + "' && '" + bin + "' '" + model + "' -t -r -a > /dev/null 2>&1";
    chk = "cd '" + work.string() + "' && '" + vt + "' '" + model + "' volume.tet" + strictflag +
          " > vt.out 2>&1";
#endif
    std::system(gen.c_str());
    std::system(chk.c_str());

    std::string out, line;
    std::ifstream f((work / "vt.out").string());
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    fs::remove_all(work, ec);
    return out;
}

// Like run_tracking_check, but also inserts extra edges/points (-e/-p, which imply -a)
// and checks them with verify_tracking's --edges/--points exact provenance tests.
std::string run_tracking_check_ep(
    const std::string& model, const std::string& edges, const std::string& points) {
    static int counter = 2000000;
    fs::path work = fs::temp_directory_path() / ("vrtrackep_" + std::to_string(counter++));
    std::error_code ec;
    fs::remove_all(work, ec);
    fs::create_directories(work, ec);

    const std::string bin = VRTEST_MESH_GENERATOR;
    const std::string vt = VRTEST_VERIFY_TRACKING;
    std::string gen, chk;
#ifdef _WIN32
    gen = "cd /d \"" + work.string() + "\" && \"" + bin + "\" \"" + model + "\" -t -r -a -e \"" +
          edges + "\" -p \"" + points + "\" > nul 2>&1";
    chk = "cd /d \"" + work.string() + "\" && \"" + vt + "\" \"" + model + "\" volume.tet --strict"
          " --edges \"" + edges + "\" --points \"" + points + "\" > vt.out 2>&1";
#else
    gen = "cd '" + work.string() + "' && '" + bin + "' '" + model + "' -t -r -a -e '" + edges +
          "' -p '" + points + "' > /dev/null 2>&1";
    chk = "cd '" + work.string() + "' && '" + vt + "' '" + model + "' volume.tet --strict --edges '" +
          edges + "' --points '" + points + "' > vt.out 2>&1";
#endif
    std::system(gen.c_str());
    std::system(chk.c_str());

    std::ifstream f((work / "vt.out").string());
    std::ostringstream ss;
    ss << f.rdbuf();
    std::string out = ss.str();
    fs::remove_all(work, ec);
    return out;
}

// Deterministic random coordinate in [lo,hi). Uses mt19937 (fixed, standard-specified
// sequence) with a manual mantissa build -- NOT uniform_real_distribution, whose output is
// implementation-defined -- so the generated inputs are identical on every platform.
double rnd(std::mt19937_64& g, double lo, double hi) {
    const double u = (double)(g() >> 11) / (double)(1ull << 53); // 53-bit fraction in [0,1)
    return lo + (hi - lo) * u;
}

void write_random_points(const fs::path& path, uint32_t n, uint64_t seed) {
    std::mt19937_64 g(seed);
    std::ofstream f(path, std::ios::binary);
    f.precision(17);
    f << "OFF\n" << n << " 0 0\n";
    for (uint32_t i = 0; i < n; i++)
        f << rnd(g, 0.05, 0.95) << " " << rnd(g, 0.05, 0.95) << " " << rnd(g, 0.05, 0.95) << "\n";
}

// n short scattered segments (a random start + a small random offset), so the edges do not
// densely cross each other.
void write_random_edges(const fs::path& path, uint32_t n, uint64_t seed) {
    std::mt19937_64 g(seed);
    std::ofstream f(path, std::ios::binary);
    f.precision(17);
    f << "OFF\n" << (2 * n) << " " << n << " 0\n";
    for (uint32_t i = 0; i < n; i++) {
        const double ax = rnd(g, 0.05, 0.95), ay = rnd(g, 0.05, 0.95), az = rnd(g, 0.05, 0.95);
        f << ax << " " << ay << " " << az << "\n";
        f << (ax + rnd(g, -0.03, 0.03)) << " " << (ay + rnd(g, -0.03, 0.03)) << " "
          << (az + rnd(g, -0.03, 0.03)) << "\n";
    }
    for (uint32_t i = 0; i < n; i++) f << "2 " << (2 * i) << " " << (2 * i + 1) << "\n";
}

// Extract the unique undirected edges of an OFF triangle mesh into an edge OFF (vertices +
// "2 i j" records). Returns false if the model can't be read.
bool extract_edges_off(const fs::path& model, const fs::path& out) {
    std::ifstream in(model);
    if (!in) return false;
    std::string tok;
    in >> tok;
    if (tok != "OFF") return false;
    uint32_t nv, nf, ne;
    in >> nv >> nf >> ne;
    std::vector<std::array<double, 3>> V(nv);
    for (uint32_t i = 0; i < nv; i++) in >> V[i][0] >> V[i][1] >> V[i][2];
    std::set<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t i = 0; i < nf; i++) {
        uint32_t k, a, b, c;
        in >> k >> a >> b >> c;
        auto add = [&](uint32_t x, uint32_t y) {
            edges.insert(x < y ? std::make_pair(x, y) : std::make_pair(y, x));
        };
        add(a, b);
        add(b, c);
        add(c, a);
    }
    if (!in) return false;
    std::ofstream f(out, std::ios::binary);
    f.precision(17);
    f << "OFF\n" << nv << " " << edges.size() << " 0\n";
    for (const auto& v : V) f << v[0] << " " << v[1] << " " << v[2] << "\n";
    for (const auto& e : edges) f << "2 " << e.first << " " << e.second << "\n";
    return true;
}

// Run mesh_generator (-t -r -a, plus any of surface / -e / -p) then verify_tracking with the
// matching checks, and return verify_tracking's output. Empty paths omit that input.
std::string run_combo(
    const std::string& surface, const std::string& edges, const std::string& points) {
    static int counter = 3000000;
    fs::path work = fs::temp_directory_path() / ("vrcombo_" + std::to_string(counter++));
    std::error_code ec;
    fs::remove_all(work, ec);
    fs::create_directories(work, ec);

    const std::string bin = VRTEST_MESH_GENERATOR;
    const std::string vt = VRTEST_VERIFY_TRACKING;
    auto q = [](const std::string& s) { return "\"" + s + "\""; };
    std::string gflags = " -t -r -a";
    if (!surface.empty()) gflags += " " + q(surface);
    if (!edges.empty()) gflags += " -e " + q(edges);
    if (!points.empty()) gflags += " -p " + q(points);
    std::string cflags = " " + q(surface.empty() ? std::string("none") : surface) + " volume.tet --strict";
    if (!edges.empty()) cflags += " --edges " + q(edges);
    if (!points.empty()) cflags += " --points " + q(points);

    std::string gen, chk;
#ifdef _WIN32
    gen = "cd /d \"" + work.string() + "\" && \"" + bin + "\"" + gflags + " > nul 2>&1";
    chk = "cd /d \"" + work.string() + "\" && \"" + vt + "\"" + cflags + " > vt.out 2>&1";
#else
    gen = "cd '" + work.string() + "' && '" + bin + "'" + gflags + " > /dev/null 2>&1";
    chk = "cd '" + work.string() + "' && '" + vt + "'" + cflags + " > vt.out 2>&1";
#endif
    std::system(gen.c_str());
    std::system(chk.c_str());

    std::ifstream f((work / "vt.out").string());
    std::ostringstream ss;
    ss << f.rdbuf();
    std::string out = ss.str();
    fs::remove_all(work, ec);
    return out;
}

} // namespace

// A small, representative set covering every tracking behavior: exactly-coplanar
// clusters (cube_subdiv, upsample_box), disconnected coplanar planes (two_cubes),
// coplanar overlap (cube_on_cube), curved surfaces (sphere, double_sphere), an OPEN
// boundary (open_box, the regression guard for open-surface preservation), and
// dropped-degenerate-triangle indexing (112856simplified). Because the tracking runs
// with -a (keep_all_cells), every well-formed input -- flat, curved, or open -- is
// reconstructed EXACTLY, so all of these are checked strictly (under=0).
//
// 112856simplified is the lone exception: it contains overlapping coplanar input
// triangles (a defective "simplified" mesh -- one 241-triangle group whose summed
// area exceeds the area of their union). The mesher correctly realizes the union, so
// the sum-of-areas coverage check under-counts by construction; it is checked
// non-strictly (soundness + flag only). Large meshes (bunny/Octocat/...) are excluded:
// the exact-rational verifier is slow on them with the no-GMP bignum backend, and they
// add no new behavior -- verify them locally with tests/verify_tracking.
struct TrackModel { const char* file; bool strict; };
static const TrackModel kTrackModels[] = {
    {"cube_subdiv.off", true}, {"two_cubes.off", true}, {"cube_on_cube.off", true},
    {"upsample_box.off", true}, {"sphere.off", true}, {"double_sphere.off", true},
    {"open_box.off", true}, {"112856simplified.off", false},
};

TEST_CASE("integration: -t face-provenance tracking is exact", "[integration][tracking]") {
    for (const auto& tm : kTrackModels) {
        DYNAMIC_SECTION(tm.file << (tm.strict ? " [strict]" : " [lenient]")) {
            const fs::path model = fs::path(VRTEST_MODELS_DIR) / tm.file;
            if (!fs::exists(model)) {
                WARN("model not present, skipping: " << tm.file);
            } else {
                const std::string out = run_tracking_check(model.string(), tm.strict);
                INFO(out);
                REQUIRE(out.find("TRACKING OK") != std::string::npos);
            }
        }
    }
}

// Random 3D polylines + points inserted into a flat and a curved mesh. The forcing
// triangles pin each edge/point into the (keep-all-cells) output, and verify_tracking
// checks in exact rational that every input edge is tiled by output tet edges and every
// input point equals an output vertex. The <model>_edges.off / <model>_points.off files
// were generated once with a fixed RNG seed (see the commit) so the test is reproducible.
struct EdgePointModel { const char* file; const char* edges; const char* points; };
static const EdgePointModel kEdgePointModels[] = {
    {"cube_subdiv.off", "cube_subdiv_edges.off", "cube_subdiv_points.off"},
    {"sphere.off", "sphere_edges.off", "sphere_points.off"},
};

TEST_CASE("integration: -t inserted edges and points are tracked exactly",
    "[integration][tracking]") {
    for (const auto& em : kEdgePointModels) {
        DYNAMIC_SECTION(em.file) {
            const fs::path model = fs::path(VRTEST_MODELS_DIR) / em.file;
            const fs::path edges = fs::path(VRTEST_MODELS_DIR) / em.edges;
            const fs::path points = fs::path(VRTEST_MODELS_DIR) / em.points;
            if (!fs::exists(model) || !fs::exists(edges) || !fs::exists(points)) {
                WARN("model/edges/points not present, skipping: " << em.file);
            } else {
                const std::string out =
                    run_tracking_check_ep(model.string(), edges.string(), points.string());
                INFO(out);
                REQUIRE(out.find("TRACKING OK") != std::string::npos);
            }
        }
    }
}

// Every non-empty combination of the three optional inputs -- surface triangles, edges,
// points -- with large random edge/point sets (10k each). Confirms both that each
// combination runs and that the tracking is EXACT (verify_tracking with --strict, --edges,
// --points): every surface group tiles exactly, every input edge is tiled by output tet
// edges, and every input point is an output vertex. Inputs are generated with a fixed
// deterministic RNG (see write_random_*), so the test is reproducible.
//
// Release-only (#ifndef NDEBUG skips it): the exact-rational + BSP work on ~millions of
// tets is far too slow in an unoptimized Debug build. The smaller tracking tests above
// still exercise the same code paths in Debug.
#ifdef NDEBUG
TEST_CASE("integration: -t all input combinations (surface/edges/points) tracked exactly",
    "[integration][tracking][large]") {
    const uint32_t N = 10000;
    const fs::path surface = fs::path(VRTEST_MODELS_DIR) / "cube_subdiv.off";
    const fs::path dir = fs::temp_directory_path() / "vrcombo_inputs";
    std::error_code ec;
    fs::create_directories(dir, ec);
    const fs::path edges = dir / "edges10k.off";
    const fs::path points = dir / "points10k.off";
    write_random_edges(edges, N, 12345);
    write_random_points(points, N, 67890);

    struct Combo { const char* name; bool tri, edge, point; };
    static const Combo combos[] = {
        {"tri", true, false, false},
        {"edges", false, true, false},
        {"points", false, false, true},
        {"tri+edges", true, true, false},
        {"tri+points", true, false, true},
        {"edges+points", false, true, true},
        {"tri+edges+points", true, true, true},
    };
    for (const auto& c : combos) {
        DYNAMIC_SECTION(c.name) {
            if (c.tri && !fs::exists(surface)) {
                WARN("cube_subdiv.off not present, skipping: " << c.name);
            } else {
                const std::string out = run_combo(
                    c.tri ? surface.string() : std::string(),
                    c.edge ? edges.string() : std::string(),
                    c.point ? points.string() : std::string());
                INFO(out);
                REQUIRE(out.find("TRACKING OK") != std::string::npos);
            }
        }
    }
    fs::remove_all(dir, ec);
}

// Insert ALL edges of a real triangle mesh (112856simplified) as an edges-only input (no
// surface, no points). Unlike the random scattered edges above, these form a connected
// edge network sharing many vertices -- a good stress of the dedup of coincident forcing-
// triangle apexes -- and each input edge must still be tiled exactly by output tet edges.
TEST_CASE("integration: -t all edges of a real mesh inserted as edges", "[integration][tracking][large]") {
    const fs::path model = fs::path(VRTEST_MODELS_DIR) / "112856simplified.off";
    if (!fs::exists(model)) {
        WARN("112856simplified.off not present, skipping");
    } else {
        const fs::path dir = fs::temp_directory_path() / "vredges_real";
        std::error_code ec;
        fs::create_directories(dir, ec);
        const fs::path edges = dir / "all_edges.off";
        REQUIRE(extract_edges_off(model, edges));
        const std::string out = run_combo(std::string(), edges.string(), std::string());
        INFO(out);
        REQUIRE(out.find("TRACKING OK") != std::string::npos);
        fs::remove_all(dir, ec);
    }
}
#endif // NDEBUG

TEST_CASE("integration: model outputs are byte-stable across platforms", "[integration]") {
    const auto entries = load_manifest();
    REQUIRE_FALSE(entries.empty()); // reference_hashes.txt must not be empty

    for (const auto& e : entries) {
        DYNAMIC_SECTION(e.file << " [" << e.mode << "]") {
            const fs::path model = fs::path(VRTEST_MODELS_DIR) / e.file;
            if (!fs::exists(model)) {
                WARN("model not present, skipping: " << e.file
                     << " (commit it under models/ to exercise it in CI)");
            } else {
                const std::string got = run_and_hash(model.string(), e.mode);
                INFO("model:    " << e.file);
                INFO("mode:     " << e.mode);
                INFO("expected: " << e.sha256);
                INFO("actual:   " << (got.empty() ? std::string("<no output produced>") : got));
                REQUIRE(got == e.sha256);
            }
        }
    }
}
