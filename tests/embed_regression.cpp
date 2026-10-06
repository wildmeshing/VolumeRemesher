// Regression: embed_tri_in_poly_mesh emits inverted tets on a near-planar input.
//
// The input is a replay capture -- the exact four arrays handed to
// embed_tri_in_poly_mesh by wildmeshing-toolkit's tetwild for Thingi10K model 100727,
// dumped verbatim so the call can be reproduced without tetwild in the loop.
//
// What went wrong (with the tetrahedralization of the time, before the Delaunay3D
// backend): six of the returned tets had NEGATIVE signed volume in exact rational
// arithmetic. Two of them shared a face, and because both were stored with the wrong
// winding they presented that face with the SAME orientation as their neighbour rather
// than the opposite one -- so downstream the mesh looked like two tets overlapping on a
// shared face, which is what tetwild's own consistency check reported:
//
//     Face [3863, 3880, 1267799] appears more than once in the tet list
//
// The vertices involved were not coincident (nothing within 1e-3) and the tets were
// genuinely distinct and correctly placed on opposite sides of the face. Only the stored
// winding was wrong.
//
// Cause (MarcoAttene/Indirect_Predicates#15): nothing is decided in floating point -- the
// exact tier is simply never reached, because the filter wrongly reports that it does not
// need to be. makeTetrahedra fixes each tet's winding with genericPoint::orient3D, which
// returns sgn(det) of a determinant scaled by the operands' homogeneous denominators and
// is therefore correct only if every denominator is positive. implicitPoint3D_TBC -- the
// barycenter apex of a cell that cannot be tetrahedralized from one of its vertices --
// neither normalized its denominator nor reported when the interval filter could not
// decide its sign, so orient3D answered (true orientation) * sgn(d) and the flip inverted
// the tet. A TBC's denominator is the product of its generators', so it inherits an
// undecided sign from any generator whose own filter gave up: here a TPI, in a near-planar
// region (the z coordinates around the failure are 0, 0, -6.5e-16 and -6.1e-4). All six
// bad tets have a TBC apex, and |vol| between 6.7e-07 and 1.2e-05.
//
// NOT RUN BY DEFAULT: tagged [.] and run with
//     ./embed_regression "[embed_regression]"
//
// What it checks: the exact orientation of EVERY output tet (about 12M), plus the
// combinatorial consistency of the whole complex. For a while it checked only the six tets
// above, looked up by vertex index, because the all-tets check took over an hour; the
// Delaunay3D backend then renumbered the output, and the lookup found none of them. Checking
// every tet does not depend on numbering, and it still catches this bug: with the
// Indirect_Predicates#15 fix reverted, the current backend emits one inverted tet here
// (6*vol ~ -1.1e-4, measured 2026-10-06), and this test fails on it -- and on the four
// faces it then presents with the same winding as its neighbours.
//
// About 4 minutes on an M3 Max with 16 threads (measured 2026-10-06), nearly all of it
// inside embed_tri_in_poly_mesh (the arrangement and the tetrahedralization); the
// orientation check is about 4 CPU-minutes of that, spread over all threads. It was ~37
// minutes while NFG's GCD was Euclid's and the exact output coordinates were computed
// serially.

#include <catch2/catch_test_macros.hpp>

#include "VolumeRemesher/embed.h"
#include "VolumeRemesher/parallel.h"
#include "tet_orientation.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct EmbedInput
{
    std::vector<double> tri_vrt_coords;
    std::vector<uint32_t> triangle_indexes;
    std::vector<double> tet_vrt_coords;
    std::vector<uint32_t> tet_indexes;
};

// Layout: 8-byte magic "VREMBD01", then four length-prefixed arrays in the order the
// embed call takes them. Lengths are uint64, payloads are native-endian double/uint32.
EmbedInput read_embed_input(const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    REQUIRE(is.good());

    char magic[8] = {};
    is.read(magic, 8);
    REQUIRE(std::memcmp(magic, "VREMBD01", 8) == 0);

    const auto read_u64 = [&is]() {
        uint64_t n = 0;
        is.read(reinterpret_cast<char*>(&n), sizeof(n));
        return n;
    };
    const auto read_d = [&is, &read_u64]() {
        std::vector<double> v(read_u64());
        is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(double));
        return v;
    };
    const auto read_u = [&is, &read_u64]() {
        std::vector<uint32_t> v(read_u64());
        is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(uint32_t));
        return v;
    };

    EmbedInput in;
    in.tri_vrt_coords = read_d();
    in.triangle_indexes = read_u();
    in.tet_vrt_coords = read_d();
    in.tet_indexes = read_u();
    REQUIRE(is.good());
    return in;
}

// Signed volume of (a,b,c,d) as ((b-a) x (c-a)) . (d-a), in exact rational arithmetic.
// Positive means positively oriented under the convention embed_tri_in_poly_mesh is
// documented to emit.
vol_rem::bigrational signed_volume_x6(
    const std::vector<vol_rem::bigrational>& c,
    const std::array<uint32_t, 4>& t)
{
    const auto co = [&c](uint32_t v, int k) -> const vol_rem::bigrational& {
        return c[3 * static_cast<size_t>(v) + k];
    };
    vol_rem::bigrational e1[3], e2[3], e3[3];
    for (int k = 0; k < 3; ++k) {
        e1[k] = co(t[1], k) - co(t[0], k);
        e2[k] = co(t[2], k) - co(t[0], k);
        e3[k] = co(t[3], k) - co(t[0], k);
    }
    // (e1 x e2) . e3
    return (e1[1] * e2[2] - e1[2] * e2[1]) * e3[0] + (e1[2] * e2[0] - e1[0] * e2[2]) * e3[1] +
           (e1[0] * e2[1] - e1[1] * e2[0]) * e3[2];
}

#ifndef USE_GNU_GMP_CLASSES
// The same sign, without a single GCD. signed_volume_x6 subtracts rationals, and NFG's
// bigrational sum runs a GCD on every call -- about 300 us per tet here, an hour of CPU over
// all of them. Instead each vertex is scaled to integer homogeneous coordinates
// (X, Y, Z, W) = W * (x, y, z, 1) with W = Dx * Dy * Dz > 0, and the orientation is the sign
// of the 4x4 determinant of those rows: det4 = -(Wa * Wb * Wc * Wd) * vol6, so
// sgn(vol6) = -sgn(det4). The integers are bigrationals with denominator 1, built with the
// (non-reducing) (num, den, sign) constructor, so every sum is a GCD of 1 and 1.
using Homogeneous = std::array<vol_rem::bigrational, 4>;

Homogeneous homogeneous(const std::vector<vol_rem::bigrational>& c, uint32_t v)
{
    const vol_rem::bigrational* x = &c[3 * static_cast<size_t>(v)];
    const vol_rem::bignatural one(uint32_t(1));
    // A zero coordinate has an empty denominator; it scales nothing.
    const auto den = [&](int k) -> const vol_rem::bignatural& {
        return sgn(x[k]) == 0 ? one : x[k].get_den();
    };
    const vol_rem::bignatural w = den(0) * den(1) * den(2);
    Homogeneous h;
    for (int k = 0; k < 3; k++) {
        if (sgn(x[k]) == 0) continue; // h[k] stays zero
        h[k] = vol_rem::bigrational(
            x[k].get_num() * den((k + 1) % 3) * den((k + 2) % 3),
            one,
            sgn(x[k]));
    }
    h[3] = vol_rem::bigrational(w, one, 1);
    return h;
}

// sgn(vol6) of the tet with these four homogeneous vertices: the determinant by its 2x2
// minors (rows 0-1 against rows 2-3).
int orientation_sign(
    const Homogeneous& p0,
    const Homogeneous& p1,
    const Homogeneous& p2,
    const Homogeneous& p3)
{
    const auto m = [](const Homogeneous& a, const Homogeneous& b, int i, int j) {
        return a[i] * b[j] - a[j] * b[i];
    };
    const vol_rem::bigrational det4 =
        m(p0, p1, 0, 1) * m(p2, p3, 2, 3) - m(p0, p1, 0, 2) * m(p2, p3, 1, 3) +
        m(p0, p1, 0, 3) * m(p2, p3, 1, 2) + m(p0, p1, 1, 2) * m(p2, p3, 0, 3) -
        m(p0, p1, 1, 3) * m(p2, p3, 0, 2) + m(p0, p1, 2, 3) * m(p2, p3, 0, 1);
    return -sgn(det4);
}
#endif

} // namespace

TEST_CASE("embed_tri_in_poly_mesh emits positively oriented tets", "[embed_regression][.]")
{
    const EmbedInput in = read_embed_input(std::string(VRTEST_DATA_DIR) + "/embed_thingi100727.bin");

    INFO(
        "input: " << in.tri_vrt_coords.size() / 3 << " surface vertices, "
                  << in.triangle_indexes.size() / 3 << " triangles, " << in.tet_vrt_coords.size() / 3
                  << " background vertices, " << in.tet_indexes.size() / 4 << " background tets");

    std::vector<vol_rem::bigrational> out_vrt_coords;
    std::vector<uint32_t> out_poly_vindexes, out_cell_findexes, final_tets_parent, facets_on_input;
    std::vector<std::array<uint32_t, 4>> out_tets;
    std::vector<bool> cells_with_faces_on_input;
    std::vector<std::vector<uint32_t>> final_tets_parent_faces;

    // tetwild embeds triangles only: no extra edges or points, and the provenance
    // outputs go unused.
    const std::vector<double> no_edge_coords, no_point_coords;
    const std::vector<uint32_t> no_edge_indexes;
    std::vector<std::vector<std::array<uint32_t, 4>>> tri_provenance;
    std::vector<uint32_t> tri_group;
    std::vector<std::vector<std::array<uint32_t, 3>>> edge_provenance;
    std::vector<std::array<uint32_t, 2>> point_provenance;

    vol_rem::embed_tri_in_poly_mesh(
        in.tri_vrt_coords,
        in.triangle_indexes,
        in.tet_vrt_coords,
        in.tet_indexes,
        out_vrt_coords,
        out_poly_vindexes,
        out_cell_findexes,
        out_tets,
        final_tets_parent,
        facets_on_input,
        cells_with_faces_on_input,
        final_tets_parent_faces,
        no_edge_coords,
        no_edge_indexes,
        no_point_coords,
        tri_provenance,
        tri_group,
        edge_provenance,
        point_provenance,
        /*verbose=*/false);

    REQUIRE(!out_tets.empty());

    // The coplanar-group map is the key for tri_provenance: one entry per input triangle,
    // each either a valid index into it or UINT32_MAX for a triangle dropped as degenerate.
    REQUIRE(tri_group.size() == in.triangle_indexes.size() / 3);
    for (const uint32_t g : tri_group) {
        REQUIRE((g == UINT32_MAX || g < tri_provenance.size()));
    }

    // 1. Geometric: every returned tet has strictly positive exact signed volume.
    //
    // All of them, not a hand-picked list (the header says why). The volumes are computed
    // on parallel_blocks: each task reads the returned coordinates and makes its own
    // rationals, and only tet indices leave it (NFG numbers cannot cross threads -- see
    // parallel.h).
    //
    // Each sign comes from orientation_sign (no GCDs, see above). On every 997th tet it is
    // also computed the direct way, with signed_volume_x6, and the two must agree -- so a
    // slip in the determinant cannot make this check pass vacuously.
    std::vector<std::pair<size_t, int>> bad; // (tet, sign of its volume), sorted below
    std::atomic<uint64_t> cross_checked{0}, disagreements{0};
    std::mutex bad_mutex;
    vol_rem::parallel_blocks(out_tets.size(), [&](uint64_t lo, uint64_t hi) {
#ifndef USE_GNU_GMP_CLASSES
        // Consecutive tets are one cell's fan and share vertices.
        std::unordered_map<uint32_t, Homogeneous> cache;
        const auto h = [&](uint32_t v) -> const Homogeneous& {
            auto it = cache.find(v);
            if (it == cache.end()) it = cache.emplace(v, homogeneous(out_vrt_coords, v)).first;
            return it->second;
        };
#endif
        for (uint64_t i = lo; i < hi; ++i) {
            const auto& t = out_tets[i];
#ifndef USE_GNU_GMP_CLASSES
            const int s = orientation_sign(h(t[0]), h(t[1]), h(t[2]), h(t[3]));
            if (i % 997 == 0) {
                cross_checked++;
                if (sgn(signed_volume_x6(out_vrt_coords, t)) != s) disagreements++;
            }
#else
            // GMP's rationals are fast enough to use directly.
            const int s = sgn(signed_volume_x6(out_vrt_coords, t));
#endif
            if (s > 0) continue;
            std::lock_guard<std::mutex> lock(bad_mutex);
            bad.emplace_back(i, s);
        }
    });
    std::sort(bad.begin(), bad.end());
    INFO("orientation_sign cross-checked on " << cross_checked << " tets");
    CHECK(disagreements == 0);

    // A scoped INFO, not UNSCOPED_INFO: the latter would be consumed by the first CHECK
    // below even when that one passes.
    size_t inverted = 0, degenerate = 0;
    std::ostringstream listing;
    for (const auto& [i, s] : bad) {
        (s < 0 ? inverted : degenerate)++;
        if (inverted + degenerate > 10) continue;
        const auto& t = out_tets[i];
        listing << (s < 0 ? "inverted" : "degenerate") << " tet #" << i << " = [" << t[0] << ", "
                << t[1] << ", " << t[2] << ", " << t[3] << "], 6*vol ~ "
                << signed_volume_x6(out_vrt_coords, t).get_d() << "\n";
    }
    INFO(listing.str());
    CHECK(degenerate == 0);
    CHECK(inverted == 0);

    // 2. Combinatorial: a valid complex gives every internal face opposite windings from
    // its two tets. same_winding > 0 means two tets sit on the same side of a face, which
    // is the downstream symptom of (1).
    const vrtest::OrientationResult r = vrtest::check_tet_orientation(out_tets);
    UNSCOPED_INFO(
        "tets " << r.tets << ", boundary " << r.boundary << ", internal " << r.internal
                << ", same_winding " << r.same_winding << ", over_shared " << r.over_shared);
    CHECK(r.same_winding == 0);
    CHECK(r.over_shared == 0);
    CHECK(r.ok);
}
