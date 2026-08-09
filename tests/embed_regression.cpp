// Regression: embed_tri_in_poly_mesh emits inverted tets on a near-planar input.
//
// The input is a replay capture -- the exact four arrays handed to
// embed_tri_in_poly_mesh by wildmeshing-toolkit's tetwild for Thingi10K model 100727,
// dumped verbatim so the call can be reproduced without tetwild in the loop.
//
// What goes wrong: six of the returned tets have NEGATIVE signed volume in exact
// rational arithmetic. Two of them share a face, and because both are stored with the
// wrong winding they present that face with the SAME orientation as their neighbour
// rather than the opposite one -- so downstream the mesh looks like two tets overlapping
// on a shared face, which is what tetwild's own consistency check reports:
//
//     Face [3863, 3880, 1267799] appears more than once in the tet list
//
// The vertices involved are not coincident (nothing within 1e-3) and the tets are
// genuinely distinct and correctly placed on opposite sides of the face. Only the stored
// winding is wrong. The region is near-planar -- the z coordinates around the failure are
// 0, 0, -6.5e-16 and -6.1e-4 -- and every bad tet has a tiny volume, |vol| between 6.7e-07
// and 1.2e-05, which is where a sign decided in floating point would go astray.
//
// NOT RUN BY DEFAULT: tagged [.] because the embedding takes ~80 minutes. Run it with
//     ./embed_regression "[embed_regression]"

#include <catch2/catch_test_macros.hpp>

#include "VolumeRemesher/embed.h"
#include "tet_orientation.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
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
    const auto co = [&c](uint32_t v, int k) { return c[3 * static_cast<size_t>(v) + k]; };
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
        edge_provenance,
        point_provenance,
        /*verbose=*/false);

    REQUIRE(!out_tets.empty());

    // 1. Geometric: every tet must have strictly positive exact signed volume.
    std::vector<size_t> inverted, degenerate;
    for (size_t i = 0; i < out_tets.size(); ++i) {
        const vol_rem::bigrational v = signed_volume_x6(out_vrt_coords, out_tets[i]);
        if (v.sgn() < 0) {
            inverted.push_back(i);
        } else if (v.sgn() == 0) {
            degenerate.push_back(i);
        }
    }
    for (size_t i = 0; i < inverted.size() && i < 10; ++i) {
        const auto& t = out_tets[inverted[i]];
        UNSCOPED_INFO(
            "inverted tet #" << inverted[i] << " = [" << t[0] << ", " << t[1] << ", " << t[2]
                             << ", " << t[3] << "]");
    }
    CHECK(degenerate.empty());
    CHECK(inverted.empty());

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
