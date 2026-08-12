// Unit tests for the 3D embedding pipeline's public API.
//
// Like unit_tests_2d.cpp, this target links the library and calls the API directly rather than
// shelling out to mesh_generator. The inputs are deliberately tiny (a cube of six tets, a
// handful of triangles) so the whole file runs in well under a second and can stay in the
// default ctest set -- embed_regression.cpp covers a real model but is hidden behind [.]
// because it takes over an hour.

#include <catch2/catch_test_macros.hpp>

#include <VolumeRemesher/embed.h>

#include <array>
#include <cstdint>
#include <vector>

namespace {

// Kuhn's triangulation of the axis-aligned box [lo,hi]^3 into six tets, all positively
// oriented under ((v1-v0)x(v2-v0)).(v3-v0) > 0. Cube vertex i has bit k set iff its k-th
// coordinate is hi.
struct Box
{
    std::vector<double> coords;
    std::vector<uint32_t> tets;
};

Box unit_box(double lo, double hi)
{
    Box b;
    for (int i = 0; i < 8; i++) {
        b.coords.push_back((i & 1) ? hi : lo);
        b.coords.push_back((i & 2) ? hi : lo);
        b.coords.push_back((i & 4) ? hi : lo);
    }
    // One tet per permutation of the axes: walk 0 -> e_i -> e_i+e_j -> 7.
    const int axis[3] = {1, 2, 4};
    const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    // The corner coordinates are exact in binary, so this determinant is exact in double and
    // is only ever compared against zero -- no tolerance needed.
    auto positive = [&](uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3) {
        auto d = [&](uint32_t a, uint32_t o, int k) {
            return b.coords[3 * a + k] - b.coords[3 * o + k];
        };
        const double a0 = d(i1, i0, 0), a1 = d(i1, i0, 1), a2 = d(i1, i0, 2);
        const double b0 = d(i2, i0, 0), b1 = d(i2, i0, 1), b2 = d(i2, i0, 2);
        const double c0 = d(i3, i0, 0), c1 = d(i3, i0, 1), c2 = d(i3, i0, 2);
        return (a1 * b2 - a2 * b1) * c0 + (a2 * b0 - a0 * b2) * c1 + (a0 * b1 - a1 * b0) * c2 > 0.0;
    };
    for (const auto& p : perm) {
        const uint32_t a = uint32_t(axis[p[0]]);
        const uint32_t c = uint32_t(axis[p[0]] | axis[p[1]]);
        if (positive(0u, a, c, 7u)) {
            b.tets.insert(b.tets.end(), {0u, a, c, 7u});
        } else {
            b.tets.insert(b.tets.end(), {0u, c, a, 7u});
        }
    }
    return b;
}

struct EmbedOut
{
    std::vector<vol_rem::bigrational> vertices;
    std::vector<uint32_t> facets, cells, final_tets_parent, facets_on_input;
    std::vector<std::array<uint32_t, 4>> tets;
    std::vector<bool> cells_with_faces_on_input;
    std::vector<std::vector<uint32_t>> final_tets_parent_faces;
    std::vector<std::vector<std::array<uint32_t, 4>>> tri_provenance;
    std::vector<uint32_t> tri_group;
    std::vector<std::vector<std::array<uint32_t, 3>>> edge_provenance;
    std::vector<std::array<uint32_t, 2>> point_provenance;
};

EmbedOut embed(const std::vector<double>& tri_coords, const std::vector<uint32_t>& tri_idx)
{
    const Box box = unit_box(-1.0, 2.0);
    EmbedOut o;
    const std::vector<double> no_edge_coords, no_point_coords;
    const std::vector<uint32_t> no_edge_indexes;
    vol_rem::embed_tri_in_poly_mesh(
        tri_coords,
        tri_idx,
        box.coords,
        box.tets,
        o.vertices,
        o.facets,
        o.cells,
        o.tets,
        o.final_tets_parent,
        o.facets_on_input,
        o.cells_with_faces_on_input,
        o.final_tets_parent_faces,
        no_edge_coords,
        no_edge_indexes,
        no_point_coords,
        o.tri_provenance,
        o.tri_group,
        o.edge_provenance,
        o.point_provenance,
        /*verbose=*/false);
    return o;
}

} // namespace

// out_triangle_group is what makes out_triangle_provenance usable: without it a caller holds
// per-group face lists and no way back to the input triangle, hence no way to carry a
// per-triangle attribute (a material tag, which input file a surface came from) through the
// arrangement.
TEST_CASE("embed_tri_in_poly_mesh reports the coplanar group of every input triangle", "[embed3d]")
{
    // Two coplanar triangles tiling the square z=0.5 (sharing the diagonal s0-s2), one
    // triangle folded out of that plane across the edge s0-s1, and one degenerate
    // (collinear) triangle that the reader drops before the arrangement.
    // clang-format off
    const std::vector<double> tri_coords = {
        0.0,  0.0, 0.5,  // s0
        1.0,  0.0, 0.5,  // s1
        1.0,  1.0, 0.5,  // s2
        0.0,  1.0, 0.5,  // s3
        0.5, -0.5, 1.0,  // s4, the fold's apex
        0.2,  0.2, 0.2,  // s5 \.
        0.4,  0.4, 0.4,  // s6  | collinear
        0.6,  0.6, 0.6,  // s7 /
    };
    // The degenerate one sits in the MIDDLE on purpose: it is dropped as the constraints are
    // read, so every later triangle has a constraint index one below its input index. A map
    // that returned the internal numbering unchanged would pass with it last and fail here.
    const std::vector<uint32_t> tri_idx = {
        0, 1, 2,  // t0 \_ coplanar, share edge s0-s2
        5, 6, 7,  // t1, degenerate
        0, 2, 3,  // t2 /
        1, 0, 4,  // t3, shares edge s0-s1 with t0 but is not coplanar with it
    };
    // clang-format on

    const EmbedOut o = embed(tri_coords, tri_idx);

    REQUIRE(!o.tets.empty());
    // One entry per input triangle, in the caller's numbering.
    REQUIRE(o.tri_group.size() == tri_idx.size() / 3);

    // Every group named is a valid index into the provenance, and names a non-empty face list:
    // each of these triangles is inside the box, so the arrangement must tile it.
    for (const size_t t : {0u, 2u, 3u}) {
        INFO("input triangle " << t);
        REQUIRE(o.tri_group[t] < o.tri_provenance.size());
        REQUIRE(!o.tri_provenance[o.tri_group[t]].empty());
    }

    // Edge-adjacent AND exactly coplanar: one group.
    CHECK(o.tri_group[0] == o.tri_group[2]);
    // Edge-adjacent but not coplanar: a different group.
    CHECK(o.tri_group[3] != o.tri_group[0]);
    // Dropped as degenerate, so it has no group at all.
    CHECK(o.tri_group[1] == UINT32_MAX);

    // The faces attributed to a group are faces of the tets they name.
    for (const auto& faces : o.tri_provenance) {
        for (const auto& f : faces) {
            REQUIRE(f[0] < o.tets.size());
            const auto& tet = o.tets[f[0]];
            for (int k = 1; k < 4; k++) {
                const bool is_corner =
                    f[k] == tet[0] || f[k] == tet[1] || f[k] == tet[2] || f[k] == tet[3];
                REQUIRE(is_corner);
            }
        }
    }
}

// A triangle with no coplanar neighbour is its own singleton group, so with a fan of
// mutually non-coplanar triangles the map must be injective.
TEST_CASE("non-coplanar input triangles get distinct groups", "[embed3d]")
{
    // clang-format off
    const std::vector<double> tri_coords = {
        0.0, 0.0, 0.0,  // apex, shared by all three
        1.0, 0.0, 0.5,
        0.0, 1.0, 0.5,
        0.0, 0.0, 1.0,
    };
    const std::vector<uint32_t> tri_idx = {
        0, 1, 2,
        0, 2, 3,
        0, 3, 1,
    };
    // clang-format on

    const EmbedOut o = embed(tri_coords, tri_idx);

    REQUIRE(o.tri_group.size() == 3);
    CHECK(o.tri_group[0] != o.tri_group[1]);
    CHECK(o.tri_group[1] != o.tri_group[2]);
    CHECK(o.tri_group[0] != o.tri_group[2]);
}
