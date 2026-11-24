// Fixtures shared by feature_vertex_regression.cpp and feature_vertex_errors.cpp: small
// background tet meshes, a cube surface, and a wrapper around embed_tri_in_poly_mesh that keeps
// the outputs the tests look at.
#pragma once

#include "VolumeRemesher/embed.h"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace feature_vertex_fixtures {

struct Background
{
    std::vector<double> coords;
    std::vector<uint32_t> tets;
};

// Appends the 6 Kuhn tets of the cell with corner ids vid(i..i+1, j..j+1, k..k+1), every tet
// oriented so that (v1-v0)x(v2-v0).(v3-v0) > 0.
template <class Vid>
inline void kuhn_cell(Background& b, const Vid& vid, int i, int j, int k)
{
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& p : perms) {
        int c[3] = {i, j, k};
        std::array<uint32_t, 4> t;
        t[0] = vid(c[0], c[1], c[2]);
        for (int s = 0; s < 3; s++) {
            c[p[s]]++;
            t[s + 1] = vid(c[0], c[1], c[2]);
        }
        const double* a = &b.coords[3 * t[0]];
        double e[3][3];
        for (int s = 0; s < 3; s++)
            for (int d = 0; d < 3; d++) e[s][d] = b.coords[3 * t[s + 1] + d] - a[d];
        const double vol = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                           e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
                           e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
        if (vol < 0) std::swap(t[2], t[3]);
        b.tets.insert(b.tets.end(), t.begin(), t.end());
    }
}

// A grid with lines at -0.2, 0, 2, 2.2 on each axis (through the corners of the cube [0,2]^3,
// padded by 0.2, like wildmeshing-toolkit's tetwild background at a small scale), each cell split
// into 6 Kuhn tets.
inline Background grid_background()
{
    const double xs[4] = {-0.2, 0.0, 2.0, 2.2};
    Background b;
    const auto vid = [](int i, int j, int k) { return uint32_t(i + 4 * j + 16 * k); };
    for (int k = 0; k < 4; k++)
        for (int j = 0; j < 4; j++)
            for (int i = 0; i < 4; i++) {
                b.coords.push_back(xs[i]);
                b.coords.push_back(xs[j]);
                b.coords.push_back(xs[k]);
            }
    for (int k = 0; k < 3; k++)
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++) kuhn_cell(b, vid, i, j, k);
    return b;
}

// The box [lo, hi]^3 as 6 Kuhn tets.
inline Background box_background(double lo, double hi)
{
    Background b;
    const auto vid = [](int i, int j, int k) { return uint32_t(i + 2 * j + 4 * k); };
    for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
                b.coords.push_back(i ? hi : lo);
                b.coords.push_back(j ? hi : lo);
                b.coords.push_back(k ? hi : lo);
            }
    kuhn_cell(b, vid, 0, 0, 0);
    return b;
}

// The cube [lo, hi]^3, outward triangles.
inline std::vector<double> cube_coords(double lo = 0.0, double hi = 2.0)
{
    std::vector<double> c = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    for (double& x : c) x = (x == 0.0) ? lo : hi;
    return c;
}
inline const std::vector<uint32_t> cube_tris = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
                                                2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3};

struct Result
{
    std::vector<vol_rem::bigrational> coords;
    std::vector<std::array<uint32_t, 4>> tets;
    std::vector<std::vector<std::array<uint32_t, 3>>> edge_provenance;
    std::vector<std::array<uint32_t, 2>> point_provenance;
};

inline Result embed(
    const Background& bg,
    const std::vector<double>& surface_coords,
    const std::vector<double>& edge_coords,
    const std::vector<uint32_t>& edge_indexes,
    const std::vector<double>& point_coords,
    const std::vector<uint32_t>& surface_tris = cube_tris)
{
    Result r;
    std::vector<uint32_t> poly_vindexes, cell_findexes, tets_parent, facets_on_input;
    std::vector<bool> cells_with_faces_on_input;
    std::vector<std::vector<uint32_t>> tets_parent_faces;
    std::vector<std::vector<std::array<uint32_t, 4>>> tri_provenance;
    std::vector<uint32_t> tri_group;
    vol_rem::embed_tri_in_poly_mesh(
        surface_coords,
        surface_tris,
        bg.coords,
        bg.tets,
        r.coords,
        poly_vindexes,
        cell_findexes,
        r.tets,
        tets_parent,
        facets_on_input,
        cells_with_faces_on_input,
        tets_parent_faces,
        edge_coords,
        edge_indexes,
        point_coords,
        tri_provenance,
        tri_group,
        r.edge_provenance,
        r.point_provenance,
        /*verbose=*/false);
    return r;
}

} // namespace feature_vertex_fixtures
