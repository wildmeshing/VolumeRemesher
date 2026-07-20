/*
 *  Copyright (c) 2000-2022 Inria
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contact: Bruno Levy
 *
 *     https://www.inria.fr/fr/bruno-levy
 *
 *     Inria,
 *     Domaine de Voluceau,
 *     78150 Le Chesnay - Rocquencourt
 *     FRANCE
 *
 */

// ---------------------------------------------------------------------------------------------
// ADAPTED FROM geogram/src/lib/geogram/delaunay/delaunay_2d.{h,cpp} (Delaunay2d).
//
// The incremental Bowyer-Watson algorithm, the triangle store, the vertex-at-infinity handling,
// the conflict-zone search and the stellation are kept as close to the original as practical, so
// that the two can still be diffed. Modifications, all deliberate:
//
//   1. PREDICATES. PCK::orient_2d / in_circle_2d_SOS / points_are_identical_2d now resolve to
//      predicates2d.h, which is built on this repository's exact kernel (numerics.h) instead of
//      geogram's expansion_nt. in_circle_2d_SOS is a new implementation; see predicates2d.cpp.
//
//   2. DETERMINISM. Geogram calls Numeric::random_int32() in two places: to pick the starting
//      edge of the visibility walk in locate(), and to pick a starting triangle when no hint is
//      given. Both are replaced with deterministic choices. This repository guarantees
//      byte-identical output on Linux/macOS/Windows, and a platform RNG would break that.
//      For the same reason compute_BRIO_order (which uses random rounds) is replaced by a pure
//      deterministic Hilbert sort.
//
//   3. SCOPE. The weighted / regular-triangulation path (heights_, orient_2dlifted_SOS,
//      has_empty_cells_), nearest_vertex() and its kd-tree, the Delaunay base class, and the
//      GEO logging/benchmarking scaffolding are all dropped. GEO::vector -> std::vector,
//      geo_debug_assert -> assert.
//
//   4. INFINITE TRIANGLES ARE KEPT in the final output (geogram's keep_infinite_ mode). The 2D
//      arrangement that consumes this triangulation walks the mesh, and keeping the infinite
//      triangles makes triangle_adjacent() total, so the walk needs no boundary special case.
// ---------------------------------------------------------------------------------------------

#ifndef VOLUMEREMESHER_2D_DELAUNAY2D_H
#define VOLUMEREMESHER_2D_DELAUNAY2D_H

#include <cassert>
#include <cstdint>
#include <vector>

namespace vol_rem {
namespace vr2d {

typedef uint32_t index_t;

class Delaunay2d
{
public:
    static constexpr index_t NO_INDEX = ~index_t(0);
    static constexpr index_t NO_TRIANGLE = NO_INDEX;

    // A triangle incident to this symbolic vertex is "virtual": it covers the region outside the
    // convex hull of the inserted points. Its edge opposite the infinite vertex is a hull edge.
    static constexpr index_t VERTEX_AT_INFINITY = NO_INDEX;

    Delaunay2d() = default;

    // points holds 2*nb_vertices doubles, (x,y) interleaved, and must be free of exact
    // duplicates (the caller deduplicates; see arrangement2d). The pointer must outlive this
    // object: coordinates are not copied.
    //
    // Returns false (leaving the triangulation empty) if all the points are collinear.
    bool set_vertices(index_t nb_vertices, const double* points);

    index_t nb_vertices() const { return nb_vertices_; }
    const double* vertex_ptr(index_t v) const { return points_ + 2 * v; }

    // Total triangle count, INCLUDING the virtual ones. After set_vertices the triangles are
    // ordered so that the finite ones occupy [0, nb_finite_triangles).
    index_t nb_triangles() const { return index_t(cell_to_v_store_.size() / 3); }
    index_t nb_finite_triangles() const { return nb_finite_triangles_; }

    index_t triangle_vertex(index_t t, index_t lv) const
    {
        assert(t < nb_triangles() && lv < 3);
        return cell_to_v_store_[3 * t + lv];
    }
    index_t triangle_adjacent(index_t t, index_t le) const
    {
        assert(t < nb_triangles() && le < 3);
        return cell_to_cell_store_[3 * t + le];
    }
    bool triangle_is_finite(index_t t) const
    {
        return cell_to_v_store_[3 * t] != NO_INDEX && cell_to_v_store_[3 * t + 1] != NO_INDEX &&
               cell_to_v_store_[3 * t + 2] != NO_INDEX;
    }

    // Hand the connectivity arrays to the caller (they are moved out; the object is left empty).
    // The 2D arrangement takes ownership and then modifies the triangulation in place.
    void steal_arrays(std::vector<index_t>& tri_node, std::vector<index_t>& tri_neigh)
    {
        tri_node = std::move(cell_to_v_store_);
        tri_neigh = std::move(cell_to_cell_store_);
        cell_to_v_store_.clear();
        cell_to_cell_store_.clear();
    }

    // Deterministic replacement for geogram's compute_BRIO_order: sorts point indices along a
    // Hilbert curve over a 21-bit-per-axis quantization of the bounding box. Exposed because the
    // unit tests check that it is a permutation and that it is platform-stable.
    static void hilbert_order(index_t nb_vertices, const double* points,
                              std::vector<index_t>& order);

    // For debugging / tests: checks combinatorial consistency of the triangulation.
    bool check_combinatorics() const;

protected:
    // --- triangle list / free list bookkeeping (geogram, unchanged) ---------------------------
    static constexpr index_t NOT_IN_LIST = ~index_t(0);
    static constexpr index_t NOT_IN_LIST_BIT = index_t(1) << (sizeof(index_t) * 8 - 1);
    static constexpr index_t END_OF_LIST = ~NOT_IN_LIST_BIT;

    index_t max_t() const { return index_t(cell_to_v_store_.size() / 3); }

    bool triangle_is_in_list(index_t t) const { return (cell_next_[t] & NOT_IN_LIST_BIT) == 0; }
    index_t triangle_next(index_t t) const { return cell_next_[t]; }
    void add_triangle_to_list(index_t t, index_t& first, index_t& last)
    {
        if (last == END_OF_LIST) {
            first = last = t;
            cell_next_[t] = END_OF_LIST;
        } else {
            cell_next_[t] = first;
            first = t;
        }
    }
    void remove_triangle_from_list(index_t t) { cell_next_[t] = NOT_IN_LIST; }
    bool triangle_is_free(index_t t) const { return triangle_is_in_list(t); }
    bool triangle_is_virtual(index_t t) const
    {
        return !triangle_is_free(t) && (cell_to_v_store_[3 * t] == VERTEX_AT_INFINITY ||
                                        cell_to_v_store_[3 * t + 1] == VERTEX_AT_INFINITY ||
                                        cell_to_v_store_[3 * t + 2] == VERTEX_AT_INFINITY);
    }
    bool triangle_is_real(index_t t) const { return !triangle_is_free(t) && triangle_is_finite(t); }

    void set_triangle_mark_stamp(index_t stamp) { cur_stamp_ = (stamp | NOT_IN_LIST_BIT); }
    bool triangle_is_marked(index_t t) const { return cell_next_[t] == cur_stamp_; }
    void mark_triangle(index_t t) { cell_next_[t] = cur_stamp_; }

    index_t new_triangle();
    index_t new_triangle(index_t v1, index_t v2, index_t v3);

    // triangle_edge_vertex_[le][lv]: the triangle formed by vertex lv, edge_vertex(lv,0) and
    // edge_vertex(lv,1) has the same orientation as the original triangle, for any lv.
    static index_t triangle_edge_vertex(index_t e, index_t v)
    {
        static const char tev[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        return index_t(tev[e][v]);
    }

    index_t finite_triangle_vertex(index_t t, index_t lv) const
    {
        assert(cell_to_v_store_[3 * t + lv] != NO_INDEX);
        return cell_to_v_store_[3 * t + lv];
    }
    void set_triangle_adjacent(index_t t1, index_t le1, index_t t2)
    {
        assert(t1 != t2);
        cell_to_cell_store_[3 * t1 + le1] = t2;
    }
    static index_t find_3(const index_t* T, index_t v)
    {
        const index_t result = index_t((T[1] == v) | ((T[2] == v) * 2));
        assert(T[result] == v);
        return result;
    }
    index_t find_triangle_vertex(index_t t, index_t v) const
    {
        return find_3(&cell_to_v_store_[3 * t], v);
    }
    index_t find_triangle_adjacent(index_t t1, index_t t2) const
    {
        return find_3(&cell_to_cell_store_[3 * t1], t2);
    }

    // --- algorithm ---------------------------------------------------------------------------
    bool create_first_triangle(index_t& iv0, index_t& iv1, index_t& iv2);
    index_t locate(const double* p, index_t hint, int* orient) const;
    index_t locate_inexact(const double* p, index_t hint, index_t max_iter) const;
    index_t insert(index_t v, index_t hint);
    void find_conflict_zone(index_t v, index_t t, const int* orient, index_t& t_bndry,
                            index_t& e_bndry, index_t& first, index_t& last);
    void find_conflict_zone_iterative(const double* p, index_t t, index_t& t_bndry,
                                      index_t& e_bndry, index_t& first, index_t& last);
    index_t stellate_conflict_zone(index_t v, index_t t_bndry, index_t e_bndry);
    bool triangle_is_conflict(index_t t, const double* p) const;

private:
    const double* points_ = nullptr;
    index_t nb_vertices_ = 0;
    index_t nb_finite_triangles_ = 0;

    std::vector<index_t> cell_to_v_store_;
    std::vector<index_t> cell_to_cell_store_;
    std::vector<index_t> cell_next_;
    std::vector<index_t> reorder_;
    std::vector<index_t> S_; // find_conflict_zone_iterative work stack
    index_t cur_stamp_ = 0;
    index_t first_free_ = END_OF_LIST;
};

} // namespace vr2d
} // namespace vol_rem

#endif
