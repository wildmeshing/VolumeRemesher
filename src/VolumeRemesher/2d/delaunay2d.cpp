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

// Adapted from geogram/src/lib/geogram/delaunay/delaunay_2d.cpp. See delaunay2d.h for the list
// of deliberate modifications.

#include "delaunay2d.h"

#include "predicates2d.h"

#include <algorithm>
#include <cmath>

namespace vol_rem {
namespace vr2d {

namespace {

// Inexact orientation, used only by locate_inexact() for structural filtering. Its result is a
// hint; every decision that matters is re-taken with the exact PCK::orient_2d.
inline Sign orient_2d_inexact(const double* p0, const double* p1, const double* p2)
{
    const double a11 = p1[0] - p0[0];
    const double a12 = p1[1] - p0[1];
    const double a21 = p2[0] - p0[0];
    const double a22 = p2[1] - p0[1];
    const double Delta = a11 * a22 - a12 * a21;
    return geo_sgn(Delta);
}

// --- deterministic Hilbert ordering -----------------------------------------------------------
// Replaces geogram's compute_BRIO_order, whose random rounds depend on Numeric::random_int32().

constexpr uint32_t HILBERT_BITS = 21;              // 21 bits/axis -> index fits in 42 bits
constexpr uint32_t HILBERT_N = 1u << HILBERT_BITS; // side length of the quantization grid

inline void hilbert_rot(uint32_t n, uint32_t& x, uint32_t& y, uint32_t rx, uint32_t ry)
{
    if (ry == 0) {
        if (rx == 1) {
            x = n - 1 - x;
            y = n - 1 - y;
        }
        const uint32_t t = x;
        x = y;
        y = t;
    }
}

inline uint64_t hilbert_xy2d(uint32_t x, uint32_t y)
{
    uint64_t d = 0;
    for (uint32_t s = HILBERT_N / 2; s > 0; s /= 2) {
        const uint32_t rx = (x & s) ? 1u : 0u;
        const uint32_t ry = (y & s) ? 1u : 0u;
        d += uint64_t(s) * uint64_t(s) * uint64_t((3 * rx) ^ ry);
        hilbert_rot(HILBERT_N, x, y, rx, ry);
    }
    return d;
}

} // namespace

void Delaunay2d::hilbert_order(index_t nb_vertices, const double* points,
                               std::vector<index_t>& order)
{
    order.resize(nb_vertices);
    for (index_t i = 0; i < nb_vertices; i++) order[i] = i;
    if (nb_vertices < 2) return;

    double xmin = points[0], xmax = points[0], ymin = points[1], ymax = points[1];
    for (index_t i = 1; i < nb_vertices; i++) {
        xmin = std::min(xmin, points[2 * i]);
        xmax = std::max(xmax, points[2 * i]);
        ymin = std::min(ymin, points[2 * i + 1]);
        ymax = std::max(ymax, points[2 * i + 1]);
    }
    // Plain IEEE double arithmetic, so the quantization is bit-identical on every platform
    // (the build forces /fp:strict on MSVC and -ffp-contract=off elsewhere).
    const double dx = xmax - xmin;
    const double dy = ymax - ymin;
    const double sx = (dx > 0.0) ? (double(HILBERT_N - 1) / dx) : 0.0;
    const double sy = (dy > 0.0) ? (double(HILBERT_N - 1) / dy) : 0.0;

    std::vector<uint64_t> key(nb_vertices);
    for (index_t i = 0; i < nb_vertices; i++) {
        const double fx = (points[2 * i] - xmin) * sx;
        const double fy = (points[2 * i + 1] - ymin) * sy;
        const uint32_t qx = uint32_t(std::min(double(HILBERT_N - 1), std::max(0.0, fx)));
        const uint32_t qy = uint32_t(std::min(double(HILBERT_N - 1), std::max(0.0, fy)));
        key[i] = hilbert_xy2d(qx, qy);
    }

    // Tie-break on the index so the comparator is a STRICT TOTAL ORDER. Without it, points
    // sharing a Hilbert cell would be ordered by std::sort's unspecified tie handling, which
    // differs between libstdc++, libc++ and MSVC -- exactly the class of bug that made a previous
    // MSVC build diverge (see conforming_mesh.cpp's half-edge sort).
    std::sort(order.begin(), order.end(), [&](index_t a, index_t b) {
        return key[a] != key[b] ? key[a] < key[b] : a < b;
    });
}

// --- combinatorics ----------------------------------------------------------------------------

index_t Delaunay2d::new_triangle()
{
    index_t result;
    if (first_free_ == END_OF_LIST) {
        cell_to_v_store_.resize(cell_to_v_store_.size() + 3, NO_INDEX);
        cell_to_cell_store_.resize(cell_to_cell_store_.size() + 3, NO_INDEX);
        cell_next_.push_back(NOT_IN_LIST);
        result = max_t() - 1;
    } else {
        result = first_free_;
        first_free_ = triangle_next(first_free_);
        remove_triangle_from_list(result);
    }
    cell_to_cell_store_[3 * result] = NO_INDEX;
    cell_to_cell_store_[3 * result + 1] = NO_INDEX;
    cell_to_cell_store_[3 * result + 2] = NO_INDEX;
    return result;
}

index_t Delaunay2d::new_triangle(index_t v1, index_t v2, index_t v3)
{
    const index_t result = new_triangle();
    cell_to_v_store_[3 * result] = v1;
    cell_to_v_store_[3 * result + 1] = v2;
    cell_to_v_store_[3 * result + 2] = v3;
    return result;
}

bool Delaunay2d::triangle_is_conflict(index_t t, const double* p) const
{
    const double* pv[3];
    for (index_t i = 0; i < 3; ++i) {
        const index_t v = triangle_vertex(t, i);
        pv[i] = (v == NO_INDEX) ? nullptr : vertex_ptr(v);
    }

    // Virtual triangle: in_circle is replaced by orient_2d on the hull edge.
    for (index_t le = 0; le < 3; ++le) {
        if (pv[le] == nullptr) {
            pv[le] = p;
            const Sign sign = PCK::orient_2d(pv[0], pv[1], pv[2]);
            if (sign > 0) return true;
            if (sign < 0) return false;

            // p is exactly on the hull edge: defer to the real triangle behind it.
            const index_t t2 = triangle_adjacent(t, le);
            assert(t2 != NO_INDEX);
            if (triangle_is_in_list(t2)) return true;
            if (triangle_is_marked(t2)) return false;
            return triangle_is_conflict(t2, p);
        }
    }

    return PCK::in_circle_2d_SOS(pv[0], pv[1], pv[2], p) > 0;
}

// --- point location ---------------------------------------------------------------------------

index_t Delaunay2d::locate_inexact(const double* p, index_t hint, index_t max_iter) const
{
    // DEVIATION (determinism): geogram picks a random triangle here via Numeric::random_int32().
    // A deterministic scan for the first live triangle gives the same guarantee (a valid,
    // non-free, real starting point) without a platform-dependent RNG.
    if (hint == NO_TRIANGLE || hint >= max_t() || triangle_is_free(hint)) {
        hint = NO_TRIANGLE;
        for (index_t t = 0; t < max_t(); ++t) {
            if (!triangle_is_free(t)) {
                hint = t;
                break;
            }
        }
        if (hint == NO_TRIANGLE) return NO_TRIANGLE;
    }

    // Always start from a real triangle. If virtual, step to its real neighbour (always the one
    // opposite the infinite vertex).
    if (triangle_is_virtual(hint)) {
        for (index_t lf = 0; lf < 3; ++lf) {
            if (triangle_vertex(hint, lf) == VERTEX_AT_INFINITY) {
                hint = triangle_adjacent(hint, lf);
                assert(hint != NO_TRIANGLE);
                break;
            }
        }
    }

    index_t t = hint;
    index_t t_pred = NO_TRIANGLE;

still_walking : {
    const double* pv[3];
    pv[0] = vertex_ptr(finite_triangle_vertex(t, 0));
    pv[1] = vertex_ptr(finite_triangle_vertex(t, 1));
    pv[2] = vertex_ptr(finite_triangle_vertex(t, 2));

    for (index_t le = 0; le < 3; ++le) {
        const index_t t_next = triangle_adjacent(t, le);
        if (t_next == NO_INDEX) return NO_TRIANGLE;
        if (t_next == t_pred) continue;

        // Orientation of p w.r.t. edge le of t: replace vertex le with p (CGAL's convention).
        const double* pv_bkp = pv[le];
        pv[le] = p;
        const Sign ori = orient_2d_inexact(pv[0], pv[1], pv[2]);
        if (ori != NEGATIVE) {
            pv[le] = pv_bkp;
            continue;
        }
        if (triangle_is_virtual(t_next)) return t_next;

        t_pred = t;
        t = t_next;
        if (--max_iter != 0) goto still_walking;
    }
}
    return t;
}

index_t Delaunay2d::locate(const double* p, index_t hint, int* orient) const
{
    // Structural filtering: improve the hint with the inexact walk first. Bounded at 2500 steps
    // because there exist configurations in which the inexact walk cycles forever.
    hint = locate_inexact(p, hint, 2500);

    if (hint == NO_TRIANGLE || hint >= max_t() || triangle_is_free(hint)) {
        hint = NO_TRIANGLE;
        for (index_t t = 0; t < max_t(); ++t) {
            if (!triangle_is_free(t)) {
                hint = t;
                break;
            }
        }
        if (hint == NO_TRIANGLE) return NO_TRIANGLE;
    }

    if (triangle_is_virtual(hint)) {
        for (index_t le = 0; le < 3; ++le) {
            if (triangle_vertex(hint, le) == VERTEX_AT_INFINITY) {
                hint = triangle_adjacent(hint, le);
                assert(hint != NO_TRIANGLE);
                break;
            }
        }
    }
    assert(!triangle_is_free(hint) && !triangle_is_virtual(hint));

    index_t t = hint;
    index_t t_pred = NO_TRIANGLE;
    int orient_local[3];
    if (orient == nullptr) orient = orient_local;

still_walking : {
    const double* pv[3];
    pv[0] = vertex_ptr(finite_triangle_vertex(t, 0));
    pv[1] = vertex_ptr(finite_triangle_vertex(t, 1));
    pv[2] = vertex_ptr(finite_triangle_vertex(t, 2));

    // DEVIATION (determinism): geogram starts from a random edge, Numeric::random_int32() % 3.
    // The purpose is to avoid systematically favouring one edge, which can make the walk cycle on
    // adversarial inputs. Deriving the start from the current triangle index keeps that variation
    // while being a pure function of the mesh state, hence reproducible everywhere.
    const index_t e0 = t % 3;
    for (index_t de = 0; de < 3; ++de) {
        const index_t le = (e0 + de) % 3;
        const index_t t_next = triangle_adjacent(t, le);
        if (t_next == NO_INDEX) return NO_TRIANGLE;

        if (t_next == t_pred) {
            orient[le] = POSITIVE;
            continue;
        }

        const double* pv_bkp = pv[le];
        pv[le] = p;
        orient[le] = PCK::orient_2d(pv[0], pv[1], pv[2]);

        if (orient[le] != NEGATIVE) {
            pv[le] = pv_bkp;
            continue;
        }

        if (triangle_is_virtual(t_next)) {
            for (index_t tle = 0; tle < 3; ++tle) orient[tle] = POSITIVE;
            return t_next;
        }

        t_pred = t;
        t = t_next;
        goto still_walking;
    }
}
    return t;
}

// --- conflict zone ----------------------------------------------------------------------------

void Delaunay2d::find_conflict_zone(index_t v, index_t t, const int* orient, index_t& t_bndry,
                                    index_t& e_bndry, index_t& first, index_t& last)
{
    first = last = END_OF_LIST;
    set_triangle_mark_stamp(v);

    const double* p = vertex_ptr(v);
    assert(t != NO_TRIANGLE);

    // The point already exists in the triangulation if it lies on two edges of t (i.e. it is a
    // vertex of t).
    const int nb_zero = (orient[0] == ZERO) + (orient[1] == ZERO) + (orient[2] == ZERO);
    if (nb_zero >= 2) return;

    add_triangle_to_list(t, first, last);

    // If p lies exactly on an edge of t, the neighbour across it is necessarily in conflict; add
    // it directly rather than paying for the predicate.
    if (nb_zero != 0) {
        for (index_t le = 0; le < 3; ++le) {
            if (orient[le] == ZERO) add_triangle_to_list(triangle_adjacent(t, le), first, last);
        }
        for (index_t le = 0; le < 3; ++le) {
            if (orient[le] == ZERO) {
                find_conflict_zone_iterative(p, triangle_adjacent(t, le), t_bndry, e_bndry, first,
                                             last);
            }
        }
    }

    find_conflict_zone_iterative(p, t, t_bndry, e_bndry, first, last);
}

void Delaunay2d::find_conflict_zone_iterative(const double* p, index_t t_in, index_t& t_bndry,
                                              index_t& e_bndry, index_t& first, index_t& last)
{
    S_.push_back(t_in);
    while (!S_.empty()) {
        const index_t t = S_.back();
        S_.pop_back();

        for (index_t le = 0; le < 3; ++le) {
            const index_t t2 = triangle_adjacent(t, le);
            if (triangle_is_in_list(t2) || triangle_is_marked(t2)) continue;

            if (triangle_is_conflict(t2, p)) {
                add_triangle_to_list(t2, first, last);
                S_.push_back(t2);
                continue;
            }

            // t is in conflict and t2 is not: this edge is on the conflict-zone boundary.
            t_bndry = t;
            e_bndry = le;
            mark_triangle(t2);
        }
    }
}

index_t Delaunay2d::stellate_conflict_zone(index_t v_in, index_t t1, index_t t1ebord)
{
    index_t t = t1;
    index_t e = t1ebord;
    index_t t_adj = triangle_adjacent(t, e);

    assert(t_adj != NO_INDEX);
    assert(triangle_is_in_list(t));
    assert(!triangle_is_in_list(t_adj));

    index_t new_t_first = NO_INDEX;
    index_t new_t_prev = NO_INDEX;

    do {
        const index_t v1 = triangle_vertex(t, (e + 1) % 3);
        const index_t v2 = triangle_vertex(t, (e + 2) % 3);

        const index_t new_t = new_triangle(v_in, v1, v2);

        // Connect the new triangle across the conflict-zone boundary.
        set_triangle_adjacent(new_t, 0, t_adj);
        const index_t adj_e = find_triangle_adjacent(t_adj, t);
        set_triangle_adjacent(t_adj, adj_e, new_t);

        // Move to the next boundary edge, turning around v2.
        e = (e + 1) % 3;
        t_adj = triangle_adjacent(t, e);
        while (triangle_is_in_list(t_adj)) {
            t = t_adj;
            e = (find_triangle_vertex(t, v2) + 2) % 3;
            t_adj = triangle_adjacent(t, e);
            assert(t_adj != NO_INDEX);
        }

        if (new_t_prev == NO_INDEX) {
            new_t_first = new_t;
        } else {
            set_triangle_adjacent(new_t_prev, 1, new_t);
            set_triangle_adjacent(new_t, 2, new_t_prev);
        }
        new_t_prev = new_t;

    } while ((t != t1) || (e != t1ebord));

    set_triangle_adjacent(new_t_prev, 1, new_t_first);
    set_triangle_adjacent(new_t_first, 2, new_t_prev);

    return new_t_prev;
}

index_t Delaunay2d::insert(index_t v, index_t hint)
{
    index_t t_bndry = NO_TRIANGLE;
    index_t e_bndry = NO_INDEX;
    index_t first_conflict = NO_TRIANGLE;
    index_t last_conflict = NO_TRIANGLE;

    const double* p = vertex_ptr(v);

    int orient[3];
    const index_t t = locate(p, hint, orient);
    find_conflict_zone(v, t, orient, t_bndry, e_bndry, first_conflict, last_conflict);

    // Empty conflict list: vertex v already exists in the triangulation.
    if (first_conflict == END_OF_LIST) return NO_TRIANGLE;

    const index_t new_t = stellate_conflict_zone(v, t_bndry, e_bndry);

    // Recycle the conflict-zone triangles.
    cell_next_[last_conflict] = first_free_;
    first_free_ = first_conflict;

    return new_t;
}

bool Delaunay2d::create_first_triangle(index_t& iv0, index_t& iv1, index_t& iv2)
{
    if (nb_vertices() < 3) return false;

    iv0 = 0;

    iv1 = 1;
    while (iv1 < nb_vertices() &&
           PCK::points_are_identical_2d(vertex_ptr(iv0), vertex_ptr(iv1))) {
        ++iv1;
    }
    if (iv1 == nb_vertices()) return false;

    iv2 = iv1 + 1;
    Sign s = ZERO;
    while (iv2 < nb_vertices() &&
           (s = PCK::orient_2d(vertex_ptr(iv0), vertex_ptr(iv1), vertex_ptr(iv2))) == ZERO) {
        ++iv2;
    }
    if (iv2 == nb_vertices()) return false;
    if (s == NEGATIVE) std::swap(iv1, iv2);

    const index_t t0 = new_triangle(iv0, iv1, iv2);

    // The three virtual triangles surrounding it.
    index_t t[3];
    for (index_t e = 0; e < 3; ++e) {
        // Reverse order, since this is an adjacent triangle.
        const index_t v1 = triangle_vertex(t0, triangle_edge_vertex(e, 1));
        const index_t v2 = triangle_vertex(t0, triangle_edge_vertex(e, 0));
        t[e] = new_triangle(VERTEX_AT_INFINITY, v1, v2);
    }
    for (index_t e = 0; e < 3; ++e) {
        set_triangle_adjacent(t[e], 0, t0);
        set_triangle_adjacent(t0, e, t[e]);
    }
    for (index_t e = 0; e < 3; ++e) {
        const index_t lv1 = triangle_edge_vertex(e, 1);
        const index_t lv2 = triangle_edge_vertex(e, 0);
        set_triangle_adjacent(t[e], 1, t[lv1]);
        set_triangle_adjacent(t[e], 2, t[lv2]);
    }

    return true;
}

// --- driver -----------------------------------------------------------------------------------

bool Delaunay2d::set_vertices(index_t nb_vertices, const double* points)
{
    points_ = points;
    nb_vertices_ = nb_vertices;
    cur_stamp_ = 0;

    const index_t expected_triangles = nb_vertices * 2;
    cell_to_v_store_.clear();
    cell_to_cell_store_.clear();
    cell_next_.clear();
    cell_to_v_store_.reserve(expected_triangles * 3);
    cell_to_cell_store_.reserve(expected_triangles * 3);
    cell_next_.reserve(expected_triangles);
    first_free_ = END_OF_LIST;

    // Spatial sort, so the incremental location walks stay short.
    hilbert_order(nb_vertices, points, reorder_);

    index_t v0, v1, v2;
    if (!create_first_triangle(v0, v1, v2)) return false; // all points collinear

    index_t hint = NO_TRIANGLE;
    for (index_t i = 0; i < nb_vertices; ++i) {
        const index_t v = reorder_[i];
        if (v != v0 && v != v1 && v != v2) {
            const index_t new_hint = insert(v, hint);
            if (new_hint != NO_TRIANGLE) hint = new_hint;
        }
    }

    // Compress: drop the free list, keep the virtual triangles (see delaunay2d.h, note 4).
    //
    // DEVIATION: geogram aliases old2new onto cell_next_ to save memory, relying on
    // old2new[t] <= t so nothing is overwritten before use. A separate array costs 4 bytes per
    // triangle and removes that subtlety.
    std::vector<index_t> old2new(max_t(), NO_INDEX);
    index_t nb_tri = 0;
    for (index_t t = 0; t < max_t(); ++t) {
        if (!triangle_is_free(t)) {
            if (t != nb_tri) {
                for (index_t k = 0; k < 3; ++k) {
                    cell_to_v_store_[nb_tri * 3 + k] = cell_to_v_store_[t * 3 + k];
                    cell_to_cell_store_[nb_tri * 3 + k] = cell_to_cell_store_[t * 3 + k];
                }
            }
            old2new[t] = nb_tri;
            ++nb_tri;
        }
    }
    cell_to_v_store_.resize(3 * nb_tri);
    cell_to_cell_store_.resize(3 * nb_tri);
    for (index_t i = 0; i < 3 * nb_tri; ++i) {
        const index_t t = cell_to_cell_store_[i];
        assert(t != NO_INDEX);
        cell_to_cell_store_[i] = old2new[t];
    }

    // Reorder so that finite triangles occupy [0, nb_finite_triangles_) and virtual ones follow.
    {
        old2new.assign(nb_tri, NO_INDEX);
        nb_finite_triangles_ = 0;
        index_t finite_ptr = 0;
        index_t infinite_ptr = nb_tri - 1;
        for (;;) {
            while (finite_ptr < nb_tri && triangle_is_finite(finite_ptr)) {
                old2new[finite_ptr] = finite_ptr;
                ++finite_ptr;
                ++nb_finite_triangles_;
            }
            while (infinite_ptr != NO_INDEX && !triangle_is_finite(infinite_ptr)) {
                old2new[infinite_ptr] = infinite_ptr;
                --infinite_ptr;
            }
            if (infinite_ptr == NO_INDEX || finite_ptr > infinite_ptr) break;
            old2new[finite_ptr] = infinite_ptr;
            old2new[infinite_ptr] = finite_ptr;
            ++nb_finite_triangles_;
            for (index_t lf = 0; lf < 3; ++lf) {
                std::swap(cell_to_cell_store_[3 * finite_ptr + lf],
                          cell_to_cell_store_[3 * infinite_ptr + lf]);
                std::swap(cell_to_v_store_[3 * finite_ptr + lf],
                          cell_to_v_store_[3 * infinite_ptr + lf]);
            }
            ++finite_ptr;
            --infinite_ptr;
        }
        for (index_t i = 0; i < 3 * nb_tri; ++i) {
            const index_t t = cell_to_cell_store_[i];
            assert(t != NO_INDEX && old2new[t] != NO_INDEX);
            cell_to_cell_store_[i] = old2new[t];
        }
    }

    cell_next_.assign(cell_next_.size(), NO_INDEX);
    return true;
}

bool Delaunay2d::check_combinatorics() const
{
    const index_t nt = nb_triangles();
    for (index_t t = 0; t < nt; ++t) {
        for (index_t le = 0; le < 3; ++le) {
            const index_t t2 = triangle_adjacent(t, le);
            if (t2 == NO_INDEX || t2 >= nt || t2 == t) return false;
            // Reciprocity: t2 must point back at t across exactly one edge, and the two triangles
            // must agree on the shared edge's endpoints (in opposite order).
            index_t back = NO_INDEX;
            for (index_t le2 = 0; le2 < 3; ++le2) {
                if (triangle_adjacent(t2, le2) == t) {
                    if (back != NO_INDEX) return false; // adjacent twice
                    back = le2;
                }
            }
            if (back == NO_INDEX) return false;
            if (triangle_vertex(t, (le + 1) % 3) != triangle_vertex(t2, (back + 2) % 3)) return false;
            if (triangle_vertex(t, (le + 2) % 3) != triangle_vertex(t2, (back + 1) % 3)) return false;
        }
        // Finite triangles must be strictly counterclockwise.
        if (triangle_is_finite(t)) {
            if (PCK::orient_2d(vertex_ptr(triangle_vertex(t, 0)), vertex_ptr(triangle_vertex(t, 1)),
                               vertex_ptr(triangle_vertex(t, 2))) != POSITIVE) {
                return false;
            }
        }
    }
    return true;
}

} // namespace vr2d
} // namespace vol_rem
