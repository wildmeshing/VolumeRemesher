#include "arrangement2d.h"

#include "delaunay2d.h"
#include "predicates2d.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

namespace vol_rem {
namespace vr2d {

namespace {

// =============================================================================================
// Triangle store
// =============================================================================================

uint32_t alloc_tri(Arrangement2D& A)
{
    if (!A.free_tris.empty()) {
        // LIFO. The alloc/free sequence is itself deterministic, so this is too.
        const uint32_t t = A.free_tris.back();
        A.free_tris.pop_back();
        A.tri_dead[t] = 0;
        for (uint32_t k = 0; k < 3; k++) {
            A.tri_node[3 * t + k] = INVALID;
            A.tri_neigh[3 * t + k] = INVALID;
            A.tri_cons[3 * t + k] = 0;
        }
        return t;
    }
    const uint32_t t = A.num_triangles();
    A.tri_node.resize(A.tri_node.size() + 3, INVALID);
    A.tri_neigh.resize(A.tri_neigh.size() + 3, INVALID);
    A.tri_cons.resize(A.tri_cons.size() + 3, 0);
    A.tri_dead.push_back(0);
    return t;
}

void free_tri(Arrangement2D& A, uint32_t t)
{
    A.tri_dead[t] = 1;
    A.free_tris.push_back(t);
}

// Rotate counterclockwise around vertex v. Given (t, i) with tri_node[3t+i] == v, the CCW-next
// incident triangle is reached by crossing edge (i+1)%3, which is the edge (v_{i+2}, v_i) -- one
// of the two edges incident to v. Returns the corresponding local index in the new triangle.
inline void rot_ccw(const Arrangement2D& A, uint32_t& t, uint32_t& i)
{
    const uint32_t tn = A.tri_neigh[3 * t + (i + 1) % 3];
    assert(tn != INVALID);
    i = A.local_index(tn, A.tri_node[3 * t + i]);
    t = tn;
}

// Finds the edge (u,w) in the triangulation. Returns false if u and w are not adjacent.
bool find_edge(const Arrangement2D& A, uint32_t u, uint32_t w, uint32_t& t_out, uint32_t& le_out)
{
    uint32_t t = A.vert_tri[u];
    if (t == INVALID) return false;
    uint32_t i = A.local_index(t, u);
    const uint32_t t0 = t, i0 = i;
    do {
        if (A.tri_node[3 * t + (i + 1) % 3] == w) {
            t_out = t;
            le_out = (i + 2) % 3; // edge (v_i, v_{i+1}) = (u, w)
            return true;
        }
        if (A.tri_node[3 * t + (i + 2) % 3] == w) {
            t_out = t;
            le_out = (i + 1) % 3; // edge (v_{i+2}, v_i) = (w, u)
            return true;
        }
        rot_ccw(A, t, i);
    } while (t != t0 || i != i0);
    return false;
}

void set_constrained(Arrangement2D& A, uint32_t t, uint32_t le)
{
    A.tri_cons[3 * t + le] = 1;
    const uint32_t t2 = A.tri_neigh[3 * t + le];
    assert(t2 != INVALID);
    A.tri_cons[3 * t2 + A.local_adj(t2, t)] = 1;
}

// =============================================================================================
// Cavity replacement
//
// Deletes a set of triangles and rebuilds the region from a list of new vertex triples. The new
// adjacency is derived by matching directed edges rather than being threaded through the
// construction, which keeps the two callers (constrained-edge insertion and edge splitting) free
// of stitching logic.
// =============================================================================================

void replace_cavity(Arrangement2D& A, const std::vector<uint32_t>& cav,
                    const std::vector<std::array<uint32_t, 3>>& new_tris)
{
    // Boundary of the cavity: directed edge (as seen from inside) -> (outside triangle, its edge).
    // Lookup-only map, never iterated.
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> bnd;
    const auto dkey = [](uint32_t a, uint32_t b) { return (uint64_t(a) << 32) | b; };

    for (const uint32_t t : cav) {
        for (uint32_t le = 0; le < 3; le++) {
            const uint32_t t2 = A.tri_neigh[3 * t + le];
            if (std::find(cav.begin(), cav.end(), t2) != cav.end()) continue;
            const uint32_t a = A.tri_node[3 * t + (le + 1) % 3];
            const uint32_t b = A.tri_node[3 * t + (le + 2) % 3];
            bnd[dkey(a, b)] = {t2, A.local_adj(t2, t)};
        }
    }

    for (const uint32_t t : cav) free_tri(A, t);

    std::vector<uint32_t> made;
    made.reserve(new_tris.size());
    for (const auto& nt : new_tris) {
        const uint32_t t = alloc_tri(A);
        A.tri_node[3 * t] = nt[0];
        A.tri_node[3 * t + 1] = nt[1];
        A.tri_node[3 * t + 2] = nt[2];
        made.push_back(t);
    }

    // Directed edges of the new triangles, for internal matching.
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> fresh;
    for (const uint32_t t : made) {
        for (uint32_t le = 0; le < 3; le++) {
            const uint32_t a = A.tri_node[3 * t + (le + 1) % 3];
            const uint32_t b = A.tri_node[3 * t + (le + 2) % 3];
            fresh[dkey(a, b)] = {t, le};
        }
    }

    for (const uint32_t t : made) {
        for (uint32_t le = 0; le < 3; le++) {
            const uint32_t a = A.tri_node[3 * t + (le + 1) % 3];
            const uint32_t b = A.tri_node[3 * t + (le + 2) % 3];
            // A neighbour inside the cavity traverses the same edge in the opposite direction.
            const auto it = fresh.find(dkey(b, a));
            if (it != fresh.end()) {
                A.tri_neigh[3 * t + le] = it->second.first;
                continue;
            }
            const auto ob = bnd.find(dkey(a, b));
            assert(ob != bnd.end() && "cavity boundary is not closed");
            const uint32_t t2 = ob->second.first;
            const uint32_t e2 = ob->second.second;
            A.tri_neigh[3 * t + le] = t2;
            A.tri_neigh[3 * t2 + e2] = t;
            // Boundary edges may be constrained; the flag survives on the outside triangle.
            A.tri_cons[3 * t + le] = A.tri_cons[3 * t2 + e2];
        }
        for (uint32_t k = 0; k < 3; k++) A.vert_tri[A.tri_node[3 * t + k]] = t;
    }
}

// =============================================================================================
// Pseudo-polygon triangulation (Anglada 1997)
//
// Fills the polygon a -> P[lo..hi) -> b -> a, which is counterclockwise and whose base edge (b,a)
// is the constrained segment. At each step it picks the chain vertex whose circumcircle with the
// base contains no other chain vertex, then recurses on the two sub-chains.
//
// NOTE ON MESH QUALITY. The requirement is only "a valid triangulation", which plain ear clipping
// would satisfy. Anglada is used anyway because it costs essentially the same (it replaces an
// O(m) containment test with an O(1) incircle) and produces the constrained-Delaunay fill. That
// matters here for a reason specific to this codebase: sliver triangles feed near-degenerate
// inputs to orient2D/incircle, the interval filter fails, and the predicate falls through to the
// bigfloat stage, which allocates and is orders of magnitude slower. Keeping the fill Delaunay
// keeps the filters succeeding. To fall back to ear clipping, drop the incircle comparison below
// and take the first candidate.
//
// DETERMINISM. Ties (cocircular candidates) give incircle == 0, and the strict > keeps the
// earlier index. The chain order comes from the walk, which is deterministic, so the choice is
// total. Deliberately not a std::sort with a comparator that can report "equal".
// =============================================================================================

void triangulate_pseudopolygon(const Arrangement2D& A, uint32_t a, const std::vector<uint32_t>& P,
                               size_t lo, size_t hi, uint32_t b,
                               std::vector<std::array<uint32_t, 3>>& out)
{
    if (lo >= hi) return;
    if (hi - lo == 1) {
        out.push_back({a, P[lo], b});
        return;
    }

    // Restrict to candidates strictly off the base line, so a degenerate triangle can never be
    // emitted. At the top level pass 1 already guarantees this; inside the recursion the base
    // changes, and Anglada's guarantee is thinner there, so it is checked rather than assumed.
    size_t c = SIZE_MAX;
    for (size_t i = lo; i < hi; i++) {
        if (genericPoint::orient2D(*A.V[a], *A.V[P[i]], *A.V[b]) <= 0) continue;
        if (c == SIZE_MAX) {
            c = i;
            continue;
        }
        if (genericPoint::incircle(*A.V[a], *A.V[P[c]], *A.V[b], *A.V[P[i]]) > 0) c = i;
    }
    assert(c != SIZE_MAX && "pseudo-polygon has empty interior");

    triangulate_pseudopolygon(A, a, P, lo, c, P[c], out);
    triangulate_pseudopolygon(A, P[c], P, c + 1, hi, b, out);
    out.push_back({a, P[c], b});
}

// =============================================================================================
// Provenance
// =============================================================================================

void record_subedge(Arrangement2D& A, uint32_t s, uint32_t v0, uint32_t v1)
{
    const uint32_t k = uint32_t(A.subedges.size());
    A.subedges.push_back({s, v0, v1, INVALID});
    if (A.seg_head[s] == INVALID) {
        A.seg_head[s] = k;
    } else {
        A.subedges[A.seg_tail[s]].next = k;
    }
    A.seg_tail[s] = k;

    auto& v = A.edge_to_subedges[Arrangement2D::ekey(v0, v1)];
    v.push_back(k); // ids only ever increase, so the vector stays sorted ascending
}

// Edge (u,w) has been split at q. Every input segment covering it inherits both halves, in order.
void split_provenance(Arrangement2D& A, uint32_t u, uint32_t w, uint32_t q)
{
    const auto it = A.edge_to_subedges.find(Arrangement2D::ekey(u, w));
    if (it == A.edge_to_subedges.end()) return;

    std::vector<uint32_t> lo, hi; // built in ascending id order because it->second is sorted
    for (const uint32_t k : it->second) {
        const uint32_t nk = uint32_t(A.subedges.size());
        A.subedges.push_back({A.subedges[k].seg, q, A.subedges[k].v1, A.subedges[k].next});
        A.subedges[k].v1 = q;
        A.subedges[k].next = nk;
        if (A.seg_tail[A.subedges[k].seg] == k) A.seg_tail[A.subedges[k].seg] = nk;
        if (A.subedges[k].v0 == u) {
            lo.push_back(k);
            hi.push_back(nk);
        } else {
            hi.push_back(k);
            lo.push_back(nk);
        }
    }
    A.edge_to_subedges.erase(it);
    std::sort(lo.begin(), lo.end());
    std::sort(hi.begin(), hi.end());
    A.edge_to_subedges.emplace(Arrangement2D::ekey(u, q), std::move(lo));
    A.edge_to_subedges.emplace(Arrangement2D::ekey(q, w), std::move(hi));
}

// =============================================================================================
// Splitting a constrained edge at a new intersection vertex
// =============================================================================================

// Which input segment supports constrained edge (u,w)? The smallest covering subedge id, which is
// deterministic because edge_to_subedges values are kept sorted.
uint32_t edge_owner_segment(const Arrangement2D& A, uint32_t u, uint32_t w)
{
    const auto it = A.edge_to_subedges.find(Arrangement2D::ekey(u, w));
    assert(it != A.edge_to_subedges.end() && !it->second.empty());
    return A.subedges[it->second.front()].seg;
}

// Creates the intersection of segment s with the constrained edge (u,w), inserts it as a new
// vertex, and returns its id.
//
// The SSI is built from the ORIGINAL endpoints of the two input segments, never from the split
// sub-edges. Those endpoints are explicit points, so every SSI in this pipeline has four explicit
// parents: depth 1, never cascaded. That keeps orient2D off its all-implicit path and lets
// getExpansionLambda take its all-explicit fast path.
//
// The denominator cannot vanish: this is reached only when (u,w) STRICTLY straddles the line of
// s, so the two lines are not parallel. The degenerate case is unreachable by construction rather
// than defended against.
uint32_t split_constrained_edge(Arrangement2D& A, uint32_t t, uint32_t le, uint32_t s)
{
    const uint32_t u = A.tri_node[3 * t + (le + 1) % 3];
    const uint32_t w = A.tri_node[3 * t + (le + 2) % 3];
    const uint32_t r = edge_owner_segment(A, u, w);

    genericPoint* q = A.arena.new_ssi(*A.V[A.seg[r][0]], *A.V[A.seg[r][1]], *A.V[A.seg[s][0]],
                                      *A.V[A.seg[s][1]]);
    const uint32_t qi = uint32_t(A.V.size());
    A.V.push_back(q);
    A.vert_tri.push_back(INVALID);

    // 2 -> 4 split of the two triangles sharing (u,w).
    const uint32_t t2 = A.tri_neigh[3 * t + le];
    const uint32_t e2 = A.local_adj(t2, t);
    const uint32_t x1 = A.tri_node[3 * t + le];   // apex of t,  so t  is (x1, u, w) CCW
    const uint32_t x2 = A.tri_node[3 * t2 + e2];  // apex of t2, so t2 is (x2, w, u) CCW

    const std::vector<uint32_t> cav = {t, t2};
    const std::vector<std::array<uint32_t, 3>> made = {
        {x1, u, qi}, {x1, qi, w}, {x2, w, qi}, {x2, qi, u}};
    replace_cavity(A, cav, made);

    split_provenance(A, u, w, qi);

    uint32_t ft, fe;
    if (find_edge(A, u, qi, ft, fe)) set_constrained(A, ft, fe);
    if (find_edge(A, qi, w, ft, fe)) set_constrained(A, ft, fe);

    return qi;
}

// =============================================================================================
// Pass 1: walk the segment, collecting the mesh vertices on it
// =============================================================================================

// Is x strictly forward of cur, toward b? Sign of (x - cur) . (b - cur).
inline bool forward(const Arrangement2D& A, uint32_t x, uint32_t cur, uint32_t b)
{
    return genericPoint::dotProductSign2D(*A.V[x], *A.V[b], *A.V[cur]) > 0;
}

// Returns the next mesh vertex along segment s after `cur`, creating an intersection vertex if
// the segment crosses an already-constrained edge on the way.
uint32_t advance(Arrangement2D& A, uint32_t cur, uint32_t s)
{
    const uint32_t a = A.seg[s][0], b = A.seg[s][1];
    const genericPoint& PA = *A.V[a];
    const genericPoint& PB = *A.V[b];

    // --- leave the star of cur -------------------------------------------------------------
    uint32_t t = A.vert_tri[cur];
    assert(t != INVALID);
    uint32_t i = A.local_index(t, cur);
    const uint32_t t0 = t, i0 = i;

    uint32_t start_t = INVALID, start_e = INVALID;
    do {
        const uint32_t u = A.tri_node[3 * t + (i + 1) % 3];
        const uint32_t w = A.tri_node[3 * t + (i + 2) % 3];

        // b adjacent: done. Checked by index before any predicate.
        if (u == b || w == b) return b;

        if (u != INVALID && w != INVALID) {
            // The wedge at cur in the CCW triangle (cur,u,w) contains the direction to b iff
            // orient2D(cur,u,b) > 0 and orient2D(cur,w,b) < 0.
            const int su = genericPoint::orient2D(*A.V[cur], *A.V[u], PB);
            const int sw = genericPoint::orient2D(*A.V[cur], *A.V[w], PB);

            // Collinear with an incident edge, in the forward direction: the segment runs along
            // that edge. The backward case must be rejected explicitly or the walk reverses.
            if (su == 0 && forward(A, u, cur, b)) return u;
            if (sw == 0 && forward(A, w, cur, b)) return w;

            if (su > 0 && sw < 0) {
                start_t = t;
                start_e = i; // edge i of t is (v_{i+1}, v_{i+2}) = (u, w)
                break;
            }
        }
        rot_ccw(A, t, i);
    } while (t != t0 || i != i0);

    assert(start_t != INVALID && "segment direction left no wedge at its own endpoint");

    // --- cross the strip ---------------------------------------------------------------------
    t = start_t;
    uint32_t e = start_e;
    // The wedge test selected (cur,u,w) with orient2D(cur,u,b) > 0, and orient2D(cur,u,b) =
    // -orient2D(cur,b,u), so u = v_{e+1} lies on the NEGATIVE side of the directed line a->b and
    // w = v_{e+2} on the positive side. (cur is itself on that line, so measuring from cur and
    // measuring from a give the same sign.)
    uint32_t neg = A.tri_node[3 * t + (e + 1) % 3];
    uint32_t pos = A.tri_node[3 * t + (e + 2) % 3];

    for (;;) {
        if (A.tri_cons[3 * t + e]) return split_constrained_edge(A, t, e, s);

        const uint32_t t2 = A.tri_neigh[3 * t + e];
        assert(t2 != INVALID);
        const uint32_t bk = A.local_adj(t2, t);
        const uint32_t x = A.tri_node[3 * t2 + bk]; // apex of the triangle we are entering

        if (x == b) return b;
        assert(x != INVALID && "walk left the convex hull");

        const int ox = genericPoint::orient2D(PA, PB, *A.V[x]);

        // An apex exactly on the line blocks the exit, so the segment passes through it. It is
        // necessarily strictly between cur and b: the entry edge strictly straddles the line, so
        // b is not one of its endpoints, and if b were in this triangle it would have to be x.
        if (ox == 0) return x;

        const uint32_t keep = (ox > 0) ? neg : pos;
        const uint32_t s1 = A.tri_node[3 * t2 + (bk + 1) % 3];
        const uint32_t s2 = A.tri_node[3 * t2 + (bk + 2) % 3];
        // Edge (bk+1)%3 of t2 is (s2, x); edge (bk+2)%3 is (x, s1).
        e = (s2 == keep) ? (bk + 1) % 3 : (bk + 2) % 3;
        t = t2;
        if (ox > 0) {
            pos = x;
        } else {
            neg = x;
        }
    }
}

// =============================================================================================
// Pass 2: force one edge of the chain into the triangulation
// =============================================================================================

void insert_constrained_edge(Arrangement2D& A, uint32_t u, uint32_t w, uint32_t s)
{
    // Already present. This branch absorbs duplicate segments, full overlap, partial collinear
    // overlap, and every edge pass 1 stepped along.
    uint32_t ft, fe;
    if (find_edge(A, u, w, ft, fe)) {
        set_constrained(A, ft, fe);
        record_subedge(A, s, u, w);
        return;
    }

    const genericPoint& PA = *A.V[u];
    const genericPoint& PB = *A.V[w];

    // Find the triangle at u that the segment enters.
    uint32_t t = A.vert_tri[u];
    uint32_t i = A.local_index(t, u);
    const uint32_t t0 = t, i0 = i;
    uint32_t e = INVALID;
    do {
        const uint32_t p = A.tri_node[3 * t + (i + 1) % 3];
        const uint32_t q = A.tri_node[3 * t + (i + 2) % 3];
        if (p != INVALID && q != INVALID) {
            if (genericPoint::orient2D(*A.V[u], *A.V[p], PB) > 0 &&
                genericPoint::orient2D(*A.V[u], *A.V[q], PB) < 0) {
                e = i;
                break;
            }
        }
        rot_ccw(A, t, i);
    } while (t != t0 || i != i0);
    assert(e != INVALID);

    std::vector<uint32_t> cav;
    std::vector<uint32_t> chain_pos, chain_neg;

    uint32_t pos = A.tri_node[3 * t + (e + 2) % 3];
    uint32_t neg = A.tri_node[3 * t + (e + 1) % 3];
    cav.push_back(t);
    chain_pos.push_back(pos);
    chain_neg.push_back(neg);

    for (;;) {
        assert(!A.tri_cons[3 * t + e] && "pass 1 should have removed every constrained crossing");
        const uint32_t t2 = A.tri_neigh[3 * t + e];
        const uint32_t bk = A.local_adj(t2, t);
        const uint32_t x = A.tri_node[3 * t2 + bk];
        cav.push_back(t2);
        if (x == w) break;

        const int ox = genericPoint::orient2D(PA, PB, *A.V[x]);
        assert(ox != 0 && "pass 1 should have removed every vertex on the segment");

        const uint32_t keep = (ox > 0) ? neg : pos;
        const uint32_t s2v = A.tri_node[3 * t2 + (bk + 2) % 3];
        e = (s2v == keep) ? (bk + 1) % 3 : (bk + 2) % 3;
        t = t2;
        if (ox > 0) {
            pos = x;
            chain_pos.push_back(x);
        } else {
            neg = x;
            chain_neg.push_back(x);
        }
    }

    // Negative side: the polygon [u] + chain_neg + [w] is counterclockwise, base edge (u,w).
    // Positive side: [w] + reverse(chain_pos) + [u] is counterclockwise, base edge (w,u).
    std::vector<std::array<uint32_t, 3>> made;
    triangulate_pseudopolygon(A, u, chain_neg, 0, chain_neg.size(), w, made);
    std::reverse(chain_pos.begin(), chain_pos.end());
    triangulate_pseudopolygon(A, w, chain_pos, 0, chain_pos.size(), u, made);

    replace_cavity(A, cav, made);

    if (find_edge(A, u, w, ft, fe)) {
        set_constrained(A, ft, fe);
    } else {
        assert(false && "constrained edge missing after cavity fill");
    }
    record_subedge(A, s, u, w);
}

} // namespace

// =============================================================================================
// Driver
// =============================================================================================

bool build_arrangement(const std::vector<double>& seg_coords,
                       const std::vector<uint32_t>& seg_indexes, Arrangement2D& A, bool verbose)
{
    const uint32_t n_in = uint32_t(seg_coords.size() / 2);
    if (n_in == 0) return false;

    // --- Phase 0a: deduplicate points ----------------------------------------------------------
    // Exact double comparison is complete here: two different doubles are never geometrically
    // equal, and equal doubles are equal. No predicate needed. The index tie-break makes the
    // comparator a strict total order, so the result does not depend on sort stability.
    std::vector<uint32_t> order(n_in);
    for (uint32_t i = 0; i < n_in; i++) order[i] = i;
    const auto lex = [&](uint32_t p, uint32_t q) {
        const double px = seg_coords[2 * p], py = seg_coords[2 * p + 1];
        const double qx = seg_coords[2 * q], qy = seg_coords[2 * q + 1];
        if (px != qx) return px < qx;
        if (py != qy) return py < qy;
        return p < q;
    };
    std::sort(order.begin(), order.end(), lex);

    A.input_point_vertex.assign(n_in, INVALID);
    A.coords.clear();
    for (uint32_t k = 0; k < n_in; k++) {
        const uint32_t p = order[k];
        // -0.0 and +0.0 compare equal, so they merge; normalise so the stored coordinate is
        // reproducible regardless of which of the two came first.
        double x = seg_coords[2 * p], y = seg_coords[2 * p + 1];
        if (x == 0.0) x = 0.0;
        if (y == 0.0) y = 0.0;
        if (k > 0) {
            const uint32_t prev = order[k - 1];
            if (seg_coords[2 * prev] == seg_coords[2 * p] &&
                seg_coords[2 * prev + 1] == seg_coords[2 * p + 1]) {
                A.input_point_vertex[p] = A.input_point_vertex[prev];
                continue;
            }
        }
        A.input_point_vertex[p] = uint32_t(A.coords.size() / 2);
        A.coords.push_back(x);
        A.coords.push_back(y);
    }
    A.num_input_pts = uint32_t(A.coords.size() / 2);

    // --- Phase 0b: the four corners of the 10%-expanded bounding box ----------------------------
    // Geogram's Delaunay2d surrounds the convex hull with virtual triangles. Placing the four
    // corners first makes every input point strictly interior, so the segment walks never reach a
    // virtual triangle and need no boundary special case.
    {
        double xmin = A.coords[0], xmax = A.coords[0];
        double ymin = A.coords[1], ymax = A.coords[1];
        for (uint32_t i = 1; i < A.num_input_pts; i++) {
            xmin = std::min(xmin, A.coords[2 * i]);
            xmax = std::max(xmax, A.coords[2 * i]);
            ymin = std::min(ymin, A.coords[2 * i + 1]);
            ymax = std::max(ymax, A.coords[2 * i + 1]);
        }
        double dx = (xmax - xmin) * 0.1;
        double dy = (ymax - ymin) * 0.1;
        // A degenerate box (all points on a line, or a single point) still needs a real rectangle.
        if (!(dx > 0.0)) dx = (std::abs(xmax) > 0.0) ? std::abs(xmax) * 0.1 : 1.0;
        if (!(dy > 0.0)) dy = (std::abs(ymax) > 0.0) ? std::abs(ymax) * 0.1 : 1.0;
        const double cx[4] = {xmin - dx, xmax + dx, xmax + dx, xmin - dx};
        const double cy[4] = {ymin - dy, ymin - dy, ymax + dy, ymax + dy};
        for (int k = 0; k < 4; k++) {
            A.coords.push_back(cx[k]);
            A.coords.push_back(cy[k]);
        }
    }
    A.num_explicit_pts = uint32_t(A.coords.size() / 2);

    A.V.clear();
    A.V.reserve(A.num_explicit_pts);
    for (uint32_t i = 0; i < A.num_explicit_pts; i++)
        A.V.push_back(A.arena.new_explicit(A.coords[2 * i], A.coords[2 * i + 1]));

    // --- Phase 0c: segments ---------------------------------------------------------------------
    const uint32_t n_seg_in = uint32_t(seg_indexes.size() / 2);
    A.input_to_seg.assign(n_seg_in, INVALID);
    {
        // Group duplicates. Sorting by (min,max,id) is a strict total order.
        std::vector<uint32_t> sorder;
        sorder.reserve(n_seg_in);
        std::vector<std::array<uint32_t, 2>> key(n_seg_in);
        for (uint32_t i = 0; i < n_seg_in; i++) {
            const uint32_t v0 = A.input_point_vertex[seg_indexes[2 * i]];
            const uint32_t v1 = A.input_point_vertex[seg_indexes[2 * i + 1]];
            key[i] = {std::min(v0, v1), std::max(v0, v1)};
            if (v0 != v1) sorder.push_back(i); // drop zero-length segments
        }
        std::sort(sorder.begin(), sorder.end(), [&](uint32_t p, uint32_t q) {
            if (key[p] != key[q]) return key[p] < key[q];
            return p < q;
        });
        for (size_t k = 0; k < sorder.size(); k++) {
            const uint32_t i = sorder[k];
            if (k > 0 && key[sorder[k - 1]] == key[i]) {
                const uint32_t si = A.input_to_seg[sorder[k - 1]];
                A.input_to_seg[i] = si;
                A.seg_inputs[si].push_back(i);
                continue;
            }
            A.input_to_seg[i] = uint32_t(A.seg.size());
            // Orient the stored segment by its first input occurrence, so provenance parameters
            // run from the input's own first endpoint.
            A.seg.push_back({A.input_point_vertex[seg_indexes[2 * i]],
                             A.input_point_vertex[seg_indexes[2 * i + 1]]});
            A.seg_inputs.push_back({i});
        }
    }

    if (verbose) {
        printf("2D arrangement: %u input points -> %u unique, %u segments -> %u unique\n", n_in,
               A.num_input_pts, n_seg_in, uint32_t(A.seg.size()));
    }

    // --- Phase A: Delaunay ----------------------------------------------------------------------
    {
        Delaunay2d D;
        if (!D.set_vertices(A.num_explicit_pts, A.coords.data())) return false;
        D.steal_arrays(A.tri_node, A.tri_neigh);
    }
    A.tri_cons.assign(A.tri_node.size(), 0);
    A.tri_dead.assign(A.tri_node.size() / 3, 0);
    A.vert_tri.assign(A.num_explicit_pts, INVALID);
    for (uint32_t t = 0; t < A.num_triangles(); t++)
        for (uint32_t k = 0; k < 3; k++) {
            const uint32_t v = A.tri_node[3 * t + k];
            if (v != INVALID) A.vert_tri[v] = t;
        }

    // --- Phase B: insert the segments -----------------------------------------------------------
    A.seg_head.assign(A.seg.size(), INVALID);
    A.seg_tail.assign(A.seg.size(), INVALID);

    std::vector<uint32_t> chain;
    for (uint32_t s = 0; s < A.seg.size(); s++) {
        const uint32_t a = A.seg[s][0], b = A.seg[s][1];

        chain.clear();
        chain.push_back(a);
        uint32_t cur = a;
        while (cur != b) {
            cur = advance(A, cur, s);
            chain.push_back(cur);
        }

        for (size_t j = 0; j + 1 < chain.size(); j++)
            insert_constrained_edge(A, chain[j], chain[j + 1], s);

        if (verbose && (s % 20000 == 0) && s) printf("  ... %u / %zu segments\n", s, A.seg.size());
    }

    return true;
}

// =============================================================================================
// Validation
// =============================================================================================

bool check_arrangement(const Arrangement2D& A, std::string* error)
{
    const auto fail = [&](const std::string& m) {
        if (error) *error = m;
        return false;
    };

    const uint32_t nt = A.num_triangles();
    for (uint32_t t = 0; t < nt; t++) {
        if (A.tri_is_dead(t)) continue;
        for (uint32_t le = 0; le < 3; le++) {
            const uint32_t t2 = A.tri_neigh[3 * t + le];
            if (t2 == INVALID || t2 >= nt || t2 == t || A.tri_is_dead(t2))
                return fail("bad adjacency at triangle " + std::to_string(t));
            uint32_t back = INVALID;
            for (uint32_t k = 0; k < 3; k++)
                if (A.tri_neigh[3 * t2 + k] == t) {
                    if (back != INVALID) return fail("doubly adjacent at " + std::to_string(t));
                    back = k;
                }
            if (back == INVALID) return fail("non-reciprocal adjacency at " + std::to_string(t));
            if (A.tri_node[3 * t + (le + 1) % 3] != A.tri_node[3 * t2 + (back + 2) % 3] ||
                A.tri_node[3 * t + (le + 2) % 3] != A.tri_node[3 * t2 + (back + 1) % 3])
                return fail("shared edge mismatch at triangle " + std::to_string(t));
            if (A.tri_cons[3 * t + le] != A.tri_cons[3 * t2 + back])
                return fail("constrained flag not mirrored at triangle " + std::to_string(t));
        }
        if (A.tri_is_finite(t)) {
            if (genericPoint::orient2D(*A.V[A.tri_node[3 * t]], *A.V[A.tri_node[3 * t + 1]],
                                       *A.V[A.tri_node[3 * t + 2]]) <= 0)
                return fail("triangle " + std::to_string(t) + " is not strictly counterclockwise");
        }
    }
    return true;
}

} // namespace vr2d
} // namespace vol_rem
