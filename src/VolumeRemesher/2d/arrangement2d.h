// 2D arrangement: a soup of double-precision segments -> a triangulation in which every input
// segment is exactly a union of output triangle edges, with provenance.
//
// This is the 2D counterpart of BSP.cpp, but it does NOT follow the 3D strategy. The 3D pipeline
// classifies constraints against tetrahedra and then splits BSP cells; here we insert each
// segment directly by walking the triangulation. That is much simpler in 2D and avoids
// materialising the intersection set up front (which would be quadratic in the worst case).
//
// PIPELINE
//   Phase 0  preprocess: dedup points, drop degenerate segments, group duplicate segments,
//            add the four corners of the 10%-expanded bounding box
//   Phase A  Delaunay triangulation of every point (delaunay2d.h)
//   Phase B  for each segment, in id order:
//              pass 1  walk A -> B collecting the mesh vertices that lie on it, creating a
//                      new SSI vertex wherever it crosses an already-constrained edge
//              pass 2  for each consecutive pair of chain vertices, force that edge into the
//                      triangulation
//
// WHY TWO PASSES. A single-pass "collect the cavity and split constrained edges as you go" loop
// has to reconcile a half-built cavity with a topology change inside it. Splitting gives clean
// invariants instead: after pass 1, no mesh vertex lies strictly inside any (chain[j],chain[j+1])
// and no constrained edge crosses it, so pass 2's walk has exactly ONE case to handle -- it
// crosses the interior of an unconstrained edge. Every degeneracy is absorbed by pass 1.
//
// Pass 2 cannot invalidate that for later pairs: it adds no vertices, and it cannot destroy a
// constrained edge (one crossing (v_j,v_j+1) would contradict pass 1; one with an endpoint
// strictly inside the cavity is impossible because a cavity has no interior vertices; and a
// constrained chord across the cavity is impossible because every edge of a crossed triangle is
// either crossed -- hence unconstrained -- or on the cavity boundary).
//
// DETERMINISM. This repository guarantees byte-identical output on Linux/macOS/Windows. Every
// ordering decision here is therefore a strict total order on integer ids. The two hash maps are
// LOOKUP-ONLY and their value vectors are kept in ascending id order; nothing ever iterates an
// unordered container in an order that reaches the output.

#ifndef VOLUMEREMESHER_2D_ARRANGEMENT2D_H
#define VOLUMEREMESHER_2D_ARRANGEMENT2D_H

#include <VolumeRemesher/implicit_point.h>

#include <array>
#include <cstdint>
#include <deque>
#include <unordered_map>
#include <vector>

namespace vol_rem {
namespace vr2d {

static constexpr uint32_t INVALID = UINT32_MAX;

// Address-stable storage for the points.
//
// implicitPoint2D_SSI holds `const genericPoint&` references to its four parents
// (implicit_point.h). If the points lived in a std::vector that reallocated, every SSI ever built
// would silently dangle -- undefined behaviour that a hash comparison could not detect. std::deque
// never invalidates references to existing elements on push_back, which is exactly the guarantee
// needed here.
class PointArena
{
public:
    explicitPoint2D* new_explicit(double x, double y)
    {
        epts_.emplace_back(x, y);
        return &epts_.back();
    }
    // The four parents must outlive the result. In this pipeline they are always the original
    // explicit endpoints of two input segments, so SSI points are never cascaded (see
    // arrangement2d.cpp, split_constrained_edge).
    implicitPoint2D_SSI* new_ssi(const genericPoint& a, const genericPoint& b,
                                 const genericPoint& c, const genericPoint& d)
    {
        ipts_.emplace_back(a, b, c, d);
        return &ipts_.back();
    }
    size_t num_explicit() const { return epts_.size(); }
    size_t num_ssi() const { return ipts_.size(); }

private:
    std::deque<explicitPoint2D> epts_;
    std::deque<implicitPoint2D_SSI> ipts_;
};

// One output edge of one input segment. Split-stable: when a later segment cuts this edge, the
// record is split in two and the halves stay linked in order along the segment.
struct SubEdge
{
    uint32_t seg;      // index into Arrangement2D::seg
    uint32_t v0, v1;   // oriented: v0 is the endpoint nearer seg's first endpoint
    uint32_t next;     // next SubEdge of the same segment, INVALID at the tail
};

struct Arrangement2D
{
    // --- geometry ---------------------------------------------------------------------------
    PointArena arena;
    std::vector<genericPoint*> V;   // vertex id -> point (EXPLICIT2D, or SSI for intersections)
    std::vector<double> coords;     // 2 * num_explicit_pts, the Delaunay input

    uint32_t num_input_pts = 0;     // V[0 .. num_input_pts) are the deduplicated input points
    uint32_t num_explicit_pts = 0;  // + the four bounding-box corners; SSI vertices follow

    // --- input segments, after preprocessing ------------------------------------------------
    std::vector<std::array<uint32_t, 2>> seg;        // unique segments, as vertex-id pairs
    std::vector<std::vector<uint32_t>> seg_inputs;   // -> the input segment ids that map to each
    std::vector<uint32_t> input_to_seg;              // input segment id -> index in seg, or INVALID
    // Per INPUT segment, the vertex id of ITS OWN first endpoint. seg[] is oriented by the first
    // input occurrence, so a duplicate given the other way round needs its provenance reversed;
    // this is what lets the writers do that. INVALID for dropped (zero-length) segments.
    std::vector<uint32_t> input_seg_v0;
    std::vector<uint32_t> input_point_vertex;        // input point id -> vertex id

    // --- triangulation ------------------------------------------------------------------------
    // Geogram's convention, inherited from Delaunay2d: triangle t has vertices tri_node[3t+k] in
    // counterclockwise order, and edge le is the one OPPOSITE vertex le, i.e. the edge
    // (tri_node[3t+(le+1)%3], tri_node[3t+(le+2)%3]). tri_neigh[3t+le] is the adjacent TRIANGLE
    // index (not a corner); the shared edge appears in the opposite order in the neighbour.
    std::vector<uint32_t> tri_node;
    std::vector<uint32_t> tri_neigh;
    std::vector<uint8_t> tri_cons;   // 3 per triangle: 1 if that edge is constrained
    std::vector<uint8_t> tri_dead;   // 1 if the slot is on the free list
    std::vector<uint32_t> vert_tri;  // per vertex: some incident triangle
    std::vector<uint32_t> free_tris; // recycled triangle slots

    // --- provenance ---------------------------------------------------------------------------
    std::vector<SubEdge> subedges;      // append-only, so ids are stable and deterministic
    std::vector<uint32_t> seg_head, seg_tail;

    // Lookup-only index, keyed by unordered vertex pair. Keyed by VERTEX PAIR rather than by
    // (triangle, local edge) because pass 2 destroys and recreates triangles, so corner ids of a
    // surviving constrained edge change -- vertex ids never do.
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_subedges;

    static uint64_t ekey(uint32_t a, uint32_t b)
    {
        return (a < b) ? ((uint64_t(a) << 32) | b) : ((uint64_t(b) << 32) | a);
    }

    // --- accessors ----------------------------------------------------------------------------
    uint32_t num_triangles() const { return uint32_t(tri_node.size() / 3); }
    // NOTE: a *virtual* triangle (one covering the region outside the convex hull) carries
    // Delaunay2d::VERTEX_AT_INFINITY == INVALID as its vertex 0, so deadness cannot be inferred
    // from the vertex array; it needs its own flag.
    bool tri_is_dead(uint32_t t) const { return tri_dead[t] != 0; }
    bool tri_is_finite(uint32_t t) const
    {
        return !tri_is_dead(t) && tri_node[3 * t] != INVALID && tri_node[3 * t + 1] != INVALID &&
               tri_node[3 * t + 2] != INVALID;
    }
    uint32_t local_index(uint32_t t, uint32_t v) const
    {
        const uint32_t* T = &tri_node[3 * t];
        return uint32_t((T[1] == v) | ((T[2] == v) * 2));
    }
    uint32_t local_adj(uint32_t t, uint32_t t2) const
    {
        const uint32_t* T = &tri_neigh[3 * t];
        return uint32_t((T[1] == t2) | ((T[2] == t2) * 2));
    }
};

// Build the arrangement.
//
//   seg_coords    2 doubles per input point (x,y)
//   seg_indexes   2 indices per input segment, into seg_coords
//
// Returns false only if the input has fewer than one usable point.
bool build_arrangement(const std::vector<double>& seg_coords,
                       const std::vector<uint32_t>& seg_indexes, Arrangement2D& A, bool verbose);

// Validation used by the tests and by the debug builds: every finite triangle strictly positive,
// adjacency reciprocal, constrained flags mirrored, and no vertex strictly inside an edge.
bool check_arrangement(const Arrangement2D& A, std::string* error);

} // namespace vr2d
} // namespace vol_rem

#endif
