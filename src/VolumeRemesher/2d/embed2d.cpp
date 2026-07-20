#include "embed2d.h"

#include "arrangement2d.h"

#include <algorithm>
#include <unordered_map>

namespace vol_rem {

bool embed_seg_in_tri_mesh(const std::vector<double>& seg_vrt_coords,
                           const std::vector<uint32_t>& segment_indexes,
                           std::vector<bigrational>& vertices,
                           std::vector<std::array<uint32_t, 3>>& out_tris,
                           std::vector<std::vector<std::array<uint32_t, 3>>>& out_segment_provenance,
                           std::vector<std::array<uint32_t, 2>>& out_point_provenance,
                           bool verbose)
{
    using namespace vol_rem::vr2d;

    Arrangement2D A;
    if (!build_arrangement(seg_vrt_coords, segment_indexes, A, verbose)) return false;

    // Compact to the vertices actually used by a finite triangle. Only the virtual triangles
    // (outside the convex hull, carrying the symbolic vertex at infinity) are dropped; the four
    // bounding-box corners are kept, since they bound the domain the arrangement covers.
    std::vector<uint32_t> vmap(A.V.size(), INVALID);
    std::vector<uint32_t> used;
    std::vector<uint32_t> tri_id(A.num_triangles(), INVALID);

    out_tris.clear();
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (!A.tri_is_finite(t)) continue;
        std::array<uint32_t, 3> tv{};
        for (uint32_t k = 0; k < 3; k++) {
            const uint32_t v = A.tri_node[3 * t + k];
            if (vmap[v] == INVALID) {
                vmap[v] = uint32_t(used.size());
                used.push_back(v);
            }
            tv[k] = vmap[v];
        }
        tri_id[t] = uint32_t(out_tris.size());
        out_tris.push_back(tv);
    }

    vertices.clear();
    vertices.reserve(2 * used.size());
    for (const uint32_t v : used) {
        bigrational x, y;
        if (!A.V[v]->getExactXYCoordinates(x, y)) return false;
        vertices.push_back(x);
        vertices.push_back(y);
    }

    // Index every emitted edge once; doing this per sub-edge instead would be quadratic.
    std::unordered_map<uint64_t, uint32_t> edge_tri;
    edge_tri.reserve(out_tris.size() * 3);
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (tri_id[t] == INVALID) continue;
        for (uint32_t le = 0; le < 3; le++) {
            const uint64_t k = Arrangement2D::ekey(A.tri_node[3 * t + (le + 1) % 3],
                                                   A.tri_node[3 * t + (le + 2) % 3]);
            const auto it = edge_tri.find(k);
            if (it == edge_tri.end() || tri_id[t] < it->second) edge_tri[k] = tri_id[t];
        }
    }

    out_segment_provenance.assign(A.input_to_seg.size(), {});
    for (uint32_t i = 0; i < A.input_to_seg.size(); i++) {
        const uint32_t s = A.input_to_seg[i];
        if (s == INVALID) continue; // dropped as zero-length; empty list is the answer
        auto& dst = out_segment_provenance[i];
        for (uint32_t k = A.seg_head[s]; k != INVALID; k = A.subedges[k].next) {
            const uint32_t v0 = A.subedges[k].v0, v1 = A.subedges[k].v1;
            const auto it = edge_tri.find(Arrangement2D::ekey(v0, v1));
            if (it == edge_tri.end()) continue;
            dst.push_back({it->second, vmap[v0], vmap[v1]});
        }
        // The stored segment is oriented by its first input occurrence; if THIS input segment
        // was given the other way round, reverse so the list runs from its own first endpoint.
        const uint32_t first = A.input_point_vertex[segment_indexes[2 * i]];
        if (first != A.seg[s][0]) {
            std::reverse(dst.begin(), dst.end());
            for (auto& e : dst) std::swap(e[1], e[2]);
        }
    }

    out_point_provenance.assign(A.input_point_vertex.size(), {UINT32_MAX, UINT32_MAX});
    for (uint32_t i = 0; i < A.input_point_vertex.size(); i++) {
        const uint32_t v = A.input_point_vertex[i];
        if (v == INVALID) continue;
        const uint32_t t = A.vert_tri[v];
        if (t == INVALID || tri_id[t] == INVALID || vmap[v] == INVALID) continue;
        out_point_provenance[i] = {tri_id[t], vmap[v]};
    }

    return true;
}

} // namespace vol_rem
