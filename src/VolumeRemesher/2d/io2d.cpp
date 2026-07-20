#include "io2d.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <array>
#include <fstream>
#include <unordered_map>

namespace vol_rem {
namespace vr2d {

namespace {

// Portable decimal "[-]num[/den]", independent of the bignum backend. Mirrors
// rational_to_string in BSP.cpp.
inline std::string rational_to_string(const bigrational& r)
{
#ifdef USE_GNU_GMP_CLASSES
    return r.get_str();
#else
    return r.get_dec_str();
#endif
}

} // namespace

bool read_OBJ_segments(const std::string& path, std::vector<double>& coords,
                       std::vector<uint32_t>& indexes)
{
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;

    std::vector<double> verts; // (x,y) pairs, in file order, 1-based when referenced
    std::vector<std::array<uint32_t, 2>> lines;

    char buf[8192];
    while (fgets(buf, sizeof(buf), f)) {
        const char* p = buf;
        while (*p == ' ' || *p == '\t') p++;
        if (p[0] == 'v' && (p[1] == ' ' || p[1] == '\t')) {
            double x = 0, y = 0, z = 0;
            if (sscanf(p + 1, "%lf %lf %lf", &x, &y, &z) < 2) {
                fclose(f);
                return false;
            }
            verts.push_back(x);
            verts.push_back(y); // z is dropped: the models are planar curves in 3D
        } else if (p[0] == 'l' && (p[1] == ' ' || p[1] == '\t')) {
            // "l" may carry a polyline of more than two indices; emit consecutive pairs.
            std::vector<long> idx;
            const char* q = p + 1;
            while (*q) {
                while (*q == ' ' || *q == '\t') q++;
                if (!*q || *q == '\n' || *q == '\r') break;
                char* end = nullptr;
                const long v = strtol(q, &end, 10);
                if (end == q) break;
                idx.push_back(v);
                q = end;
                while (*q && *q != ' ' && *q != '\t' && *q != '\n' && *q != '\r') q++; // skip v/vt
            }
            for (size_t k = 0; k + 1 < idx.size(); k++) {
                long a = idx[k], b = idx[k + 1];
                const long nv = long(verts.size() / 2);
                if (a < 0) a += nv + 1; // OBJ allows negative (relative) indices
                if (b < 0) b += nv + 1;
                if (a < 1 || b < 1 || a > nv || b > nv) {
                    fclose(f);
                    return false;
                }
                lines.push_back({uint32_t(a - 1), uint32_t(b - 1)});
            }
        }
        // Everything else (f, vn, vt, comments, groups) is ignored.
    }
    fclose(f);

    if (verts.empty()) return false;

    coords = std::move(verts);
    indexes.clear();
    indexes.reserve(lines.size() * 2);
    for (const auto& l : lines) {
        indexes.push_back(l[0]);
        indexes.push_back(l[1]);
    }
    return true;
}

bool write_arrangement(const std::string& base, const Arrangement2D& A, bool export_rational)
{
    // Emit every finite triangle. That includes the ones incident to the four bounding-box
    // corners: the corners define the domain the arrangement was computed over (the 10%-expanded
    // box), so keeping them yields a clean triangulation of that whole rectangle. Only the
    // virtual triangles -- the ones outside the convex hull, carrying the symbolic vertex at
    // infinity -- are dropped.
    std::vector<uint32_t> vmap(A.V.size(), INVALID);
    std::vector<uint32_t> used;
    std::vector<std::array<uint32_t, 3>> tris;
    std::vector<uint32_t> tri_id(A.num_triangles(), INVALID);

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
        tri_id[t] = uint32_t(tris.size());
        tris.push_back(tv);
    }

    // Binary mode so line endings stay LF on every OS (byte-identical output).
    std::ofstream f(base, std::ios::binary);
    if (!f) return false;
    // Full double precision: the default 6 significant figures collapses fine meshes.
    f.precision(17);
    f << "OFF\n" << used.size() << " " << tris.size() << " 0\n";
    for (const uint32_t v : used) {
        double x = 0, y = 0;
        A.V[v]->getApproxXYCoordinates(x, y);
        f << x << " " << y << " 0\n";
    }
    for (const auto& t : tris) f << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    f.close();

    if (export_rational) {
        std::ofstream rf(base + ".rational", std::ios::binary);
        if (!rf) return false;
        rf << used.size() << "\n";
        bigrational rx, ry;
        for (const uint32_t v : used) {
            A.V[v]->getExactXYCoordinates(rx, ry);
            rf << rational_to_string(rx) << " " << rational_to_string(ry) << "\n";
        }
    }

    // One pass to index every emitted edge by its (unordered) vertex pair. Doing this per
    // sub-edge instead would be quadratic in the mesh size. Lookup-only map; the smallest
    // triangle id wins so the reported representative is deterministic.
    std::unordered_map<uint64_t, uint32_t> edge_tri;
    edge_tri.reserve(tris.size() * 3);
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (tri_id[t] == INVALID) continue;
        for (uint32_t le = 0; le < 3; le++) {
            const uint64_t k = Arrangement2D::ekey(A.tri_node[3 * t + (le + 1) % 3],
                                                   A.tri_node[3 * t + (le + 2) % 3]);
            const auto it = edge_tri.find(k);
            if (it == edge_tri.end() || tri_id[t] < it->second) edge_tri[k] = tri_id[t];
        }
    }

    // Segment provenance, keyed by INPUT segment id (not by the deduplicated segment), so the
    // caller can index it with the id it supplied. Duplicated segments therefore each get the
    // full edge list, and a segment dropped as zero-length gets an empty one.
    {
        std::ofstream pf(base + ".segmentprov", std::ios::binary);
        if (!pf) return false;
        pf << A.input_to_seg.size() << "\n";
        for (uint32_t i = 0; i < A.input_to_seg.size(); i++) {
            const uint32_t s = A.input_to_seg[i];
            pf << i;
            if (s == INVALID) {
                pf << " 0\n";
                continue;
            }
            std::vector<std::array<uint32_t, 3>> out;
            for (uint32_t k = A.seg_head[s]; k != INVALID; k = A.subedges[k].next) {
                const uint32_t v0 = A.subedges[k].v0, v1 = A.subedges[k].v1;
                const auto it = edge_tri.find(Arrangement2D::ekey(v0, v1));
                if (it == edge_tri.end()) continue;
                out.push_back({it->second, vmap[v0], vmap[v1]});
            }
            pf << " " << out.size();
            for (const auto& e : out) pf << " " << e[0] << " " << e[1] << " " << e[2];
            pf << "\n";
        }
    }

    // Point provenance: one line per INPUT point, giving an output triangle it belongs to and its
    // output vertex index, or "-1 -1" if it did not survive.
    {
        std::ofstream pf(base + ".pointprov", std::ios::binary);
        if (!pf) return false;
        pf << A.input_point_vertex.size() << "\n";
        for (uint32_t i = 0; i < A.input_point_vertex.size(); i++) {
            const uint32_t v = A.input_point_vertex[i];
            const uint32_t ov = (v != INVALID) ? vmap[v] : INVALID;
            const uint32_t t = (v != INVALID) ? A.vert_tri[v] : INVALID;
            const uint32_t ot = (t != INVALID) ? tri_id[t] : INVALID;
            if (ov == INVALID || ot == INVALID) pf << i << " -1 -1\n";
            else pf << i << " " << ot << " " << ov << "\n";
        }
    }

    return true;
}

} // namespace vr2d
} // namespace vol_rem
