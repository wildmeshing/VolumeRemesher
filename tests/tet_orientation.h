#pragma once
// Combinatorial orientation / manifoldness check for a tetrahedral mesh.
//
// A valid tetrahedralization is a simplicial complex: every triangular face is shared by
// exactly one tetrahedron (a boundary face) or two (an internal face), and an internal face
// must appear with OPPOSITE winding as seen from its two tets. If two tets give a shared face
// the SAME winding they lie on the same side of it -- i.e. they overlap -- which is an invalid
// mesh even if every tet has positive volume on its own. This check catches exactly that (and
// non-manifold faces shared by >2 tets). It is purely combinatorial: it uses only the tet
// vertex indices, no coordinates.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace vrtest {

struct OrientationResult {
    bool ok = false;
    uint64_t tets = 0;
    uint64_t boundary = 0;     // triangular faces shared by exactly 1 tet
    uint64_t internal = 0;     // faces shared by exactly 2 tets (opposite winding)
    uint64_t over_shared = 0;  // faces shared by >2 tets (non-manifold)
    uint64_t same_winding = 0; // internal faces given the SAME winding by both tets (overlap)
    std::string error;         // non-empty if the mesh could not be read
};

namespace detail {
// Canonical cyclic rotation (smallest vertex first) -- PRESERVES winding, so (a,b,c) and its
// reverse (a,c,b) get different keys.
inline std::array<uint32_t, 3> rot_min(uint32_t a, uint32_t b, uint32_t c) {
    if (a <= b && a <= c) return {a, b, c};
    if (b <= c) return {b, c, a};
    return {c, a, b};
}
struct Tri3Hash {
    size_t operator()(const std::array<uint32_t, 3>& t) const {
        size_t h = 1469598103934665603ull;
        for (uint32_t x : t) h = (h ^ x) * 1099511628211ull;
        return h;
    }
};
} // namespace detail

// Check the orientation consistency of an in-memory tet list.
inline OrientationResult check_tet_orientation(const std::vector<std::array<uint32_t, 4>>& tets) {
    OrientationResult r;
    r.tets = tets.size();
    // The four outward-oriented faces of a positively-oriented tet (v0,v1,v2,v3).
    static const int F[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    std::unordered_map<std::array<uint32_t, 3>, uint32_t, detail::Tri3Hash> oriented; // winding key -> count
    std::unordered_map<std::array<uint32_t, 3>, uint32_t, detail::Tri3Hash> incident; // {i,j,k} -> #tets
    oriented.reserve(tets.size() * 4);
    incident.reserve(tets.size() * 4);
    for (const auto& t : tets)
        for (const auto& f : F) {
            const uint32_t a = t[f[0]], b = t[f[1]], c = t[f[2]];
            oriented[detail::rot_min(a, b, c)]++;
            std::array<uint32_t, 3> u = {a, b, c};
            std::sort(u.begin(), u.end());
            incident[u]++;
        }
    for (const auto& kv : oriented)
        if (kv.second != 1) r.same_winding++; // an oriented face emitted >1x == same winding on both sides
    for (const auto& kv : incident) {
        if (kv.second == 1) r.boundary++;
        else if (kv.second == 2) r.internal++;
        else r.over_shared++;
    }
    r.ok = (r.same_winding == 0 && r.over_shared == 0);
    return r;
}

// Read a .tet file ("<N> vertices" / "<M> tets" headers, N coordinate lines, then M
// "4 v0 v1 v2 v3" lines) and check its orientation.
inline OrientationResult check_tet_orientation_file(const std::string& tet_path) {
    OrientationResult r;
    std::ifstream in(tet_path);
    if (!in) {
        r.error = "cannot open " + tet_path;
        return r;
    }
    auto next = [&](std::string& line) {
        while (std::getline(in, line)) {
            const size_t p = line.find_first_not_of(" \t\r\n");
            if (p != std::string::npos && line[p] != '#') return true;
        }
        return false;
    };
    std::string line;
    uint64_t nv = 0, nt = 0;
    if (!next(line)) { r.error = "bad .tet header"; return r; }
    std::istringstream(line) >> nv;
    if (!next(line)) { r.error = "bad .tet header"; return r; }
    std::istringstream(line) >> nt;
    for (uint64_t i = 0; i < nv; i++)
        if (!next(line)) { r.error = "truncated .tet vertices"; return r; }
    std::vector<std::array<uint32_t, 4>> tets;
    tets.reserve(nt);
    for (uint64_t i = 0; i < nt; i++) {
        if (!next(line)) { r.error = "truncated .tet tets"; return r; }
        uint32_t k, a, b, c, d;
        std::istringstream(line) >> k >> a >> b >> c >> d;
        tets.push_back({a, b, c, d});
    }
    return check_tet_orientation(tets);
}

} // namespace vrtest
