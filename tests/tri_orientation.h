// Orientation / manifoldness check for a 2D triangle mesh.
//
// This is the 2D counterpart of tet_orientation.h, and like it the check is purely
// COMBINATORIAL: it uses only the triangle vertex indices, no coordinates. That makes it
// independent of the exact-arithmetic backend and cheap enough to run on every output.
//
// For each triangle it emits the three directed edges (v0,v1), (v1,v2), (v2,v0).
//
//   - A directed edge must appear AT MOST ONCE. Two triangles giving the same edge the same
//     direction sit on the same side of it, i.e. they overlap. In 3D this is what caught the
//     mis-placed cell barycentre; the 2D failure mode is identical.
//   - An undirected edge must be incident to 1 triangle (boundary) or 2 (internal). More than 2
//     is non-manifold.

#ifndef VRTEST_TRI_ORIENTATION_H
#define VRTEST_TRI_ORIENTATION_H

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace vrtest {

struct TriOrientationResult
{
    bool ok = false;
    uint64_t tris = 0;
    uint64_t boundary = 0;     // undirected edges used by exactly 1 triangle
    uint64_t internal = 0;     // undirected edges used by exactly 2 triangles
    uint64_t over_shared = 0;  // undirected edges used by more than 2 (non-manifold)
    uint64_t same_winding = 0; // directed edges emitted more than once (overlapping triangles)
    std::string error;         // non-empty if the mesh could not be read
};

namespace detail {

inline uint64_t dir_key(uint32_t a, uint32_t b) { return (uint64_t(a) << 32) | uint64_t(b); }
inline uint64_t und_key(uint32_t a, uint32_t b)
{
    return (a < b) ? dir_key(a, b) : dir_key(b, a);
}

} // namespace detail

inline TriOrientationResult check_tri_orientation(const std::vector<std::array<uint32_t, 3>>& tris)
{
    TriOrientationResult r;
    r.tris = tris.size();

    std::unordered_map<uint64_t, uint32_t> directed;
    std::unordered_map<uint64_t, uint32_t> incident;
    directed.reserve(tris.size() * 3);
    incident.reserve(tris.size() * 3);

    for (const auto& t : tris) {
        for (int k = 0; k < 3; k++) {
            const uint32_t a = t[k];
            const uint32_t b = t[(k + 1) % 3];
            directed[detail::dir_key(a, b)]++;
            incident[detail::und_key(a, b)]++;
        }
    }

    for (const auto& kv : directed)
        if (kv.second != 1) r.same_winding += kv.second - 1;

    for (const auto& kv : incident) {
        if (kv.second == 1) r.boundary++;
        else if (kv.second == 2) r.internal++;
        else r.over_shared++;
    }

    r.ok = (r.same_winding == 0 && r.over_shared == 0);
    return r;
}

// Reads an OFF file (as written by the 2D pipeline) and checks it. Tolerant of blank lines and
// '#' comments, and requires every face to be a triangle.
inline TriOrientationResult check_tri_orientation_file(const std::string& off_path)
{
    TriOrientationResult r;
    FILE* f = fopen(off_path.c_str(), "r");
    if (!f) {
        r.error = "cannot open " + off_path;
        return r;
    }

    char line[4096];
    const auto next = [&]() -> bool {
        while (fgets(line, sizeof(line), f)) {
            char* p = line;
            while (*p == ' ' || *p == '\t') p++;
            if (*p && *p != '\n' && *p != '\r' && *p != '#') return true;
        }
        return false;
    };

    if (!next()) {
        r.error = "empty file";
        fclose(f);
        return r;
    }
    // The magic line may be "OFF" alone or "OFF" followed by the counts.
    long nv = 0, nf = 0, ne = 0;
    if (sscanf(line, "OFF %ld %ld %ld", &nv, &nf, &ne) != 3) {
        if (!next() || sscanf(line, "%ld %ld %ld", &nv, &nf, &ne) != 3) {
            r.error = "bad OFF header";
            fclose(f);
            return r;
        }
    }

    for (long i = 0; i < nv; i++) {
        if (!next()) {
            r.error = "truncated vertex list";
            fclose(f);
            return r;
        }
    }

    std::vector<std::array<uint32_t, 3>> tris;
    tris.reserve(size_t(nf));
    for (long i = 0; i < nf; i++) {
        if (!next()) {
            r.error = "truncated face list";
            fclose(f);
            return r;
        }
        unsigned n = 0, a = 0, b = 0, c = 0;
        if (sscanf(line, "%u %u %u %u", &n, &a, &b, &c) != 4 || n != 3) {
            r.error = "face " + std::to_string(i) + " is not a triangle";
            fclose(f);
            return r;
        }
        tris.push_back({a, b, c});
    }
    fclose(f);

    return check_tri_orientation(tris);
}

} // namespace vrtest

#endif
