// Independent verifier for the 2D pipeline's provenance, mirroring verify_tracking.cpp.
//
// Like its 3D counterpart it shares NO code with the mesh generator (only the bigrational type,
// for exact arithmetic) and re-derives everything from geometry. The checks are the 2D forms of
// 3D's Check 3 (edge provenance), Check 4 (point provenance) and Check 5 (orientation).
//
//   usage: verify_tracking_2d <input.obj> <triangulation.off>
//
// Reads the sidecars <triangulation.off>.rational / .segmentprov / .pointprov, which the
// generator writes with -2d -r. Exits 0 and prints "TRACKING OK" on success.

#include <VolumeRemesher/numerics.h>

#include "tri_orientation.h"

#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using vol_rem::bigrational;

namespace {

struct R2
{
    bigrational x, y;
};

bool read_obj(const std::string& path, std::vector<R2>& pts,
              std::vector<std::array<uint32_t, 2>>& segs)
{
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    char buf[8192];
    std::vector<long> ids;
    while (fgets(buf, sizeof(buf), f)) {
        const char* p = buf;
        while (*p == ' ' || *p == '\t') p++;
        if (p[0] == 'v' && (p[1] == ' ' || p[1] == '\t')) {
            double x = 0, y = 0;
            if (sscanf(p + 1, "%lf %lf", &x, &y) < 2) {
                fclose(f);
                return false;
            }
            pts.push_back({bigrational(x), bigrational(y)});
        } else if (p[0] == 'l' && (p[1] == ' ' || p[1] == '\t')) {
            ids.clear();
            const char* q = p + 1;
            while (*q) {
                while (*q == ' ' || *q == '\t') q++;
                if (!*q || *q == '\n' || *q == '\r') break;
                char* e = nullptr;
                const long v = strtol(q, &e, 10);
                if (e == q) break;
                ids.push_back(v);
                q = e;
                while (*q && *q != ' ' && *q != '\t' && *q != '\n' && *q != '\r') q++;
            }
            const long nv = long(pts.size());
            for (size_t k = 0; k + 1 < ids.size(); k++) {
                long a = ids[k], b = ids[k + 1];
                if (a < 0) a += nv + 1;
                if (b < 0) b += nv + 1;
                if (a < 1 || b < 1 || a > nv || b > nv) {
                    fclose(f);
                    return false;
                }
                segs.push_back({uint32_t(a - 1), uint32_t(b - 1)});
            }
        }
    }
    fclose(f);
    return true;
}

// Parses "[-]num[/den]" as written by rational_to_string. bigrational has no string constructor
// that is common to both backends, so build it digit by digit -- same approach as parse_rat in
// verify_tracking.cpp.
bigrational parse_rational(const std::string& s)
{
    size_t i = 0;
    bool neg = false;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) {
        neg = (s[i] == '-');
        i++;
    }
    bigrational num(0.0), ten(10.0);
    for (; i < s.size() && std::isdigit((unsigned char)s[i]); i++)
        num = num * ten + bigrational((double)(s[i] - '0'));
    bigrational den(1.0);
    if (i < s.size() && s[i] == '/') {
        i++;
        den = bigrational(0.0);
        for (; i < s.size() && std::isdigit((unsigned char)s[i]); i++)
            den = den * ten + bigrational((double)(s[i] - '0'));
    }
    bigrational r = num / den;
    // Split rather than `neg ? -r : r`: with gmpxx, `-r` is a lazy __gmp_expr and MSVC finds the
    // ternary's common type ambiguous. Same workaround as verify_tracking.cpp.
    if (neg) return -r;
    return r;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s <input.obj> <triangulation.off>\n", argv[0]);
        return 1;
    }
    const std::string obj = argv[1];
    const std::string off = argv[2];

    std::vector<R2> in_pts;
    std::vector<std::array<uint32_t, 2>> in_segs;
    if (!read_obj(obj, in_pts, in_segs)) {
        fprintf(stderr, "cannot read %s\n", obj.c_str());
        return 1;
    }

    // --- output triangles (from the OFF) and exact vertex coordinates (from .rational) --------
    std::vector<std::array<uint32_t, 3>> tris;
    {
        std::ifstream f(off);
        if (!f) {
            fprintf(stderr, "cannot read %s\n", off.c_str());
            return 1;
        }
        std::string magic;
        long nv = 0, nf = 0, ne = 0;
        f >> magic >> nv >> nf >> ne;
        double d;
        for (long i = 0; i < 3 * nv; i++) f >> d;
        for (long i = 0; i < nf; i++) {
            long n, a, b, c;
            f >> n >> a >> b >> c;
            if (n != 3) {
                fprintf(stderr, "face %ld is not a triangle\n", i);
                return 1;
            }
            tris.push_back({uint32_t(a), uint32_t(b), uint32_t(c)});
        }
    }

    std::vector<R2> out_pts;
    {
        std::ifstream f(off + ".rational");
        if (!f) {
            fprintf(stderr, "cannot read %s.rational (run with -r)\n", off.c_str());
            return 1;
        }
        size_t n = 0;
        f >> n;
        out_pts.resize(n);
        for (size_t i = 0; i < n; i++) {
            std::string sx, sy;
            f >> sx >> sy;
            out_pts[i] = {parse_rational(sx), parse_rational(sy)};
        }
    }

    // --- Check 1: segment provenance ----------------------------------------------------------
    // The 2D form of verify_tracking.cpp's Check 3. For input segment [A,B], every output edge
    // endpoint must be EXACTLY collinear with AB, its parameter t must lie in [0,1], and the
    // sub-intervals must form an exact partition of [0,1]: both endpoints reached, no gap, no
    // overlap. All in exact rational arithmetic; no tolerance anywhere.
    size_t seg_fail = 0;
    {
        std::ifstream f(off + ".segmentprov");
        if (!f) {
            fprintf(stderr, "cannot read %s.segmentprov\n", off.c_str());
            return 1;
        }
        size_t n = 0;
        f >> n;
        if (n != in_segs.size()) {
            fprintf(stderr, "segmentprov has %zu entries, input has %zu segments\n", n,
                    in_segs.size());
            return 1;
        }
        for (size_t i = 0; i < n; i++) {
            size_t id = 0, m = 0;
            f >> id >> m;
            std::vector<std::array<uint32_t, 3>> edges(m);
            for (size_t k = 0; k < m; k++) f >> edges[k][0] >> edges[k][1] >> edges[k][2];

            const R2& A = in_pts[in_segs[i][0]];
            const R2& B = in_pts[in_segs[i][1]];
            const bigrational dx = B.x - A.x, dy = B.y - A.y;
            const bigrational len2 = dx * dx + dy * dy;
            if (len2.sgn() == 0) continue; // degenerate input segment, legitimately has no edges

            const char* why = nullptr;
            std::vector<std::array<bigrational, 2>> iv;
            if (edges.empty()) why = "no output edges";
            for (const auto& e : edges) {
                bigrational tt[2];
                bool ok = true;
                for (int k = 0; k < 2; k++) {
                    const R2& P = out_pts[e[1 + k]];
                    const bigrational qx = P.x - A.x, qy = P.y - A.y;
                    if ((qx * dy - qy * dx).sgn() != 0) {
                        why = "output edge not collinear with segment";
                        ok = false;
                        break;
                    }
                    tt[k] = (qx * dx + qy * dy) / len2;
                }
                if (!ok) break;
                if (tt[0].sgn() < 0 || tt[1].sgn() < 0 || tt[0] > bigrational(1.0) ||
                    tt[1] > bigrational(1.0) || !(tt[0] < tt[1])) {
                    why = "output edge outside segment or degenerate";
                    break;
                }
                iv.push_back({tt[0], tt[1]});
            }
            if (!why) {
                std::sort(iv.begin(), iv.end(),
                          [](const std::array<bigrational, 2>& a,
                             const std::array<bigrational, 2>& b) { return a[0] < b[0]; });
                if (iv.front()[0].sgn() != 0 || !(iv.back()[1] == bigrational(1.0)))
                    why = "output edges do not reach both endpoints";
                else
                    for (size_t k = 1; k < iv.size() && !why; k++)
                        if (!(iv[k][0] == iv[k - 1][1])) why = "gap or overlap between output edges";
            }
            if (why) {
                if (seg_fail < 10 || getenv("VT_DEBUG"))
                    fprintf(stderr, "segment %zu: %s\n", i, why);
                seg_fail++;
            }
        }
    }

    // --- Check 2: point provenance (3D's Check 4) ---------------------------------------------
    size_t pt_fail = 0;
    {
        std::ifstream f(off + ".pointprov");
        if (!f) {
            fprintf(stderr, "cannot read %s.pointprov\n", off.c_str());
            return 1;
        }
        size_t n = 0;
        f >> n;
        for (size_t i = 0; i < n; i++) {
            long id = 0, t = 0, v = 0;
            f >> id >> t >> v;
            if (v < 0) continue; // did not survive
            if (size_t(v) >= out_pts.size() || !(out_pts[size_t(v)].x == in_pts[i].x) ||
                !(out_pts[size_t(v)].y == in_pts[i].y)) {
                if (pt_fail < 10 || getenv("VT_DEBUG"))
                    fprintf(stderr, "point %zu: output vertex does not match exactly\n", i);
                pt_fail++;
            }
        }
    }

    // --- Check 3: orientation / manifoldness (3D's Check 5) -----------------------------------
    const vrtest::TriOrientationResult orient = vrtest::check_tri_orientation(tris);

    printf("triangles     %zu\n", tris.size());
    printf("vertices      %zu\n", out_pts.size());
    printf("segments      %zu  (failures %zu)\n", in_segs.size(), seg_fail);
    printf("points        %zu  (failures %zu)\n", in_pts.size(), pt_fail);
    printf("orient        ok=%d boundary=%llu internal=%llu over_shared=%llu same_winding=%llu\n",
           int(orient.ok), (unsigned long long)orient.boundary, (unsigned long long)orient.internal,
           (unsigned long long)orient.over_shared, (unsigned long long)orient.same_winding);

    const bool ok = (seg_fail == 0) && (pt_fail == 0) && orient.ok;
    printf("%s\n", ok ? "TRACKING OK" : "TRACKING FAILED");
    return ok ? 0 : 1;
}
