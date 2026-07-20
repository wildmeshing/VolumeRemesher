// Unit tests for the 2D arrangement pipeline.
//
// Unlike integration_tests.cpp (which shells out to the mesh_generator binary), this target links
// the library and exercises the internal 2D APIs directly.

#include <catch2/catch_test_macros.hpp>

#include "tri_orientation.h"

#include <VolumeRemesher/2d/arrangement2d.h>
#include <VolumeRemesher/2d/delaunay2d.h>
#include <VolumeRemesher/2d/predicates2d.h>
#include <VolumeRemesher/numerics.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <unordered_set>
#include <fstream>
#include <string>
#include <vector>

using vol_rem::bigrational;
using vol_rem::explicitPoint2D;
using vol_rem::genericPoint;
using vol_rem::implicitPoint2D_SSI;
using namespace vol_rem::vr2d;

namespace {

// Deterministic RNG. Deliberately NOT std::uniform_real_distribution, whose output is
// implementation-defined -- the existing integration suite avoids it for the same reason.
struct Rnd
{
    uint64_t s;
    explicit Rnd(uint64_t seed)
        : s(seed ? seed : 0x9E3779B97F4A7C15ull)
    {}
    uint64_t next()
    {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return s;
    }
    // Small integers keep the exact predicates in their cheap filter path while still producing
    // plenty of exactly-degenerate configurations (collinear, cocircular).
    double coord(int range) { return double(int64_t(next() % uint64_t(2 * range)) - range); }
    // Uniform in [0,1) with full double precision: almost surely general position.
    double unit() { return double(next() >> 11) * (1.0 / 9007199254740992.0); }
};

// --- exact reference determinants, computed independently in bigrational -------------------

bigrational det3(
    const bigrational& a,
    const bigrational& b,
    const bigrational& c,
    const bigrational& d,
    const bigrational& e,
    const bigrational& f,
    const bigrational& g,
    const bigrational& h,
    const bigrational& i)
{
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

// The lifted 4x4 determinant
//     | x0 y0 w0 1 ; x1 y1 w1 1 ; x2 y2 w2 1 ; x3 y3 w3 1 |,  w = x^2 + y^2
// expanded along the last column. in_circle_2d_SOS's symbolic-perturbation cofactors are derived
// from exactly this matrix, so this reference is what pins the derivation.
bigrational det4_lifted(const double* p0, const double* p1, const double* p2, const double* p3)
{
    const double* p[4] = {p0, p1, p2, p3};
    bigrational x[4], y[4], w[4];
    for (int k = 0; k < 4; k++) {
        x[k] = bigrational(p[k][0]);
        y[k] = bigrational(p[k][1]);
        w[k] = x[k] * x[k] + y[k] * y[k];
    }
    bigrational acc(0.0);
    for (int i = 0; i < 4; i++) {
        int r[3], n = 0;
        for (int k = 0; k < 4; k++)
            if (k != i) r[n++] = k;
        const bigrational m = det3(
            x[r[0]], y[r[0]], w[r[0]],
            x[r[1]], y[r[1]], w[r[1]],
            x[r[2]], y[r[2]], w[r[2]]);
        // cofactor sign of entry (i, 3)
        acc = ((i + 3) % 2 == 0) ? (acc + m) : (acc - m);
    }
    return acc;
}

int sgn_rat(const bigrational& r) { return r.sgn(); }

} // namespace

TEST_CASE("2d predicates: orient_2d is positive for counterclockwise", "[2d][predicates]")
{
    const double a[2] = {0.0, 0.0};
    const double b[2] = {1.0, 0.0};
    const double c[2] = {0.0, 1.0};

    CHECK(PCK::orient_2d(a, b, c) == POSITIVE);
    CHECK(PCK::orient_2d(a, c, b) == NEGATIVE);

    const double collinear[2] = {2.0, 0.0};
    CHECK(PCK::orient_2d(a, b, collinear) == ZERO);
}

TEST_CASE("2d predicates: points_are_identical_2d", "[2d][predicates]")
{
    const double a[2] = {1.5, -2.5};
    const double b[2] = {1.5, -2.5};
    const double c[2] = {1.5, -2.4};
    CHECK(PCK::points_are_identical_2d(a, b));
    CHECK_FALSE(PCK::points_are_identical_2d(a, c));
}

// The SOS cofactors C_i are read off the lifted 4x4 determinant assuming incircle() computes that
// determinant (not its negation, and not some other scaling). If that ever stopped holding, the
// symbolic perturbation would silently invert and the Delaunay would tear only on cocircular
// input. This test is the guard.
TEST_CASE("2d predicates: incircle equals the lifted 4x4 determinant", "[2d][predicates]")
{
    // A hand-checked case first: right triangle CCW, query point inside the circumcircle.
    {
        const double a[2] = {0.0, 0.0};
        const double b[2] = {1.0, 0.0};
        const double c[2] = {0.0, 1.0};
        const double inside[2] = {0.25, 0.25};
        const double outside[2] = {5.0, 5.0};

        REQUIRE(PCK::orient_2d(a, b, c) == POSITIVE);
        CHECK(vol_rem::incircle(a[0], a[1], b[0], b[1], c[0], c[1], inside[0], inside[1]) > 0);
        CHECK(vol_rem::incircle(a[0], a[1], b[0], b[1], c[0], c[1], outside[0], outside[1]) < 0);
        CHECK(sgn_rat(det4_lifted(a, b, c, inside)) > 0);
        CHECK(sgn_rat(det4_lifted(a, b, c, outside)) < 0);
    }

    // Then agreement on many random configurations, including degenerate ones.
    Rnd rnd(0xC0FFEEu);
    for (int it = 0; it < 4000; it++) {
        double p[4][2];
        for (int k = 0; k < 4; k++) {
            p[k][0] = rnd.coord(6);
            p[k][1] = rnd.coord(6);
        }
        const int got =
            vol_rem::incircle(p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);
        const int want = sgn_rat(det4_lifted(p[0], p[1], p[2], p[3]));
        REQUIRE(got == want);
    }
}

TEST_CASE("2d predicates: in_circle_2d_SOS agrees with incircle when non-degenerate",
          "[2d][predicates]")
{
    Rnd rnd(0xBEEF01u);
    int degenerate = 0;
    for (int it = 0; it < 4000; it++) {
        double p[4][2];
        for (int k = 0; k < 4; k++) {
            p[k][0] = rnd.coord(6);
            p[k][1] = rnd.coord(6);
        }
        // precondition: (p0,p1,p2) must be a real triangle
        if (PCK::orient_2d(p[0], p[1], p[2]) == ZERO) continue;

        const int plain =
            vol_rem::incircle(p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);
        const Sign sos = PCK::in_circle_2d_SOS(p[0], p[1], p[2], p[3]);

        REQUIRE(sos != ZERO); // SOS must never be undecided
        if (plain != 0) {
            REQUIRE(sos == to_sign(plain));
        } else {
            degenerate++;
        }
    }
    // The small-integer coordinates are chosen so cocircular quadruples actually occur; if this
    // ever hits zero the test below is the only thing exercising the perturbation path.
    INFO("cocircular quadruples encountered: " << degenerate);
    CHECK(degenerate > 0);
}

// REGRESSION TEST for a bug in the generated indirect predicates (indirect_predicates.h).
//
// genericPoint::dotProductSign2D returned the WRONG SIGN for roughly a third of inputs in the two
// configurations that route to dotProductSign2D_IEI_t -- (explicit, implicit, implicit) and
// (implicit, explicit, implicit) -- because the homogeneous scaling of the first term carried a
// spurious extra factor of the first point's denominator. That is correct only when the first
// argument is explicit, which is never the case on those two paths.
//
// It went unnoticed because nothing in the 3D pipeline calls the 2D dot-product predicate. The
// upstream project (github.com/MarcoAttene/Indirect_Predicates, HEAD 2026-05-21) still carries
// the same defect, so this cannot be fixed by updating the vendored headers.
//
// The 2D segment walk relies on this predicate to decide whether a collinear neighbour lies
// forward or backward along the segment, so the bug made the walk fail outright.
TEST_CASE("2d predicates: dotProductSign2D agrees with exact arithmetic", "[2d][predicates]")
{
    Rnd rnd(0xD07u);

    // Address-stable storage: implicitPoint2D_SSI holds references to its parents.
    std::deque<explicitPoint2D> E;
    std::deque<implicitPoint2D_SSI> I;
    for (int i = 0; i < 400; i++) E.emplace_back(rnd.coord(20), rnd.coord(20));
    for (size_t i = 0; i + 3 < E.size(); i += 4) {
        I.emplace_back(E[i], E[i + 1], E[i + 2], E[i + 3]);
        double x, y;
        if (!I.back().getApproxXYCoordinates(x, y)) I.pop_back(); // parallel lines
    }
    REQUIRE(I.size() > 20);

    const auto exact = [](const genericPoint& a, const genericPoint& b, const genericPoint& c) {
        bigrational ax, ay, bx, by, cx, cy;
        REQUIRE(a.getExactXYCoordinates(ax, ay));
        REQUIRE(b.getExactXYCoordinates(bx, by));
        REQUIRE(c.getExactXYCoordinates(cx, cy));
        return ((ax - cx) * (bx - cx) + (ay - cy) * (by - cy)).sgn();
    };

    // Every combination of explicit/implicit in each of the three argument slots.
    int checked = 0;
    for (int mask = 0; mask < 8; mask++) {
        int per_config = 0;
        for (int t = 0; t < 600; t++) {
            const auto pick = [&](int implicit) -> const genericPoint* {
                return implicit ? (const genericPoint*)&I[rnd.next() % I.size()]
                                : (const genericPoint*)&E[rnd.next() % E.size()];
            };
            const genericPoint* a = pick(mask & 4);
            const genericPoint* b = pick(mask & 2);
            const genericPoint* c = pick(mask & 1);
            if (a == b || a == c || b == c) continue;
            INFO("config mask " << mask);
            REQUIRE(genericPoint::dotProductSign2D(*a, *b, *c) == exact(*a, *b, *c));
            per_config++;
            checked++;
        }
        REQUIRE(per_config > 100); // each configuration actually exercised
    }
    INFO("comparisons: " << checked);
    CHECK(checked > 4000);
}

// The load-bearing property of the symbolic perturbation: because eps is attached to the POINT
// (via lexicographic rank) and not to the argument position, the predicate must behave exactly
// like the determinant it perturbs -- antisymmetric under transposition of any two arguments.
//
// Bowyer-Watson depends on this directly. Two triangles (a,b,c) and (b,a,d) sharing edge ab must
// agree on whether the shared edge is Delaunay; that is the even permutation
// (a,b,c,d) -> (b,a,d,c), which must preserve the sign. If it did not, the two triangles would
// disagree and the cavity would stop being star-shaped.
TEST_CASE("2d predicates: in_circle_2d_SOS is antisymmetric under argument permutation",
          "[2d][predicates]")
{
    // Four exactly cocircular point sets (integer coordinates, so cocircularity is exact).
    const std::vector<std::array<std::array<double, 2>, 4>> sets = {
        {{{{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}}},          // unit square
        {{{{5.0, 0.0}}, {{3.0, 4.0}}, {{-5.0, 0.0}}, {{0.0, -5.0}}}},        // radius-5 circle
        {{{{-1.0, 0.0}}, {{0.0, 1.0}}, {{1.0, 0.0}}, {{0.0, -1.0}}}},        // unit circle
        {{{{0.0, 0.0}}, {{4.0, 0.0}}, {{4.0, 2.0}}, {{0.0, 2.0}}}},          // rectangle
    };

    const int perm[24][4] = {
        {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
        {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
        {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
        {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0},
    };

    for (const auto& S : sets) {
        // Confirm the set really is cocircular, otherwise the test proves nothing.
        REQUIRE(sgn_rat(det4_lifted(S[0].data(), S[1].data(), S[2].data(), S[3].data())) == 0);

        const Sign base = PCK::in_circle_2d_SOS(S[0].data(), S[1].data(), S[2].data(), S[3].data());
        REQUIRE(base != ZERO);

        for (const auto& q : perm) {
            // parity of the permutation, by counting inversions
            int inv = 0;
            for (int i = 0; i < 4; i++)
                for (int j = i + 1; j < 4; j++)
                    if (q[i] > q[j]) inv++;
            const bool even = (inv % 2) == 0;

            const double* a = S[q[0]].data();
            const double* b = S[q[1]].data();
            const double* c = S[q[2]].data();
            const double* d = S[q[3]].data();
            if (PCK::orient_2d(a, b, c) == ZERO) continue; // precondition violated, skip

            const Sign got = PCK::in_circle_2d_SOS(a, b, c, d);
            REQUIRE(got != ZERO);
            INFO("permutation " << q[0] << q[1] << q[2] << q[3] << (even ? " even" : " odd"));
            REQUIRE(got == (even ? base : Sign(-base)));
        }
    }
}

// ==============================================================================================
// Delaunay2d
// ==============================================================================================

namespace {

using vol_rem::vr2d::Delaunay2d;
using vol_rem::vr2d::index_t;

struct DelaunayReport
{
    bool combinatorics_ok = false;
    bool empty_circle_ok = false;
    bool euler_ok = false;
    bool all_vertices_used = false;
    index_t nb_finite = 0;
    index_t nb_virtual = 0;
    std::string detail;

    bool ok() const
    {
        return combinatorics_ok && empty_circle_ok && euler_ok && all_vertices_used;
    }
};

// Full exact validation of a computed triangulation. This is the "self-check" half of test (a):
// it is strictly stronger than matching geogram, because for points in general position the
// Delaunay triangulation is UNIQUE -- so anything passing these checks IS the triangulation
// geogram would produce, whatever route it took to get there.
DelaunayReport validate_delaunay(const Delaunay2d& D, index_t nb_pts)
{
    DelaunayReport r;
    r.combinatorics_ok = D.check_combinatorics();
    if (!r.combinatorics_ok) r.detail = "combinatorics";

    const index_t nt = D.nb_triangles();
    r.nb_finite = D.nb_finite_triangles();
    r.nb_virtual = nt - r.nb_finite;

    // --- exact empty-circumcircle test ---
    // A triangulation is globally Delaunay iff it is LOCALLY Delaunay across every internal edge:
    // for adjacent finite triangles t and t2, the vertex of t2 opposite the shared edge must not
    // lie strictly inside t's circumcircle. Checked with the exact predicate, no tolerance.
    r.empty_circle_ok = true;
    for (index_t t = 0; t < nt && r.empty_circle_ok; t++) {
        if (!D.triangle_is_finite(t)) continue;
        const double* a = D.vertex_ptr(D.triangle_vertex(t, 0));
        const double* b = D.vertex_ptr(D.triangle_vertex(t, 1));
        const double* c = D.vertex_ptr(D.triangle_vertex(t, 2));
        for (index_t le = 0; le < 3; le++) {
            const index_t t2 = D.triangle_adjacent(t, le);
            if (t2 == Delaunay2d::NO_INDEX || !D.triangle_is_finite(t2)) continue;
            // vertex of t2 not shared with t
            index_t opp = Delaunay2d::NO_INDEX;
            for (index_t k = 0; k < 3; k++) {
                const index_t v = D.triangle_vertex(t2, k);
                if (v != D.triangle_vertex(t, 0) && v != D.triangle_vertex(t, 1) &&
                    v != D.triangle_vertex(t, 2)) {
                    opp = v;
                    break;
                }
            }
            if (opp == Delaunay2d::NO_INDEX) continue;
            if (vol_rem::incircle(a[0], a[1], b[0], b[1], c[0], c[1],
                                  D.vertex_ptr(opp)[0], D.vertex_ptr(opp)[1]) > 0) {
                r.empty_circle_ok = false;
                r.detail = "not locally Delaunay at triangle " + std::to_string(t);
                break;
            }
        }
    }

    // --- Euler ---
    // A triangulation of n points with h of them on the convex hull has exactly 2n - 2 - h
    // triangles, and there is one virtual triangle per hull edge, i.e. per hull vertex.
    const index_t h = r.nb_virtual;
    r.euler_ok = (h >= 3) && (r.nb_finite + 2 + h == 2 * nb_pts);
    if (!r.euler_ok) {
        r.detail = "euler: finite=" + std::to_string(r.nb_finite) + " hull=" + std::to_string(h) +
                   " n=" + std::to_string(nb_pts);
    }

    // --- every input point is a vertex of some triangle ---
    std::vector<bool> used(nb_pts, false);
    for (index_t t = 0; t < nt; t++) {
        for (index_t k = 0; k < 3; k++) {
            const index_t v = D.triangle_vertex(t, k);
            if (v != Delaunay2d::NO_INDEX) used[v] = true;
        }
    }
    r.all_vertices_used = std::all_of(used.begin(), used.end(), [](bool b) { return b; });
    if (!r.all_vertices_used) r.detail = "some input vertex is in no triangle";

    return r;
}

} // namespace

TEST_CASE("2d delaunay: hilbert_order is a permutation", "[2d][delaunay]")
{
    Rnd rnd(0x51D5077u);
    for (const index_t n : {1u, 2u, 3u, 17u, 500u}) {
        std::vector<double> pts(2 * n);
        for (auto& c : pts) c = rnd.unit();
        std::vector<index_t> order;
        Delaunay2d::hilbert_order(n, pts.data(), order);
        REQUIRE(order.size() == n);
        std::vector<index_t> sorted = order;
        std::sort(sorted.begin(), sorted.end());
        for (index_t i = 0; i < n; i++) REQUIRE(sorted[i] == i);
    }

    // Degenerate bounding boxes must not divide by zero or collapse the order.
    {
        const index_t n = 8;
        std::vector<double> pts(2 * n);
        for (index_t i = 0; i < n; i++) {
            pts[2 * i] = 3.0; // all points on a vertical line
            pts[2 * i + 1] = double(i);
        }
        std::vector<index_t> order;
        Delaunay2d::hilbert_order(n, pts.data(), order);
        std::vector<index_t> sorted = order;
        std::sort(sorted.begin(), sorted.end());
        for (index_t i = 0; i < n; i++) REQUIRE(sorted[i] == i);
    }
}

TEST_CASE("2d delaunay: degenerate inputs", "[2d][delaunay]")
{
    SECTION("fewer than three points")
    {
        const double pts[4] = {0.0, 0.0, 1.0, 1.0};
        Delaunay2d D;
        CHECK_FALSE(D.set_vertices(2, pts));
    }

    SECTION("all points collinear")
    {
        std::vector<double> pts;
        for (int i = 0; i < 10; i++) {
            pts.push_back(double(i));
            pts.push_back(2.0 * double(i)); // exactly on y = 2x
        }
        Delaunay2d D;
        CHECK_FALSE(D.set_vertices(10, pts.data()));
    }

    SECTION("minimal triangle")
    {
        const double pts[6] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
        Delaunay2d D;
        REQUIRE(D.set_vertices(3, pts));
        const DelaunayReport r = validate_delaunay(D, 3);
        INFO(r.detail);
        CHECK(r.ok());
        CHECK(r.nb_finite == 1);
        CHECK(r.nb_virtual == 3);
    }
}

TEST_CASE("2d delaunay: random point sets are exactly Delaunay", "[2d][delaunay]")
{
    Rnd rnd(0xD3134Eu);
    for (const index_t n : {4u, 5u, 10u, 50u, 200u, 1000u}) {
        std::vector<double> pts(2 * n);
        for (auto& c : pts) c = rnd.unit();

        Delaunay2d D;
        REQUIRE(D.set_vertices(n, pts.data()));
        const DelaunayReport r = validate_delaunay(D, n);
        INFO("n=" << n << " " << r.detail);
        CHECK(r.combinatorics_ok);
        CHECK(r.empty_circle_ok);
        CHECK(r.euler_ok);
        CHECK(r.all_vertices_used);
    }
}

// Cocircular and collinear configurations are where the SOS layer earns its keep: a plain exact
// incircle returns 0 here and Bowyer-Watson would have no consistent answer.
TEST_CASE("2d delaunay: heavily degenerate point sets", "[2d][delaunay]")
{
    SECTION("regular grid (every unit cell is a cocircular quadruple)")
    {
        const index_t k = 12;
        std::vector<double> pts;
        for (index_t i = 0; i < k; i++)
            for (index_t j = 0; j < k; j++) {
                pts.push_back(double(i));
                pts.push_back(double(j));
            }
        const index_t n = k * k;
        Delaunay2d D;
        REQUIRE(D.set_vertices(n, pts.data()));
        const DelaunayReport r = validate_delaunay(D, n);
        INFO(r.detail);
        CHECK(r.combinatorics_ok);
        CHECK(r.empty_circle_ok);
        CHECK(r.euler_ok);
        CHECK(r.all_vertices_used);
    }

    SECTION("points exactly on a circle (all n cocircular)")
    {
        // Pythagorean points on the radius-1105 circle, plus multiples, all exact integers.
        const int pyth[][2] = {{1105, 0},  {1104, 47},  {1073, 264}, {1052, 340}, {1020, 425},
                               {943, 576}, {884, 663},  {816, 745},  {780, 782},  {744, 817}};
        std::vector<double> pts;
        for (const auto& p : pyth) {
            for (const int sx : {1, -1})
                for (const int sy : {1, -1}) {
                    pts.push_back(double(sx * p[0]));
                    pts.push_back(double(sy * p[1]));
                }
        }
        // dedup exact duplicates (points with a zero coordinate appear twice)
        std::vector<std::array<double, 2>> uniq;
        for (size_t i = 0; i + 1 < pts.size(); i += 2) {
            const std::array<double, 2> q = {pts[i], pts[i + 1]};
            if (std::find(uniq.begin(), uniq.end(), q) == uniq.end()) uniq.push_back(q);
        }
        std::vector<double> flat;
        for (const auto& q : uniq) {
            flat.push_back(q[0]);
            flat.push_back(q[1]);
        }
        const index_t n = index_t(uniq.size());

        Delaunay2d D;
        REQUIRE(D.set_vertices(n, flat.data()));
        const DelaunayReport r = validate_delaunay(D, n);
        INFO(r.detail);
        CHECK(r.combinatorics_ok);
        CHECK(r.empty_circle_ok);
        CHECK(r.euler_ok);
        CHECK(r.all_vertices_used);
        // Every point is on the hull, so every triangle is virtual except the fan.
        CHECK(r.nb_virtual == n);
    }

    SECTION("many collinear points plus two off-line points")
    {
        std::vector<double> pts;
        for (int i = 0; i < 40; i++) {
            pts.push_back(double(i));
            pts.push_back(0.0);
        }
        pts.push_back(20.0);
        pts.push_back(7.0);
        pts.push_back(20.0);
        pts.push_back(-7.0);
        const index_t n = 42;
        Delaunay2d D;
        REQUIRE(D.set_vertices(n, pts.data()));
        const DelaunayReport r = validate_delaunay(D, n);
        INFO(r.detail);
        CHECK(r.combinatorics_ok);
        CHECK(r.empty_circle_ok);
        CHECK(r.euler_ok);
        CHECK(r.all_vertices_used);
    }
}

// Cross-check against real geogram. See tests/data/2d/delaunay_ref/README.md for how the
// reference files were produced and why they are stored rather than generated in CI.
TEST_CASE("2d delaunay: matches geogram", "[2d][delaunay]")
{
    const char* names[] = {"random_00010", "random_00050", "random_00200", "random_01000",
                           "random_05000"};

    for (const char* name : names) {
        DYNAMIC_SECTION(name)
        {
            const std::string path =
                std::string(VRTEST_DATA_DIR) + "/2d/delaunay_ref/" + name + ".txt";
            std::ifstream in(path);
            REQUIRE(in.good());

            std::vector<double> pts;
            std::vector<std::array<uint32_t, 3>> want;
            std::string tok;
            while (in >> tok) {
                if (tok[0] == '#') {
                    std::getline(in, tok);
                } else if (tok == "POINTS") {
                    size_t n = 0;
                    in >> n;
                    pts.resize(2 * n);
                    for (size_t i = 0; i < 2 * n; i++) {
                        std::string s;
                        in >> s;
                        pts[i] = std::strtod(s.c_str(), nullptr); // %a round-trips exactly
                    }
                } else if (tok == "TRIANGLES") {
                    size_t m = 0;
                    in >> m;
                    want.resize(m);
                    for (size_t i = 0; i < m; i++) in >> want[i][0] >> want[i][1] >> want[i][2];
                }
            }
            REQUIRE_FALSE(pts.empty());
            REQUIRE_FALSE(want.empty());

            const index_t n = index_t(pts.size() / 2);
            Delaunay2d D;
            REQUIRE(D.set_vertices(n, pts.data()));

            // Same canonical form as the generator: smallest vertex first (winding preserved),
            // then sort. Makes the comparison independent of triangle numbering.
            std::vector<std::array<uint32_t, 3>> got;
            got.reserve(D.nb_finite_triangles());
            for (index_t t = 0; t < D.nb_triangles(); t++) {
                if (!D.triangle_is_finite(t)) continue;
                uint32_t v[3];
                for (int k = 0; k < 3; k++) v[k] = D.triangle_vertex(t, index_t(k));
                const int m = (v[0] < v[1]) ? ((v[0] < v[2]) ? 0 : 2) : ((v[1] < v[2]) ? 1 : 2);
                got.push_back({v[m], v[(m + 1) % 3], v[(m + 2) % 3]});
            }
            std::sort(got.begin(), got.end());

            INFO("geogram triangles: " << want.size() << ", ours: " << got.size());
            REQUIRE(got.size() == want.size());
            for (size_t i = 0; i < got.size(); i++) {
                INFO("triangle " << i << ": got (" << got[i][0] << "," << got[i][1] << ","
                                 << got[i][2] << ") want (" << want[i][0] << "," << want[i][1]
                                 << "," << want[i][2] << ")");
                REQUIRE(got[i] == want[i]);
            }
        }
    }
}

// ==============================================================================================
// Arrangement2D
// ==============================================================================================

namespace {

using vol_rem::vr2d::Arrangement2D;
using vol_rem::vr2d::build_arrangement;
using vol_rem::vr2d::check_arrangement;
using vol_rem::vr2d::INVALID;

struct R2
{
    bigrational x, y;
};

R2 exact_of(const Arrangement2D& A, uint32_t v)
{
    R2 p;
    const bool ok = A.V[v]->getExactXYCoordinates(p.x, p.y);
    REQUIRE(ok);
    return p;
}

// The 2D form of verify_tracking.cpp's edge-provenance check (Check 3), in exact rationals:
// every output edge endpoint must be exactly collinear with the input segment, its parameter t
// must lie in [0,1], and after sorting, the sub-intervals must form an EXACT partition of [0,1]
// -- both endpoints reached, no gap, no overlap.
std::string check_provenance(const Arrangement2D& A)
{
    // Every recorded sub-edge must still BE an edge of the output triangulation. Checking only
    // the geometry (below) is not enough: a pair can stay geometrically valid on the segment
    // while no longer bounding any triangle, in which case the writer silently drops it and the
    // emitted provenance has a hole.
    std::unordered_set<uint64_t> mesh_edges;
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (!A.tri_is_finite(t)) continue;
        for (uint32_t le = 0; le < 3; le++) {
            mesh_edges.insert(Arrangement2D::ekey(A.tri_node[3 * t + (le + 1) % 3],
                                                  A.tri_node[3 * t + (le + 2) % 3]));
        }
    }
    for (uint32_t s = 0; s < A.seg.size(); s++) {
        for (uint32_t k = A.seg_head[s]; k != INVALID; k = A.subedges[k].next) {
            if (!mesh_edges.count(Arrangement2D::ekey(A.subedges[k].v0, A.subedges[k].v1)))
                return "segment " + std::to_string(s) +
                       ": recorded sub-edge is not an edge of the output triangulation";
        }
    }

    for (uint32_t s = 0; s < A.seg.size(); s++) {
        const R2 P0 = exact_of(A, A.seg[s][0]);
        const R2 P1 = exact_of(A, A.seg[s][1]);
        const bigrational dx = P1.x - P0.x;
        const bigrational dy = P1.y - P0.y;
        const bigrational len2 = dx * dx + dy * dy;
        if (len2.sgn() == 0) return "segment " + std::to_string(s) + ": zero length";

        std::vector<std::array<bigrational, 2>> iv;
        for (uint32_t k = A.seg_head[s]; k != INVALID; k = A.subedges[k].next) {
            bigrational tt[2];
            const uint32_t vv[2] = {A.subedges[k].v0, A.subedges[k].v1};
            for (int e = 0; e < 2; e++) {
                const R2 Q = exact_of(A, vv[e]);
                const bigrational qx = Q.x - P0.x;
                const bigrational qy = Q.y - P0.y;
                if ((qx * dy - qy * dx).sgn() != 0)
                    return "segment " + std::to_string(s) + ": output edge not collinear";
                tt[e] = (qx * dx + qy * dy) / len2;
            }
            if (tt[0].sgn() < 0 || tt[1].sgn() < 0) return "segment " + std::to_string(s) + ": t<0";
            if (tt[0] > bigrational(1.0) || tt[1] > bigrational(1.0))
                return "segment " + std::to_string(s) + ": t>1";
            if (!(tt[0] < tt[1]))
                return "segment " + std::to_string(s) + ": output edge reversed or degenerate";
            iv.push_back({tt[0], tt[1]});
        }
        if (iv.empty()) return "segment " + std::to_string(s) + ": no output edges";

        std::sort(iv.begin(), iv.end(),
                  [](const std::array<bigrational, 2>& a, const std::array<bigrational, 2>& b) {
                      return a[0] < b[0];
                  });
        if (iv.front()[0].sgn() != 0 || !(iv.back()[1] == bigrational(1.0)))
            return "segment " + std::to_string(s) + ": output edges do not reach both endpoints";
        for (size_t i = 1; i < iv.size(); i++)
            if (!(iv[i][0] == iv[i - 1][1]))
                return "segment " + std::to_string(s) + ": gap or overlap between output edges";
    }
    return {};
}

std::vector<std::array<uint32_t, 3>> finite_triangles(const Arrangement2D& A)
{
    std::vector<std::array<uint32_t, 3>> out;
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (!A.tri_is_finite(t)) continue;
        out.push_back({A.tri_node[3 * t], A.tri_node[3 * t + 1], A.tri_node[3 * t + 2]});
    }
    return out;
}

// Runs the whole battery on one segment soup: structural validity, exact positive orientation,
// combinatorial orientability, and exact provenance.
void check_all(const std::vector<double>& coords, const std::vector<uint32_t>& idx,
               const char* what)
{
    Arrangement2D A;
    INFO("case: " << what);
    REQUIRE(build_arrangement(coords, idx, A, false));

    std::string err;
    INFO("structure: " << err);
    REQUIRE(check_arrangement(A, &err));

    // (c) every output triangle strictly positively oriented, in exact rational arithmetic --
    // independent of the orient2D predicate the pipeline itself used.
    for (uint32_t t = 0; t < A.num_triangles(); t++) {
        if (!A.tri_is_finite(t)) continue;
        const R2 a = exact_of(A, A.tri_node[3 * t]);
        const R2 b = exact_of(A, A.tri_node[3 * t + 1]);
        const R2 c = exact_of(A, A.tri_node[3 * t + 2]);
        const bigrational area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        INFO("triangle " << t << " has non-positive exact area");
        REQUIRE(area2.sgn() > 0);
    }

    // (c) orientable and manifold
    const vrtest::TriOrientationResult o = vrtest::check_tri_orientation(finite_triangles(A));
    INFO("orientation: same_winding=" << o.same_winding << " over_shared=" << o.over_shared);
    REQUIRE(o.ok);

    // (d) provenance
    const std::string perr = check_provenance(A);
    INFO("provenance: " << perr);
    REQUIRE(perr.empty());
}

} // namespace

TEST_CASE("2d arrangement: hand-built degenerate cases", "[2d][arrangement]")
{
    SECTION("two crossing segments")
    {
        check_all({0, 0, 10, 10, 0, 10, 10, 0}, {0, 1, 2, 3}, "X");
    }
    SECTION("three segments concurrent at one point")
    {
        // The same geometric point arises from three different segment pairs. It must become ONE
        // vertex: the first crossing creates it, and the third segment's walk then meets it as an
        // existing vertex rather than constructing a coincident duplicate.
        Arrangement2D A;
        REQUIRE(build_arrangement({-10, 0, 10, 0, 0, -10, 0, 10, -10, -10, 10, 10},
                                  {0, 1, 2, 3, 4, 5}, A, false));
        CHECK(A.arena.num_ssi() == 1);
        check_all({-10, 0, 10, 0, 0, -10, 0, 10, -10, -10, 10, 10}, {0, 1, 2, 3, 4, 5},
                  "three concurrent");
    }
    SECTION("duplicate segments, both orientations")
    {
        check_all({0, 0, 10, 0}, {0, 1, 0, 1, 1, 0}, "duplicates");
    }
    SECTION("zero-length segments are dropped")
    {
        Arrangement2D A;
        REQUIRE(build_arrangement({0, 0, 10, 0, 5, 5}, {0, 0, 0, 1, 2, 2}, A, false));
        CHECK(A.seg.size() == 1);
        check_all({0, 0, 10, 0, 5, 5}, {0, 0, 0, 1, 2, 2}, "zero-length");
    }
    SECTION("duplicate input points are merged")
    {
        Arrangement2D A;
        REQUIRE(build_arrangement({0, 0, 10, 0, 0, 0, 5, 5}, {0, 1, 2, 3}, A, false));
        CHECK(A.num_input_pts == 3);
        check_all({0, 0, 10, 0, 0, 0, 5, 5}, {0, 1, 2, 3}, "duplicate points");
    }
    SECTION("partially overlapping collinear chain")
    {
        // (0,0)-(4,0) and (2,0)-(6,0): the shared span [2,4] must be reported for BOTH segments.
        check_all({0, 0, 4, 0, 2, 0, 6, 0}, {0, 1, 2, 3}, "collinear overlap");
    }
    SECTION("fully overlapping collinear chain")
    {
        check_all({0, 0, 6, 0, 2, 0, 4, 0}, {0, 1, 2, 3}, "collinear containment");
    }
    SECTION("segments sharing an endpoint")
    {
        check_all({0, 0, 5, 5, 5, -5}, {0, 1, 0, 2}, "shared endpoint");
    }
    SECTION("closed triangle outline")
    {
        check_all({0, 0, 10, 0, 5, 8}, {0, 1, 1, 2, 2, 0}, "triangle");
    }
    SECTION("T-junction: endpoint in the interior of another segment")
    {
        check_all({0, 0, 10, 0, 5, 0, 5, 7}, {0, 1, 2, 3}, "T-junction");
    }
}

TEST_CASE("2d arrangement: grids of crossing segments", "[2d][arrangement]")
{
    for (const int k : {3, 5, 8}) {
        DYNAMIC_SECTION("grid " << k << "x" << k)
        {
            std::vector<double> c;
            std::vector<uint32_t> s;
            for (int i = 0; i < k; i++) {
                c.insert(c.end(), {double(i), -1.0});
                c.insert(c.end(), {double(i), double(k)});
                c.insert(c.end(), {-1.0, double(i)});
                c.insert(c.end(), {double(k), double(i)});
            }
            for (uint32_t i = 0; i < uint32_t(4 * k); i += 2) s.insert(s.end(), {i, i + 1});

            Arrangement2D A;
            REQUIRE(build_arrangement(c, s, A, false));
            // Every vertical meets every horizontal exactly once, and no three are concurrent.
            CHECK(A.arena.num_ssi() == size_t(k) * size_t(k));
            check_all(c, s, "grid");
        }
    }
}

TEST_CASE("2d arrangement: random segment soups", "[2d][arrangement]")
{
    for (const int n : {5, 20, 60}) {
        DYNAMIC_SECTION("n=" << n)
        {
            Rnd rnd(0xA11CEu + uint64_t(n));
            std::vector<double> c;
            std::vector<uint32_t> s;
            for (int i = 0; i < n; i++) {
                c.push_back(rnd.unit());
                c.push_back(rnd.unit());
                c.push_back(rnd.unit());
                c.push_back(rnd.unit());
                s.push_back(uint32_t(4 * i) / 2);
                s.push_back(uint32_t(4 * i) / 2 + 1);
            }
            check_all(c, s, "random soup");
        }
    }
}

// Small integer coordinates make exact degeneracies (collinear, cocircular, concurrent) common
// rather than measure-zero, so this is where the walk's degenerate branches actually get hit.
TEST_CASE("2d arrangement: random segments on a small integer lattice", "[2d][arrangement]")
{
    for (const int n : {10, 40, 100}) {
        DYNAMIC_SECTION("n=" << n)
        {
            Rnd rnd(0x1A771CEu + uint64_t(n));
            std::vector<double> c;
            std::vector<uint32_t> s;
            for (int i = 0; i < n; i++) {
                c.push_back(rnd.coord(8));
                c.push_back(rnd.coord(8));
                c.push_back(rnd.coord(8));
                c.push_back(rnd.coord(8));
                s.push_back(uint32_t(2 * i));
                s.push_back(uint32_t(2 * i) + 1);
            }
            check_all(c, s, "lattice soup");
        }
    }
}

// ==============================================================================================
// KNOWN-FAILING regression cases
//
// These were found by a stress sweep over inputs with far more overlap and concurrency than the
// hand-built cases above, and they expose two genuine defects in the provenance tracking. The
// TRIANGULATION is correct in both -- structurally valid, exactly positively oriented, orientable
// and non-overlapping -- so the bug is confined to which output edges get attributed to which
// input segment.
//
//   1. collinear overlap  -> "output edge outside segment or degenerate": a sub-edge is recorded
//      for a segment that does not actually contain it.
//   2. dense concurrency  -> "output edges do not reach both endpoints": the recorded sub-edges
//      leave part of the segment uncovered.
//
// The smaller hand-built cases pass, so the defects only surface once many degenerate
// configurations interact; the sweep needed ~160 random segments on a coarse integer lattice, or
// ~480 mutually overlapping collinear segments, before one tripped.
//
// Both are cheap (well under 0.1 s each). They are expected to FAIL until the defects are fixed.
// ==============================================================================================

TEST_CASE("2d arrangement: heavy collinear overlap (KNOWN FAILURE)", "[2d][arrangement][known-bug]")
{
    // 480 collinear segments spread over 20 horizontal lines: 24 mutually overlapping segments
    // per line, no crossings at all (zero intersection points are constructed).
    Rnd rnd(0xFEEDu);
    std::vector<double> c;
    std::vector<uint32_t> idx;
    for (int i = 0; i < 480; i++) {
        const double y = double(i % 20);
        const double s = rnd.coord(20);
        const double e = s + 1 + rnd.coord(19);
        c.push_back(s);
        c.push_back(y);
        c.push_back(e);
        c.push_back(y);
        idx.push_back(uint32_t(2 * i));
        idx.push_back(uint32_t(2 * i) + 1);
    }
    check_all(c, idx, "heavy collinear overlap");
}

TEST_CASE("2d arrangement: dense lattice concurrency (KNOWN FAILURE)",
          "[2d][arrangement][known-bug]")
{
    // 160 random segments on a coarse integer lattice: heavy concurrency and many exactly
    // collinear configurations, ~2600 intersection points.
    Rnd rnd(0x1A77u + 500u);
    std::vector<double> c;
    std::vector<uint32_t> idx;
    for (int i = 0; i < 160; i++) {
        for (int k = 0; k < 4; k++) c.push_back(rnd.coord(64));
        idx.push_back(uint32_t(2 * i));
        idx.push_back(uint32_t(2 * i) + 1);
    }
    check_all(c, idx, "dense lattice concurrency");
}
