// Unit tests for the 2D arrangement pipeline.
//
// Unlike integration_tests.cpp (which shells out to the mesh_generator binary), this target links
// the library and exercises the internal 2D APIs directly.

#include <catch2/catch_test_macros.hpp>

#include <VolumeRemesher/2d/delaunay2d.h>
#include <VolumeRemesher/2d/predicates2d.h>
#include <VolumeRemesher/numerics.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

using vol_rem::bigrational;
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
