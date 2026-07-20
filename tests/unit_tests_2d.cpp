// Unit tests for the 2D arrangement pipeline.
//
// Unlike integration_tests.cpp (which shells out to the mesh_generator binary), this target links
// the library and exercises the internal 2D APIs directly.

#include <catch2/catch_test_macros.hpp>

#include <VolumeRemesher/2d/predicates2d.h>
#include <VolumeRemesher/numerics.h>

#include <array>
#include <cstdint>
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
