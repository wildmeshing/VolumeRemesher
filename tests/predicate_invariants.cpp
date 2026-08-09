// Invariants this library relies on from Indirect_Predicates, checked directly so a repin
// that loses one fails here in seconds instead of surfacing as a bad mesh hours later.
//
// The indirect predicates work on homogeneous coordinates (lx, ly, lz, d) and return
// sgn(det), where det is the true determinant scaled by the operands' denominators --
// orient3d_indirect_IEEE_t and friends multiply the explicit coordinates through by d and
// never look at its sign. So the predicates are correct only if every d > 0, which two
// things maintain:
//
//   1. the constructor negates (lx, ly, lz, d) when d is RELIABLY negative;
//   2. getIntervalLambda returns false when d's sign is NOT reliable, so the predicate
//      abandons the interval filter and redoes the work in exact bigfloat.
//
// implicitPoint3D_TBC is what makeTetrahedra uses for the barycenter apex of a cell that
// cannot be tetrahedralized from one of its vertices, and it used to do NEITHER: never
// normalized, and getIntervalLambda returned an unconditional true. Because a TBC's
// denominator is the PRODUCT of its generators' denominators, it inherits an undecided
// sign from any generator whose own filter gave up -- routine for an LPI whose line is
// nearly parallel to its plane, i.e. anywhere near-planar. The predicate then answered
// (true orientation) * sgn(d), and makeTetrahedra's winding-flip turned that into tets
// with negative volume. See MarcoAttene/Indirect_Predicates#15.

#include <catch2/catch_test_macros.hpp>

#include <VolumeRemesher/implicit_point.h>

#include <cstdint>
#include <memory>
#include <vector>

using namespace vol_rem;

namespace {

// Deterministic: a failure here has to be reproducible.
struct Rng
{
    uint64_t s = 0x9E3779B97F4A7C15ull;
    double operator()(double lo, double hi)
    {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return lo + (hi - lo) * ((double)(s >> 11) / (double)(1ull << 53));
    }
};

// Owns every point for the lifetime of one test: the implicit types hold their generators
// by reference, so the generators must outlive them.
struct Arena
{
    std::vector<std::unique_ptr<genericPoint>> owned;

    const explicitPoint3D& expl(double x, double y, double z)
    {
        owned.push_back(std::make_unique<explicitPoint3D>(x, y, z));
        return static_cast<const explicitPoint3D&>(*owned.back());
    }
    template <class T, class... A>
    const T& make(A&&... a)
    {
        owned.push_back(std::make_unique<T>(std::forward<A>(a)...));
        return static_cast<const T&>(*owned.back());
    }
};

bool denominator_sign_is_decided(const genericPoint& p)
{
    if (p.isExplicit3D()) return true;
    interval_number lx, ly, lz, d;
    p.getIntervalLambda(lx, ly, lz, d);
    return d.signIsReliable();
}

// An LPI whose line is nearly parallel to its plane, so its denominator suffers near-total
// cancellation and its interval filter gives up. `tilt` is how far from parallel.
const genericPoint& near_parallel_lpi(Arena& A, Rng& rnd, double tilt)
{
    const double r[3] = {rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)};
    const double s[3] = {rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)};
    const double t[3] = {rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)};
    const double u[3] = {s[0] - r[0], s[1] - r[1], s[2] - r[2]};
    const double v[3] = {t[0] - r[0], t[1] - r[1], t[2] - r[2]};
    const double n[3] = {
        u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};

    const double a = rnd(-1, 1), b = rnd(-1, 1);
    double dir[3];
    for (int k = 0; k < 3; k++) dir[k] = a * u[k] + b * v[k] + tilt * n[k];

    const double p[3] = {rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)};
    return A.make<implicitPoint3D_LPI>(
        A.expl(p[0], p[1], p[2]), A.expl(p[0] + dir[0], p[1] + dir[1], p[2] + dir[2]),
        A.expl(r[0], r[1], r[2]), A.expl(s[0], s[1], s[2]), A.expl(t[0], t[1], t[2]));
}

// Sign of the exact rational signed volume of (a,b,c,d), the ground truth orient3D must
// agree with. Same convention makeTetrahedra's winding-flip uses.
int exact_orientation(
    const genericPoint& a,
    const genericPoint& b,
    const genericPoint& c,
    const genericPoint& d)
{
    const genericPoint* g[4] = {&a, &b, &c, &d};
    bigrational p[4][3];
    for (int i = 0; i < 4; i++) g[i]->getExactXYZCoordinates(p[i][0], p[i][1], p[i][2]);
    bigrational e1[3], e2[3], e3[3];
    for (int j = 0; j < 3; j++) {
        e1[j] = p[1][j] - p[0][j];
        e2[j] = p[2][j] - p[0][j];
        e3[j] = p[3][j] - p[0][j];
    }
    const bigrational v = (e1[1] * e2[2] - e1[2] * e2[1]) * e3[0] +
        (e1[2] * e2[0] - e1[0] * e2[2]) * e3[1] + (e1[0] * e2[1] - e1[1] * e2[0]) * e3[2];
    return v.sgn();
}

} // namespace

TEST_CASE("predicates: an undecided denominator sign is never reported as usable", "[predicates]")
{
    Rng rnd;
    size_t built = 0, undecided = 0, reported_usable_anyway = 0;

    for (int trial = 0; trial < 4000 && undecided < 200; trial++) {
        Arena A;
        const genericPoint& g0 = near_parallel_lpi(A, rnd, 1e-18 * (1 + (trial % 64)));
        if (denominator_sign_is_decided(g0)) continue; // not a useful trial

        const implicitPoint3D_TBC& tbc = A.make<implicitPoint3D_TBC>(
            g0, A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)),
            A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)),
            A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)));
        built++;

        interval_number lx, ly, lz, d;
        const bool usable = tbc.getIntervalLambda(lx, ly, lz, d);
        if (d.signIsReliable() && !d.isNegative()) continue; // filter is fine here
        undecided++;
        // The contract: `true` promises the predicate that d is positive, so the scaled
        // determinant's sign IS the orientation. Anything else must return false.
        if (usable) reported_usable_anyway++;
    }

    INFO(
        "TBCs built from a generator whose filter gave up: "
        << built << ", of which the TBC's own denominator sign is undecided: " << undecided);
    REQUIRE(undecided > 0); // otherwise the test proves nothing
    CHECK(reported_usable_anyway == 0);
}

TEST_CASE("predicates: orient3D on a TBC agrees with exact rational coordinates", "[predicates]")
{
    Rng rnd;
    size_t queries = 0, wrong = 0;

    for (int trial = 0; trial < 4000 && queries < 4000; trial++) {
        Arena A;
        const genericPoint& g0 = near_parallel_lpi(A, rnd, 1e-18 * (1 + (trial % 64)));
        if (denominator_sign_is_decided(g0)) continue;

        const implicitPoint3D_TBC& tbc = A.make<implicitPoint3D_TBC>(
            g0, A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)),
            A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)),
            A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1)));

        bigrational bx, by, bz;
        if (!tbc.getExactXYZCoordinates(bx, by, bz)) continue;

        for (int q = 0; q < 40; q++) {
            const explicitPoint3D& a = A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
            const explicitPoint3D& b = A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
            const explicitPoint3D& c = A.expl(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
            const int exact = exact_orientation(a, b, c, tbc);
            if (exact == 0) continue;
            queries++;
            if (genericPoint::orient3D(a, b, c, tbc) != exact) wrong++;
        }
    }

    INFO("orientation queries against a TBC with an undecided denominator: " << queries);
    REQUIRE(queries > 0); // otherwise the test proves nothing
    CHECK(wrong == 0);
}
