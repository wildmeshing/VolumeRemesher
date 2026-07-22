#include "predicates2d.h"

#include <cassert>

namespace vol_rem {
namespace vr2d {
namespace PCK {

// ---------------------------------------------------------------------------------------------
// in_circle_2d_SOS
//
// This is the only predicate the 2D pipeline needs that numerics.h does not already provide.
// vol_rem::incircle() is exact but returns 0 for four cocircular points; geogram's Delaunay2d
// requires a predicate that never returns 0, because triangle_is_conflict() must partition
// triangles into "in conflict" / "not in conflict" with no third outcome. Geogram gets that from
// side3_exact_SOS on its own expansion_nt stack; importing that would defeat the point of using
// this repository's kernel, so we build the SOS layer on top of the exact predicates we have.
//
// DERIVATION
//
// Write w_i = x_i^2 + y_i^2 and consider the lifted determinant
//
//         | x0 y0 w0 1 |
//     D = | x1 y1 w1 1 |
//         | x2 y2 w2 1 |
//         | x3 y3 w3 1 |
//
// D is exactly what vol_rem::incircle(p0,p1,p2,p3) computes (verified numerically and pinned by
// PredicateConventions in the unit tests): D > 0 iff p3 is inside circumcircle(p0,p1,p2) when
// (p0,p1,p2) is counterclockwise.
//
// Simulation of Simplicity perturbs each point's LIFTED coordinate, w_i -> w_i + eps_i, with
// eps_i > 0 infinitesimal and all eps_i of pairwise different magnitude. Geometrically this lifts
// p_i off the paraboloid, which turns the Delaunay triangulation into a regular (weighted)
// triangulation -- always a valid triangulation, which is what makes this a legitimate
// perturbation rather than an arbitrary tie-break.
//
// Expanding D along the w column:
//
//     D(eps) = D + sum_i eps_i * C_i
//
// where C_i is the cofactor of entry (i, 2). Deleting row i and the w column leaves the 3x3
// matrix of the other three points over columns (x, y, 1), whose determinant is precisely
// orient_2d of those three points in increasing row order. With the (-1)^(i+2) cofactor sign:
//
//     C_0 = +orient_2d(p1, p2, p3)
//     C_1 = -orient_2d(p0, p2, p3)
//     C_2 = +orient_2d(p0, p1, p3)
//     C_3 = -orient_2d(p0, p1, p2)
//
// When D == 0 the sign of D(eps) is the sign of the surviving term with the LARGEST eps. To make
// that choice independent of the order in which the caller passed the points -- without which the
// predicate would be inconsistent under permutation and the triangulation would tear -- eps must
// be a function of the POINT, not of its argument position. We therefore order the four points
// lexicographically and give the lexicographically smallest point the largest eps. Lexicographic
// rank restricted to a subset preserves the relative order of a global ranking, so this is
// consistent with a single fixed assignment of infinitesimals over the whole point set.
//
// The loop below is that rule: walk the points from lexicographically smallest, return the first
// non-zero cofactor.
//
// TERMINATION: C_3 = -orient_2d(p0,p1,p2), and the precondition says (p0,p1,p2) is a real
// triangle, so C_3 != 0 and the loop always returns.
// ---------------------------------------------------------------------------------------------

Sign in_circle_2d_SOS(const double* p0, const double* p1, const double* p2, const double* p3)
{
    const int s = incircle(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
    if (s != 0) return to_sign(s);

    // The four points are exactly cocircular. Disambiguate symbolically.
    const double* p[4] = {p0, p1, p2, p3};

    // Lexicographic order of the four points, smallest first. Insertion sort: 4 elements, and it
    // is a total order because the points are pairwise distinct (precondition), so the result is
    // deterministic on every platform.
    int order[4] = {0, 1, 2, 3};
    for (int i = 1; i < 4; i++) {
        const int v = order[i];
        int j = i - 1;
        while (j >= 0 && (p[order[j]][0] > p[v][0] ||
                          (p[order[j]][0] == p[v][0] && p[order[j]][1] > p[v][1]))) {
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = v;
    }

    for (int k = 0; k < 4; k++) {
        int c = 0;
        switch (order[k]) {
        case 0: c = +orient2d(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]); break;
        case 1: c = -orient2d(p0[0], p0[1], p2[0], p2[1], p3[0], p3[1]); break;
        case 2: c = +orient2d(p0[0], p0[1], p1[0], p1[1], p3[0], p3[1]); break;
        default: c = -orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]); break;
        }
        if (c != 0) return to_sign(c);
    }

    // Unreachable: would require all four points collinear, contradicting the precondition that
    // (p0,p1,p2) is a real triangle.
    assert(false && "in_circle_2d_SOS: degenerate input (collinear p0,p1,p2)");
    return POSITIVE;
}

} // namespace PCK
} // namespace vr2d
} // namespace vol_rem
