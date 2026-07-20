// Exact 2D predicate layer for the 2D arrangement pipeline.
//
// This header exposes the three predicates that geogram's Delaunay2d needs, under geogram's own
// names and signatures, so that delaunay2d.cpp can stay a near-verbatim port of
// geogram/delaunay/delaunay_2d.cpp. Underneath, every one of them is implemented with THIS
// repository's exact kernel (numerics.h / hand_optimized_predicates.hpp / indirect_predicates.h)
// rather than geogram's expansion_nt.
//
//   PCK::orient_2d              -> vol_rem::orient2d  (hand_optimized_predicates.hpp)
//   PCK::points_are_identical_2d -> exact double equality
//   PCK::in_circle_2d_SOS       -> vol_rem::incircle  (indirect_predicates.h) + the SOS layer below
//
// The SOS (Simulation of Simplicity) layer is the only genuinely new predicate. See the comment
// on in_circle_2d_SOS in predicates2d.cpp for the derivation.

#ifndef VOLUMEREMESHER_2D_PREDICATES2D_H
#define VOLUMEREMESHER_2D_PREDICATES2D_H

#include <VolumeRemesher/implicit_point.h>
#include <VolumeRemesher/indirect_predicates.h>

namespace vol_rem {
namespace vr2d {

// geogram's Sign type, reproduced so ported code reads unchanged.
enum Sign { NEGATIVE = -1, ZERO = 0, POSITIVE = 1 };

inline Sign geo_sgn(double x) { return (x > 0) ? POSITIVE : ((x < 0) ? NEGATIVE : ZERO); }
inline Sign to_sign(int s) { return (s > 0) ? POSITIVE : ((s < 0) ? NEGATIVE : ZERO); }

namespace PCK {

// Sign of det | x1 y1 1 ; x2 y2 1 ; x3 y3 1 |, i.e. POSITIVE when (p0,p1,p2) is counterclockwise.
// Exact, never approximate.
inline Sign orient_2d(const double* p0, const double* p1, const double* p2)
{
    return to_sign(orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]));
}

inline bool points_are_identical_2d(const double* p0, const double* p1)
{
    return p0[0] == p1[0] && p0[1] == p1[1];
}

// POSITIVE iff p3 lies strictly inside the circle through (p0,p1,p2), assuming (p0,p1,p2) is
// counterclockwise. NEVER returns ZERO: four cocircular points are disambiguated by symbolic
// perturbation, so the caller always gets a strict, consistent answer.
//
// Precondition: the four points are pairwise distinct and (p0,p1,p2) is not collinear.
Sign in_circle_2d_SOS(const double* p0, const double* p1, const double* p2, const double* p3);

} // namespace PCK

} // namespace vr2d
} // namespace vol_rem

#endif
