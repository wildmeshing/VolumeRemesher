#!/usr/bin/env python3
"""Playground for the degenerate face polygons captured from `mesh_generator -t`.

mesh_generator's ear-clipping (BSPcomplex::triangulateFace) fails to triangulate
these BSP faces into all-positive-area triangles; it force-clips a zero-area
triangle instead (and warns). Each face was dumped as its ORIGINAL boundary
polygon in *exact rational* 3D coordinates (VR_DUMP_DEGENERATE_FACES), so we can
reason about them exactly here with fractions.Fraction and figure out how a
correct triangulation should work.

Data format (one block per face):
    POLYGON <face_id> <n>
    <x> <y> <z>        # n lines, each an exact rational "num/den" (or integer)

Usage:
    triangulation_playground.py degenerate_faces.txt
"""

import sys
from fractions import Fraction


# ---- exact rational 3D vector helpers ---------------------------------------

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def is_zero(v):
    return v[0] == 0 and v[1] == 0 and v[2] == 0


def collinear(a, b, c):
    """Exact: are the three 3D points collinear (zero-area triangle)?"""
    return is_zero(cross(sub(b, a), sub(c, a)))


def newell_normal(pts):
    """Exact polygon normal via Newell's method (robust for any planar polygon)."""
    nx = ny = nz = Fraction(0)
    m = len(pts)
    for i in range(m):
        a, b = pts[i], pts[(i + 1) % m]
        nx += (a[1] - b[1]) * (a[2] + b[2])
        ny += (a[2] - b[2]) * (a[0] + b[0])
        nz += (a[0] - b[0]) * (a[1] + b[1])
    return (nx, ny, nz)


def dominant_axis(normal):
    return max(range(3), key=lambda k: abs(normal[k]))


def project(pts, drop_axis):
    keep = [k for k in range(3) if k != drop_axis]
    return [(p[keep[0]], p[keep[1]]) for p in pts]


def twice_signed_area(poly2d):
    """Exact 2x signed area (shoelace) of a 2D polygon."""
    s = Fraction(0)
    m = len(poly2d)
    for i in range(m):
        x1, y1 = poly2d[i]
        x2, y2 = poly2d[(i + 1) % m]
        s += x1 * y2 - x2 * y1
    return s


# ---- parsing ----------------------------------------------------------------

def parse(path):
    polys = []
    with open(path) as f:
        toks = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    i = 0
    while i < len(toks):
        assert toks[i][0] == "POLYGON", toks[i]
        fid, n = int(toks[i][1]), int(toks[i][2])
        i += 1
        pts = []
        for _ in range(n):
            pts.append(tuple(Fraction(t) for t in toks[i]))
            i += 1
        polys.append((fid, pts))
    return polys


# ---- analysis ---------------------------------------------------------------

def distinct_points(pts):
    seen, out = set(), []
    for p in pts:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def all_collinear(pts):
    d = distinct_points(pts)
    if len(d) < 3:
        return True
    a, b = d[0], d[1]
    return all(collinear(a, b, c) for c in d[2:])


def orient2d(a, b, c):
    """Exact 2D orientation = 2x signed area of triangle a,b,c."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def is_convex(poly2d):
    """Convex (allowing collinear/flat vertices): all turns share one sign."""
    m = len(poly2d)
    sign = 0
    for i in range(m):
        o = orient2d(poly2d[i], poly2d[(i + 1) % m], poly2d[(i + 2) % m])
        if o > 0:
            if sign < 0:
                return False
            sign = 1
        elif o < 0:
            if sign > 0:
                return False
            sign = -1
    return True


def fan_triangulate(poly2d):
    """Fan from the first strictly-convex vertex. Returns index triples."""
    m = len(poly2d)
    apex = next((i for i in range(m)
                 if orient2d(poly2d[(i - 1) % m], poly2d[i], poly2d[(i + 1) % m]) != 0), None)
    if apex is None:
        return None
    return [(apex, (apex + j) % m, (apex + j + 1) % m) for j in range(1, m - 1)]


def earclip_strictly_convex(poly2d, want_sign):
    """Triangulate by repeatedly clipping any STRICTLY-convex ear (nonzero,
    correctly-signed ear triangle). Valid for convex polygons even with flat
    (collinear) vertices: a positive-area convex polygon always has a strictly-
    convex ear, and clipping it keeps the remainder convex + positive-area, so the
    flat vertices get absorbed without ever forming a collinear triangle."""
    idx = list(range(len(poly2d)))
    tris = []
    while len(idx) > 3:
        m = len(idx)
        for k in range(m):
            a, b, c = idx[(k - 1) % m], idx[k], idx[(k + 1) % m]
            o = orient2d(poly2d[a], poly2d[b], poly2d[c])
            if o != 0 and (o > 0) == (want_sign > 0):
                tris.append((a, b, c))
                idx.pop(k)
                break
        else:
            return None  # no strictly-convex ear (should not happen here)
    tris.append(tuple(idx))
    return tris


def incircle(a, b, c, d):
    """Exact in-circle test. For a,b,c in CCW order returns >0 iff d is strictly
    inside the circumcircle of a,b,c, ==0 if on it, <0 if outside."""
    ax, ay = a[0] - d[0], a[1] - d[1]
    bx, by = b[0] - d[0], b[1] - d[1]
    cx, cy = c[0] - d[0], c[1] - d[1]
    return ((ax * ax + ay * ay) * (bx * cy - cx * by)
            - (bx * bx + by * by) * (ax * cy - cx * ay)
            + (cx * cx + cy * cy) * (ax * by - bx * ay))


def triangulate_convex(poly2d):
    """Robust, orient2D-only triangulation of a convex polygon that may carry
    collinear (flat / Steiner) boundary vertices -- exactly our BSP faces.

    Provably all-positive:
      1. Fan the strictly-convex "corners" (orient2D != 0). Corners of a convex
         polygon are never 3-collinear, so every fan triangle is positive.
      2. Insert each flat vertex f by splitting the one triangle whose *boundary*
         edge (u,v) contains f: (u,v,w) -> (u,f,w)+(f,v,w). Since w is off the line
         u-v, both are positive. (f is a boundary vertex, so it can only lie on a
         boundary edge, never strictly inside an interior diagonal -> unambiguous.)

    Uses ONLY orient2D and a between-test -- no in-circle, no super-triangle, no
    tie-breaks. In C++ this maps to genericPoint::orient2D{xy,yz,zx} (dominant
    plane) + dotProductSign3D for "f strictly between u and v"."""
    n = len(poly2d)
    corners = [i for i in range(n)
               if orient2d(poly2d[(i - 1) % n], poly2d[i], poly2d[(i + 1) % n]) != 0]
    if len(corners) < 3:
        return None  # degenerate (zero-area) polygon -- should not happen for a BSP face
    apex = corners[0]
    tris = [[apex, corners[j], corners[j + 1]] for j in range(1, len(corners) - 1)]
    for f in range(n):
        if orient2d(poly2d[(f - 1) % n], poly2d[f], poly2d[(f + 1) % n]) != 0:
            continue  # a corner, already placed
        for ti in range(len(tris)):
            x, y, z = tris[ti]
            for u, v, w in ((x, y, z), (y, z, x), (z, x, y)):
                if orient2d(poly2d[u], poly2d[v], poly2d[f]) != 0:
                    continue
                fu = (poly2d[u][0] - poly2d[f][0], poly2d[u][1] - poly2d[f][1])
                fv = (poly2d[v][0] - poly2d[f][0], poly2d[v][1] - poly2d[f][1])
                if fu[0] * fv[0] + fu[1] * fv[1] < 0:  # f strictly between u and v
                    tris[ti] = [u, f, w]
                    tris.append([f, v, w])
                    break
            else:
                continue
            break
    return tris


def incircle_sos(a, b, c, d, P):
    """In-circle test with Simulation-of-Simplicity tie-breaking, so it is never
    zero even for cocircular / collinear inputs (which our faces have, from
    Steiner points). a,b,c,d are indices into P with a,b,c CCW; returns +1 if d is
    (virtually) inside the circumcircle, -1 otherwise.

    SoS = perturb each lifted coordinate w_i -> w_i + eps_i with eps ordered by
    global index (smaller index = larger eps). When the exact determinant is 0,
    its sign is that of the first nonzero perturbation term; the coefficient of
    eps_i is the cofactor of point i, i.e. +/- orient2d of the other three. So the
    tie collapses to an orient2d on the smallest-index point that gives a nonzero
    one -- purely integer/rational, no actual epsilon."""
    v = incircle(P[a], P[b], P[c], P[d])
    if v != 0:
        return 1 if v > 0 else -1
    rows = (a, b, c, d)
    for r in sorted(range(4), key=lambda k: rows[k]):
        others = [rows[k] for k in range(4) if k != r]
        o = orient2d(P[others[0]], P[others[1]], P[others[2]])
        if o != 0:
            s = 1 if o > 0 else -1
            return s if r % 2 == 0 else -s
    return -1  # all four collinear (never happens for a non-degenerate triangle)


def delaunay_triangulate(pts):
    """From-scratch 2D Delaunay triangulation (Bowyer-Watson incremental) using
    only exact predicates (orient2d, incircle) -- no external libraries, so it
    ports directly to C++. Returns triangles as CCW index triples into `pts`.

    For points in convex position (our BSP faces) the Delaunay triangulation of
    the vertices already contains every boundary edge and covers the polygon, so
    no constraint-enforcement step is needed."""
    n = len(pts)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    D = max(max(xs) - min(xs), max(ys) - min(ys))
    if D == 0:
        D = 1
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    # The super-triangle must strictly contain every input point AND lie outside
    # every input triangle's circumcircle. Near-collinear vertices (our Steiner
    # points) form very "thin" triangles with enormous circumcircles, so a small
    # super-triangle gets swallowed and the cavity corrupts (dropped triangles).
    # Size it from the input: a triangle's circumradius R = (prod of sides)/(4 area)
    # <= (2D)^3 / (4 * A_min/2) = 4 D^3 / A_min, where A_min is the smallest positive
    # |orient2d| over vertex triples (all rational -- no sqrt). The circumcircle lies
    # within D + 2R of the centre, so reach = D + 8 D^3 / A_min is safe; double it for
    # margin. (Fully robust alternative for the C++ port: ghost vertices at infinity.)
    a_min = None
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                o = abs(orient2d(pts[i], pts[j], pts[k]))
                if o != 0 and (a_min is None or o < a_min):
                    a_min = o
    if a_min is None:
        a_min = 1  # degenerate polygon (all collinear); never happens for our faces
    M = 2 * (D + 8 * D ** 3 / a_min)
    P = list(pts) + [(cx - 3 * M, cy - 2 * M), (cx + 3 * M, cy - 2 * M), (cx, cy + 4 * M)]
    s0, s1, s2 = n, n + 1, n + 2
    if orient2d(P[s0], P[s1], P[s2]) < 0:
        s1, s2 = s2, s1
    tris = [(s0, s1, s2)]  # kept CCW throughout

    for i in range(n):
        bad = [t for t in tris if incircle_sos(t[0], t[1], t[2], i, P) > 0]
        # Cavity boundary = directed edges of bad triangles whose reverse is not
        # also a bad-triangle edge. Being CCW edges of the removed cavity, (u,v,p)
        # comes out CCW.
        directed = []
        for t in bad:
            directed += [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
        dset = set(directed)
        tris = [t for t in tris if t not in bad]
        for (u, v) in directed:
            if (v, u) not in dset:  # boundary edge
                tris.append((u, v, i))

    return [t for t in tris if t[0] < n and t[1] < n and t[2] < n]


def all_triangles_positive(poly2d, tris, want_sign):
    """Do all triangles have the target (nonzero) orientation sign?"""
    for (a, b, c) in tris:
        o = orient2d(poly2d[a], poly2d[b], poly2d[c])
        if o == 0 or (o > 0) != (want_sign > 0):
            return False
    return True


def analyze(polys):
    zero_area = positive_area = 0
    convex = fan_ok = earclip_ok = delaunay_ok = cornersflats_ok = 0
    for fid, pts in polys:
        d = distinct_points(pts)
        normal = newell_normal(pts)
        poly2d = project(pts, dominant_axis(normal)) if not is_zero(normal) else project(pts, 2)
        area2 = twice_signed_area(poly2d)
        degenerate = all_collinear(pts) or area2 == 0
        if degenerate:
            zero_area += 1
            print(f"face {fid:>8}: n={len(pts)} -> ZERO-AREA (collinear)")
            continue
        positive_area += 1
        cvx = is_convex(poly2d)
        convex += cvx
        fan = fan_triangulate(poly2d)
        fan_pos = fan is not None and len(fan) == len(pts) - 2 and \
            all_triangles_positive(poly2d, fan, area2)
        fan_ok += fan_pos
        ec = earclip_strictly_convex(poly2d, area2)
        ec_pos = ec is not None and len(ec) == len(pts) - 2 and \
            all_triangles_positive(poly2d, ec, area2)
        earclip_ok += ec_pos
        # Delaunay: all triangles CCW (orient2d > 0), count n-2, every vertex used.
        dt = delaunay_triangulate(poly2d)
        dt_pos = (len(dt) == len(pts) - 2
                  and all(orient2d(poly2d[a], poly2d[b], poly2d[c]) > 0 for a, b, c in dt)
                  and {v for t in dt for v in t} == set(range(len(pts))))
        delaunay_ok += dt_pos
        cf = triangulate_convex(poly2d)
        cf_pos = cf is not None and len(cf) == len(pts) - 2 and \
            all_triangles_positive(poly2d, cf, area2) and \
            {v for t in cf for v in t} == set(range(len(pts)))
        cornersflats_ok += cf_pos
        print(f"face {fid:>8}: n={len(pts)} 2*area={float(area2):+.4g} convex={cvx} "
              f"earclip_ok={ec_pos} delaunay_ok={dt_pos} cornersflats_ok={cf_pos}")

    print()
    print(f"total polygons                         : {len(polys)}")
    print(f"  genuinely zero-area                  : {zero_area}")
    print(f"  positive-area                        : {positive_area}")
    print(f"  ... convex                           : {convex}")
    print(f"  ... naive fan -> all positive        : {fan_ok}")
    print(f"  ... strictly-convex ear-clip -> all +: {earclip_ok}")
    print(f"  ... from-scratch Delaunay -> all +   : {delaunay_ok}")
    print(f"  ... corners+flats (orient2D only)    : {cornersflats_ok}")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "degenerate_faces.txt"
    polys = parse(path)
    analyze(polys)


if __name__ == "__main__":
    main()
