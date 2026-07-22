# Wrong-sign bug in `dotProductSign2D` / `dotProductSign3D` for implicit points

**Affects:** [`MarcoAttene/Indirect_Predicates`](https://github.com/MarcoAttene/Indirect_Predicates),
verified present at HEAD `0b71b2f` (2026-05-21), and in the copy vendored in
[`wildmeshing/VolumeRemesher`](https://github.com/wildmeshing/VolumeRemesher).

**Severity:** the predicate returns the **wrong sign** — not an unreliable zero — for roughly a
third of inputs in two of its eight argument configurations. There is no filter fallback that can
recover it: the error is in the algebra, so the `interval_number` stage and the `bigfloat` stage
agree on the same wrong answer.

**Affected functions**

| function | upstream location |
| --- | --- |
| `dotProductSign2D_IEI_t` | `include/indirect_predicates.h:456` |
| `dotProductSign3D_IEI_t` | `include/indirect_predicates.h:711` |

Reached through `genericPoint::dotProductSign2D` (`include/implicit_point.hpp:289`) in exactly two
configurations, and the corresponding two for the 3D version:

```cpp
if (a.isExplicit2D() && b.isSSI()        && c.isSSI())  return dotproductSign2D_IEI(b, a, c);  // implicit_point.hpp:296
if (a.isSSI()        && b.isExplicit2D() && c.isSSI())  return dotproductSign2D_IEI(a, b, c);  // implicit_point.hpp:297
```

The other six configurations (`EEE`, `EEI`, `EIE`, `IEE`, `IIE`, `III`) are correct.

---

## The defect

`dotProductSign2D` computes the sign of `(a − c) · (b − c)`. The `IEI` path binds `p = a`,
`q = c`, `r = b`, so it needs the sign of `(p − q) · (r − q)`, where `p` and `q` are implicit
points supplied in homogeneous form,

```
p = (lpx/dp, lpy/dp)      q = (lqx/dq, lqy/dq)      r = (rx, ry)   exact
```

Clearing denominators — legitimate here because both are normalised non-negative (the
`implicitPoint2D_SSI` constructor negates when the denominator is negative, and explicit points
report `d == 1`) — gives

```
dp*dq * (p − q).x  =  lpx*dq − lqx*dp                (A)
     dq * (r − q).x  =  rx*dq  − lqx                   (B)
```

so `sign( (p−q)·(r−q) ) = sign( (A_x·B_x + A_y·B_y) )`.

The generated code computes `A` as:

```cpp
const T dqp  = (dq * dp);
const T pxq  = (lpx * dqp);     // = lpx * dq * dp     <-- extra factor of dp
const T pyq  = (lpy * dqp);
const T rxq  = (rx * dq);
const T ryq  = (ry * dq);
const T lqxd = (lqx * dp);
const T lqyd = (lqy * dp);
const T lx   = (pxq - lqxd);    // = lpx*dq*dp − lqx*dp   instead of   lpx*dq − lqx*dp
const T ly   = (pyq - lqyd);
const T gx   = (rxq - lqx);     // correct
const T gy   = (ryq - lqy);     // correct
const T d    = (lx * gx) + (ly * gy);
return sgn(d);
```

The first term of `lx` carries an extra factor of `dp` that the second term does not. Factoring:

```
lx_code = dp * (lpx*dq − lqx)  =  dp*dq * (lpx − qx)
lx_true =      lpx*dq − lqx*dp =  dp*dq * (px  − qx)
```

The code has substituted `lpx` for `px = lpx/dp`. The two agree **exactly when `dp == 1`**, i.e.
when the first argument is an *explicit* point.

That is precisely the case the dispatcher never produces: both call sites at
`implicit_point.hpp:296-297` pass an implicit point as `p`, so `dp != 1` and the computed vector
is not parallel to `p − q`. The resulting dot product has an essentially arbitrary sign.

`dotProductSign3D_IEI_t` (`indirect_predicates.h:711`) has the identical structure and the
identical defect, with `pzq = lpz * dqp` as the third component.

---

## Evidence

Measured against exact rational arithmetic (`bigrational`, via `getExactXYCoordinates`) over 400
random explicit points and 100 `implicitPoint2D_SSI` points built from them, ~3000 comparisons per
configuration, integer coordinates in [−20, 20]:

| config | comparisons | wrong sign | spurious zero |
| --- | --- | --- | --- |
| EEE | 2976 | 0 | 0 |
| EEI | 2990 | 0 | 0 |
| EIE | 2993 | 0 | 0 |
| **EII** | 2958 | **858 (29.0%)** | 0 |
| IEE | 2990 | 0 | 0 |
| **IEI** | 2966 | **869 (29.3%)** | 0 |
| IIE | 2969 | 0 | 0 |
| III | 2917 | 0 | 0 |

After the fix below, all eight configurations report 0 wrong signs.

Note the failure is a *wrong sign*, not a *zero*: the predicate reports a confident, incorrect
answer. Callers that treat 0 as "degenerate, fall back" get no signal at all.

### Reproduction

```cpp
#include <implicit_point.h>
#include <cstdio>

static int exact_sign(const genericPoint& A, const genericPoint& B, const genericPoint& C) {
    bigrational ax, ay, bx, by, cx, cy;
    A.getExactXYCoordinates(ax, ay);
    B.getExactXYCoordinates(bx, by);
    C.getExactXYCoordinates(cx, cy);
    return ((ax - cx) * (bx - cx) + (ay - cy) * (by - cy)).sgn();
}

int main() {
    // Two SSI points on the line y = x, and one explicit point.
    explicitPoint2D p00(0, 0), p22(2, 2), p20(2, 0), p02(0, 2), p03(0, 3), p30(3, 0);
    implicitPoint2D_SSI m1(p00, p22, p20, p02);   // (1.0, 1.0)
    implicitPoint2D_SSI m2(p00, p22, p03, p30);   // (1.5, 1.5)
    explicitPoint2D r(10, 0);

    // config IEI: (implicit, explicit, implicit)
    printf("predicate = %d, exact = %d\n",
           genericPoint::dotProductSign2D(m2, r, m1), exact_sign(m2, r, m1));
    return 0;
}
```

(Any pair of distinct SSI points and one explicit point exercises the path; the table above is
the statistical picture over many such triples.)

---

## Fix

In `dotProductSign2D_IEI_t`, drop the extra factor:

```diff
-    const T dqp = (dq * dp);
-    const T pxq = (lpx * dqp);
-    const T pyq = (lpy * dqp);
+    const T pxq = (lpx * dq);
+    const T pyq = (lpy * dq);
     const T rxq = (rx * dq);
     const T ryq = (ry * dq);
     const T lqxd = (lqx * dp);
     const T lqyd = (lqy * dp);
```

and the same in `dotProductSign3D_IEI_t`, additionally replacing `pzq = (lpz * dqp)` with
`pzq = (lpz * dq)`. `dqp` then becomes unused and can be removed.

The overall expression is scaled by `dp * dq²` relative to the true dot product. Both denominators
are normalised non-negative and `dq² > 0`, so the scaling is sign-preserving and no further change
is needed.

This is a two-line change per function and does not alter the filtering structure, the rounding-
mode handling, or the `interval_number → bigfloat` cascade.

If these files are generated, the defect is presumably in the generator's homogeneous-coordinate
handling for the case where two *different* implicit points appear in the same expression — note
that the sibling `IIE` and `III` variants, which also take multiple implicit points, are correct,
so it looks specific to the `IEI` argument pattern rather than to implicit points in general.

---

## Why this was not noticed

`dotProductSign2D` and `dotProductSign3D` appear to have no callers in the projects that vendor
these headers. In `wildmeshing/VolumeRemesher` the entire 3D pipeline never calls either, so the
functions were dead code until a new 2D feature used `dotProductSign2D` to decide whether a
collinear neighbour lies forward or backward along a segment. That made a mesh traversal fail
outright, which is how the bug surfaced.

Every *other* arithmetic statement in `indirect_predicates.h` is byte-identical between upstream
HEAD and the vendored copy, so this is the only known divergence.

---

## Suggested regression test

Compare all eight explicit/implicit configurations against exact rational arithmetic, rather than
against a reference implementation, so the test is independent of the filtering cascade. The
version used downstream is `"2d predicates: dotProductSign2D agrees with exact arithmetic"` in
[`tests/unit_tests_2d.cpp`](../tests/unit_tests_2d.cpp); it asserts each configuration is actually
exercised (>100 comparisons) so a dispatch change cannot silently skip a path.
