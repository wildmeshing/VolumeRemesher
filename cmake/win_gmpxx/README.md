# Bundled `gmpxx.h` for Windows

`gmpxx.h` here is a verbatim copy from **GMP 5.0.5** (`https://gmplib.org/`, LGPL — see the
license header inside the file).

## Why

The precompiled GMP used on Windows (the CGAL auxiliary-libraries package,
`CGAL-5.2.1-win64-auxiliary-libraries-gmp-mpfr.zip`, which is **GMP 5.0.1**) ships `gmp.h`
and `libgmp` but **not** the C++ wrapper `gmpxx.h`. Without it the build's `HAS_GMPXX` probe
fails, `USE_GNU_GMP_CLASSES` is left off, and the exact-arithmetic code falls back to a much
slower in-house bignum backend.

`gmpxx.h` is a header-only wrapper over GMP's C library — every class operation is inline and
calls the C functions we already link (`libgmp`); no separate `libgmp**xx**` library is
needed (the code never calls a `libgmpxx` function, and streaming goes through a custom
`operator<<` in `numerics.h`). So `CMakeLists.txt` adds this directory to `gmp::gmp`'s include
path **on Windows only**, letting the probe succeed and enabling the fast GMP-backed numbers.

## Version choice

GMP 5.0.5's `gmpxx.h` is API-compatible with the 5.0.1 `gmp.h` (newer gmpxx versions call
functions such as `mpq_cmp_z` / `mpz_primorial_ui` that were added after 5.0.1 and are absent
from that `gmp.h`). It also uses no GCC-only builtins, which helps under MSVC.

If a future toolchain cannot build it, the `HAS_GMPXX` probe (`cmake/test_gmpxx.cpp`) simply
fails and the build reverts to the in-house bignum — this never breaks the build, it only
enables the faster path when it compiles.
