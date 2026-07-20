# geogram Delaunay2d reference triangulations

These files pin `vol_rem::vr2d::Delaunay2d` (a port of geogram's `Delaunay2d`, see
`src/VolumeRemesher/2d/delaunay2d.h`) against the output of **real geogram**.

Each file is self-contained: it carries both the input points and the triangulation geogram
produced for them, so the unit test never has to reproduce the generator's RNG.

    # comment lines
    POINTS <n>
    <x> <y>            # n lines, printf %a (exact binary), so they round-trip bit-for-bit
    TRIANGLES <m>
    <a> <b> <c>        # m lines, canonical form: smallest vertex first (winding preserved),
                       # then the whole list sorted lexicographically

Consumed by the `2d delaunay: matches geogram` test case in `tests/unit_tests_2d.cpp`.

## Why a stored reference rather than a live dependency

Building geogram in CI would add several minutes to all six jobs and pull a large third-party
tree into the build for one test. The stored references give the same signal at no build cost.

They are also not the primary correctness check. For points in **general position** the Delaunay
triangulation is *unique*, so the exact empty-circumcircle self-check in the same test file
(`validate_delaunay`) already proves any passing triangulation is the one geogram would produce.
These files add an independent-implementation cross-check on top of that.

Note this equivalence only holds in general position. For cocircular input the triangulation is
not unique and the choice comes down to the tie-breaking rule; ours is a different SOS
implementation from geogram's (`side3_exact_SOS` vs. our lifted-cofactor perturbation), so a
byte-match on degenerate input is not expected and is deliberately not tested. The degenerate
cases are covered instead by the self-check, which is valid either way.

## Regenerating

    git clone https://github.com/BrunoLevy/geogram.git
    cd geogram && git submodule update --init --recursive
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGEOGRAM_LIB_ONLY=ON \
          -DGEOGRAM_WITH_GRAPHICS=OFF -DGEOGRAM_WITH_LUA=OFF -DGEOGRAM_WITH_HLBFGS=OFF \
          -DGEOGRAM_WITH_TETGEN=OFF -DGEOGRAM_WITH_TRIANGLE=OFF -DGEOGRAM_WITH_EXPLORAGRAM=OFF
    cmake --build build -j

    c++ -std=c++17 -O2 gen_delaunay_ref.cpp \
        -I <geogram>/src/lib -I <geogram>/build/src/lib \
        -L <geogram>/build/lib -lgeogram -o gen_delaunay_ref
    ./gen_delaunay_ref <this directory>

`gen_delaunay_ref.cpp` is kept here for that purpose. It is not part of the VolumeRemesher build.
