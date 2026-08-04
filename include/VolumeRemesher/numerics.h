#pragma once

// ---------------------------------------------------------------------------
// NFG's number types, defined INSIDE namespace vol_rem.
//
// They live upstream in NFG (Numbers For Geometry, https://github.com/MarcoAttene/nfg),
// a header-only LGPL-3.0 library by Marco Attene (IMATI-GE / CNR), fetched at
// configure time and pinned (see the root CMakeLists.txt).
//
// The upstream header is included *inside* namespace vol_rem rather than at global
// scope with using-declarations. That is deliberate and load-bearing.
//
// These are header-only inline/template definitions, so every library that includes
// them emits its own weak copies and the linker keeps one for the whole binary. This
// library compiles with -mavx2 -mfma (see the root CMakeLists.txt); a downstream that
// links both this and another user of the same upstream headers -- wildmeshing-toolkit
// links FastEnvelope, which does -- compiles its copies without them. -mfma contracts
// a*b+c into a single fused multiply-add with one rounding instead of two, so the two
// copies are NOT numerically identical, and whichever the linker keeps is then used by
// both. That is an ODR violation with observable consequences: it silently changed the
// arrangement's output (three vertices of a 6104-vertex output landed elsewhere) on a
// wildmeshing-toolkit model, while mesh_generator on its own -- with nothing to collide
// with -- stayed byte-identical.
//
// Keeping the definitions in vol_rem gives them distinct mangled names, so this library
// always calls the copy it was compiled with. That is how it behaved before the
// externalization, when these types were vendored inside the namespace.
//
// The system headers upstream needs are included at global scope first; their include
// guards then make the nested includes no-ops, so nothing from the standard library is
// dragged into vol_rem.
// ---------------------------------------------------------------------------

#include <assert.h>
#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <bit>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef USE_GNU_GMP_CLASSES
#include <gmpxx.h>
#endif

// The same SIMD selection upstream makes, hoisted to global scope -- conditions copied
// verbatim from NFG's numerics.h so the two agree. SIMDE_ENABLE_NATIVE_ALIASES is what
// gives SIMDe's headers the __m128d spelling NFG uses.
#if INTPTR_MAX == INT64_MAX
#ifdef __ARM_NEON
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <x86/avx2.h>
#include <x86/fma.h>
#else
#if defined(__SSE2__) || defined(__AVX2__)
#ifdef __AVX2__
#include <immintrin.h>
#else
#include <emmintrin.h>
#endif
#endif
#endif
#endif

namespace vol_rem {

// NFG declares sqrt/fabs overloads for its own number types. Inside a namespace those
// hide the global ones, so a plain sqrt(double) in its own code would find only the
// class overloads and be ambiguous. Pull the global math functions in first: for a
// double argument they are an exact match and win. (The copy this library vendored
// before the externalization had the same problem and solved it by writing std::sqrt.)
using ::fabs;
using ::sqrt;

#include <numerics.h> // NFG
} // namespace vol_rem
