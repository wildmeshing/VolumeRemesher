// Configure-time probe: does this target have SIMD the library can use?
//
// The exact predicates use a SIMD interval_number when available (numerics.h). On x86/x64 that
// is SSE2/AVX2 (<emmintrin.h> etc.); on ARM the baseline is NEON, on top of which numerics.h
// emulates the x86 SSE2/AVX2 intrinsics via SIMDe. This probe must therefore accept NEON too --
// otherwise Apple Silicon (arm64) fails it and CMake wrongly reports "no SIMD" even though the
// AVX2-via-NEON path is what actually gets compiled.
//
// It fails only on a target with neither SSE2 nor NEON, where CMake defines VR_DISABLE_SIMD and
// numerics.h falls back to scalar interval arithmetic.
#include <climits>

#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(_M_ARM64)
// ARM NEON present -> SIMD available (numerics.h uses SIMDe to emulate x86 AVX2/SSE2 on it).
int main() {}
#else
// x86/x64 (and MSVC): require the SSE2 intrinsics header. A target with no NEON and no SSE2
// header fails to compile here, which is the intended "no SIMD" signal.
#include <emmintrin.h>
int main() {}
#endif
