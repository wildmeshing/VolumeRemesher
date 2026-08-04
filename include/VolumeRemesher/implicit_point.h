#pragma once

// ---------------------------------------------------------------------------
// Compatibility shim.
//
// The implicit/explicit point classes and the fast double-precision predicates
// used to be vendored here. They now live upstream in Indirect_Predicates
// (https://github.com/MarcoAttene/Indirect_Predicates), a header-only LGPL-3.0
// library by Marco Attene (IMATI-GE / CNR), fetched at configure time and pinned
// to a specific commit (see the root CMakeLists.txt). Its implicit_point.h pulls
// in hand_optimized_predicates.hpp and implicit_point.hpp, and #includes
// "numerics.h" from NFG (kept on the include path).
//
// This header keeps the historical include path <VolumeRemesher/implicit_point.h>
// and re-exports the global-namespace point types (and the double-precision
// orient predicates) into namespace vol_rem, the public API this library relies
// on (e.g. vol_rem::genericPoint, vol_rem::implicitPoint2D_SSI).
// ---------------------------------------------------------------------------

#include <implicit_point.h>          // Indirect_Predicates, global namespace
#include <VolumeRemesher/numerics.h> // re-exports the vol_rem number types

namespace vol_rem {
using ::genericPoint;

using ::explicitPoint2D;
using ::implicitPoint2D_SSI;

using ::explicitPoint3D;
using ::implicitPoint3D_LPI;
using ::implicitPoint3D_TPI;

// Fast double-precision predicates from hand_optimized_predicates.hpp.
using ::orient2d;
using ::orient3d;
} // namespace vol_rem
