#pragma once

// ---------------------------------------------------------------------------
// Compatibility shim.
//
// The exact and indirect predicates (orient2d/orient3D/incircle and their
// implicit-point variants) used to be vendored here. They now live upstream in
// Indirect_Predicates (https://github.com/MarcoAttene/Indirect_Predicates), a
// header-only LGPL-3.0 library by Marco Attene (IMATI-GE / CNR), fetched at
// configure time and pinned to a specific commit (see the root CMakeLists.txt).
//
// This header keeps the historical include path
// <VolumeRemesher/indirect_predicates.h> and re-exports the exact incircle into
// namespace vol_rem (referenced as vol_rem::incircle in the tests).
// ---------------------------------------------------------------------------

#include <indirect_predicates.h> // Indirect_Predicates, global namespace

namespace vol_rem {
using ::incircle;
} // namespace vol_rem
