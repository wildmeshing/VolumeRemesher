#pragma once

// ---------------------------------------------------------------------------
// Compatibility shim.
//
// The number types used to be vendored here verbatim. They now live upstream in
// NFG (Numbers For Geometry, https://github.com/MarcoAttene/nfg), a header-only
// LGPL-3.0 library by Marco Attene (IMATI-GE / CNR), fetched at configure time
// and pinned to a specific commit (see the root CMakeLists.txt).
//
// This header keeps the historical include path <VolumeRemesher/numerics.h> and,
// crucially, re-exports NFG's global-namespace number types into namespace
// vol_rem, which is the public API this library and its downstreams rely on
// (e.g. vol_rem::bigrational, referenced across wildmeshing-toolkit).
// ---------------------------------------------------------------------------

#include <numerics.h> // NFG, global namespace

namespace vol_rem {
using ::bigfloat;
using ::bignatural;
using ::bigrational;
using ::interval_number;
} // namespace vol_rem
