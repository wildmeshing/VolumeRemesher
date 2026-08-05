#pragma once

// ---------------------------------------------------------------------------
// Indirect_Predicates' point classes and fast predicates, defined INSIDE
// namespace vol_rem. See VolumeRemesher/numerics.h for why the upstream headers are
// included inside the namespace rather than re-exported from global scope.
//
// Upstream: https://github.com/MarcoAttene/Indirect_Predicates (header-only, LGPL-3.0,
// Marco Attene, IMATI-GE / CNR), fetched at configure time and pinned.
//
// implicit_point.h pulls in hand_optimized_predicates.hpp, implicit_point.hpp and
// indirect_predicates.h, and includes "numerics.h"; the numerics shim below has already
// placed NFG inside vol_rem, so that nested include resolves within the namespace too.
// ---------------------------------------------------------------------------

#include <VolumeRemesher/numerics.h> // NFG, already inside vol_rem

#include <iostream>

namespace vol_rem {
#include <implicit_point.h> // Indirect_Predicates
} // namespace vol_rem
