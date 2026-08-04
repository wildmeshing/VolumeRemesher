# Make Indirect_Predicates' implicit_point.hpp includable from inside a namespace.
#
# genericPoint::inSphere reaches the free inSphere() by qualifying the call as
# ::inSphere, because the member of the same name would otherwise hide it. That is
# correct at global scope but unresolvable when the header is included inside a
# namespace, which VolumeRemesher does (see include/VolumeRemesher/numerics.h).
#
# The rewrite must keep working at global scope too: this same fetched header is shared
# with other consumers in a downstream build -- wildmeshing-toolkit's FastEnvelope
# includes it globally -- so it cannot be rewritten to name vol_rem. Instead add
# inSphere_EEEEE, the all-explicit member of the existing inSphere_* family (_IEEEE,
# _IIEEE, _IIIEE, _IIIIE) which was simply missing, and call that. A distinct name is not
# hidden by the member, so it needs no qualification and resolves in either scope.
#
# Idempotent: FetchContent does not re-download on reconfigure, so a second application
# must be a no-op.
#
# Proposed upstream as MarcoAttene/Indirect_Predicates#13; drop this once that lands.

if(NOT DEFINED IP_FILE)
	message(FATAL_ERROR "patch_indirect_predicates.cmake: IP_FILE not set")
endif()
if(NOT EXISTS "${IP_FILE}")
	message(FATAL_ERROR "patch_indirect_predicates.cmake: ${IP_FILE} does not exist")
endif()

file(READ "${IP_FILE}" _ip)
string(FIND "${_ip}" "inSphere_EEEEE" _already)
if(NOT _already EQUAL -1)
	message(STATUS "Indirect_Predicates: already namespace-safe, nothing to patch")
	return()
endif()

set(_helper
"inline int inSphere_EEEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere(
		a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(),
		b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IEEEE(const genericPoint& a,")

string(REPLACE "inline int inSphere_IEEEE(const genericPoint& a," "${_helper}" _ip "${_ip}")

# Rewrite both ::inSphere call sites. They span several lines, so match the whole
# argument list non-greedily up to the terminating ");".
string(REGEX REPLACE
	"if \\(num_explicit == 5\\) return ::inSphere\\([^;]*\\);"
	"if (num_explicit == 5) return inSphere_EEEEE(a, b, c, d, e);"
	_ip "${_ip}")

if(_ip MATCHES "::inSphere\\(")
	message(FATAL_ERROR "patch_indirect_predicates.cmake: a ::inSphere call site was not rewritten")
endif()

file(WRITE "${IP_FILE}" "${_ip}")
message(STATUS "Indirect_Predicates: added inSphere_EEEEE so the header resolves at global scope and inside a namespace")
