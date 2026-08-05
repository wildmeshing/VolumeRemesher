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

# Insert a forwarder taking the same 15 doubles the call sites already pass, so only the
# callee name changes and the (multi-line) argument lists are left exactly as they are.
set(_helper
"inline int inSphere_EEEEE(const double& pax, const double& pay, const double& paz, const double& pbx, const double& pby, const double& pbz, const double& pcx, const double& pcy, const double& pcz, const double& pdx, const double& pdy, const double& pdz, const double& pex, const double& pey, const double& pez) {
	return inSphere(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}

inline int inSphere_IEEEE(const genericPoint& a,")
string(REPLACE "inline int inSphere_IEEEE(const genericPoint& a," "${_helper}" _ip "${_ip}")

string(REPLACE "return ::inSphere(" "return inSphere_EEEEE(" _ip "${_ip}")

# Only the CALL sites, not the "genericPoint::inSphere(" definitions, which also
# contain the "::inSphere(" spelling.
string(FIND "${_ip}" "return ::inSphere(" _leftover)
if(NOT _leftover EQUAL -1)
	message(FATAL_ERROR "patch_indirect_predicates.cmake: a ::inSphere call site was not rewritten")
endif()
string(FIND "${_ip}" "inSphere_EEEEE" _ok)
if(_ok EQUAL -1)
	message(FATAL_ERROR "patch_indirect_predicates.cmake: helper was not inserted")
endif()

file(WRITE "${IP_FILE}" "${_ip}")
message(STATUS "Indirect_Predicates: added inSphere_EEEEE so the header resolves at global scope and inside a namespace")
