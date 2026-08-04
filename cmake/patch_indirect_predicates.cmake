# Make Indirect_Predicates' implicit_point.hpp includable from inside a namespace.
#
# Two calls in it are explicitly qualified to the global scope:
#
#     if (num_explicit == 5) return ::inSphere(...);
#
# because genericPoint::inSphere would otherwise find the member and hide the free
# function of the same name. That is correct while the header lives at global scope,
# but VolumeRemesher includes it inside namespace vol_rem (see
# include/VolumeRemesher/numerics.h for why), where ::inSphere does not exist.
#
# Rewriting the qualification to vol_rem::inSphere keeps the member/free distinction --
# which is the whole point of the qualification -- while naming the scope the header is
# actually instantiated in.
#
# Run as a FetchContent PATCH_COMMAND, so it must be idempotent: re-configuring does not
# re-download, and a second application must be a no-op.
#
# Proposed upstream at MarcoAttene/Indirect_Predicates; drop this once that lands.

if(NOT DEFINED IP_FILE)
	message(FATAL_ERROR "patch_indirect_predicates.cmake: IP_FILE not set")
endif()
if(NOT EXISTS "${IP_FILE}")
	message(FATAL_ERROR "patch_indirect_predicates.cmake: ${IP_FILE} does not exist")
endif()

file(READ "${IP_FILE}" _ip_contents)
string(FIND "${_ip_contents}" "return ::inSphere(" _needs_patch)
if(_needs_patch EQUAL -1)
	message(STATUS "Indirect_Predicates: already namespace-safe, nothing to patch")
	return()
endif()

string(REPLACE "return ::inSphere(" "return vol_rem::inSphere(" _ip_patched "${_ip_contents}")
file(WRITE "${IP_FILE}" "${_ip_patched}")
message(STATUS "Indirect_Predicates: qualified inSphere to vol_rem so the header can be included inside the namespace")
