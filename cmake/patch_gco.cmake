# Applied to the GCO 3.0 sources right after FetchContent downloads them (CMakeLists.txt,
# VOLUMEREMESHER_WITH_GRAPHCUT). Run as `cmake -P patch_gco.cmake <source dir>`.
#
# Two changes:
# - block.h allocates its blocks with `new char[...]` and frees them with `delete first`: a
#   mismatched new[]/delete that valgrind reports and that is undefined behaviour. Free them with
#   delete[] on the char array they are.
# - energy.h declares three locals `register`, a storage class C++17 removed; drop the keyword.
# Nothing else is changed: the energy types are selected through GCO's own GCO_ENERGYTYPE /
# GCO_ENERGYTERMTYPE macros (compile definitions).
set(dir "${CMAKE_ARGV3}")
set(block_h "${dir}/block.h")
if(NOT EXISTS "${block_h}")
	message(FATAL_ERROR "patch_gco.cmake: ${block_h} not found")
endif()
file(READ "${block_h}" contents)
string(REPLACE "block *next = first -> next; delete first; first = next;"
	"block *next = first -> next; delete[] (char*)first; first = next;" patched "${contents}")
if(patched STREQUAL contents)
	message(FATAL_ERROR "patch_gco.cmake: block.h did not contain the expected text; the archive changed")
endif()
file(WRITE "${block_h}" "${patched}")

set(energy_h "${dir}/energy.h")
file(READ "${energy_h}" contents)
string(REGEX REPLACE "([\n\t ])register " "\\1" patched "${contents}")
if(patched STREQUAL contents)
	message(FATAL_ERROR "patch_gco.cmake: energy.h did not contain the expected text; the archive changed")
endif()
file(WRITE "${energy_h}" "${patched}")
