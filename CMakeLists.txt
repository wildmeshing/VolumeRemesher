cmake_minimum_required(VERSION 3.25)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
option(VOLREM_ENABLE_SSE2 "Enable SSE2 instructions" ON)
option(VOLREM_ENABLE_AVX2 "Enable AVX2 instructions" ON)

#####################

# set the project name
project(mesh_generator)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (SOURCES
    src/VolumeRemesher/makePolyhedralMesh.cpp
    src/VolumeRemesher/embed.h
    src/VolumeRemesher/embed.cpp
    src/VolumeRemesher/conforming_mesh.h
    src/VolumeRemesher/conforming_mesh.cpp
    src/VolumeRemesher/BSP.h
    src/VolumeRemesher/BSP.cpp
    src/VolumeRemesher/inOutPartition.cpp
    src/VolumeRemesher/winding_number.h
    src/VolumeRemesher/winding_number.cpp
    src/VolumeRemesher/parallel.h
    src/VolumeRemesher/parallel.cpp
    src/VolumeRemesher/exact_coords.h
    src/VolumeRemesher/io.h
    src/VolumeRemesher/io.cpp
    src/VolumeRemesher/mesh_output.h
    src/VolumeRemesher/mesh_output.cpp
    src/VolumeRemesher/driver.h
    src/VolumeRemesher/driver.cpp
    # 2D pipeline: segment arrangement -> triangles
    src/VolumeRemesher/2d/predicates2d.h
    src/VolumeRemesher/2d/predicates2d.cpp
    src/VolumeRemesher/2d/delaunay2d.h
    src/VolumeRemesher/2d/delaunay2d.cpp
    src/VolumeRemesher/2d/arrangement2d.h
    src/VolumeRemesher/2d/arrangement2d.cpp
    src/VolumeRemesher/2d/io2d.h
    src/VolumeRemesher/2d/io2d.cpp
    src/VolumeRemesher/2d/embed2d.h
    src/VolumeRemesher/2d/embed2d.cpp
)

# Does this target have SIMD (x86 SSE2/AVX2, or ARM NEON on which numerics.h emulates them via
# SIMDe)? cmake/test_simd.cpp accepts both; it fails only when neither is present.
try_compile(HAS_SIMD_INSTRUCTIONS ${CMAKE_BINARY_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_simd.cpp")


# make a static library for use within external programs
add_library(${PROJECT_NAME}_lib
    ${SOURCES}
)


if(NOT HAS_SIMD_INSTRUCTIONS)
	message(WARNING "No SIMD (SSE2/NEON) available on this target; using the scalar (slower) "
		"interval arithmetic. Output is unchanged.")
	target_compile_definitions(${PROJECT_NAME}_lib PUBLIC VR_DISABLE_SIMD)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME}_lib PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${PROJECT_NAME}_lib PUBLIC "/Oi")
# turn off annoying warnings
	target_compile_options(${PROJECT_NAME}_lib PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE "_CRT_NONSTDC_NO_WARNINGS")
	# workaround for known bug in C++17 on MSVC
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE "_HAS_STD_BYTE=0")
	# set target architecture: /arch:AVX2 selects the SIMD interval_number path
	# (as -mavx2 does on Linux and the NEON/SIMDe path on Apple Silicon). Verified
	# to produce byte-identical, deterministic output vs the scalar path and vs
	# Linux/macOS (the integration test suite passes on all three). /fp:strict +
	# fenv_access keep the SIMD interval arithmetic rounding outward under
	# FE_UPWARD. (This was disabled while an unrelated half-edge ordering
	# nondeterminism made Windows diverge; that bug is now fixed in conforming_mesh.cpp.)
	if(VOLREM_ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME}_lib PUBLIC "/arch:AVX2")
	endif()
# reserve enough stack size
	target_link_options(${PROJECT_NAME}_lib PUBLIC "/STACK:8421376")
else()
# set standard optimization level
	target_compile_options(${PROJECT_NAME}_lib PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME}_lib PUBLIC -frounding-math)
# the exact predicates rely on error-free transformations (e.g. Dekker
# two-product): forbid the compiler from contracting a*b+c into an FMA, which
# rounds once instead of twice and silently breaks those identities. This
# mirrors the /fp:strict used on the MSVC path.
	target_compile_options(${PROJECT_NAME}_lib PUBLIC -ffp-contract=off)
	# set target architecture
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
		# ARM (e.g. Apple Silicon): the x86 -mavx2/-mfma/-msse2 flags are
		# invalid here. The code emulates the x86 AVX2/FMA intrinsics on top
		# of NEON via SIMDe, which is fetched below.
	else()
		if(VOLREM_ENABLE_AVX2)
			target_compile_options(${PROJECT_NAME}_lib PUBLIC "-mavx2")
			target_compile_options(${PROJECT_NAME}_lib PUBLIC "-mfma")
		elseif(VOLREM_ENABLE_SSE2)
			target_compile_options(${PROJECT_NAME}_lib PUBLIC "-msse2")
		endif()
	endif()
# reserve enough stack size (GNU ld only; the macOS linker uses the 8 MB default)
	if(NOT APPLE)
		target_compile_options(${PROJECT_NAME}_lib PUBLIC -Wl,-z,stacksize=8421376)
	endif()
endif()

# --- Optional graph-cut classifier (NOT LGPL) ---------------------------------
# The historical inside/outside classification is a min-cut over the dual graph of the cells,
# using GCO 3.0 (O. Veksler, A. Delong) and the maxflow library (Y. Boykov, V. Kolmogorov) it
# bundles. GCO's license allows research use only (and cites a US patent on the algorithm) and
# reserves redistribution to its author; maxflow is GPL. Neither is compatible with this
# library's LGPL, so this repository does not contain that code: when this option is ON the
# archive is downloaded from the author's site at configure time (pinned by hash, one fix
# applied by cmake/patch_gco.cmake), compiled into the library, and `mesh_generator -c
# graphcut` (Python classifier="graphcut") becomes available. The resulting binary is then NOT
# distributable under the LGPL, and by downloading the archive you accept GCO's terms
# (GCO_README.TXT in it). The default classifier (winding number, inOutPartition.cpp) has no
# such restriction, and the Python wheels and sdist are built without this option.
option(VOLUMEREMESHER_WITH_GRAPHCUT
	"Download and compile the graph-cut cell classifier (GCO/maxflow: research-only/GPL code)" OFF)
if(VOLUMEREMESHER_WITH_GRAPHCUT)
	message(WARNING "VOLUMEREMESHER_WITH_GRAPHCUT is ON: GCO 3.0 and maxflow are downloaded and "
		"compiled in; their licenses (GCO: research use only; maxflow: GPL) are not compatible "
		"with the LGPL. Do not redistribute this binary under the LGPL.")
	include(FetchContent)
	FetchContent_Declare(gco
		URL https://vision.cs.uwaterloo.ca/files/gco-v3.0.zip
		URL_HASH SHA256=7562be2fa7ef19404a24a29d15414b949bdcab785370479ecc3a57a8f1cfb538
		PATCH_COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_gco.cmake" <SOURCE_DIR>
		SOURCE_SUBDIR do-not-configure
	)
	FetchContent_MakeAvailable(gco)
	message(STATUS "GCO is available in " ${gco_SOURCE_DIR})
	target_sources(${PROJECT_NAME}_lib PRIVATE src/VolumeRemesher/inOutPartition_graphcut.cpp)
	target_include_directories(${PROJECT_NAME}_lib PRIVATE ${gco_SOURCE_DIR})
	# The cell costs and face areas are doubles; GCO defaults to integer energies.
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE GCO_ENERGYTYPE=double GCO_ENERGYTERMTYPE=double)
	target_compile_definitions(${PROJECT_NAME}_lib PUBLIC VOLUMEREMESHER_HAVE_GRAPHCUT)
endif()

# --- Exact-arithmetic dependencies -------------------------------------------
# The number types (numerics.h) and exact/indirect predicates used to be vendored
# under include/VolumeRemesher/. They now live upstream as two header-only
# libraries and are fetched here, pinned to a specific commit for reproducible,
# byte-identical output:
#   NFG (Numbers For Geometry): intervals, expansions, bignatural, bigfloat,
#       bigrational -- i.e. numerics.h.
#   Indirect_Predicates: indirect_predicates.h, implicit_point.{h,hpp},
#       hand_optimized_predicates.hpp. Its implicit_point.h #includes "numerics.h",
#       so NFG must be on the include path.
# The thin shims under include/VolumeRemesher/ pull these in and re-export the
# public types into namespace vol_rem (the API this library and its downstreams
# depend on, e.g. vol_rem::bigrational).
include(FetchContent)

# Local-checkout override for upstream development. If a file named .local.cmake
# exists one level above this repository, NFG and Indirect_Predicates are taken
# from sibling working copies (../NFG and ../Indirect_Predicates) instead of being
# fetched. This lets a developer who has all three repos checked out side by side
# (as Marco Attene does) edit the number/predicate libraries and rebuild here
# immediately, with no re-fetch. The marker file is untracked (see .gitignore), so
# ordinary builds always take the FetchContent path below and are unaffected.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.local.cmake")
	set(nfg_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../NFG")
	set(indirect_predicates_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Indirect_Predicates")
	set(delaunay3d_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Delaunay3D")
	message(STATUS "VolumeRemesher: .local.cmake present -- using local NFG (${nfg_SOURCE_DIR}) "
		"and Indirect_Predicates (${indirect_predicates_SOURCE_DIR}) "
		"and Delaunay3D (${delaunay3d_SOURCE_DIR})")
else()
	# Pinned to exact commits for reproducible, byte-identical output across platforms
	# (moving branch heads would let different CI runners fetch different sources).
	# NFG: upstream main aeec975. Among what it carries: MarcoAttene/NFG#4 counts _M_X64 as
	# SSE2, so on MSVC x64 every translation unit gets the same (SIMD) interval_number layout
	# whether or not it is built with /arch:AVX2 (downstreams that link this library without
	# that flag share NFG:: with it, and a layout mismatch swaps interval bounds across the
	# link); #5 adds bignatural::trimMemoryPool(), which lets a caller return the memory the
	# thread-local bignatural pool grew into during an arrangement; #6 makes ip_error exit with
	# a failure status instead of 0, so ctest sees a test that trips it as failed; and #7 makes
	# it throw a std::runtime_error instead of exiting when NFG_IP_ERROR_THROWS is defined,
	# which the Python module needs so a bad input becomes an exception rather than the end of
	# the interpreter. Same numerics.h as the wildmeshing fork commit pinned before (2330267).
	# fast-envelope pins NFG too: see its cmake/recipes/ipred.cmake.
	FetchContent_Declare(nfg
		GIT_REPOSITORY https://github.com/MarcoAttene/NFG.git
		GIT_TAG aeec97547a67dc6904ab84b39b15d20c6d62f65a
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_Declare(indirect_predicates
		GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
		GIT_TAG a5cb2c4e6a2e8ae137f23ffcee91d113e8f3a1c0
		SOURCE_SUBDIR do-not-configure
	)

	# Delaunay3D: upstream main 362502e. Besides the 3D tetrahedrization itself:
	#   - 2d0279a inserts the vertices of tetrahedrize() in Morton order above 1150 of them, and
	#     #3 breaks ties between points in the same grid cell by vertex index (std::sort would
	#     leave them in a library-dependent order, and the tets in a different order on each
	#     platform). That changes the order tetrahedrize() creates tets in, not which tets are
	#     created: symbolic perturbation makes the set order-independent. embed_tri_in_poly_mesh
	#     never calls tetrahedrize() -- it inserts into the caller's tet mesh.
	#   - #2 adds searchTetrahedron for a point that is not a vertex (locate_point relies on it).
	#     Inserting a point at an existing vertex's position is unsupported; remakePolyhedralMesh
	#     checks for that itself.
	#   - 60d00a8: tetrahedrize() returns no tets, instead of stopping in ip_error, when the
	#     vertices do not span a volume; makePolyhedralMesh reports that itself.
	FetchContent_Declare(Delaunay3D
		GIT_REPOSITORY https://github.com/MarcoAttene/Delaunay3D.git
		GIT_TAG 362502e63abb7adb3bf5a7a3b1e1a303640f8812
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_MakeAvailable(nfg indirect_predicates delaunay3d)
	message(STATUS "NFG is available in " ${nfg_SOURCE_DIR})
	message(STATUS "Indirect_Predicates is available in " ${indirect_predicates_SOURCE_DIR})
	message(STATUS "Delaunay3D is available in " ${delaunay3d_SOURCE_DIR})
endif()

# On ARM (e.g. Apple Silicon) the exact predicates rely on SIMDe to emulate the
# x86 AVX2/FMA intrinsics (<x86/avx2.h>, <x86/fma.h>) on top of NEON. NFG's
# numerics.h includes <x86/avx2.h>/<x86/fma.h> directly under __ARM_NEON.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
	target_include_directories(${PROJECT_NAME}_lib PUBLIC ${simde_SOURCE_DIR}/simde)
endif()


# Public include directory. Order: our headers (the shims) first, then NFG and
# Indirect_Predicates. The shims are addressed as <VolumeRemesher/...> and the
# upstream headers as bare <numerics.h>/<implicit_point.h>, so the two never clash.
target_include_directories(${PROJECT_NAME}_lib PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
	${delaunay3d_SOURCE_DIR}/include
)
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_20)

# Parallelize the tetrahedralization (makeTetrahedra: the per-cell decomposition and the
# winding-flip) and the exact output coordinates of embed_tri_in_poly_mesh and
# embed_seg_in_tri_mesh with std::thread. On by default; turn OFF for a fully serial build.
# The output is byte-identical either way.
option(VOLUMEREMESHER_PARALLEL_TETRAHEDRALIZATION
	"Parallelize makeTetrahedra with std::thread" ON)
if(VOLUMEREMESHER_PARALLEL_TETRAHEDRALIZATION)
	find_package(Threads REQUIRED) # std::thread needs the platform threads library (-pthread)
	target_link_libraries(${PROJECT_NAME}_lib PUBLIC Threads::Threads)
else()
	target_compile_definitions(${PROJECT_NAME}_lib PUBLIC VOLUMEREMESHER_SERIAL_TET)
endif()

# Exact-arithmetic backend. By default the library uses its built-in arbitrary-
# precision bignum (the hand-written `bigfloat`): the exact predicates only add,
# subtract and multiply -- they never divide -- so they never pay for GMP
# mpq_class's per-operation GCD canonicalization. That makes the built-in type
# ~2-2.4x faster for mesh generation (measured across the model set, up to a 1.6M
# triangle mesh) with byte-identical output. GMP is only faster for the rational-
# heavy verify_tracking test tool. Turn this option ON to use GMP (mpq_class)
# instead; it pulls in the GMP dependency (downloaded on Windows, plus its runtime
# DLL) and enables USE_GNU_GMP_CLASSES.
option(VOLUMEREMESHER_WITH_GMP "Use GMP (mpq_class) for exact arithmetic instead of the built-in bignum" OFF)

if(VOLUMEREMESHER_WITH_GMP)
	include(gmp)
	# Enable the GMP classes only if the C++ wrapper (gmpxx.h) actually compiles and
	# links here; otherwise warn and fall back to the built-in bignum rather than break.
	try_compile(HAS_GMPXX ${CMAKE_BINARY_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_gmpxx.cpp" LINK_LIBRARIES gmp::gmp)
	if(HAS_GMPXX)
		target_link_libraries(${PROJECT_NAME}_lib PUBLIC gmp::gmp)
		target_compile_definitions(${PROJECT_NAME}_lib PUBLIC USE_GNU_GMP_CLASSES)
	else()
		message(WARNING "VOLUMEREMESHER_WITH_GMP is ON but gmpxx is not usable here; using the built-in bignum")
	endif()
endif()


add_library(VolumeRemesher::VolumeRemesher ALIAS ${PROJECT_NAME}_lib)

# --- Python module ------------------------------------------------------------
# volume_remesher._core (nanobind), see python/. pip builds it through scikit-build-core, which
# turns this on and the CLI/tests off (pyproject.toml); a plain CMake configure with it ON also
# builds it, into python/volume_remesher/, so `PYTHONPATH=python pytest` works from the tree.
option(VOLUMEREMESHER_BUILD_PYTHON "Build the volume_remesher Python module" OFF)
if(VOLUMEREMESHER_BUILD_PYTHON)
	# The static library is linked into a shared module.
	set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
	# A bad input must become a Python exception, not the end of the interpreter: ip_error
	# throws std::runtime_error instead of exiting. PUBLIC so every translation unit sees the
	# same inline function (mesh_generator, built in the same tree, catches it in main()).
	target_compile_definitions(${PROJECT_NAME}_lib PUBLIC NFG_IP_ERROR_THROWS)
endif()


# add the executable
option(VOLUMEREMESHER_BUILD_CLI "Build the mesh_generator executable" ON)
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(NOT HAS_PARENT AND VOLUMEREMESHER_BUILD_CLI)
	add_executable(${PROJECT_NAME}
    	src/VolumeRemesher/main.cpp
	)
	target_link_libraries(${PROJECT_NAME} PUBLIC VolumeRemesher::VolumeRemesher)
	target_compile_definitions(${PROJECT_NAME} PRIVATE VOL_REM_MODELS_PATH=\"${CMAKE_SOURCE_DIR}/models\")
	# When the optional GMP backend is used on Windows it is a downloaded DLL the
	# executable depends on at runtime, so copy it next to the .exe (its own directory
	# is the first place Windows looks). Skipped entirely without GMP.
	if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
		add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				$<TARGET_RUNTIME_DLLS:${PROJECT_NAME}> $<TARGET_FILE_DIR:${PROJECT_NAME}>
			COMMAND_EXPAND_LISTS)
	endif()
endif()

# Integration tests (only when this is the top-level project; they run mesh_generator)
option(VOLUMEREMESHER_BUILD_TESTS "Build the integration test suite" ON)
if(NOT HAS_PARENT AND VOLUMEREMESHER_BUILD_TESTS AND VOLUMEREMESHER_BUILD_CLI)
	enable_testing()
	add_subdirectory(tests)
endif()

if(VOLUMEREMESHER_BUILD_PYTHON)
	add_subdirectory(python)
endif()
