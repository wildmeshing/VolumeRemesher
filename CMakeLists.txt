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
    src/VolumeRemesher/delaunay.h
    src/VolumeRemesher/delaunay.cpp
    src/VolumeRemesher/embed.h
    src/VolumeRemesher/embed.cpp
    src/VolumeRemesher/conforming_mesh.h
    src/VolumeRemesher/conforming_mesh.cpp
    src/VolumeRemesher/extended_predicates.h
    src/VolumeRemesher/extended_predicates.cpp
    src/VolumeRemesher/BSP.h
    src/VolumeRemesher/BSP.cpp
    src/VolumeRemesher/inOutPartition.cpp
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
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME}_lib PUBLIC "/D USE_MAROTS_METHOD")
	endif()
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
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME}_lib PUBLIC "-DUSE_MAROTS_METHOD")
	endif()
endif()

# On ARM (e.g. Apple Silicon) the exact predicates rely on SIMDe to emulate the
# x86 AVX2/FMA intrinsics (<x86/avx2.h>, <x86/fma.h>) on top of NEON.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
	include(FetchContent)
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
	target_include_directories(${PROJECT_NAME}_lib PUBLIC ${simde_SOURCE_DIR}/simde)
endif()


# Public include directory
target_include_directories(${PROJECT_NAME}_lib PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_20)

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


# add the executable
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(NOT HAS_PARENT)
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

# Integration tests (only when this is the top-level project)
option(VOLUMEREMESHER_BUILD_TESTS "Build the integration test suite" ON)
if(NOT HAS_PARENT AND VOLUMEREMESHER_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
