# Integration test suite: runs mesh_generator on the models and checks the
# SHA-256 of the generated files against tests/reference_hashes.txt.

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
)
FetchContent_MakeAvailable(Catch2)

# Independent exact-rational verifier for the -t face-provenance tracking. Shares
# no mesh-generation code with mesh_generator (only bigrational, via the lib).
add_executable(verify_tracking verify_tracking.cpp)
target_link_libraries(verify_tracking PRIVATE VolumeRemesher::VolumeRemesher)
target_compile_features(verify_tracking PRIVATE cxx_std_20)
# When the optional GMP backend is used on Windows, copy its DLL next to
# verify_tracking so it launches (bigrational == mpq_class). See CMakeLists.txt.
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET verify_tracking POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:verify_tracking> $<TARGET_FILE_DIR:verify_tracking>
		COMMAND_EXPAND_LISTS)
endif()

add_executable(integration_tests integration_tests.cpp sha256.h)
target_link_libraries(integration_tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(integration_tests PRIVATE cxx_std_20) # std::filesystem

# The tests invoke the built mesh_generator binary as a subprocess; make sure it
# is built first and hand its path (plus the models dir and the hash manifest)
# to the test via compile definitions.
add_dependencies(integration_tests mesh_generator verify_tracking)
target_compile_definitions(integration_tests PRIVATE
    VRTEST_MESH_GENERATOR="$<TARGET_FILE:mesh_generator>"
    VRTEST_VERIFY_TRACKING="$<TARGET_FILE:verify_tracking>"
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
    VRTEST_HASHES_FILE="${CMAKE_CURRENT_SOURCE_DIR}/reference_hashes.txt"
)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(integration_tests)

# Unit tests for the 2D pipeline. Unlike integration_tests (which shells out to the
# mesh_generator binary and links nothing), these link the library and exercise the internal 2D
# APIs directly -- the lib exports both include/ and src/ as public include dirs.
add_executable(unit_tests_2d
    unit_tests_2d.cpp
)
target_link_libraries(unit_tests_2d PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(unit_tests_2d PRIVATE cxx_std_20)
target_compile_definitions(unit_tests_2d PRIVATE
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
    VRTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET unit_tests_2d POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:unit_tests_2d> $<TARGET_FILE_DIR:unit_tests_2d>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(unit_tests_2d)
