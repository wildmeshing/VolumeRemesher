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

# Independent exact-rational verifier for the 2D pipeline's segment provenance. Like
# verify_tracking, it shares no mesh-generation code (only bigrational, via the lib).
add_executable(verify_tracking_2d verify_tracking_2d.cpp tri_orientation.h)
target_link_libraries(verify_tracking_2d PRIVATE VolumeRemesher::VolumeRemesher)
target_compile_features(verify_tracking_2d PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET verify_tracking_2d POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:verify_tracking_2d> $<TARGET_FILE_DIR:verify_tracking_2d>
		COMMAND_EXPAND_LISTS)
endif()

add_executable(integration_tests integration_tests.cpp sha256.h)
target_link_libraries(integration_tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(integration_tests PRIVATE cxx_std_20) # std::filesystem

# The tests invoke the built mesh_generator binary as a subprocess; make sure it
# is built first and hand its path (plus the models dir and the hash manifest)
# to the test via compile definitions.
add_dependencies(integration_tests mesh_generator verify_tracking verify_tracking_2d)
target_compile_definitions(integration_tests PRIVATE
    VRTEST_MESH_GENERATOR="$<TARGET_FILE:mesh_generator>"
    VRTEST_VERIFY_TRACKING="$<TARGET_FILE:verify_tracking>"
    VRTEST_VERIFY_TRACKING_2D="$<TARGET_FILE:verify_tracking_2d>"
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
    VRTEST_HASHES_FILE="${CMAKE_CURRENT_SOURCE_DIR}/reference_hashes.txt"
)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(integration_tests)

# The same suite against the graph-cut classifier (`-c graphcut`), checked against the manifest
# the suite used before the winding-number classifier became the default. The max-flow changed
# from GCO/maxflow to graphcut.h since; the hashes did not.
add_executable(integration_tests_graphcut integration_tests.cpp sha256.h)
target_link_libraries(integration_tests_graphcut PRIVATE Catch2::Catch2WithMain)
target_compile_features(integration_tests_graphcut PRIVATE cxx_std_20)
add_dependencies(integration_tests_graphcut mesh_generator verify_tracking verify_tracking_2d)
target_compile_definitions(integration_tests_graphcut PRIVATE
    VRTEST_MESH_GENERATOR="$<TARGET_FILE:mesh_generator>"
    VRTEST_VERIFY_TRACKING="$<TARGET_FILE:verify_tracking>"
    VRTEST_VERIFY_TRACKING_2D="$<TARGET_FILE:verify_tracking_2d>"
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
    VRTEST_HASHES_FILE="${CMAKE_CURRENT_SOURCE_DIR}/reference_hashes.graphcut.txt"
    VRTEST_EXTRA_FLAGS="-c graphcut"
)
catch_discover_tests(integration_tests_graphcut TEST_PREFIX "graphcut: ")

# The max-flow behind the graph-cut classifier (src/VolumeRemesher/graphcut.h) against an
# exhaustive search for the minimum cut on small random graphs. Header only: no library.
add_executable(graphcut_tests graphcut_tests.cpp)
target_link_libraries(graphcut_tests PRIVATE Catch2::Catch2WithMain)
target_include_directories(graphcut_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_features(graphcut_tests PRIVATE cxx_std_20)
catch_discover_tests(graphcut_tests)

# The winding-number machinery behind the default cell classifier (winding_number_tests.cpp):
# the deterministic atan2 against the system one, the fast winding number against the exact
# sum, and the winding and graph-cut classifiers against each other on the closed, consistently
# oriented models, where they must agree cell by cell.
add_executable(winding_number_tests winding_number_tests.cpp)
target_link_libraries(winding_number_tests PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(winding_number_tests PRIVATE cxx_std_20)
target_compile_definitions(winding_number_tests PRIVATE
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET winding_number_tests POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:winding_number_tests> $<TARGET_FILE_DIR:winding_number_tests>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(winding_number_tests)

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

# Regression for embed_tri_in_poly_mesh emitting inverted tets on a near-planar input.
# Replays a captured call (tests/data/embed_thingi100727.bin) so no mesh generation or
# external tool is involved. NOT run by default: the embedding takes minutes (~4 on an
# M3 Max), so the test case carries Catch2's [.] hidden tag and must be asked for by name.
add_executable(embed_regression embed_regression.cpp tet_orientation.h)
target_link_libraries(embed_regression PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(embed_regression PRIVATE cxx_std_20)
target_compile_definitions(embed_regression PRIVATE
    VRTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET embed_regression POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:embed_regression> $<TARGET_FILE_DIR:embed_regression>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(embed_regression)

# Invariants this library relies on from Indirect_Predicates. Fast (seconds) and
# self-contained: a repin that loses one fails here instead of surfacing as a bad mesh.
add_executable(predicate_invariants predicate_invariants.cpp)
target_link_libraries(predicate_invariants PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(predicate_invariants PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET predicate_invariants POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:predicate_invariants> $<TARGET_FILE_DIR:predicate_invariants>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(predicate_invariants)

# Unit tests for the 3D embedding API. Tiny synthetic inputs (a six-tet cube, a few
# triangles), so unlike embed_regression these run in the default ctest set.
add_executable(unit_tests_3d unit_tests_3d.cpp)
target_link_libraries(unit_tests_3d PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(unit_tests_3d PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET unit_tests_3d POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:unit_tests_3d> $<TARGET_FILE_DIR:unit_tests_3d>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(unit_tests_3d)

# Forced points and edges on a background that already has the surface's vertices, whose
# forcing-triangle corners would leave the background, or at a boundary vertex of the background
# with no inward axis direction (feature_vertex_regression.cpp).
add_executable(feature_vertex_regression feature_vertex_regression.cpp feature_vertex_fixtures.h)
target_link_libraries(feature_vertex_regression PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(feature_vertex_regression PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET feature_vertex_regression POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:feature_vertex_regression> $<TARGET_FILE_DIR:feature_vertex_regression>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(feature_vertex_regression)

# Every boundary edge of the input gets a virtual constraint, the last one in sorted order
# included, so no on-input facet crosses it (virtual_constraint_regression.cpp). One triangle in
# tet meshes of a dozen vertices: milliseconds.
add_executable(virtual_constraint_regression virtual_constraint_regression.cpp)
target_link_libraries(virtual_constraint_regression PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(virtual_constraint_regression PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET virtual_constraint_regression POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:virtual_constraint_regression> $<TARGET_FILE_DIR:virtual_constraint_regression>
		COMMAND_EXPAND_LISTS)
endif()
catch_discover_tests(virtual_constraint_regression)

# Error paths of embedding into a given tet mesh (feature_vertex_errors.cpp). Every case ends in
# ip_error, which prints its message and exits, so each case is registered here by its Catch2
# name together with the message it must print: it passes iff the message appears, and a case
# that runs to completion prints Catch's summary instead and fails. Not through
# catch_discover_tests: its build-time discovery step passes the test properties through cmd.exe
# on Windows, where the `|` of a regex alternation is a pipe.
add_executable(feature_vertex_errors feature_vertex_errors.cpp feature_vertex_fixtures.h)
target_link_libraries(feature_vertex_errors PRIVATE Catch2::Catch2WithMain VolumeRemesher::VolumeRemesher)
target_compile_features(feature_vertex_errors PRIVATE cxx_std_20)
if(WIN32 AND VOLUMEREMESHER_WITH_GMP AND HAS_GMPXX)
	add_custom_command(TARGET feature_vertex_errors POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_RUNTIME_DLLS:feature_vertex_errors> $<TARGET_FILE_DIR:feature_vertex_errors>
		COMMAND_EXPAND_LISTS)
endif()
function(feature_vertex_error_test case message)
	add_test(NAME "feature_vertex_errors: ${case}" COMMAND feature_vertex_errors "${case}")
	set_tests_properties("feature_vertex_errors: ${case}" PROPERTIES PASS_REGULAR_EXPRESSION "${message}")
endfunction()
feature_vertex_error_test("a forced point outside the tet mesh is an error" "lies outside the tet mesh")
feature_vertex_error_test("a forced point on a boundary face of the tet mesh is an error" "lies on the boundary of the tet mesh")
feature_vertex_error_test("a forced edge endpoint on a boundary edge of the tet mesh is an error" "lies on the boundary of the tet mesh")
feature_vertex_error_test("a forced edge with coincident endpoints is an error" "has coincident endpoints")
feature_vertex_error_test("a surface vertex outside the tet mesh is an error" "lies outside the tet mesh")
feature_vertex_error_test("a surface vertex on the boundary of the tet mesh is an error" "lies on the boundary of the tet mesh")
