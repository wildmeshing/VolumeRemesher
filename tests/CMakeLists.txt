# Integration test suite: runs mesh_generator on the models and checks the
# SHA-256 of the generated files against tests/reference_hashes.txt.

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(integration_tests integration_tests.cpp sha256.h)
target_link_libraries(integration_tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(integration_tests PRIVATE cxx_std_20) # std::filesystem

# The tests invoke the built mesh_generator binary as a subprocess; make sure it
# is built first and hand its path (plus the models dir and the hash manifest)
# to the test via compile definitions.
add_dependencies(integration_tests mesh_generator)
target_compile_definitions(integration_tests PRIVATE
    VRTEST_MESH_GENERATOR="$<TARGET_FILE:mesh_generator>"
    VRTEST_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
    VRTEST_HASHES_FILE="${CMAKE_CURRENT_SOURCE_DIR}/reference_hashes.txt"
)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(integration_tests)
