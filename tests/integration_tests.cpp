// Integration tests: run mesh_generator on each model listed in
// reference_hashes.txt and check the SHA-256 of the generated file against a
// stored reference. Because the exact predicates make the pipeline fully
// deterministic, these hashes must be identical on every OS/compiler -- so a
// mismatch means a portability/nondeterminism regression.
//
// The reference hashes are produced by tests/compute_hashes.py (which uses the
// same SHA-256), so adding a new model is: drop the .off in models/, run the
// python helper, paste the line into reference_hashes.txt.
//
// Models listed in the manifest but not present on disk are skipped (with a
// warning) rather than failing -- this lets CI run whatever subset of models is
// committed to the repo.

#include <catch2/catch_test_macros.hpp>

#include "sha256.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Entry {
    std::string file;   // e.g. "sphere.off"
    std::string mode;   // "default" | "blackfaces" | "skin"
    std::string sha256; // expected 64-hex digest of the generated file
};

std::vector<Entry> load_manifest() {
    std::vector<Entry> entries;
    std::ifstream in(VRTEST_HASHES_FILE);
    std::string line;
    while (std::getline(in, line)) {
        const auto p = line.find_first_not_of(" \t\r\n");
        if (p == std::string::npos || line[p] == '#')
            continue;
        std::istringstream ss(line);
        Entry e;
        ss >> e.file >> e.mode >> e.sha256;
        if (!e.file.empty() && !e.mode.empty() && !e.sha256.empty())
            entries.push_back(e);
    }
    return entries;
}

std::string flags_for(const std::string& mode) {
    if (mode == "blackfaces") return "-b";
    if (mode == "skin") return "-s";
    return ""; // default -> volume.msh
}

std::string output_for(const std::string& mode) {
    if (mode == "blackfaces") return "black_faces.off";
    if (mode == "skin") return "skin.off";
    return "volume.msh";
}

// Run mesh_generator on 'model' in a private temp dir (the tool writes its
// output to the current directory with a fixed name) and return the SHA-256 of
// the generated file. Empty string if no output was produced.
std::string run_and_hash(const std::string& model, const std::string& mode) {
    static int counter = 0;
    fs::path work = fs::temp_directory_path() / ("vrtest_" + std::to_string(counter++));
    std::error_code ec;
    fs::remove_all(work, ec);
    fs::create_directories(work, ec);

    const std::string bin = VRTEST_MESH_GENERATOR;
    const std::string flags = flags_for(mode);
    std::string cmd;
#ifdef _WIN32
    cmd = "cd /d \"" + work.string() + "\" && \"" + bin + "\" \"" + model + "\" " + flags + " > nul 2>&1";
#else
    cmd = "cd '" + work.string() + "' && '" + bin + "' '" + model + "' " + flags + " > /dev/null 2>&1";
#endif
    std::system(cmd.c_str());

    const fs::path out = work / output_for(mode);
    std::string digest = fs::exists(out) ? vrtest::sha256_file(out.string()) : std::string();
    fs::remove_all(work, ec);
    return digest;
}

} // namespace

TEST_CASE("integration: model outputs are byte-stable across platforms", "[integration]") {
    const auto entries = load_manifest();
    REQUIRE_FALSE(entries.empty()); // reference_hashes.txt must not be empty

    for (const auto& e : entries) {
        DYNAMIC_SECTION(e.file << " [" << e.mode << "]") {
            const fs::path model = fs::path(VRTEST_MODELS_DIR) / e.file;
            if (!fs::exists(model)) {
                WARN("model not present, skipping: " << e.file
                     << " (commit it under models/ to exercise it in CI)");
            } else {
                const std::string got = run_and_hash(model.string(), e.mode);
                INFO("model:    " << e.file);
                INFO("mode:     " << e.mode);
                INFO("expected: " << e.sha256);
                INFO("actual:   " << (got.empty() ? std::string("<no output produced>") : got));
                REQUIRE(got == e.sha256);
            }
        }
    }
}
