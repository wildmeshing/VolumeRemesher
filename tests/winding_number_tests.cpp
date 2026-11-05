// Unit tests for the winding-number machinery behind the default cell classifier
// (src/VolumeRemesher/winding_number.{h,cpp}, inOutPartition.cpp):
//   - det_atan2 / det_atan agree with the system libm to within a couple of ulps everywhere
//     (they are what makes the generalized-winding-number fallback bit-identical across
//     platforms, so they must also be correct);
//   - the fast winding number is within its published accuracy of the exact sum, and the exact
//     sum is the winding number (1 inside a closed surface, 0 outside);
//   - the winding and graph-cut classifiers label the cells of the closed, consistently oriented
//     models identically.

#include <catch2/catch_test_macros.hpp>

#include "VolumeRemesher/BSP.h"
#include "VolumeRemesher/io.h"
#include "VolumeRemesher/winding_number.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace vol_rem;

namespace {

std::string model(const char* name)
{
    return std::string(VRTEST_MODELS_DIR) + "/" + name;
}

// Distance in ulps between two doubles of the same sign (large when they differ in sign).
uint64_t ulp_distance(double a, double b)
{
    if (a == b) return 0;
    if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<uint64_t>::max();
    if ((a < 0) != (b < 0)) return std::numeric_limits<uint64_t>::max();
    uint64_t ua, ub;
    std::memcpy(&ua, &a, sizeof ua);
    std::memcpy(&ub, &b, sizeof ub);
    return ua > ub ? ua - ub : ub - ua;
}

} // namespace

TEST_CASE("det_atan2 matches the system atan2 to within 2 ulps", "[winding]")
{
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> mag(-30.0, 30.0); // exponents: 1e-30 .. 1e30
    std::uniform_real_distribution<double> sign01(0.0, 1.0);
    uint64_t worst = 0;
    for (int i = 0; i < 200000; i++) {
        const double y = (sign01(rng) < 0.5 ? -1 : 1) * std::pow(10.0, mag(rng)) * sign01(rng);
        const double x = (sign01(rng) < 0.5 ? -1 : 1) * std::pow(10.0, mag(rng)) * sign01(rng);
        const double d = det_atan2(y, x), s = std::atan2(y, x);
        const uint64_t u = ulp_distance(d, s);
        worst = std::max(worst, u);
        if (u > 2) {
            INFO("atan2(" << y << ", " << x << "): det=" << d << " libm=" << s);
            REQUIRE(u <= 2);
        }
    }
    INFO("worst ulp distance " << worst);
    CHECK(worst <= 2);

    // Special values, as fdlibm defines them
    CHECK(det_atan2(0.0, 1.0) == 0.0);
    CHECK(det_atan2(1.0, 0.0) == std::atan2(1.0, 0.0));
    CHECK(det_atan2(-1.0, 0.0) == std::atan2(-1.0, 0.0));
    CHECK(det_atan2(0.0, -1.0) == std::atan2(0.0, -1.0));
    CHECK(det_atan2(-0.0, -1.0) == std::atan2(-0.0, -1.0));
    CHECK(ulp_distance(det_atan2(1.0, 1.0), std::atan2(1.0, 1.0)) <= 1);
    CHECK(std::isnan(det_atan2(std::nan(""), 1.0)));
    CHECK(
        det_atan2(std::numeric_limits<double>::infinity(), 1.0) ==
        std::atan2(std::numeric_limits<double>::infinity(), 1.0));
    CHECK(det_atan(0.0) == 0.0);
    CHECK(ulp_distance(det_atan(1.0), std::atan(1.0)) <= 1);
    CHECK(ulp_distance(det_atan(1e-10), std::atan(1e-10)) <= 1);
    CHECK(ulp_distance(det_atan(1e300), std::atan(1e300)) <= 1);
}

TEST_CASE("det_atan2 is a pure function of its bits", "[winding]")
{
    // The same inputs give the same bits on every call (no state), and the function is odd in y.
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-5.0, 5.0);
    for (int i = 0; i < 1000; i++) {
        const double y = u(rng), x = u(rng);
        CHECK(det_atan2(y, x) == det_atan2(y, x));
        if (y != 0.0) CHECK(det_atan2(-y, x) == -det_atan2(y, x));
    }
}

TEST_CASE("the winding number of a closed surface is 1 inside and 0 outside", "[winding]")
{
    const TriMeshData sphere = read_OFF(model("sphere.off").c_str(), false);
    const FastWindingNumber fwn(sphere.coords.data(), sphere.tris.data(), sphere.num_triangles());

    // Bounding box and center
    double lo[3], hi[3], c[3];
    for (int k = 0; k < 3; k++) lo[k] = hi[k] = sphere.coords[k];
    for (uint32_t v = 0; v < sphere.num_vertices(); v++)
        for (int k = 0; k < 3; k++) {
            lo[k] = std::min(lo[k], sphere.coords[3 * v + k]);
            hi[k] = std::max(hi[k], sphere.coords[3 * v + k]);
        }
    for (int k = 0; k < 3; k++) c[k] = 0.5 * (lo[k] + hi[k]);
    const double r = 0.5 * (hi[0] - lo[0]);

    // Deep inside, far outside, exact and fast. The sphere is oriented outward (or inward: the
    // winding number is then -1 inside), so compare magnitudes.
    const double w_center = fwn.exact(c);
    CHECK(std::fabs(std::fabs(w_center) - 1.0) < 1e-6);
    CHECK(std::fabs(fwn(c) - w_center) < 2e-2);
    const double far[3] = {c[0] + 10 * r, c[1] + 10 * r, c[2] + 10 * r};
    CHECK(std::fabs(fwn.exact(far)) < 1e-6);
    CHECK(std::fabs(fwn(far)) < 1e-3);

    // Random points: fast vs exact within the expansion's accuracy on this coarse sphere (760
    // triangles, so the accepted far-field nodes are large patches; measured below 1e-2), and
    // the inside/outside decision identical wherever the exact value is not near 1/2.
    std::mt19937_64 rng(99);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    int checked = 0;
    for (int i = 0; i < 2000; i++) {
        const double q[3] = {c[0] + u(rng) * r, c[1] + u(rng) * r, c[2] + u(rng) * r};
        const double e = fwn.exact(q), f = fwn(q);
        INFO("q=(" << q[0] << "," << q[1] << "," << q[2] << ") exact=" << e << " fast=" << f);
        CHECK(std::fabs(e - f) < 2e-2);
        if (std::fabs(std::fabs(e) - 0.5) > 0.05) {
            CHECK((std::fabs(e) > 0.5) == (std::fabs(f) > 0.5));
            checked++;
        }
    }
    CHECK(checked > 1500);
}

TEST_CASE("the fast winding number is deterministic", "[winding]")
{
    const TriMeshData bunny = read_OFF(model("bunny.off").c_str(), false);
    const FastWindingNumber a(bunny.coords.data(), bunny.tris.data(), bunny.num_triangles());
    const FastWindingNumber b(bunny.coords.data(), bunny.tris.data(), bunny.num_triangles());
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> u(-0.2, 0.2);
    for (int i = 0; i < 200; i++) {
        const double q[3] = {u(rng), u(rng), u(rng)};
        CHECK(a(q) == b(q)); // bitwise
    }
}

namespace {
std::vector<uint32_t> cell_places(const char* name, CellClassifier classifier)
{
    TriMeshData A = read_OFF(model(name).c_str(), false);
    BSPcomplex* cpx = makePolyhedralMesh(
        name,
        A.coords.data(),
        A.num_vertices(),
        A.tris.data(),
        A.num_triangles(),
        NULL,
        NULL,
        0,
        NULL,
        0,
        '0',
        false,
        false,
        false,
        false,
        false,
        nullptr,
        classifier);
    std::vector<uint32_t> places;
    places.reserve(cpx->cells.size());
    for (const BSPcell& c : cpx->cells) places.push_back(c.place);
    delete cpx;
    return places;
}
} // namespace

TEST_CASE(
    "winding and graph-cut classifiers agree on closed, consistently oriented input",
    "[winding][graphcut]")
{
    // Every model whose input surface is closed and consistently oriented (so the exact winding
    // number applies). Both classifiers must then find the same solid, cell by cell.
    for (const char* name :
         {"sphere.off",
          "upsample_box.off",
          "100071_sf.off",
          "Octocat.off",
          "simwild_double_sphere_notop.off",
          "cube_subdiv.off",
          "cube_on_cube.off"}) {
        INFO(name);
        const std::vector<uint32_t> gc = cell_places(name, CellClassifier::GRAPH_CUT);
        const std::vector<uint32_t> wn = cell_places(name, CellClassifier::WINDING_NUMBER);
        REQUIRE(gc.size() == wn.size());
        size_t diff = 0;
        for (size_t i = 0; i < gc.size(); i++) diff += gc[i] != wn[i];
        INFO(diff << " of " << gc.size() << " cells differ");
        CHECK(diff == 0);
    }
}
