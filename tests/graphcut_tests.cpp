// Unit tests for the max-flow / min-cut behind the graph-cut cell classifier
// (src/VolumeRemesher/graphcut.h): on small random graphs, the flow it returns and the cut of the
// labeling it returns must both equal the minimum cut, found by trying every labeling.

#include <catch2/catch_test_macros.hpp>

#include "VolumeRemesher/graphcut.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

using vol_rem::GC_Graph;

namespace {
struct Edge
{
    int a, b;
    double capacity;
};

struct Problem
{
    int n = 0;
    std::vector<Edge> edges;
    std::vector<double> source, sink; // terminal capacities, as passed to AddCapacity
};

// The cost of a labeling: label 1 (sink side) pays the source capacity, label 0 the sink
// capacity, and every edge between different labels its capacity.
template <typename Label>
double cut_cost(const Problem& p, Label label)
{
    double c = 0;
    for (int i = 0; i < p.n; i++) c += label(i) ? p.source[i] : p.sink[i];
    for (const Edge& e : p.edges)
        if (label(e.a) != label(e.b)) c += e.capacity;
    return c;
}

double brute_force_min_cut(const Problem& p)
{
    double best = 1e300;
    for (uint32_t mask = 0; mask < (1u << p.n); mask++)
        best = std::min(best, cut_cost(p, [mask](int i) { return int((mask >> i) & 1); }));
    return best;
}

struct Result
{
    double flow, cut;
};

Result solve(const Problem& p)
{
    GC_Graph g(p.n);
    for (const Edge& e : p.edges) g.AddEdge(e.a, e.b, e.capacity);
    for (int i = 0; i < p.n; i++) g.AddCapacity(i, p.source[i], p.sink[i]);
    const double flow = g.ComputeMaxFlow();
    return {flow, cut_cost(p, [&g](int i) { return g.WhichLabel(i); })};
}
} // namespace

TEST_CASE("graph cut: a node dropped from the search trees is queued again", "[graphcut]")
{
    // The smallest graph on which graphcut.h, as first published in VolumeMesher, returned the
    // right flow (2) with a labeling (1 0 1) that cuts 6: the node it was growing left its tree
    // with its active flag still set, so it was never queued again (MarcoAttene/VolumeMesher#5).
    Problem p;
    p.n = 3;
    p.edges = {{0, 1, 2}, {0, 2, 1}, {1, 2, 1}};
    p.source = {2, 0, 0};
    p.sink = {0, 1, 2};
    const Result r = solve(p);
    CHECK(r.flow == 2);
    CHECK(r.cut == 2);
}

TEST_CASE("graph cut: flow and labeling are a minimum cut on random graphs", "[graphcut]")
{
    // Integer capacities, so every sum is exact and the comparisons can be ==. Graphs of up to 9
    // nodes keep the exhaustive search at 512 labelings.
    std::mt19937_64 rng(20261105);
    for (int trial = 0; trial < 20000; trial++) {
        Problem p;
        p.n = 2 + int(rng() % 8);
        const unsigned density = 1 + unsigned(rng() % 4); // an edge with probability density/4
        for (int i = 0; i < p.n; i++)
            for (int j = i + 1; j < p.n; j++)
                if (rng() % 4 < density) p.edges.push_back({i, j, double(rng() % 5)});
        for (int i = 0; i < p.n; i++) {
            p.source.push_back(rng() % 3 == 0 ? double(rng() % 6) : 0.0);
            p.sink.push_back(rng() % 3 == 0 ? double(rng() % 6) : 0.0);
        }
        const double min_cut = brute_force_min_cut(p);
        const Result r = solve(p);
        INFO("trial " << trial << ", " << p.n << " nodes");
        REQUIRE(r.flow == min_cut);
        REQUIRE(r.cut == min_cut);
    }
}
