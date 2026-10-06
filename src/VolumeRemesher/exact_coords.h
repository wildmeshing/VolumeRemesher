#ifndef VOLUMEREMESHER_EXACT_COORDS_H
#define VOLUMEREMESHER_EXACT_COORDS_H

// Exact output coordinates, in lowest terms, computed in parallel.
//
// genericPoint::getExactXYZCoordinates (and the 2D getExactXYCoordinates) return x = lx / d
// as a bigrational quotient. Since NFG a5ab699, bigrational's (num, den, sign) constructor and
// operator* no longer reduce, so those values come back unreduced (3/6, not 1/2). Callers rely
// on the reduced form: wildmeshing-toolkit hands the coordinates to GMP through mpq_set_str,
// and GMP requires canonical operands. So each coordinate is reduced here, at the cost of a
// GCD, and the whole pass runs on parallel_blocks.
//
// The catch is that NFG numbers cannot leave the thread that made them (see parallel.h). So
// each worker reduces its coordinates and exports them as plain 64-bit limbs; only the
// calling thread builds the returned bigrationals, by copying limbs, with no arithmetic. A
// reduced rational is unique, so the result does not depend on how the work was split.

#include <numerics.h>

#include "parallel.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <utility>
#include <vector>

namespace vol_rem {

#ifndef USE_GNU_GMP_CLASSES
namespace exact_coords_detail {

// canonicalize() is protected in NFG::bigrational. It also asserts a non-empty denominator,
// which zero does not have, so zero (already canonical) is left alone.
struct ReducibleRational : NFG::bigrational
{
    explicit ReducibleRational(NFG::bigrational&& r)
        : NFG::bigrational(std::move(r))
    {}
    void reduce()
    {
        if (sgn() != 0) canonicalize();
    }
};

// reserve() and push_back() are protected in NFG::bignatural. Limbs are most significant
// first, as bignatural::operator[] reads them.
struct LimbNatural : NFG::bignatural
{
    LimbNatural(const uint64_t* limbs, uint32_t n)
    {
        reserve(n);
        for (uint32_t i = 0; i < n; i++) push_back(limbs[i]);
    }
};

// One rational as plain words: sign + 1, #numerator limbs, #denominator limbs, then the
// numerator's limbs and the denominator's.
inline void export_rational(const NFG::bigrational& r, std::vector<uint64_t>& out)
{
    const NFG::bignatural& num = r.get_num();
    const NFG::bignatural& den = r.get_den();
    out.push_back(uint64_t(r.sgn() + 1));
    out.push_back(num.size());
    out.push_back(den.size());
    for (uint32_t i = 0; i < num.size(); i++) out.push_back(num[i]);
    for (uint32_t i = 0; i < den.size(); i++) out.push_back(den[i]);
}

// Inverse of export_rational: appends the rational at p to out and returns the word past it.
inline const uint64_t* import_rational(const uint64_t* p, std::vector<NFG::bigrational>& out)
{
    const int32_t sign = int32_t(p[0]) - 1;
    const uint32_t nn = uint32_t(p[1]), nd = uint32_t(p[2]);
    p += 3;
    if (sign == 0)
        out.emplace_back();
    else
        out.emplace_back(LimbNatural(p, nn), LimbNatural(p + nn, nd), sign);
    return p + nn + nd;
}

} // namespace exact_coords_detail
#endif

// Fills out with DIM * n exact coordinates, out[DIM * i + k] being coordinate k of point i,
// each in lowest terms. get(i, c) must store point i's exact coordinates in c[0..DIM), in any
// (possibly unreduced) form, and return false if it cannot. get runs on worker threads, so it
// must be safe to call concurrently for different i. Returns false, leaving out unspecified,
// if any call of get failed.
template <int DIM, class Get>
bool exact_coords_reduced(uint64_t n, Get&& get, std::vector<NFG::bigrational>& out)
{
    std::atomic<bool> failed{false};

#ifdef USE_GNU_GMP_CLASSES
    // mpq_class is always canonical and lives on the malloc heap, not in a thread-local
    // pool, so workers may write the output directly.
    out.assign(DIM * n, NFG::bigrational());
    parallel_blocks(n, [&](uint64_t lo, uint64_t hi) {
        for (uint64_t i = lo; i < hi; i++)
            if (!get(i, &out[DIM * i])) failed = true;
    });
#else
    using namespace exact_coords_detail;

    // Each call of the worker lambda exports [lo, hi) into one buffer of plain words.
    std::vector<std::pair<uint64_t, std::vector<uint64_t>>> chunks;
    std::mutex chunks_mutex;
    parallel_blocks(n, [&](uint64_t lo, uint64_t hi) {
        std::vector<uint64_t> words;
        for (uint64_t i = lo; i < hi; i++) {
            NFG::bigrational c[DIM];
            if (!get(i, c)) {
                failed = true;
                return;
            }
            for (int k = 0; k < DIM; k++) {
                ReducibleRational r(std::move(c[k]));
                r.reduce();
                export_rational(r, words);
            }
        }
        std::lock_guard<std::mutex> lock(chunks_mutex);
        chunks.emplace_back(lo, std::move(words));
    });
    if (failed) return false;

    // Rebuild on this thread, in point order. Each buffer is released once it is copied out,
    // so the limbs and the bigrationals are not both held in full.
    std::sort(chunks.begin(), chunks.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    out.clear();
    out.reserve(DIM * n);
    for (auto& chunk : chunks) {
        const uint64_t* p = chunk.second.data();
        const uint64_t* end = p + chunk.second.size();
        while (p != end) p = import_rational(p, out);
        std::vector<uint64_t>().swap(chunk.second);
    }
#endif

    return !failed;
}

} // namespace vol_rem

#endif // VOLUMEREMESHER_EXACT_COORDS_H
