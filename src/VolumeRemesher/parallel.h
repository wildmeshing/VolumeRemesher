#ifndef VOLUMEREMESHER_PARALLEL_H
#define VOLUMEREMESHER_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <thread>
#include <vector>

namespace vol_rem {

//
//
// Run fn over blocks of [0,n) on std::thread::hardware_concurrency() threads, dynamically
// scheduled (atomic block-fetch) so uneven per-item cost stays balanced. Falls back to a
// serial call for small n or a single hardware thread. Compiling with VOLUMEREMESHER_SERIAL_TET
// (CMake: -DVOLUMEREMESHER_PARALLEL_TETRAHEDRALIZATION=OFF) forces the serial path everywhere.
//
// CAREFUL with exact arithmetic in here: NFG's bignatural (the storage under bigfloat and
// bigrational) and expansion allocate from thread_local pools, so such a value must be born,
// used and destroyed inside ONE call of fn -- let only indices and other plain data escape.
// Two mistakes compile cleanly and assert nothing:
//   - computing values in parallel and reading them after the join is a use-after-free (the
//     worker's pool dies with the worker);
//   - destroying a value on a thread other than the one that allocated it pushes the block
//     onto the wrong thread's free stack, which is sized for its own pool and can overflow.
// So a worker must not even assign into a caller-owned bigrational: the caller destroys it.
// See exact_coords.h for how exact results are carried out of a worker as plain limbs.
template <class F>
inline void parallel_blocks(uint64_t n, F&& fn)
{
#ifdef VOLUMEREMESHER_SERIAL_TET
    fn(uint64_t(0), n);
#else
    unsigned nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    if (nthreads == 1 || n < 512) {
        fn(uint64_t(0), n);
        return;
    }
    std::atomic<uint64_t> next{0};
    const uint64_t block = 64;
    auto runner = [&]() {
        uint64_t i;
        while ((i = next.fetch_add(block)) < n) fn(i, std::min(n, i + block));
    };
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(runner);
    for (std::thread& t : pool) t.join();
#endif
}

} // namespace vol_rem

#endif // VOLUMEREMESHER_PARALLEL_H
