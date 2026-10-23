#ifndef VOLUMEREMESHER_PARALLEL_H
#define VOLUMEREMESHER_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
#ifdef __linux__
#include <sched.h>
#endif

namespace vol_rem {

// 0: no explicit setting (see num_threads).
inline std::atomic<unsigned> requested_num_threads{0};

// The number of threads parallel_blocks runs on: set_num_threads(n) when n > 0, else the
// VOLUMEREMESHER_NUM_THREADS environment variable, else the CPUs this process may run on. On
// Linux that is its affinity mask (taskset, cgroups, a batch scheduler's binding), which
// std::thread::hardware_concurrency() ignores: it counts every CPU of the machine.
inline unsigned num_threads()
{
    if (unsigned n = requested_num_threads.load(std::memory_order_relaxed)) return n;
    if (const char* env = std::getenv("VOLUMEREMESHER_NUM_THREADS")) {
        const long n = std::strtol(env, nullptr, 10);
        if (n > 0) return (unsigned)n;
    }
#ifdef __linux__
    cpu_set_t set;
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int n = CPU_COUNT(&set);
        if (n > 0) return (unsigned)n;
    }
#endif
    const unsigned n = std::thread::hardware_concurrency();
    return n ? n : 1;
}

// Set the number of threads parallel_blocks runs on; 0 restores the default (num_threads).
inline void set_num_threads(unsigned n)
{
    requested_num_threads.store(n, std::memory_order_relaxed);
}

//
//
// Run fn over blocks of [0,n) on num_threads() threads (at most one per block), dynamically
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
//
// An exception thrown by fn on a worker (ip_error throws when NFG_IP_ERROR_THROWS is defined,
// as it is in the Python module) is caught there, the remaining blocks are abandoned, and the
// first such exception is rethrown on the calling thread after every worker has joined. Left
// uncaught it would escape a std::thread and terminate the process.
template <class F>
inline void parallel_blocks(uint64_t n, F&& fn)
{
#ifdef VOLUMEREMESHER_SERIAL_TET
    fn(uint64_t(0), n);
#else
    const uint64_t block = 64;
    const unsigned nthreads = (unsigned)std::min<uint64_t>(num_threads(), (n + block - 1) / block);
    if (nthreads <= 1 || n < 512) {
        fn(uint64_t(0), n);
        return;
    }
    std::atomic<uint64_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;
    auto runner = [&]() {
        try {
            uint64_t i;
            while (!failed.load(std::memory_order_relaxed) && (i = next.fetch_add(block)) < n)
                fn(i, std::min(n, i + block));
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!error) error = std::current_exception();
            failed.store(true, std::memory_order_relaxed);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(runner);
    for (std::thread& t : pool) t.join();
    if (error) std::rethrow_exception(error);
#endif
}

} // namespace vol_rem

#endif // VOLUMEREMESHER_PARALLEL_H
