// Probe for the GMP C++ classes (gmpxx.h). Compiled AND linked by CMake against gmp::gmp;
// if it succeeds, USE_GNU_GMP_CLASSES is enabled (GMP-backed exact numbers) instead of the
// slow in-house bignum. It exercises the mpq_class / mpz_class / mpf_class operations that
// numerics.h relies on, so a header that merely parses but cannot be instantiated (e.g. under
// MSVC) is rejected and the build falls back gracefully rather than failing. Every operation
// used here is header-inline over the C library, so only libgmp -- not libgmpxx -- is needed
// at link time (mirroring the actual code, which never calls a libgmpxx function).
#include <gmpxx.h>
#include <string>

int main()
{
    // mpq_class -- bigrational / bigfloat
    mpq_class a(1.5), b(2.0), c;
    c = a + b;
    c = a - b;
    c = a * b;
    c = a / b;
    c = -a;
    bool r = (a < b) && (a != b) && (a <= b) && (b > a);
    a.canonicalize();
    double d = a.get_d();
    std::string s = a.get_str();
    mpq_class e(s, 10);
    (void)a.get_mpq_t();

    // mpz_class -- bignatural
    mpz_class z(42), z2 = z * 2 + 1;
    std::string zs = z2.get_str();

    // mpf_class -- used by the sqrt / log2 / add1ULP helpers
    mpf_class bf(a, 128);
    mpf_class sq = sqrt(bf);
    mpf_class ulp(1U, bf.get_prec());
    ulp.get_mpf_t()->_mp_exp = bf.get_mpf_t()->_mp_exp;
    bf += ulp;
    mpq_class back(sq);

    (void)c;
    (void)r;
    (void)d;
    (void)e;
    (void)zs;
    (void)back;
    return 0;
}
