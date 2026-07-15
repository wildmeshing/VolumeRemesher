if(WIN32)
    message(STATUS "Third-party: downloading gmp + mpfr")

    include(FetchContent)
    # CGAL's win64 auxiliary libraries. v6.2 ships GMP 6.3.0 -- the same version
    # Linux/macOS use, so the exact->double rounding matches and the byte-identical
    # cross-platform output is preserved -- and unlike the old v5.2.1 package (GMP
    # 5.0.1) it includes the gmpxx.h / gmpxx.lib C++ wrapper, so USE_GNU_GMP_CLASSES
    # can be enabled on Windows (much faster exact arithmetic than the in-house bignum).
    FetchContent_Declare(
        gmp_mpfr
        URL https://github.com/CGAL/cgal/releases/download/v6.2/CGAL-6.2-win64-auxiliary-libraries-gmp-mpfr.zip
        URL_MD5 3094e1aeb307c0447552a569169b5700
    )
    FetchContent_MakeAvailable(gmp_mpfr)

    # For CGAL and Cork
    set(ENV{GMP_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
    set(ENV{MPFR_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
else()
    # On Linux/macOS, gmp+mpfr should be installed system-wide
endif()