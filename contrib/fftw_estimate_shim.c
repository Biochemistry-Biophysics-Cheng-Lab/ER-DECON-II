/*
 * fftw_estimate_shim.c -- LD_PRELOAD shim that forces FFTW_ESTIMATE planning.
 *
 * The delivered code plans with flags 0 (= FFTW_MEASURE): the plan, and so
 * the last-bit rounding of every FFT, depends on timing.  On a loaded machine
 * the ORIGINAL core_decon/core2_decon are then not reproducible run to run
 * (seen 2026-09-22: 9.6e-5 rel RMS between repeated runs of core_decon
 * -sub=105 -ncycl=3, where the accept/reject steps amplify the rounding).
 * Preloading this into both the original binary and ours makes both
 * deterministic and directly comparable.
 *
 *   gcc -shared -fPIC -O2 -o fftw_estimate_shim.so fftw_estimate_shim.c -ldl
 *   LD_PRELOAD=$PWD/fftw_estimate_shim.so core2_decon ...
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fftw3.h>

#define FORCE(f) (((f) & ~(FFTW_MEASURE | FFTW_PATIENT | FFTW_EXHAUSTIVE)) | FFTW_ESTIMATE)

#define WRAP3(name, T1, T2) \
fftw_plan name(int n0, int n1, int n2, T1 in, T2 out, unsigned flags) { \
    static fftw_plan (*real)(int, int, int, T1, T2, unsigned) = 0; \
    if (!real) real = dlsym(RTLD_NEXT, #name); \
    return real(n0, n1, n2, in, out, FORCE(flags)); }
#define WRAP2(name, T1, T2) \
fftw_plan name(int n0, int n1, T1 in, T2 out, unsigned flags) { \
    static fftw_plan (*real)(int, int, T1, T2, unsigned) = 0; \
    if (!real) real = dlsym(RTLD_NEXT, #name); \
    return real(n0, n1, in, out, FORCE(flags)); }
#define WRAP1(name, T1, T2) \
fftw_plan name(int n0, T1 in, T2 out, unsigned flags) { \
    static fftw_plan (*real)(int, T1, T2, unsigned) = 0; \
    if (!real) real = dlsym(RTLD_NEXT, #name); \
    return real(n0, in, out, FORCE(flags)); }
#define WRAPN(name, T1, T2) \
fftw_plan name(int rank, const int* n, T1 in, T2 out, unsigned flags) { \
    static fftw_plan (*real)(int, const int*, T1, T2, unsigned) = 0; \
    if (!real) real = dlsym(RTLD_NEXT, #name); \
    return real(rank, n, in, out, FORCE(flags)); }

WRAP3(fftw_plan_dft_r2c_3d, double*, fftw_complex*)
WRAP3(fftw_plan_dft_c2r_3d, fftw_complex*, double*)
WRAP2(fftw_plan_dft_r2c_2d, double*, fftw_complex*)
WRAP2(fftw_plan_dft_c2r_2d, fftw_complex*, double*)
WRAP1(fftw_plan_dft_r2c_1d, double*, fftw_complex*)
WRAP1(fftw_plan_dft_c2r_1d, fftw_complex*, double*)
WRAPN(fftw_plan_dft_r2c, double*, fftw_complex*)
WRAPN(fftw_plan_dft_c2r, fftw_complex*, double*)
