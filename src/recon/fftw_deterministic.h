/*
 * fftw_deterministic.h -- force-included (gcc -include) when the build is
 * configured with FFTW_DETERMINISTIC=1 (the default).
 *
 * The delivered code plans FFTs with flags 0 = FFTW_MEASURE (fft_double.h,
 * JacobiNL4Ddu.c).  MEASURE times candidate algorithms, so the chosen plan --
 * and the last bit of every transform -- depends on machine load; on
 * ER-Decon I's accept/reject path this shows up as 1e-4 run-to-run
 * differences.  These macros add
 * FFTW_ESTIMATE to every planner call made by the delivered code, which
 * makes a build bit-reproducible, without editing the delivered files.
 * contrib/fftw_estimate_shim.c does the same, as an LD_PRELOAD library,
 * for a build (or a binary) that does not include this header.
 *
 * ERDecon reference build, 2026-09-22.
 */
#ifndef ERDECON_RECON_FFTW_DETERMINISTIC_H
#define ERDECON_RECON_FFTW_DETERMINISTIC_H

#include <fftw3.h>

#define ERDECON_FFTW_FORCE(f) \
    (((unsigned) (f) & ~(unsigned) (FFTW_MEASURE | FFTW_PATIENT | FFTW_EXHAUSTIVE)) | FFTW_ESTIMATE)

/* A function-like macro may name itself: the inner name is not re-expanded. */
#define fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, f) fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_c2r_3d(n0, n1, n2, in, out, f) fftw_plan_dft_c2r_3d(n0, n1, n2, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_r2c_2d(n0, n1, in, out, f)     fftw_plan_dft_r2c_2d(n0, n1, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_c2r_2d(n0, n1, in, out, f)     fftw_plan_dft_c2r_2d(n0, n1, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_r2c_1d(n0, in, out, f)         fftw_plan_dft_r2c_1d(n0, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_c2r_1d(n0, in, out, f)         fftw_plan_dft_c2r_1d(n0, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_r2c(r, n, in, out, f)          fftw_plan_dft_r2c(r, n, in, out, ERDECON_FFTW_FORCE(f))
#define fftw_plan_dft_c2r(r, n, in, out, f)          fftw_plan_dft_c2r(r, n, in, out, ERDECON_FFTW_FORCE(f))

#endif /* include guard */
