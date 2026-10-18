#ifndef UCSF_MSG_FFT_DOUBLE_H
#define UCSF_MSG_FFT_DOUBLE_H

#include <stdlib.h>  /* free, malloc */

/*
 * Defines macros for covering the differences between FFTW 2 and FFTW 3
 * double-precision real->complex and complex->real 3D FFTs.  Basic usage is:
 *
 * initialize_fft_double()
 *     Call once before using any of the other routines.  Returns zero if
 *     an error occurred; otherwise returns a non-zero integer value.
 *
 * set_fft_double_thread_count(nthread)
 *     Sets the number of threads to use.  With FFTW 3, the number of threads
 *     currently selected when a plan is created will be the number of threads
 *     that the plan uses.  With FFTW 2, the number of threads currently
 *     selected when an FFT is calculated will be the number of threads the
 *     used when calculating the transform.
 *
 * double* p = allocate_real_complex_3d_double(nx, ny, nz)
 * double* p = allocate_real_complex_2d_double(nx, ny)
 * double* p = allocate_real_complex_1d_double(nx)
 *     Allocates memory for array which will be used for real->complex
 *     or complex->real transforms.  Returns (double*) 0 if the allocation
 *     failed.  You could use another memory allocation mechanism, but this
 *     one may align the data appropriately for a faster transform.
 *
 * double* p;
 * deallocate_fft_double_array(p)
 *     Releases memory allocated by allocate_real_complex_3d_double() or
 *     allocate_real_complex_2d_double().
 *
 * double* in;
 * double* out;
 * FFTDoublePlan plan1 = plan_real_complex_3d_fft_double(nx, ny, nz, in, out)
 * FFTDoublePlan plan2 = plan_complex_real_3d_fft_double(nx, ny, nz, in, out)
 * FFTDoublePlan plan3 =
 *         plan_no_overwrite_real_complex_3d_fft_double(nx, ny, nz, in, out)
 * FFTDoublePlan plan4 =
 *         plan_no_overwrite_complex_real_3d_fft_double(nx, ny, nz, in, out)
 *     The first initializes the auxiliary information to compute a 3D
 *     real to complex FFT.  The second initializes the auxiliary information
 *     compute a 3D complex to real FFT.  The last two are similar to the first
 *     two but do not overwrite the contents of in and out.  In all cases, the
 *     return value is not zero if the operation succeeded and zero if the
 *     operation failed.
 *
 * double* in;
 * double* out;
 * FFTDoublePlan plan1 = plan_real_complex_2d_fft_double(nx, ny, in, out)
 * FFTDoublePlan plan2 = plan_complex_real_2d_fft_double(nx, ny, in, out)
 * FFTDoublePlan plan3 =
 *         plan_no_overwrite_real_complex_2d_fft_double(nx, ny, in, out)
 * FFTDoublePlan plan4 =
 *         plan_no_overwrite_complex_real_2d_fft_double(nx, ny, in, out)
 *     The first initializes the auxiliary information to compute a 2D
 *     real to complex FFT.  The second initializes the auxiliary information
 *     compute a 2D complex to real FFT.  The last two are similar to the first
 *     two but do not overwrite the contents of in and out.  In all cases, the
 *     return value is not zero if the operation succeeded and zero if the
 *     operation failed.
 *
 * double* in;
 * double* out;
 * FFTDoublePlan plan1 = plan_real_complex_1d_fft_double(nx, in, out)
 * FFTDoublePlan plan2 = plan_complex_real_1d_fft_double(nx, in, out)
 * FFTDoublePlan plan3 =
 *         plan_no_overwrite_real_complex_1d_fft_double(nx, in, out)
 * FFTDoublePlan plan4 =
 *         plan_no_overwrite_complex_real_1d_fft_double(nx, in, out)
 *     The first initializes the auxiliary information to compute a 2D
 *     real to complex FFT.  The second initializes the auxiliary information
 *     compute a 2D complex to real FFT.  The last two are similar to the first
 *     two but do not overwrite the contents of in and out.  In all cases, the
 *     return value is not zero if the operation succeeded and zero if the
 *     operation failed.
 *
 * FFTDoublePlan plan;
 * destroy_fft_double_plan(plan)
 *     Releases the resources allocated to a plan by a prior call to
 *     plan_complex_real_3d_fft_double(), plan_real_complex_3d_fft_double(),
 *     plan_no_overwrite_real_complex_3d_fft_double(),
 *     plan_no_overwrite_complex_real_3d_fft_double(),
 *     plan_complex_real_2d_fft_double(), plan_real_complex_2d_fft_double(),
 *     plan_no_overwrite_real_complex_2d_fft_double(), and
 *     plan_no_overwrite_complex_real_2d_fft_double().
 *
 * FFTDoublePlan plan;
 * double* in;
 * double* out;
 * compute_real_complex_3d_fft_double(plan, in, out)
 * compute_complex_real_3d_fft_double(plan, in, out)
 *     The first performs a 3D real to complex FFT.  The second performs
 *     a 3D complex to real FFT.  For compute_real_complex_3d_fft_double(),
 *     the plan must have been computed with plan_real_complex_3d_fft_double().
 *     For compute_complex_real_3d_fft_double(), the plan must have been
 *     computed with plan_complex_real_3d_fft_double().  The in and out
 *     arguments used with the plan must match those that were used when the
 *     plan was generated.
 *
 * FFTDoublePlan plan;
 * double* in;
 * double* out;
 * compute_real_complex_2d_fft_double(plan, in, out)
 * compute_complex_real_2d_fft_double(plan, in, out)
 *     The first performs a 2D real to complex FFT.  The second performs
 *     a 2D complex to real FFT.  For compute_real_complex_2d_fft_double(),
 *     the plan must have been computed with plan_real_complex_2d_fft_double().
 *     For compute_complex_real_2d_fft_double(), the plan must have been
 *     computed with plan_complex_real_2d_fft_double().  The in and out
 *     arguments used with the plan must match those that were used when the
 *     plan was generated.
 *
 * FFTDoublePlan plan;
 * double* in;
 * double* out;
 * compute_real_complex_1d_fft_double(plan, in, out)
 * compute_complex_real_1d_fft_double(plan, in, out)
 *     The first performs a 1D real to complex FFT.  The second performs
 *     a 1D complex to real FFT.  For compute_real_complex_1d_fft_double(),
 *     the plan must have been computed with plan_real_complex_1d_fft_double().
 *     For compute_complex_real_1d_fft_double(), the plan must have been
 *     computed with plan_complex_real_1d_fft_double().  The in and out
 *     arguments used with the plan must match those that were used when the
 *     plan was generated.
 *
 * Note that the result of real to complex transforms and the input to
 * complex to real transforms use a half-complex format with 2 * (nx / 2 + 1)
 * by ny by nz contiguous array (for the 3D case; for the 2D case it is
 * 2 * (nx /2 + 1) by ny) of floats where the real and imaginary parts
 * for each point are adjacent in memory with the real part first.  For
 * in-place transforms, the input to the real to complex transform or result
 * of a complex to real transform must also be structured as a 2 * (nx / 2 + 1)
 * by ny by nz array (2 * (nx /2 + t) * ny for the 2D case).  The values in
 * the additional elements in each row of the real valued array do not have an
 * effect on the in-place transform.
 */

/* If the preprocessor macro FFTW2 is set use FFTW 2; otherwise use FFTW 3. */
#ifdef FFTW2

#ifdef FFT_THREADED
#include "dfftw_threads.h"
#include "drfftw_threads.h"
#else
#include "dfftw.h"
#include "drfftw.h"
#endif

typedef rfftwnd_plan FFTDoublePlan;

#ifdef FFT_THREADED
#define initialize_fft_double() (! fftw_threads_init())
void set_fft_double_thread_count(int nthreads);
int get_fft_double_thread_count(void);
#else
#define initialize_fft_double() (1)
#define set_fft_double_thread_count(nthreads) 
#endif

#define allocate_real_complex_3d_double(nx, ny, nz) ((double*) malloc((2 * ((nx) / 2 + 1)) * (size_t) (ny) * (size_t) (nz) * sizeof(double)))
#define allocate_real_complex_2d_double(nx, ny) ((double*) malloc((2 * ((nx) / 2 + 1)) * (size_t) (ny) * sizeof(double)))
#define allocate_real_complex_1d_double(nx) ((double*) malloc((2 * ((nx) / 2 + 1)) * sizeof(double)))
#define deallocate_fft_double_array(a) free((a))
/*
 * These eight macros access in and out multiple times:  do not use in and out
 * expressions with side effects.
 */
#define plan_complex_real_3d_fft_double(nx, ny, nz, in, out) rfftw3d_create_plan_specific((nz), (ny), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_real_complex_3d_fft_double(nx, ny, nz, in, out) rfftw3d_create_plan_specific((nz), (ny), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_no_overwrite_complex_real_3d_fft_double(nx, ny, nz, in, out) rfftw3d_create_plan((nz), (ny), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define plan_no_overwrite_real_complex_3d_fft_double(nx, ny, nz, in, out) rfftw3d_create_plan((nz), (ny), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define plan_complex_real_2d_fft_double(nx, ny, in, out) rfftw2d_create_plan_specific((ny), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_real_complex_2d_fft_double(nx, ny, in, out) rfftw2d_create_plan_specific((ny), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_no_overwrite_complex_real_2d_fft_double(nx, ny, in, out) rfftw2d_create_plan((ny), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define plan_no_overwrite_real_complex_2d_fft_double(nx, ny, in, out) rfftw2d_create_plan((ny), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define plan_complex_real_1d_fft_double(nx, in, out) rfftw2d_create_plan_specific((1), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_real_complex_1d_fft_double(nx, in, out) rfftw2d_create_plan_specific((1), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM, (in), 1, (out), 1)
#define plan_no_overwrite_complex_real_1d_fft_double(nx, in, out) rfftw2d_create_plan((1), (nx), FFTW_COMPLEX_TO_REAL, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define plan_no_overwrite_real_complex_1d_fft_double(nx, in, out) rfftw2d_create_plan((1), (nx), FFTW_REAL_TO_COMPLEX, (((in) == (out)) ? FFTW_IN_PLACE : FFTW_OUT_OF_PLACE) | FFTW_USE_WISDOM)
#define destroy_fft_double_plan(plan) rfftwnd_destroy_plan((plan))
#ifdef FFT_THREADED
#define compute_complex_real_3d_fft_double(plan, in, out) rfftwnd_threads_one_complex_to_real(get_fft_double_thread_count(), (plan), (fftw_complex*) (in), (out))
#define compute_real_complex_3d_fft_double(plan, in, out) rfftwnd_threads_one_real_to_complex(get_fft_double_thread_count(), (plan), (in), (fftw_complex*) (out))
#define compute_complex_real_2d_fft_double(plan, in, out) rfftwnd_threads_one_complex_to_real(get_fft_double_thread_count(), (plan), (fftw_complex*) (in), (out))
#define compute_real_complex_2d_fft_double(plan, in, out) rfftwnd_threads_one_real_to_complex(get_fft_double_thread_count(), (plan), (in), (fftw_complex*) (out))
#define compute_complex_real_1d_fft_double(plan, in, out) rfftwnd_threads_one_complex_to_real(get_fft_double_thread_count(), (plan), (fftw_complex*) (in), (out))
#define compute_real_complex_1d_fft_double(plan, in, out) rfftwnd_threads_one_real_to_complex(get_fft_double_thread_count(), (plan), (in), (fftw_complex*) (out))
#else
#define compute_complex_real_3d_fft_double(plan, in, out) rfftwnd_one_complex_to_real((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_3d_fft_double(plan, in, out) rfftwnd_one_real_to_complex((plan), (in), (fftw_complex*) (out))
#define compute_complex_real_2d_fft_double(plan, in, out) rfftwnd_one_complex_to_real((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_2d_fft_double(plan, in, out) rfftwnd_one_real_to_complex((plan), (in), (fftw_complex*) (out))
#define compute_complex_real_1d_fft_double(plan, in, out) rfftwnd_one_complex_to_real((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_1d_fft_double(plan, in, out) rfftwnd_one_real_to_complex((plan), (in), (fftw_complex*) (out))
#endif

#else

#include "fftw3.h"

typedef fftw_plan FFTDoublePlan;

#ifdef FFT_THREADED
#define initialize_fft_double() (fftw_init_threads())
#define set_fft_double_thread_count(nthreads) fftw_plan_with_nthreads(nthreads)
#else
#define initialize_fft_double() (1)
#define set_fft_double_thread_count(nthreads) 
#endif  /* FFT_THREADED */

#define allocate_real_complex_3d_double(nx, ny, nz) ((double *) fftw_malloc((2 * ((nx) / 2 + 1)) * (size_t) (ny) * (size_t) (nz) * sizeof(double)))
#define allocate_real_complex_2d_double(nx, ny) ((double *) fftw_malloc((2 * ((nx) / 2 + 1)) * (size_t) (ny) * sizeof(double)))
#define allocate_real_complex_1d_double(nx) ((double *) fftw_malloc((2 * ((nx) / 2 + 1)) * sizeof(double)))
#define deallocate_fft_double_array(a) fftw_free((a))
#define plan_complex_real_3d_fft_double(nx, ny, nz, in, out) fftw_plan_dft_c2r_3d((nz), (ny), (nx), (fftw_complex*) (in), (out), 0)
#define plan_real_complex_3d_fft_double(nx, ny, nz, in, out) fftw_plan_dft_r2c_3d((nz), (ny), (nx), (in), (fftw_complex*) (out), 0)
#define plan_no_overwrite_complex_real_3d_fft_double(nx, ny, nz, in, out) fftw_plan_dft_c2r_3d((nz), (ny), (nx), (fftw_complex*) (in), (out), FFTW_ESTIMATE)
#define plan_no_overwrite_real_complex_3d_fft_double(nx, ny, nz, in, out) fftw_plan_dft_r2c_3d((nz), (ny), (nx), (in), (fftw_complex*) (out), FFTW_ESTIMATE)
#define plan_complex_real_2d_fft_double(nx, ny, in, out) fftw_plan_dft_c2r_2d((ny), (nx), (fftw_complex*) (in), (out), 0)
#define plan_real_complex_2d_fft_double(nx, ny, in, out) fftw_plan_dft_r2c_2d((ny), (nx), (in), (fftw_complex*) (out), 0)
#define plan_no_overwrite_complex_real_2d_fft_double(nx, ny, in, out) fftw_plan_dft_c2r_2d((ny), (nx), (fftw_complex*) (in), (out), FFTW_ESTIMATE)
#define plan_no_overwrite_real_complex_2d_fft_double(nx, ny, in, out) fftw_plan_dft_r2c_2d((ny), (nx), (in), (fftw_complex*) (out), FFTW_ESTIMATE)
#define plan_complex_real_1d_fft_double(nx, in, out) fftw_plan_dft_c2r_1d((nx), (fftw_complex*) (in), (out), 0)
#define plan_real_complex_1d_fft_double(nx, in, out) fftw_plan_dft_r2c_1d((nx), (in), (fftw_complex*) (out), 0)
#define plan_no_overwrite_complex_real_1d_fft_double(nx, in, out) fftw_plan_dft_c2r_1d((nx), (fftw_complex*) (in), (out), FFTW_ESTIMATE)
#define plan_no_overwrite_real_complex_1d_fft_double(nx, in, out) fftw_plan_dft_r2c_1d((nx), (in), (fftw_complex*) (out), FFTW_ESTIMATE)
#define destroy_fft_double_plan(plan) fftw_destroy_plan((plan))
#define compute_complex_real_3d_fft_double(plan, in, out) fftw_execute_dft_c2r((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_3d_fft_double(plan, in, out) fftw_execute_dft_r2c((plan), (in), (fftw_complex*) (out))
#define compute_complex_real_2d_fft_double(plan, in, out) fftw_execute_dft_c2r((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_2d_fft_double(plan, in, out) fftw_execute_dft_r2c((plan), (in), (fftw_complex*) (out))
#define compute_complex_real_1d_fft_double(plan, in, out) fftw_execute_dft_c2r((plan), (fftw_complex*) (in), (out))
#define compute_real_complex_1d_fft_double(plan, in, out) fftw_execute_dft_r2c((plan), (in), (fftw_complex*) (out))

#endif

#endif /* include guard */
