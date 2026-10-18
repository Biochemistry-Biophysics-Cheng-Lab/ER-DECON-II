#include "JacobiNL3Ddu.h"
#include "Res3D.h"
#include "LineSearch3D2n.h"
#include "fft_double.h"
#include <math.h>    /* fabs */
#include <string.h>  /* memcpy */


static double compute_rfactor(
    int nx,
    int ny,
    int nz,
    const double* guess,
    const double* tf,
    const double* measurement,
    FFTDoublePlan fpl,
    FFTDoublePlan rpl,
    double* work
);


int
JacobiNL3Ddu(
    int nx,
    int ny,
    int nz,
    double lamf,
    double lamratio,
    double lampos,
    double alpha,
    double laml2,
    double lammask,
    double delfrac,
    int regularization_type,
    void (*cost_reporter)(
	int iter,
	const ResReporterCosts* p_costs,
	double rfactor,
	void* closure
    ),
    void* cost_reporter_closure,
    int (*mask_updater)(
	int iter,
	int nx,
	int ny,
	int nz,
	const double* guess,
	double* mask,
	double* work,
	void* closure
    ),
    void* mask_updater_closure,
    int outputs_rfactor,
    int niter,
    const double* B,
    const double* TF,
    const double* Imn,
    double Diagd,
    const double* PCF,
    double* X,
    double* mask,
    double* p_cost,
    double* work
)
{
    double* R = work;
    FFTDoublePlan fpl = 0;
    FFTDoublePlan rpl = 0;
    int status = 1;
    int nxskip = (nx % 2 == 0) ? 2 : 1;
    int ind;

    fpl = plan_real_complex_3d_fft_double(nx, ny, nz, R, R);
    rpl = plan_complex_real_3d_fft_double(nx, ny, nz, R, R);
    if (fpl == 0 || rpl == 0) {
	goto cleanup;
    }

    for (ind = 0; ind < niter; ++ind) {
	double* D = work + (nx + nxskip) * (size_t) ny * (size_t) nz;
	double scl;
	double alphaw;
	double rfactor;
	ResReporterCosts pre_costs, post_costs;
	int iz, fresult;

	if (ind == 0 && cost_reporter != 0 && outputs_rfactor) {
	    rfactor = compute_rfactor(nx, ny, nz, X, TF, Imn, fpl, rpl, R);
	} else {
	    rfactor = 0.0;
	}

	fresult = Res3D(
	    nx,
	    ny,
	    nz,
	    lamf,
	    lamratio,
	    lampos,
	    alpha,
	    laml2,
	    lammask,
	    delfrac,
	    regularization_type,
	    X,
	    TF,
	    mask,
	    Diagd,
	    R,
	    D,
	    D + (nx + nxskip) * (size_t) ny * (size_t) nz
	);
	if (fresult != 0) {
	    goto cleanup;
	}

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) shared(nx, ny, nz, nxskip, R, B) schedule(dynamic,1)
#endif
	for (iz = 0; iz < nz; ++iz) {
	    double* rc = R + (ptrdiff_t) (nx + nxskip) * (ptrdiff_t) ny * iz;
	    const double* bc = B + (ptrdiff_t) (nx + nxskip) *
		(ptrdiff_t) ny * iz;
	    int iy;

	    for (iy = 0; iy < ny; ++iy) {
		int ix;

		for (ix = 0; ix < nx; ++ix) {
		    *rc = *bc - (*rc);
		    ++rc;
		    ++bc;
		}
		rc += nxskip;
		bc += nxskip;
	    }
	}

	compute_real_complex_3d_fft_double(fpl, R, R);

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) shared(nx, ny, nz, nxskip, R, PCF) schedule(dynamic,1)
#endif
	for (iz = 0; iz < nz; ++iz) {
	    const double* pcfc = PCF + (ptrdiff_t) (nx + nxskip) *
		(ptrdiff_t) ny * iz;
	    double* rc = R + (ptrdiff_t) (nx + nxskip) * (ptrdiff_t) ny * iz;
	    int iy;

	    for (iy = 0; iy < ny; ++iy) {
		int ix;

		for (ix = 0; ix < nx + nxskip; ix += 2) {
		    double rp = *rc;
		    double ip = *(rc + 1);

		    *rc = rp * (*pcfc) - ip * (*(pcfc + 1));
		    *(rc + 1) = rp * (*(pcfc + 1)) + ip * (*pcfc);
		    rc += 2;
		    pcfc += 2;
		}
	    }
	}
    
	compute_complex_real_3d_fft_double(rpl, R, R);

	scl = 1.0 / ((double) nx * (double) ny * (double) nz);
	scl *= scl;
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) shared(nx, ny, nz, nxskip, R, D, scl) schedule(dynamic,1)
#endif
	for (iz = 0; iz < nz; ++iz) {
	    const double* dc = D + (ptrdiff_t) (nx + nxskip) *
		(ptrdiff_t) ny * iz;
	    double* rc = R + (ptrdiff_t) (nx + nxskip) * (ptrdiff_t) ny * iz;
	    int iy;

	    for (iy = 0; iy < ny; ++iy) {
		int ix;

		for (ix = 0; ix < nx; ++ix) {
		    *rc *= scl / (*dc);
		    ++rc;
		    ++dc;
		}
		rc += nxskip;
		dc += nxskip;
	    }
	}

	compute_real_complex_3d_fft_double(fpl, R, R);

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) shared(nx, ny, nz, nxskip, R, PCF) schedule(dynamic,1)
#endif
	for (iz = 0; iz < nz; ++iz) {
	    const double* pcfc = PCF + (ptrdiff_t) (nx + nxskip) *
		(ptrdiff_t) ny * iz;
	    double* rc = R + (ptrdiff_t) (nx + nxskip) * (ptrdiff_t) ny * iz;
	    int iy;

	    for (iy = 0; iy < ny; ++iy) {
		int ix;

		for (ix = 0; ix < nx + nxskip; ix += 2) {
		    double rp = *rc;
		    double ip = *(rc + 1);

		    *rc = rp * (*pcfc) - ip * (*(pcfc + 1));
		    *(rc + 1) = rp * (*(pcfc + 1)) + ip * (*pcfc);
		    rc += 2;
		    pcfc += 2;
		}
	    }
	}

	compute_complex_real_3d_fft_double(rpl, R, R);

	fresult = LineSearch3D2n(
	    nx,
	    ny,
	    nz,
	    lamf,
	    lamratio,
	    lampos,
	    alpha,
	    laml2,
	    lammask,
	    delfrac,
	    regularization_type,
	    X,
	    R,
	    TF,
	    Imn,
	    mask,
	    &alphaw,
	    (ind == 0 && cost_reporter != 0) ?
	        &pre_costs : (ResReporterCosts*) 0,
	    &post_costs,
	    work + (nx + nxskip) * (size_t) ny * (size_t) nz
	);
	if (fresult != 0) {
	    goto cleanup;
	}
	*p_cost = post_costs.total;

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) shared(nx, ny, nz, nxskip, X, R, alphaw) schedule(dynamic,1)
#endif
	for (iz = 0; iz < nz; ++iz) {
	    const double* rc = R + (ptrdiff_t) (nx + nxskip) *
		(ptrdiff_t) ny * iz;
	    double* xc = X + (ptrdiff_t) (nx + nxskip) * (ptrdiff_t) ny * iz;
	    int iy;

	    for (iy = 0; iy < ny; ++iy) {
		int ix;

		for (ix = 0; ix < nx; ++ix) {
		    *xc += alphaw * (*rc);
		    ++xc;
		    ++rc;
		}
		rc += nxskip;
		xc += nxskip;
	    }
	}

	if (cost_reporter != 0) {
	    if (ind == 0) {
		(*cost_reporter)(
		    ind, &pre_costs, rfactor, cost_reporter_closure
		);
	    }
	    if (outputs_rfactor) {
		rfactor = compute_rfactor(nx, ny, nz, X, TF, Imn, fpl, rpl, R);
	    }
	    (*cost_reporter)(
		ind + 1, &post_costs, rfactor, cost_reporter_closure
	    );
	}

	if (mask_updater != 0) {
	    if ((*mask_updater)(
		    ind + 1,
		    nx,
		    ny,
		    nz,
		    X,
		    mask,
		    work,
		    mask_updater_closure) != 0) {
		goto cleanup;
	    }
	}
    }

    status = 0;
 cleanup:
    if (rpl != 0) {
	destroy_fft_double_plan(rpl);
    }
    if (fpl != 0) {
	destroy_fft_double_plan(fpl);
    }
    return status;
}


static
double
compute_rfactor(
    int nx,
    int ny,
    int nz,
    const double* guess,
    const double* tf,
    const double* measurement,
    FFTDoublePlan fpl,
    FFTDoublePlan rpl,
    double* work
)
{
    int nxskip = (nx % 2 == 0) ? 2 : 1;
    double num_sum = 0.0;
    double den_sum = 0.0;
    double scl = nx * (double) ny * (double) nz;
    double* cursor;
    const double* in_cursor;
    int iz;

    (void) memcpy(
	work, guess, nz * (size_t) ny * (size_t) (nx + nxskip) * sizeof(double)
    );
    compute_real_complex_3d_fft_double(fpl, work, work);

    cursor = work;
    in_cursor = tf;
    for (iz = 0; iz < nz; ++iz) {
	int iy;

	for (iy = 0; iy < ny; ++iy) {
	    int ix;

	    for (ix = 0; ix < nx + nxskip; ix += 2) {
		double rp = cursor[0];
		double ip = cursor[1];

		cursor[0] = rp * in_cursor[0] - ip * in_cursor[1];
		cursor[1] = rp * in_cursor[1] + ip * in_cursor[0];
		cursor += 2;
		in_cursor += 2;
	    }
	}
    }
    compute_complex_real_3d_fft_double(rpl, work, work);

    cursor = work;
    in_cursor = measurement;
    for (iz = 0; iz < nz; ++iz) {
	double num_sec_sum = 0.0;
	double den_sec_sum = 0.0;
	int iy;

	for (iy = 0; iy < ny; ++iy) {
	    double num_line_sum = 0.0;
	    double den_line_sum = 0.0;
	    int ix;

	    for (ix = 0; ix < nx; ++ix) {
		num_line_sum += fabs(*cursor / scl - *in_cursor);
		den_line_sum += fabs(*in_cursor);
		++cursor;
		++in_cursor;
	    }
	    cursor += nxskip;
	    in_cursor += nxskip;
	    num_sec_sum += num_line_sum;
	    den_sec_sum += den_line_sum;
	}
	num_sum += num_sec_sum;
	den_sum += den_sec_sum;
    }

    return (den_sum > 0.0) ? num_sum / den_sum : 0.0;
}
