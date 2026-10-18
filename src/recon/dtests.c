/*
 * dtests.c -- RECONSTRUCTED (not part of the delivered source).
 * See dtests.h.  ERDecon reference build, 2026-09-18.
 */
#include "dtests.h"

#include <math.h>
#include <stdio.h>

static double g_abs_tol = 1e-9;
static double g_rel_tol = 1e-9;

void
dcompare_set_tolerance(double abs_tol, double rel_tol)
{
    g_abs_tol = abs_tol;
    g_rel_tol = rel_tol;
}

void
dcompare_get_tolerance(double* p_abs_tol, double* p_rel_tol)
{
    if (p_abs_tol != 0) {
	*p_abs_tol = g_abs_tol;
    }
    if (p_rel_tol != 0) {
	*p_rel_tol = g_rel_tol;
    }
}

size_t
dcompare(
    size_t n,
    const double* expected,
    const double* actual,
    struct dcompare_result* p_res
)
{
    struct dcompare_result r;
    double sumsq = 0.0;
    size_t i;

    r.n = n;
    r.n_mismatch = 0;
    r.n_nonfinite = 0;
    r.imax = 0;
    r.max_abs = 0.0;
    r.max_rel = 0.0;
    r.rms = 0.0;
    r.expected_at_max = (n > 0) ? expected[0] : 0.0;
    r.actual_at_max = (n > 0) ? actual[0] : 0.0;

    for (i = 0; i < n; ++i) {
	double e = expected[i];
	double a = actual[i];
	double d;

	if (! isfinite(a)) {
	    ++r.n_nonfinite;
	    ++r.n_mismatch;
	    if (r.max_abs != HUGE_VAL) {
		r.max_abs = HUGE_VAL;
		r.imax = i;
		r.expected_at_max = e;
		r.actual_at_max = a;
	    }
	    continue;
	}
	d = fabs(e - a);
	sumsq += d * d;
	if (d > r.max_abs) {
	    r.max_abs = d;
	    r.imax = i;
	    r.expected_at_max = e;
	    r.actual_at_max = a;
	}
	if (e != 0.0) {
	    double rel = d / fabs(e);
	    if (rel > r.max_rel) {
		r.max_rel = rel;
	    }
	}
	if (d > g_abs_tol + g_rel_tol * fabs(e)) {
	    ++r.n_mismatch;
	}
    }
    r.rms = (n > 0) ? sqrt(sumsq / (double) n) : 0.0;

    if (p_res != 0) {
	*p_res = r;
    }
    return r.n_mismatch;
}

void
print_test_result(
    const char* what,
    int test_number,
    int nx,
    int ny,
    int nz,
    int nt,
    const struct dcompare_result* p_res,
    int* p_status
)
{
    const char* verdict = (p_res->n_mismatch == 0) ? "PASS" : "FAIL";

    (void) printf(
	"%s: %s (test %d, %dx%dx%dx%d, %lu values): "
	"max|d|=%.3e max_rel=%.3e rms=%.3e mismatches=%lu",
	verdict, what, test_number, nx, ny, nz, nt,
	(unsigned long) p_res->n, p_res->max_abs, p_res->max_rel, p_res->rms,
	(unsigned long) p_res->n_mismatch
    );
    if (p_res->n_mismatch != 0) {
	(void) printf(
	    "\n      worst at [%lu]: expected %+.16e actual %+.16e",
	    (unsigned long) p_res->imax, p_res->expected_at_max,
	    p_res->actual_at_max
	);
	if (p_res->n_nonfinite != 0) {
	    (void) printf(
		"  (%lu non-finite)", (unsigned long) p_res->n_nonfinite
	    );
	}
	if (p_status != 0) {
	    *p_status = 1;
	}
    }
    (void) printf("\n");
}
