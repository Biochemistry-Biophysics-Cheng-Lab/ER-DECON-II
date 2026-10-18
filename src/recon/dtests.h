/*
 * dtests.h -- RECONSTRUCTED (not part of the delivered source).
 *
 * Minimal comparison harness used by JacobiNL{2,3,4}DduTests.c.  The
 * interface was inferred from the call sites in those files:
 *
 *     struct dcompare_result res;
 *     dcompare(n, expected, actual, &res);
 *     print_test_result("x from JacobiNL3Ddu", test_no, nx, ny, nz, nt,
 *                       &res, &status);
 *
 * Nothing in the delivered tree depends on the *contents* of
 * struct dcompare_result, so the fields below are our own choice.
 * print_test_result() sets *p_status to nonzero when the comparison fails
 * the tolerance so that the Test*() functions return a failure code.
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_DTESTS_H
#define ERDECON_RECON_DTESTS_H

#include <stddef.h>

struct dcompare_result {
    size_t n;            /* number of elements compared */
    size_t n_mismatch;   /* elements with |e - a| > abs_tol + rel_tol*|e| */
    size_t n_nonfinite;  /* elements where actual is NaN or Inf */
    size_t imax;         /* index of the largest absolute difference */
    double max_abs;      /* largest |expected - actual| */
    double max_rel;      /* largest |expected - actual| / |expected| (expected != 0) */
    double rms;          /* root-mean-square of (expected - actual) */
    double expected_at_max;
    double actual_at_max;
};

/*
 * Tolerances used to decide pass/fail.  The golden values were produced by
 * MATLAB/NumPy on a different machine, so bit-for-bit equality is not the
 * criterion; these defaults are deliberately tight and can be relaxed via
 * dcompare_set_tolerance() by a test driver if a documented reason exists.
 */
void dcompare_set_tolerance(double abs_tol, double rel_tol);
void dcompare_get_tolerance(double* p_abs_tol, double* p_rel_tol);

/* Compare n doubles; fill *p_res.  Returns number of mismatches. */
size_t dcompare(
    size_t n,
    const double* expected,
    const double* actual,
    struct dcompare_result* p_res
);

/*
 * Print one line describing the outcome and, on failure, the worst element.
 * nx, ny, nz, nt only label the array shape in the message.  Sets
 * *p_status = 1 if the comparison failed; leaves it untouched otherwise.
 */
void print_test_result(
    const char* what,
    int test_number,
    int nx,
    int ny,
    int nz,
    int nt,
    const struct dcompare_result* p_res,
    int* p_status
);

#endif /* include guard */
