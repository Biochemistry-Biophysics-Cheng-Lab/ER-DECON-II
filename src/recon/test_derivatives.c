/*
 * test_derivatives.c -- unit test for Derivative() / Derivative2D().
 *
 *   1. adjointness   <L_f x, y> == <x, L_r y>          every slot, random x, y
 *   2. 's' mode      (L o L)^T W == sum_r' L(r',r)^2 W(r')  built explicitly
 *                    from L_f applied to unit vectors
 *   3. 'r' mode      accumulates (R += L^T X): calling twice doubles the result
 *
 * Run on a few small even/odd geometries; exits nonzero on failure.
 * ERDecon reference build, 2026-09-21.
 */
#include "volumestruct.h"
#include "derivatives.h"
#include "volumestruct2d.h"
#include "derivatives2d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double
urand(void)
{
    return 2.0 * ((double) rand() / RAND_MAX) - 1.0;
}

static void
fill_random(double* a, long nx, long ny, long nz, long nxpad)
{
    long j, ix;

    for (j = 0; j < ny * nz; ++j) {
	for (ix = 0; ix < nxpad; ++ix) {
	    a[j * nxpad + ix] = (ix < nx) ? urand() : 1e300; /* poison pads */
	}
    }
}

static double
dot(const double* a, const double* b, long nx, long ny, long nz, long nxpad)
{
    double s = 0.0;
    long j, ix;

    for (j = 0; j < ny * nz; ++j) {
	for (ix = 0; ix < nx; ++ix) {
	    s += a[j * nxpad + ix] * b[j * nxpad + ix];
	}
    }
    return s;
}

/* Apply the operator for either dimensionality through one entry point. */
typedef struct {
    int dim;
    Geometry g3;
    Geometry2D g2;
    long nx, ny, nz, nxpad, len;
} Geo;

static void
apply(const Geo* g, const double* X, double* R, const char* dir, char mode,
      double* swp)
{
    if (g->dim == 3) {
	Derivative(X, R, &g->g3, dir, mode, swp);
    } else {
	Derivative2D(X, R, &g->g2, dir, mode, swp);
    }
}

static void
setup(Geo* g, int dim, long nx, long ny, long nz)
{
    g->dim = dim;
    if (dim == 3) {
	ConstructVolumeStruct(&g->g3, nx, ny, nz, 0);
	g->nxpad = g->g3.Nx + g->g3.NxFSkip;
	g->len = g->g3.LengthF;
    } else {
	nz = 1;
	ConstructVolumeStruct2D(&g->g2, nx, ny, 0);
	g->nxpad = g->g2.Nx + g->g2.NxFSkip;
	g->len = g->g2.LengthF;
    }
    g->nx = nx; g->ny = ny; g->nz = nz;
}

static int
run_geometry(int dim, long nx, long ny, long nz)
{
    static const char* slots3[] = {"  ", "xx", "yy", "zz", "xy", "yz", "xz"};
    static const char* slots2[] = {"  ", "xx", "yy", "xy"};
    const char** slots = (dim == 3) ? slots3 : slots2;
    int nslots = (dim == 3) ? 7 : 4;
    Geo g;
    double *x, *y, *lx, *lty, *w, *s, *ref, *e, *le, *swp;
    int k, fails = 0;
    long n, j, ix, r;

    setup(&g, dim, nx, ny, nz);
    n = g.nx * g.ny * g.nz;
    x = malloc(10 * (size_t) g.len * sizeof(double));
    y = x + g.len; lx = y + g.len; lty = lx + g.len; w = lty + g.len;
    s = w + g.len; ref = s + g.len; e = ref + g.len; le = e + g.len;
    swp = le + g.len;

    fill_random(x, g.nx, g.ny, g.nz, g.nxpad);
    fill_random(y, g.nx, g.ny, g.nz, g.nxpad);
    fill_random(w, g.nx, g.ny, g.nz, g.nxpad);

    for (k = 0; k < nslots; ++k) {
	double a, b, rel, maxd = 0.0, scale = 0.0;

	/* 1. adjointness */
	apply(&g, x, lx, slots[k], 'f', 0);
	memset(lty, 0, (size_t) g.len * sizeof(double));
	apply(&g, y, lty, slots[k], 'r', 0);
	a = dot(lx, y, g.nx, g.ny, g.nz, g.nxpad);
	b = dot(x, lty, g.nx, g.ny, g.nz, g.nxpad);
	rel = fabs(a - b) / (fabs(a) + fabs(b) + 1e-300);
	if (rel > 1e-13) {
	    printf("FAIL %dD %ldx%ldx%ld slot '%s': adjointness rel %.3e\n",
		   dim, nx, ny, nz, slots[k], rel);
	    ++fails;
	}

	/* 3. 'r' accumulates */
	apply(&g, y, lty, slots[k], 'r', 0);
	b = dot(x, lty, g.nx, g.ny, g.nz, g.nxpad);
	if (fabs(b - 2.0 * a) > 1e-12 * (fabs(a) + 1.0)) {
	    printf("FAIL %dD slot '%s': 'r' does not accumulate (%.17g vs %.17g)\n",
		   dim, slots[k], b, 2.0 * a);
	    ++fails;
	}

	/* 2. 's' against explicit diag(L^T W L) */
	apply(&g, w, s, slots[k], 's', 0);
	for (r = 0; r < n; ++r) {
	    long rz = r / (g.nx * g.ny), ry = (r / g.nx) % g.ny, rx = r % g.nx;
	    long ir = rz * g.nxpad * g.ny + ry * g.nxpad + rx;
	    double acc = 0.0;

	    for (j = 0; j < g.ny * g.nz; ++j)
		for (ix = 0; ix < g.nx; ++ix) e[j * g.nxpad + ix] = 0.0;
	    e[ir] = 1.0;
	    apply(&g, e, le, slots[k], 'f', 0);   /* column r of L */
	    for (j = 0; j < g.ny * g.nz; ++j) {
		for (ix = 0; ix < g.nx; ++ix) {
		    double v = le[j * g.nxpad + ix];
		    acc += v * v * w[j * g.nxpad + ix];
		}
	    }
	    ref[ir] = acc;
	    if (fabs(acc - s[ir]) > maxd) maxd = fabs(acc - s[ir]);
	    if (fabs(acc) > scale) scale = fabs(acc);
	}
	if (maxd > 1e-13 * (scale + 1.0)) {
	    printf("FAIL %dD %ldx%ldx%ld slot '%s': 's' vs explicit max|d| %.3e\n",
		   dim, nx, ny, nz, slots[k], maxd);
	    ++fails;
	}

	/* in-place call through swpspace must equal out-of-place */
	memcpy(le, x, (size_t) g.len * sizeof(double));
	apply(&g, le, le, slots[k], 'f', swp);
	for (j = 0; j < g.ny * g.nz; ++j)
	    for (ix = 0; ix < g.nx; ++ix)
		if (le[j * g.nxpad + ix] != lx[j * g.nxpad + ix]) {
		    printf("FAIL %dD slot '%s': in-place differs\n", dim, slots[k]);
		    ++fails; j = g.ny * g.nz; break;
		}
    }

    /* pads untouched? (poisoned with 1e300 in fill_random) */
    for (j = 0; j < g.ny * g.nz; ++j)
	for (ix = g.nx; ix < g.nxpad; ++ix)
	    if (lx[j * g.nxpad + ix] != lx[j * g.nxpad + ix]) { /* NaN check */
		printf("FAIL %dD: pad column written with NaN\n", dim);
		++fails;
	    }

    printf("%s %dD %ldx%ldx%ld: %d slots (adjoint, accumulate, 's', in-place)\n",
	   fails ? "FAIL" : "PASS", dim, nx, ny, nz, nslots);
    free(x);
    return fails;
}

int
main(void)
{
    int fails = 0;

    srand(12345);
    fails += run_geometry(3, 6, 5, 8);   /* golden test 1 shape */
    fails += run_geometry(3, 5, 8, 4);   /* golden test 2 shape (odd nx) */
    fails += run_geometry(3, 4, 3, 2);   /* nz = 2: z neighbours coincide */
    fails += run_geometry(3, 3, 3, 3);
    fails += run_geometry(3, 4, 1, 3);   /* ny = 1: y slots vanish */
    fails += run_geometry(3, 2, 2, 2);
    fails += run_geometry(2, 6, 5, 1);
    fails += run_geometry(2, 5, 8, 1);
    fails += run_geometry(2, 7, 3, 1);
    printf("%s: test_derivatives\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}
