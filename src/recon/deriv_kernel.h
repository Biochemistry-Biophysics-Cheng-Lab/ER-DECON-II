/*
 * deriv_kernel.h -- RECONSTRUCTED.  Shared kernel behind Derivative() (3D)
 * and Derivative2D() (2D).
 *
 * The operators are those of the PNAS 2013 SI, eqs. S1-S6, with PERIODIC
 * (wrap-around) boundaries and the mixed differences aligned BACKWARD.  This
 * was established numerically against the four golden costs in
 * JacobiNL{2,3}DduTests.c (see docs/SPECIFICATION.md,
 * section 2): with these conventions the reported costs are matched to
 * 1e-13 .. 1e-16; every other boundary rule or alignment misses by >= 1e-5.
 *
 *   identity  "  "   (L X)(r) = X(r)
 *   pure      "aa"   (L X)(r) = X(r - e_a) - 2 X(r) + X(r + e_a)
 *   mixed     "ab"   (L X)(r) = X(r) - X(r - e_a) - X(r - e_b) + X(r - e_a - e_b)
 *
 * Modes (contract fixed by the delivered callers in regc.c and misc.c):
 *
 *   'f'  R  = L X           overwrite
 *   'r'  R += L^T X         ACCUMULATE.  RegNL() and AddDuchon-style loops
 *                           memset R to zero and then call 'r' once per
 *                           slot, so 'r' must add into R.
 *   's'  R  = (L o L)^T X   overwrite; squared stencil coefficients applied
 *                           in the adjoint direction.  This is what Regs()
 *                           needs to build diag(L^T W L) = (L o L)^T W.
 *
 * Only the first Nx entries of every padded row are data; the NxFSkip pad
 * entries are never read and never written.
 *
 * ERDecon reference build, 2026-09-21.
 */
#ifndef ERDECON_RECON_DERIV_KERNEL_H
#define ERDECON_RECON_DERIV_KERNEL_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Slot codes.  The pair of characters at dir maps to one of these. */
enum deriv_slot {
    DERIV_ID = 0,
    DERIV_XX, DERIV_YY, DERIV_ZZ,
    DERIV_XY, DERIV_YZ, DERIV_XZ,
    DERIV_BAD
};

static enum deriv_slot
deriv_parse_slot(const char* dir)
{
    char a = dir[0], b = dir[1];

    if (a == ' ' && b == ' ') return DERIV_ID;
    if (a == 'x' && b == 'x') return DERIV_XX;
    if (a == 'y' && b == 'y') return DERIV_YY;
    if (a == 'z' && b == 'z') return DERIV_ZZ;
    if (a == 'x' && b == 'y') return DERIV_XY;
    if (a == 'y' && b == 'z') return DERIV_YZ;
    if (a == 'x' && b == 'z') return DERIV_XZ;
    return DERIV_BAD;
}

/* Periodic neighbour tables for an axis of length n: prev[i], next[i]. */
static void
deriv_wrap_tables(long n, long* prev, long* next)
{
    long i;

    for (i = 0; i < n; ++i) {
	prev[i] = (i == 0) ? n - 1 : i - 1;
	next[i] = (i == n - 1) ? 0 : i + 1;
    }
}

/*
 * Apply the operator.  X and R must not alias (the public wrappers copy X
 * into swpspace when they do).  nxpad is the padded row length.
 */
static void
deriv_apply(
    const double* X,
    double* R,
    long nx, long ny, long nz, long nxpad,
    enum deriv_slot slot,
    char mode,
    const char* who
)
{
    long* xm; long* xp; long* ym; long* yp; long* zm; long* zp;
    long iz;
    const ptrdiff_t sy = nxpad;
    const ptrdiff_t sz = nxpad * ny;
    int acc = (mode == 'r');
    double s_side = 1.0, s_mid = 4.0;   /* squared coefficients of (1,-2,1) */

    if (mode != 'f' && mode != 'r' && mode != 's') {
	(void) fprintf(stderr, "%s: unknown mode '%c'\n", who, mode);
	abort();
    }
    if (slot == DERIV_BAD) {
	(void) fprintf(stderr, "%s: unknown direction\n", who);
	abort();
    }

    /* identity: no neighbours needed */
    if (slot == DERIV_ID) {
	long j;

	for (j = 0; j < ny * nz; ++j) {
	    const double* xc = X + j * sy;
	    double* rc = R + j * sy;
	    long ix;

	    if (acc) {
		for (ix = 0; ix < nx; ++ix) rc[ix] += xc[ix];
	    } else {
		for (ix = 0; ix < nx; ++ix) rc[ix] = xc[ix];
	    }
	}
	return;
    }

    /*
     * Degenerate axis lengths.  With n == 1 every difference along that
     * axis vanishes identically; with n == 2 the two neighbours of a pure
     * second difference coincide, so the stencil is (-2, +2) and its
     * squared coefficients are (4, 4) rather than (1, 4, 1).  'f' and 'r'
     * get this right automatically; 's' needs the explicit coefficients.
     */
    {
	long na = 0, nb = 0;

	switch (slot) {
	case DERIV_XX: na = nx; break;
	case DERIV_YY: na = ny; break;
	case DERIV_ZZ: na = nz; break;
	case DERIV_XY: na = nx; nb = ny; break;
	case DERIV_YZ: na = ny; nb = nz; break;
	case DERIV_XZ: na = nx; nb = nz; break;
	default: break;
	}
	if (na == 1 || nb == 1) {
	    if (!acc) {
		long j, ix;

		for (j = 0; j < ny * nz; ++j)
		    for (ix = 0; ix < nx; ++ix) R[j * sy + ix] = 0.0;
	    }
	    return;
	}
	if (nb == 0 && na == 2 && mode == 's') {
	    /* coincident neighbours: (-2,+2)^2 = (4,4); the two side terms
	     * below both address the same voxel, so 2 + 2 = 4 */
	    s_side = 2.0; s_mid = 4.0;
	}
    }

    xm = (long*) malloc(2 * (size_t) (nx + ny + nz) * sizeof(long));
    if (xm == 0) {
	(void) fprintf(stderr, "%s: out of memory\n", who);
	abort();
    }
    xp = xm + nx; ym = xp + nx; yp = ym + ny; zm = yp + ny; zp = zm + nz;
    deriv_wrap_tables(nx, xm, xp);
    deriv_wrap_tables(ny, ym, yp);
    deriv_wrap_tables(nz, zm, zp);

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(iz) schedule(static)
#endif
    for (iz = 0; iz < nz; ++iz) {
	long iy;

	for (iy = 0; iy < ny; ++iy) {
	    const double* x0 = X + iz * sz + iy * sy;          /* this row */
	    double* rc = R + iz * sz + iy * sy;
	    const double* xym = X + iz * sz + ym[iy] * sy;     /* row y-1 */
	    const double* xyp = X + iz * sz + yp[iy] * sy;     /* row y+1 */
	    const double* xzm = X + zm[iz] * sz + iy * sy;     /* row z-1 */
	    const double* xzp = X + zp[iz] * sz + iy * sy;     /* row z+1 */
	    const double* xymzm = X + zm[iz] * sz + ym[iy] * sy;
	    const double* xypzp = X + zp[iz] * sz + yp[iy] * sy;
	    long ix;

	    switch (slot) {
	    case DERIV_XX:
		if (mode == 's') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = s_side * (x0[xm[ix]] + x0[xp[ix]]) + s_mid * x0[ix];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += x0[xm[ix]] - 2.0 * x0[ix] + x0[xp[ix]];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[xm[ix]] - 2.0 * x0[ix] + x0[xp[ix]];
		}
		break;
	    case DERIV_YY:
		if (mode == 's') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = s_side * (xym[ix] + xyp[ix]) + s_mid * x0[ix];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += xym[ix] - 2.0 * x0[ix] + xyp[ix];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = xym[ix] - 2.0 * x0[ix] + xyp[ix];
		}
		break;
	    case DERIV_ZZ:
		if (mode == 's') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = s_side * (xzm[ix] + xzp[ix]) + s_mid * x0[ix];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += xzm[ix] - 2.0 * x0[ix] + xzp[ix];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = xzm[ix] - 2.0 * x0[ix] + xzp[ix];
		}
		break;
	    case DERIV_XY:
		/* f: X(r) - X(r-ex) - X(r-ey) + X(r-ex-ey)
		 * r: X(r) - X(r+ex) - X(r+ey) + X(r+ex+ey)   (accumulate)
		 * s: X(r) + X(r+ex) + X(r+ey) + X(r+ex+ey) */
		if (mode == 'f') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] - x0[xm[ix]] - xym[ix] + xym[xm[ix]];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += x0[ix] - x0[xp[ix]] - xyp[ix] + xyp[xp[ix]];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] + x0[xp[ix]] + xyp[ix] + xyp[xp[ix]];
		}
		break;
	    case DERIV_YZ:
		if (mode == 'f') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] - xym[ix] - xzm[ix] + xymzm[ix];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += x0[ix] - xyp[ix] - xzp[ix] + xypzp[ix];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] + xyp[ix] + xzp[ix] + xypzp[ix];
		}
		break;
	    case DERIV_XZ:
		if (mode == 'f') {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] - x0[xm[ix]] - xzm[ix] + xzm[xm[ix]];
		} else if (acc) {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] += x0[ix] - x0[xp[ix]] - xzp[ix] + xzp[xp[ix]];
		} else {
		    for (ix = 0; ix < nx; ++ix)
			rc[ix] = x0[ix] + x0[xp[ix]] + xzp[ix] + xzp[xp[ix]];
		}
		break;
	    default:
		break;
	    }
	}
    }

    free(xm);
}

#endif /* include guard */
