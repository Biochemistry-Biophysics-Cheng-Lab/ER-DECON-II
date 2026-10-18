/*
 * solver2014.c -- RECONSTRUCTED.  See solver2014.h.
 * ERDecon reference build, 2026-09-21.
 */
#include "solver2014.h"
#include "volumestruct.h"
#include "derivatives.h"
#include "volumestruct2d.h"
#include "derivatives2d.h"
#include "fft_double.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void (*solver2014_linesearch_trace)(int it, double a, double g, double h) = 0;

/*
 * OpenMP (2026-09-23): every voxel loop runs over rows j in parallel.  Sums
 * keep their exact serial order -- each row's partial sum is computed as
 * before and stored, then the rows are added in order j = 0, 1, ... -- and
 * maxima are order-independent, so results are bit-identical for any
 * thread count (and to the serial code used while reconstructing the
 * solver).
 */
#if defined(_OPENMP) && defined(USE_OPENMP)
#define S14_PRAGMA(x) _Pragma(#x)
#else
#define S14_PRAGMA(x)
#endif
#define S14_PARALLEL_ROWS S14_PRAGMA(omp parallel for schedule(static))

/*
 * Scratch (2026-09-24): Res and the line search each need two padded
 * volumes of their own.  They are never active at the same time, so they
 * share one buffer that is allocated once and kept for the whole run --
 * allocating and first-touching ~2 volumes per call cost ~25 % of the run
 * time at 256^3 (page faults), and the line search's three E2 volumes now
 * live in the caller's work area (see solver2014_linesearch_w).  Not
 * reentrant: the front ends run one solver at a time.
 */
static double* s14_scratch = 0;
static size_t s14_scratch_len = 0;

static void
s14_scratch_release(void)
{
    if (s14_scratch != 0) fftw_free(s14_scratch);
    s14_scratch = 0;
    s14_scratch_len = 0;
}

static double*
s14_scratch_get(size_t len)
{
    static int registered = 0;

    if (len > s14_scratch_len) {
	s14_scratch_release();
	s14_scratch = (double*) fftw_malloc(len * sizeof(double));
	if (s14_scratch == 0) return 0;
	s14_scratch_len = len;
	if (!registered) { (void) atexit(s14_scratch_release); registered = 1; }
    }
    return s14_scratch;
}

/*
 * Spatial-constraint term in the diagonal D: lamf*lammask*M, scaled like the
 * positivity term (step 6).  With the 2026-09-22 D formula this reproduces
 * golden test 3 (2D and 3D) to 1e-15; without it the match is 1e-2.
 * (The 2026-09-21 conclusion "no mask term in D" was an artefact of the
 * then-missing max(Dreg)/max(E2g) scaling.)
 */
/*
 * Spatial constraint.  The COST term is lamf*lammask*sum_r (M(r) X(r))^2 --
 * the mask enters squared ("squared product of the object guess and the
 * mask", core2_decon.hlp) -- but the GRADIENT and DIAGONAL terms use M
 * linearly (lamf*lammask*M*X and lamf*lammask*M), exactly as the help text
 * states them.  Another deliberate inconsistency of the 2014 code, like the
 * mixed weights and laml2.  Golden test 3 (binary mask) cannot tell these
 * apart; a fractional mask against the original core2_decon does: this
 * combination reproduces the 10-iteration trace, all others miss it by
 * 1-9% (2026-09-22; docs/SPECIFICATION.md section 4).
 */
#define MASKW(m)   ((m) * (m))   /* cost and line search */
#ifndef MASKW_R
#define MASKW_R(m) (m)           /* gradient R */
#endif
#ifndef MASKW_D
#define MASKW_D(m) (m)           /* diagonal D */
#endif

#ifndef SOLVER2014_MASK_IN_D
#define SOLVER2014_MASK_IN_D 1
#endif

#define S14_MAX_SLOTS 6

typedef struct Grid {
    int dim;
    long nx, ny, nz, nxpad, L, N;
    Geometry g3;
    Geometry2D g2;
    int nslots;
    const char* dirs;                 /* 2 chars per slot */
    double wgrad[S14_MAX_SLOTS];      /* 2013 weights: mixed 2, 2 d^2 */
    double wcost[S14_MAX_SLOTS];      /* cost weights: mixed 4, 4 d^2 */
} Grid;

static int
grid_setup(Grid* g, int dim, int nx, int ny, int nz, double delfrac)
{
    double d2 = delfrac * delfrac;

    memset(g, 0, sizeof(*g));
    g->dim = dim;
    g->nx = nx; g->ny = ny;
    if (dim == 3) {
	g->nz = nz;
	ConstructVolumeStruct(&g->g3, nx, ny, nz, 0);
	g->nxpad = g->g3.Nx + g->g3.NxFSkip;
	g->L = g->g3.LengthF;
	g->nslots = 6;
	g->dirs = "xxyyzzxyyzxz";
	g->wgrad[0] = 1.0; g->wgrad[1] = 1.0; g->wgrad[2] = d2 * d2;
	g->wgrad[3] = 2.0; g->wgrad[4] = 2.0 * d2; g->wgrad[5] = 2.0 * d2;
	g->wcost[0] = 1.0; g->wcost[1] = 1.0; g->wcost[2] = d2 * d2;
	g->wcost[3] = 4.0; g->wcost[4] = 4.0 * d2; g->wcost[5] = 4.0 * d2;
    } else if (dim == 2 && nz == 1) {
	g->nz = 1;
	ConstructVolumeStruct2D(&g->g2, nx, ny, 0);
	g->nxpad = g->g2.Nx + g->g2.NxFSkip;
	g->L = g->g2.LengthF;
	g->nslots = 3;
	g->dirs = "xxyyxy";
	g->wgrad[0] = 1.0; g->wgrad[1] = 1.0; g->wgrad[2] = 2.0;
	g->wcost[0] = 1.0; g->wcost[1] = 1.0; g->wcost[2] = 4.0;
    } else {
	return -1;
    }
    g->N = g->nx * g->ny * g->nz;
    return 0;
}

static void
grid_deriv(const Grid* g, const double* X, double* R, int k, char mode)
{
    if (g->dim == 3) {
	Derivative(X, R, &g->g3, g->dirs + 2 * k, mode, 0);
    } else {
	Derivative2D(X, R, &g->g2, g->dirs + 2 * k, mode, 0);
    }
}

/* ------------------------------------------------------------------------
 * FFT helpers.  Plans are created per call with FFTW_ESTIMATE (the
 * "no_overwrite" variants) on the buffer they are executed on; fftw_malloc'd
 * buffers share an alignment, so one plan serves several of them.
 * ---------------------------------------------------------------------- */
typedef struct Plans {
    FFTDoublePlan f, r;
} Plans;

static int
plans_make(Plans* p, const Grid* g, double* buf)
{
    if (g->dim == 3) {
	p->f = plan_no_overwrite_real_complex_3d_fft_double(
	    (int) g->nx, (int) g->ny, (int) g->nz, buf, buf);
	p->r = plan_no_overwrite_complex_real_3d_fft_double(
	    (int) g->nx, (int) g->ny, (int) g->nz, buf, buf);
    } else {
	p->f = plan_no_overwrite_real_complex_2d_fft_double(
	    (int) g->nx, (int) g->ny, buf, buf);
	p->r = plan_no_overwrite_complex_real_2d_fft_double(
	    (int) g->nx, (int) g->ny, buf, buf);
    }
    return (p->f == 0 || p->r == 0) ? -1 : 0;
}

static void
plans_free(Plans* p)
{
    if (p->r != 0) destroy_fft_double_plan(p->r);
    if (p->f != 0) destroy_fft_double_plan(p->f);
    p->f = p->r = 0;
}

static void
fft_fwd(const Grid* g, const Plans* p, double* buf)
{
    if (g->dim == 3) compute_real_complex_3d_fft_double(p->f, buf, buf);
    else             compute_real_complex_2d_fft_double(p->f, buf, buf);
}

static void
fft_inv(const Grid* g, const Plans* p, double* buf)
{
    if (g->dim == 3) compute_complex_real_3d_fft_double(p->r, buf, buf);
    else             compute_complex_real_2d_fft_double(p->r, buf, buf);
}

/*
 * buf <- (1/N) * IFFT( TF * FFT(buf) )         if squared == 0   (h * X)
 * buf <- (1/N) * IFFT( |TF|^2 * FFT(buf) )     if squared != 0   (H^T H X)
 * The 1/N normalises FFTW's unnormalised round trip; it is folded into the
 * spectral multiply as in Res3D_old().
 */
static void
apply_tf(const Grid* g, const Plans* p, const double* TF, double* buf, int squared)
{
    double arrsz = 1.0 / (double) g->N;
    long j;

    fft_fwd(g, p, buf);
    if (squared) {
S14_PARALLEL_ROWS
	for (j = 0; j < g->L / 2; ++j) {
	    double scl = (TF[2*j] * TF[2*j] + TF[2*j+1] * TF[2*j+1]) * arrsz;

	    buf[2*j] *= scl;
	    buf[2*j+1] *= scl;
	}
    } else {
S14_PARALLEL_ROWS
	for (j = 0; j < g->L / 2; ++j) {
	    double rp = buf[2*j], ip = buf[2*j+1];
	    double tr = TF[2*j] * arrsz, ti = TF[2*j+1] * arrsz;

	    buf[2*j] = rp * tr - ip * ti;
	    buf[2*j+1] = rp * ti + ip * tr;
	}
    }
    fft_inv(g, p, buf);
}

/* ------------------------------------------------------------------------
 * rho and derivatives
 * ---------------------------------------------------------------------- */
static double
rho(int type, double u)
{
    return (type == 1) ? log(1.0 + u) : sqrt(1.0 + u);
}

static double
drho(int type, double u)
{
    return (type == 1) ? 1.0 / (1.0 + u) : 0.5 / sqrt(1.0 + u);
}

static double
ddrho(int type, double u)
{
    if (type == 1) {
	double t = 1.0 + u;
	return -1.0 / (t * t);
    } else {
	double t = 1.0 + u;
	return -0.25 / (t * sqrt(t));
    }
}

/* ------------------------------------------------------------------------
 * Res: R = H^T H X + lamf [ lamratio W X + sum_k w_k L_k^T (W L_k X)
 *                           + lampos min(X,0) + lammask M X ]
 *      D = lamf [ lamratio W + sum_k w_k (L_k o L_k)^T W
 *                 ] + (max Dreg/max E2g) (lamf (lampos [X<0] + lammask M) + Diagd)
 *      W = rho'(alpha E2g) + laml2,  E2g with the 2013 (gradient) weights
 * work: unused (see the allocation note in the body).
 * ---------------------------------------------------------------------- */
int
solver2014_res(
    int dim, int nx, int ny, int nz,
    const Solver2014Params* p,
    const double* X,
    const double* TF,
    const double* mask,
    double Diagd,
    double* R,
    double* D,
    double* work
)
{
    Grid g;
    Plans pl = {0, 0};
    double* W;
    double* T;
    int k;
    long j;
    int use_mask;
    double maxE2 = 0.0;

    if (grid_setup(&g, dim, nx, ny, nz, p->delfrac) != 0) return -1;
    if (p->regularization_type != 1 && p->regularization_type != 2) return -1;
    use_mask = (p->lammask > 0.0 && mask != 0);
    /*
     * Two padded volumes are needed (W and a temporary), but the caller's
     * work area after R and D is only guaranteed to hold
     * max(3*nx, 4*(nx/2+1))*ny*nz - 2*(nx/2+1)*ny*nz doubles
     * (JacobiNL3Ddu.h), which is less than two padded volumes.  Allocate
     * our own rather than rely on callers being generous (the golden tests
     * and core2.c are, but the contract is not).
     */
    (void) work;
    W = s14_scratch_get(2 * (size_t) g.L);
    if (W == 0) return -1;
    T = W + g.L;

    /* 1. R <- H^T H X */
    (void) memcpy(R, X, (size_t) g.L * sizeof(double));
    if (plans_make(&pl, &g, R) != 0) { plans_free(&pl); return -1; }
    apply_tf(&g, &pl, TF, R, 1);
    plans_free(&pl);

    /*
     * 2. W <- rho'(alpha E2g) + laml2.
     * The identity slot of E2g uses max(X,0)^2, as the 2013 Regr() did
     * ("tmp = (*Rs > 0.0) ? *Rs : 0.0"), whereas the COST uses X^2 (cost
     * matched to 4e-16 with X^2, off by 1e-5 with max(X,0)^2).  On golden
     * 3D test 1 -- the only one with negative voxels -- this takes the
     * trajectory residual from 3.2e-6 to 3.7e-7.
     */
S14_PARALLEL_ROWS
    for (j = 0; j < g.ny * g.nz; ++j) {
	const double* xc = X + j * g.nxpad;
	double* wc = W + j * g.nxpad;
	long ix;

	for (ix = 0; ix < g.nx; ++ix) {
	    double xp = (xc[ix] > 0.0) ? xc[ix] : 0.0;

	    wc[ix] = p->lamratio * xp * xp;
	}
    }
    for (k = 0; k < g.nslots; ++k) {
	double w = g.wgrad[k];

	grid_deriv(&g, X, T, k, 'f');
S14_PARALLEL_ROWS
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* tc = T + j * g.nxpad;
	    double* wc = W + j * g.nxpad;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix) wc[ix] += w * tc[ix] * tc[ix];
	}
    }
S14_PRAGMA(omp parallel for schedule(static) reduction(max:maxE2))
    for (j = 0; j < g.ny * g.nz; ++j) {
	double* wc = W + j * g.nxpad;
	long ix;

	for (ix = 0; ix < g.nx; ++ix) {
	    if (wc[ix] > maxE2) maxE2 = wc[ix];
	    wc[ix] = drho(p->regularization_type, p->alpha * wc[ix]) + p->laml2;
	}
    }

    /* 3. pointwise terms */
S14_PARALLEL_ROWS
    for (j = 0; j < g.ny * g.nz; ++j) {
	const double* xc = X + j * g.nxpad;
	const double* wc = W + j * g.nxpad;
	const double* mc = use_mask ? mask + j * g.nxpad : 0;
	double* rc = R + j * g.nxpad;
	double* dc = D + j * g.nxpad;
	long ix;

	for (ix = 0; ix < g.nx; ++ix) {
	    double x = xc[ix];
	    double rr = p->lamratio * wc[ix] * x;
	    double dd = p->lamratio * wc[ix];

	    if (x < 0.0) {
		rr += p->lampos * x;
	    }
	    if (use_mask) {
		rr += p->lammask * MASKW_R(mc[ix]) * x;
	    }
	    rc[ix] += p->lamf * rr;
	    dc[ix] = p->lamf * dd;   /* positivity/mask terms added in step 6 */
	}
    }

    /* 4. R += lamf sum_k w_k L_k^T ( W L_k X )     ('r' accumulates) */
    for (k = 0; k < g.nslots; ++k) {
	double s = p->lamf * g.wgrad[k];

	grid_deriv(&g, X, T, k, 'f');
S14_PARALLEL_ROWS
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* wc = W + j * g.nxpad;
	    double* tc = T + j * g.nxpad;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix) tc[ix] *= s * wc[ix];
	}
	grid_deriv(&g, T, R, k, 'r');
    }

    /* 5. D += lamf sum_k w_k (L_k o L_k)^T W */
    for (k = 0; k < g.nslots; ++k) {
	double s = p->lamf * g.wgrad[k];

	grid_deriv(&g, W, T, k, 's');
S14_PARALLEL_ROWS
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* tc = T + j * g.nxpad;
	    double* dc = D + j * g.nxpad;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix) dc[ix] += s * tc[ix];
	}
    }

    /*
     * 6. Positivity, (mask) and data terms of D.  The 2014 solver uses
     *
     *      D = lamf Dreg + (max Dreg / max E2g) * ( lamf lampos [X<0] + Diagd )
     *
     * Dreg = lamratio W + sum_k w_k (L_k o L_k)^T W (regulariser part without
     * lamf), E2g the argument of rho' above, Diagd the caller's mean |TF|^2.
     * (Equivalently, up to the global scale the line search removes: the
     * regulariser diagonal is rescaled so that its maximum equals max E2g,
     * and the positivity term has factor 1, not the 2 of the 2013 code.)
     * Found 2026-09-22 by solving for D pointwise from the original
     * program's iterates; golden 2D 1, 2D 2, 3D 2 then match to 1e-16 and
     * 3D 1 (negative voxels) to ~1e-9.
     */
    {
	double maxDreg = 0.0;
	double ratio;

S14_PRAGMA(omp parallel for schedule(static) reduction(max:maxDreg))
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* dc = D + j * g.nxpad;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix)
		if (dc[ix] > maxDreg) maxDreg = dc[ix];
	}
	maxDreg /= p->lamf;
	ratio = (maxE2 > 0.0) ? maxDreg / maxE2 : 1.0;
S14_PARALLEL_ROWS
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* xc = X + j * g.nxpad;
	    const double* mc = use_mask ? mask + j * g.nxpad : 0;
	    double* dc = D + j * g.nxpad;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix) {
		double dd = 0.0;

		if (xc[ix] < 0.0) dd += p->lampos;
#if SOLVER2014_MASK_IN_D
		if (use_mask) dd += p->lammask * MASKW_D(mc[ix]);
#endif
		dc[ix] += ratio * (p->lamf * dd + Diagd);
	    }
	}
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * Line search: alphaw = argmin_a J(X + a U), Newton on dJ/da from a = 0.
 *
 * Along the line every ingredient of J is a quadratic in a per voxel:
 *   fidelity   sum (hX + a hU)^2  = A + 2 a Bp + a^2 C     (three scalars)
 *   E2(a)      e0 + 2 a e1 + a^2 e2                        (three volumes)
 * so after one pass of FFTs and derivatives each Newton iteration is a
 * single sweep over the voxels with no transforms.
 *
 * Memory: e0, e1, e2 (unpadded, in the caller's work area) and two padded
 * temporaries for L_k X and L_k U (the shared scratch, fftw_malloc'ed so
 * the plans see a common alignment).
 * ---------------------------------------------------------------------- */
typedef struct LineData {
    const Grid* g;
    const Solver2014Params* p;
    const double* X;
    const double* U;
    const double* mask;   /* NULL if unused */
    const double* e0;
    const double* e1;
    const double* e2;
    double A, Bp, C;      /* fidelity quadratic */
    double* part;         /* 4 row partial sums per row (ny*nz rows) */
} LineData;

/* dJ/da and d2J/da2 at a. */
static void
line_derivs(const LineData* d, double a, double* p_g, double* p_h)
{
    const Grid* g = d->g;
    const Solver2014Params* p = d->p;
    int type = p->regularization_type;
    double l2a = p->laml2 / p->alpha;
    double gsum = 0.0, hsum = 0.0;
    double* part = d->part;
    long j, nrow = g->ny * g->nz;

S14_PARALLEL_ROWS
    for (j = 0; j < nrow; ++j) {
	const double* e0 = d->e0 + j * g->nx;
	const double* e1 = d->e1 + j * g->nx;
	const double* e2 = d->e2 + j * g->nx;
	const double* xc = d->X + j * g->nxpad;
	const double* uc = d->U + j * g->nxpad;
	const double* mc = d->mask ? d->mask + j * g->nxpad : 0;
	double grow = 0.0, hrow = 0.0;
	long ix;

	for (ix = 0; ix < g->nx; ++ix) {
	    double e = e0[ix] + a * (2.0 * e1[ix] + a * e2[ix]);
	    double ep = 2.0 * (e1[ix] + a * e2[ix]);
	    double epp = 2.0 * e2[ix];
	    double u = p->alpha * e;
	    double rp = drho(type, u) + l2a;
	    double rpp = ddrho(type, u);
	    double xp = xc[ix] + a * uc[ix];

	    grow += rp * ep;
	    hrow += p->alpha * rpp * ep * ep + rp * epp;
	    if (xp < 0.0) {
		grow += 2.0 * p->lampos * xp * uc[ix];
		hrow += 2.0 * p->lampos * uc[ix] * uc[ix];
	    }
	    if (mc) {
		grow += 2.0 * p->lammask * MASKW(mc[ix]) * xp * uc[ix];
		hrow += 2.0 * p->lammask * MASKW(mc[ix]) * uc[ix] * uc[ix];
	    }
	}
	part[2 * j] = grow;
	part[2 * j + 1] = hrow;
    }
    for (j = 0; j < nrow; ++j) {
	gsum += part[2 * j];
	hsum += part[2 * j + 1];
    }
    *p_g = p->lamf * gsum + 2.0 * d->Bp + 2.0 * a * d->C;
    *p_h = p->lamf * hsum + 2.0 * d->C;
}

/* Cost decomposition at a. */
static void
line_costs(const LineData* d, double a, ResReporterCosts* c)
{
    const Grid* g = d->g;
    const Solver2014Params* p = d->p;
    int type = p->regularization_type;
    double sreg = 0.0, se = 0.0, spos = 0.0, sspa = 0.0;
    double* part = d->part;
    long j, nrow = g->ny * g->nz;

S14_PARALLEL_ROWS
    for (j = 0; j < nrow; ++j) {
	const double* e0 = d->e0 + j * g->nx;
	const double* e1 = d->e1 + j * g->nx;
	const double* e2 = d->e2 + j * g->nx;
	const double* xc = d->X + j * g->nxpad;
	const double* uc = d->U + j * g->nxpad;
	const double* mc = d->mask ? d->mask + j * g->nxpad : 0;
	double rreg = 0.0, re = 0.0, rpos = 0.0, rspa = 0.0;
	long ix;

	for (ix = 0; ix < g->nx; ++ix) {
	    double e = e0[ix] + a * (2.0 * e1[ix] + a * e2[ix]);
	    double xp = xc[ix] + a * uc[ix];

	    rreg += rho(type, p->alpha * e);
	    re += e;
	    if (xp < 0.0) rpos += xp * xp;
	    if (mc) rspa += MASKW(mc[ix]) * xp * xp;
	}
	part[4 * j] = rreg; part[4 * j + 1] = re; part[4 * j + 2] = rpos; part[4 * j + 3] = rspa;
    }
    for (j = 0; j < nrow; ++j) {
	sreg += part[4 * j]; se += part[4 * j + 1];
	spos += part[4 * j + 2]; sspa += part[4 * j + 3];
    }
    c->data_fidelity = d->A + a * (2.0 * d->Bp + a * d->C);
    c->reg = (p->lamf / p->alpha) * sreg;
    c->laml2 = (p->lamf / p->alpha) * p->laml2 * se;
    c->pos = p->lamf * p->lampos * spos;
    c->spa = d->mask ? p->lamf * p->lammask * sspa : 0.0;
    c->total = c->data_fidelity + c->reg + c->laml2 + c->pos + c->spa;
}

int
solver2014_linesearch_w(
    int dim, int nx, int ny, int nz,
    const Solver2014Params* p,
    const double* X,
    const double* U,
    const double* TF,
    const double* Imn,
    const double* mask,
    double* p_alphaw,
    ResReporterCosts* p_pre_costs,
    ResReporterCosts* p_post_costs,
    double* work
)
{
    Grid g;
    Plans pl = {0, 0};
    LineData d;
    double* buf = 0;
    double* ebuf = 0;
    double *tX, *tU, *e0, *e1, *e2;
    int k, it, status = -1;
    long j;
    double a = 0.0;

    if (grid_setup(&g, dim, nx, ny, nz, p->delfrac) != 0) return -1;
    if (p->regularization_type != 1 && p->regularization_type != 2) return -1;

    /* tX, tU: two padded volumes (FFT buffers) from the shared scratch.
       e0, e1, e2: three unpadded nx*ny*nz volumes in the caller's work
       area, which JacobiNL{2,3}Ddu.h guarantee holds max(3 nx, 4 (nx/2+1))
       doubles per row; with work == NULL they are allocated here. */
    buf = s14_scratch_get(2 * (size_t) g.L);
    if (buf == 0) return -1;
    if (work == 0) {
	ebuf = (double*) malloc(3 * (size_t) g.nx * (size_t) (g.ny * g.nz) * sizeof(double));
	if (ebuf == 0) return -1;
	work = ebuf;
    }
    d.part = (double*) malloc(4 * (size_t) (g.ny * g.nz) * sizeof(double));
    if (d.part == 0) { free(ebuf); return -1; }
    tX = buf; tU = tX + g.L;
    e0 = work; e1 = e0 + g.nx * (g.ny * g.nz); e2 = e1 + g.nx * (g.ny * g.nz);
    if (plans_make(&pl, &g, tX) != 0) goto cleanup;

    d.g = &g; d.p = p; d.X = X; d.U = U;
    d.mask = (p->lammask > 0.0 && mask != 0) ? mask : 0;
    d.e0 = e0; d.e1 = e1; d.e2 = e2;

    /* fidelity quadratic: hX = h*X - Imn, hU = h*U */
    (void) memcpy(tX, X, (size_t) g.L * sizeof(double));
    (void) memcpy(tU, U, (size_t) g.L * sizeof(double));
    apply_tf(&g, &pl, TF, tX, 0);
    apply_tf(&g, &pl, TF, tU, 0);
    d.A = d.Bp = d.C = 0.0;
S14_PARALLEL_ROWS
    for (j = 0; j < g.ny * g.nz; ++j) {
	const double* ic = Imn + j * g.nxpad;
	const double* hu = tU + j * g.nxpad;
	double* hx = tX + j * g.nxpad;
	double ra = 0.0, rb = 0.0, rc = 0.0;
	long ix;

	for (ix = 0; ix < g.nx; ++ix) {
	    double r = hx[ix] - ic[ix];

	    ra += r * r;
	    rb += r * hu[ix];
	    rc += hu[ix] * hu[ix];
	}
	d.part[4 * j] = ra; d.part[4 * j + 1] = rb; d.part[4 * j + 2] = rc;
    }
    for (j = 0; j < g.ny * g.nz; ++j) {
	d.A += d.part[4 * j]; d.Bp += d.part[4 * j + 1]; d.C += d.part[4 * j + 2];
    }

    /* E2 quadratic per voxel (cost weights) */
S14_PARALLEL_ROWS
    for (j = 0; j < g.ny * g.nz; ++j) {
	const double* xc = X + j * g.nxpad;
	const double* uc = U + j * g.nxpad;
	double* a0 = e0 + j * g.nx;
	double* a1 = e1 + j * g.nx;
	double* a2 = e2 + j * g.nx;
	long ix;

	for (ix = 0; ix < g.nx; ++ix) {
	    a0[ix] = p->lamratio * xc[ix] * xc[ix];
	    a1[ix] = p->lamratio * xc[ix] * uc[ix];
	    a2[ix] = p->lamratio * uc[ix] * uc[ix];
	}
    }
    for (k = 0; k < g.nslots; ++k) {
	double w = g.wcost[k];

	grid_deriv(&g, X, tX, k, 'f');
	grid_deriv(&g, U, tU, k, 'f');
S14_PARALLEL_ROWS
	for (j = 0; j < g.ny * g.nz; ++j) {
	    const double* lx = tX + j * g.nxpad;
	    const double* lu = tU + j * g.nxpad;
	    double* a0 = e0 + j * g.nx;
	    double* a1 = e1 + j * g.nx;
	    double* a2 = e2 + j * g.nx;
	    long ix;

	    for (ix = 0; ix < g.nx; ++ix) {
		a0[ix] += w * lx[ix] * lx[ix];
		a1[ix] += w * lx[ix] * lu[ix];
		a2[ix] += w * lu[ix] * lu[ix];
	    }
	}
    }

    if (p_pre_costs != 0) line_costs(&d, 0.0, p_pre_costs);

    /*
     * Debug dump (reference build only): with ERDECON_DUMP_DIR set, write
     * X and U of every line-search call as raw little-endian doubles in the
     * padded layout, plus a one-line text header, so that a comparison
     * with another implementation can separate direction errors from
     * step-length errors.
     */
    {
	const char* dir = getenv("ERDECON_DUMP_DIR");

	if (dir != 0) {
	    static int call_no = 0;
	    char name[4096];
	    FILE* fp;

	    (void) snprintf(name, sizeof name, "%s/ls_%03d_X.f64", dir, call_no);
	    if ((fp = fopen(name, "wb")) != 0) {
		(void) fwrite(X, sizeof(double), (size_t) g.L, fp); (void) fclose(fp);
	    }
	    (void) snprintf(name, sizeof name, "%s/ls_%03d_U.f64", dir, call_no);
	    if ((fp = fopen(name, "wb")) != 0) {
		(void) fwrite(U, sizeof(double), (size_t) g.L, fp); (void) fclose(fp);
	    }
	    if (call_no == 0) {
		(void) snprintf(name, sizeof name, "%s/TF.f64", dir);
		if ((fp = fopen(name, "wb")) != 0) {
		    (void) fwrite(TF, sizeof(double), (size_t) g.L, fp); (void) fclose(fp);
		}
		(void) snprintf(name, sizeof name, "%s/Imn.f64", dir);
		if ((fp = fopen(name, "wb")) != 0) {
		    (void) fwrite(Imn, sizeof(double), (size_t) g.L, fp); (void) fclose(fp);
		}
	    }
	    (void) snprintf(name, sizeof name, "%s/ls_%03d_hdr.txt", dir, call_no);
	    if ((fp = fopen(name, "w")) != 0) {
		(void) fprintf(fp, "dim %d nx %d ny %d nz %d nxpad %ld lamf %.17g lamratio %.17g lampos %.17g alpha %.17g laml2 %.17g lammask %.17g delfrac %.17g regtype %d\n",
			       dim, nx, ny, nz, g.nxpad, p->lamf, p->lamratio, p->lampos, p->alpha, p->laml2, p->lammask, p->delfrac, p->regularization_type);
		(void) fclose(fp);
	    }
	    ++call_no;
	}
    }

    /* Newton iteration on dJ/da from a = 0 */
    for (it = 0; it < 100; ++it) {
	double gg, hh, da;

	line_derivs(&d, a, &gg, &hh);
	if (solver2014_linesearch_trace != 0) {
	    (*solver2014_linesearch_trace)(it, a, gg, hh);
	}
	if (!(hh > 0.0)) {
	    /* not locally convex along the line: fall back to a gradient
	     * step scaled by the fidelity curvature, which is always > 0 */
	    hh = 2.0 * d.C;
	    if (!(hh > 0.0)) break;
	}
	da = -gg / hh;
	a += da;
	if (fabs(da) <= 1e-15 * (1.0 + fabs(a))) break;
    }

    *p_alphaw = a;
    if (p_post_costs != 0) line_costs(&d, a, p_post_costs);
    if (getenv("ERDECON_DUMP_DIR") != 0) {
	(void) fprintf(stderr, "erdecon linesearch: alphaw = %.17g\n", a);
    }
    status = 0;

 cleanup:
    plans_free(&pl);
    free(ebuf);
    free(d.part);
    return status;
}

int
solver2014_linesearch(
    int dim, int nx, int ny, int nz,
    const Solver2014Params* p,
    const double* X, const double* U, const double* TF, const double* Imn,
    const double* mask, double* p_alphaw,
    ResReporterCosts* p_pre_costs, ResReporterCosts* p_post_costs
)
{
    return solver2014_linesearch_w(dim, nx, ny, nz, p, X, U, TF, Imn, mask,
				   p_alphaw, p_pre_costs, p_post_costs, 0);
}
