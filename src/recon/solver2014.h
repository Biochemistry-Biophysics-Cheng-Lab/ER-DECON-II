/*
 * solver2014.h -- RECONSTRUCTED.  Shared internals of the 2014 ER-Decon II
 * solver step: the objective J, its gradient/diagonal (Res3D, Res2D) and
 * the exact line search (LineSearch3D2n, LineSearch2D2n).  The four public
 * routines are thin wrappers around the dimension-generic code here.
 *
 * What is implemented (docs/SPECIFICATION.md sections 2 and 4):
 *
 *   J(X) = sum_r (h*X - f)^2
 *        + (lamf/alpha) sum_r [ rho(alpha E2(r)) + laml2 E2(r) ]
 *        + lamf lampos  sum_r min(X(r),0)^2
 *        + lamf lammask sum_r (M(r) X(r))^2                (spatial constraint)
 *
 *   rho(u) = log(1+u)  (regularization_type 1)   sqrt(1+u)  (type 2)
 *   E2(r)  = lamratio X^2 + (Lxx X)^2 + (Lyy X)^2 + d^4 (Lzz X)^2
 *            + 4 (Lxy X)^2 + 4 d^2 [ (Lyz X)^2 + (Lxz X)^2 ],   d = delfrac
 *
 * Update direction (Res*):
 *   E2g  = E2 with lamratio max(X,0)^2 and mixed weights 2, 2 d^2
 *   W    = rho'(alpha E2g) + laml2
 *   R    = H^T H X + lamf [ lamratio W X + sum_k w_k L_k^T (W L_k X)
 *                           + lampos min(X,0) + lammask M X ]
 *   Dreg = lamratio W + sum_k w_k (L_k o L_k)^T W
 *   D    = lamf Dreg + (max Dreg / max E2g) (lamf (lampos [X<0] + lammask M) + Diagd)
 *
 * None of this is the exact gradient/Hessian diagonal of J (mixed weights,
 * laml2 without 1/alpha, max(X,0) in E2g, M vs M^2, the max-ratio scaling).
 * It is what the original does -- verified bit-identical against the
 * original Priism core2_decon on 57 Tier-1 runs + 50-iteration, mask and 2D
 * runs, and against all six golden tests -- and is reproduced deliberately.
 *
 * ERDecon reference build, 2026-09-21 / 2026-09-22.
 */
#ifndef ERDECON_RECON_SOLVER2014_H
#define ERDECON_RECON_SOLVER2014_H

#include "res_reporter.h"

typedef struct Solver2014Params {
    double lamf;
    double lamratio;
    double lampos;
    double alpha;
    double laml2;
    double lammask;
    double delfrac;      /* 1.0 in 2D */
    int regularization_type;
} Solver2014Params;

/*
 * dim is 2 or 3 (nz must be 1 for dim == 2).  Arrays are the padded
 * half-complex layout 2*(nx/2+1) x ny x nz.  Both routines share one
 * internal scratch buffer of two padded volumes (kept between calls).
 * solver2014_res ignores work.  solver2014_linesearch_w keeps its three
 * E2 volumes (unpadded, 3*nx*ny*nz doubles) in work, as JacobiNL{2,3}Ddu.h
 * guarantee; solver2014_linesearch is the same with its own allocation.
 */
int solver2014_res(
    int dim, int nx, int ny, int nz,
    const Solver2014Params* p,
    const double* X,
    const double* TF,
    const double* mask,
    double Diagd,
    double* R,
    double* D,
    double* work
);

int solver2014_linesearch(
    int dim, int nx, int ny, int nz,
    const Solver2014Params* p,
    const double* X,
    const double* U,
    const double* TF,
    const double* Imn,
    const double* mask,
    double* p_alphaw,
    ResReporterCosts* p_pre_costs,
    ResReporterCosts* p_post_costs
);

int solver2014_linesearch_w(
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
);

/*
 * Debug hook: if non-NULL, called by solver2014_linesearch with the Newton
 * iterate history (a_k, dJ/da(a_k), d2J/da2(a_k)) once per iteration.
 */
extern void (*solver2014_linesearch_trace)(int it, double a, double g, double h);

#endif /* include guard */
