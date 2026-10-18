/*
 * LineSearch3D2n.h -- RECONSTRUCTED (implemented 2026-09-21).
 *
 * LineSearch3D2n() returns in *p_alphaw the exact minimiser of the 2014 cost
 * J along X + alpha*R (Newton on dJ/dalpha from 0; the 2013 backtracking
 * rule cannot produce the golden steps, see the reconstruction note s.2).  Fills *p_pre_costs (if non-NULL) with the
 * cost decomposition at X and *p_post_costs at the accepted point.
 * Signature fixed by the single call site in JacobiNL3Ddu.c.
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_LINESEARCH3D2N_H
#define ERDECON_RECON_LINESEARCH3D2N_H

#include "res_reporter.h"

int LineSearch3D2n(
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
    const double* X,
    const double* R,
    const double* TF,
    const double* Imn,
    const double* mask,
    double* p_alphaw,
    ResReporterCosts* p_pre_costs,
    ResReporterCosts* p_post_costs,
    double* work
);

#endif /* include guard */
