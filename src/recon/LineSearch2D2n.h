/*
 * LineSearch2D2n.h -- RECONSTRUCTED (implemented 2026-09-21).
 * 2D counterpart of LineSearch3D2n.h; call site in JacobiNL2Ddu.c.
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_LINESEARCH2D2N_H
#define ERDECON_RECON_LINESEARCH2D2N_H

#include "res_reporter.h"

int LineSearch2D2n(
    int nx,
    int ny,
    double lamf,
    double lamratio,
    double lampos,
    double alpha,
    double laml2,
    double lammask,
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
