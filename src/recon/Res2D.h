/*
 * Res2D.h -- RECONSTRUCTED (implemented 2026-09-21).
 * 2D counterpart of Res3D.h; call site in JacobiNL2Ddu.c.
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_RES2D_H
#define ERDECON_RECON_RES2D_H

int Res2D(
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
    const double* TF,
    const double* mask,
    double Diagd,
    double* R,
    double* D,
    double* work
);

#endif /* include guard */
