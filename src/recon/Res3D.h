/*
 * Res3D.h -- RECONSTRUCTED (implemented 2026-09-21).
 *
 * Res3D() computes, for the 2014 ER-Decon II solver (JacobiNL3Ddu.c), the
 * gradient R of the regularised cost at X and the diagonal Hessian
 * approximation D used as the preconditioner.  Signature fixed by the single
 * call site in JacobiNL3Ddu.c.  Its 2013 ancestor is Res3D_old() in
 * jacobi2.c plus Regsum() in regc.c.
 *
 * Array layout: 2*(nx/2+1) x ny x nz padded real arrays (half-complex
 * width), TF complex (nx/2+1) x ny x nz.  work: at least
 * 2*(nx/2+1)*ny*nz doubles (that is what the caller reserves after D).
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_RES3D_H
#define ERDECON_RECON_RES3D_H

int Res3D(
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
    const double* TF,
    const double* mask,
    double Diagd,
    double* R,
    double* D,
    double* work
);

#endif /* include guard */
