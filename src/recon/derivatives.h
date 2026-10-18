/*
 * derivatives.h -- RECONSTRUCTED (implemented 2026-09-21).
 *
 * Declares Derivative(), the finite-difference operator behind every
 * regularizer term in ER-Decon I and II.  The signature is fixed by the
 * delivered call sites (applyreg.c, regc.c, misc.c):
 *
 *     Derivative(X, R, G, dir, mode, swpspace)
 *
 *   X, R      padded real arrays of G->LengthF doubles
 *             (row stride Nx + NxFSkip; only the first Nx of each row are data)
 *   dir       points at TWO characters selecting the operator, taken from
 *             "    xxyyzzxyyzxz":  "  " identity, "xx" d2/dx2, "yy", "zz",
 *             "xy", "yz", "xz" (mixed second derivatives)
 *   mode      'f' forward operator L, 'r' its adjoint L^T, 's' symmetric
 *   swpspace  scratch of at least G->LengthF doubles
 *
 * Hard constraint (from Regsum in regc.c, which forms L^T (W L g) by calling
 * 'f' then 'r'):  'r' must be the exact adjoint of 'f'.
 *
 * Settled 2026-09-21 (docs/SPECIFICATION.md, section 2):
 * SI S1-S6 stencils, backward-aligned mixed differences, periodic wrap.
 * 'r' ACCUMULATES into R (RegNL/AddDuchon zero R and call 'r' per slot);
 * 's' is the squared stencil in the adjoint direction.  See deriv_kernel.h.
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_DERIVATIVES_H
#define ERDECON_RECON_DERIVATIVES_H

#include "volumestruct.h"

void Derivative(
    const double* X,
    double* R,
    const Geometry* G,
    const char* dir,
    char mode,
    double* swpspace
);

#endif /* include guard */
