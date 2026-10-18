/*
 * derivatives2d.h -- RECONSTRUCTED (implemented 2026-09-21).
 * 2D counterpart of derivatives.h (used by regc2d.c, misc2d.c).  dir selects
 * from "  xxyyxy": "  " identity, "xx", "yy", "xy".  Same 'f'/'r'/'s'
 * modes and the same adjointness requirement.
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_DERIVATIVES2D_H
#define ERDECON_RECON_DERIVATIVES2D_H

#include "volumestruct2d.h"

void Derivative2D(
    const double* X,
    double* R,
    const Geometry2D* G,
    const char* dir,
    char mode,
    double* swpspace
);

#endif /* include guard */
