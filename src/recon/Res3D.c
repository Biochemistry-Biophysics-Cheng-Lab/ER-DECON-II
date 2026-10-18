/*
 * Res3D.c -- RECONSTRUCTED.  See Res3D.h and solver2014.h.
 * ERDecon reference build, 2026-09-21.
 */
#include "Res3D.h"
#include "solver2014.h"

int
Res3D(
    int nx, int ny, int nz,
    double lamf, double lamratio, double lampos, double alpha, double laml2,
    double lammask, double delfrac, int regularization_type,
    const double* X, const double* TF, const double* mask, double Diagd,
    double* R, double* D, double* work
)
{
    Solver2014Params p;

    p.lamf = lamf; p.lamratio = lamratio; p.lampos = lampos; p.alpha = alpha;
    p.laml2 = laml2; p.lammask = lammask; p.delfrac = delfrac;
    p.regularization_type = regularization_type;
    return solver2014_res(3, nx, ny, nz, &p, X, TF, mask, Diagd, R, D, work);
}
