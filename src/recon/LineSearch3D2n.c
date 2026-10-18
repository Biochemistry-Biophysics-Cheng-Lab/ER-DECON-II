/*
 * LineSearch3D2n.c -- RECONSTRUCTED.  See LineSearch3D2n.h and solver2014.h.
 * The step is the exact minimiser of J along X + a R (Newton on dJ/da),
 * not the 2013 backtracking rule described in the header's first draft.
 * ERDecon reference build, 2026-09-21.
 */
#include "LineSearch3D2n.h"
#include "solver2014.h"

int
LineSearch3D2n(
    int nx, int ny, int nz,
    double lamf, double lamratio, double lampos, double alpha, double laml2,
    double lammask, double delfrac, int regularization_type,
    const double* X, const double* R, const double* TF, const double* Imn,
    const double* mask, double* p_alphaw,
    ResReporterCosts* p_pre_costs, ResReporterCosts* p_post_costs,
    double* work
)
{
    Solver2014Params p;

    p.lamf = lamf; p.lamratio = lamratio; p.lampos = lampos; p.alpha = alpha;
    p.laml2 = laml2; p.lammask = lammask; p.delfrac = delfrac;
    p.regularization_type = regularization_type;
    return solver2014_linesearch_w(3, nx, ny, nz, &p, X, R, TF, Imn, mask,
				 p_alphaw, p_pre_costs, p_post_costs, work);
}
