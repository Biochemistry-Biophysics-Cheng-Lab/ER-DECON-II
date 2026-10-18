/*
 * LineSearch2D2n.c -- RECONSTRUCTED.  See LineSearch2D2n.h and solver2014.h.
 * ERDecon reference build, 2026-09-21.
 */
#include "LineSearch2D2n.h"
#include "solver2014.h"

int
LineSearch2D2n(
    int nx, int ny,
    double lamf, double lamratio, double lampos, double alpha, double laml2,
    double lammask, int regularization_type,
    const double* X, const double* R, const double* TF, const double* Imn,
    const double* mask, double* p_alphaw,
    ResReporterCosts* p_pre_costs, ResReporterCosts* p_post_costs,
    double* work
)
{
    Solver2014Params p;

    p.lamf = lamf; p.lamratio = lamratio; p.lampos = lampos; p.alpha = alpha;
    p.laml2 = laml2; p.lammask = lammask; p.delfrac = 1.0;
    p.regularization_type = regularization_type;
    return solver2014_linesearch_w(2, nx, ny, 1, &p, X, R, TF, Imn, mask,
				 p_alphaw, p_pre_costs, p_post_costs, work);
}
