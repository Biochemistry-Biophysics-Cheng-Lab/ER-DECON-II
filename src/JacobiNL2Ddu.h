#ifndef UCSF_MSG_JACOBINL2DDU_H
#define UCSF_MSG_JACOBINL2DDU_H

#include "res_reporter.h"

/*
 * Imn, X, B, and mask must be dimensioned 2 * (nx / 2 + 1) by ny (mask is only
 * used if lammask is greater than zero).  Only the first nx elements in x need
 * to be initialized.  PCF must be complex-valued and dimensions (nx / 2 + 1)
 * by ny by nz.  TF must be complex-valued and dimensioned (nx / 2 + 1) by ny
 * by nz.  work must point to
 * (2 * (nx / 2 + 1) + max(3 * nx, 4 * (nx / 2 + 1))) * ny doubles.
 */
int
JacobiNL2Ddu(
    int nx,
    int ny,
    double lamf,
    double lamratio,
    double lampos,
    double alpha,
    double laml2,
    double lammask,
    int regularization_type,
    void (*cost_reporter)(
	int iter,
	const ResReporterCosts* p_costs,
	double rfactor,
	void* closure
    ),
    void* cost_reporter_closure,
    int (*mask_updater)(
	int iter,
	int nx,
	int ny,
	const double* guess,
	double* mask,
	double* work,
	void* closure
    ),
    void* mask_updater_closure,
    int outputs_rfactor,
    int niter,
    const double* B,
    const double* TF,
    const double* Imn,
    double Diagd,
    const double* PCF,
    double* X,
    double* mask,
    double* p_cost,
    double* work
);

#endif /* include guard */
