#ifndef UCSF_MSG_JACOBINL3DDU_H
#define UCSF_MSG_JACOBINL3DDU_H

#include "res_reporter.h"

/*
 * Imn, X, B, and mask must be dimensioned 2 * (nx / 2 + 1) by ny by nz (mask
 * is only used if lammask is greater than zero).  Only the first nx elements
 * in x need to be initialized.  PCF must be complex-valued and dimensioned
 * (nx / 2 + 1) by ny by nz.  TF must be complex-valued and dimensioned
 * (nx / 2 + 1) by ny by nz.  work must point to
 * (2 * (nx / 2 + 1) + max(3 * nx, 4 * (nx / 2 + 1))) * ny * nz doubles.
 */
int
JacobiNL3Ddu(
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
	int nz,
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
