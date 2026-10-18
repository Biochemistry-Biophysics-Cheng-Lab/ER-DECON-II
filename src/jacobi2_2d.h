/*
 *  jacobi2.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 6/9/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "fft_double.h"
#include "jacobi2_common.h"

/*
 * swpspace must point to space for 11 * G->LengthF doubles.  X, TF, B, and
 * PCF must point to G->LengthF doubles.
 */
void JacobiNL2D(
    double *X,
    const double *TF,
    double Diagd,
    const double *B,
    const double *PCF,
    double lamf, 
    const double *lamratio,
    double lampos,
    const double *beta,
    const double *alpha,
    const Geometry2D *G,
    double *swpspace,
    double *omega,
    long NIter,
    double tol,
    void (*res_reporter)(
	long iiters, long iitert, double res, void* res_reporter_closure
    ),
    void *res_reporter_closure,
    double laml2
);

/*
 * swpspace must point to space for 8 * G->LengthF doubles.  X, U, TF, R, B,
 * and D must point to G->LengthF doubles.
 */
void JacobiNL2Dstep(
    double *X,
    const double *U,
    const double *TF,
    double Diagd,
    double *R,
    const double *B,
    double *D,
    double lamf, 
    const double *lamratio,
    double lampos,
    const double *beta,
    const double *alpha,
    const Geometry2D *G,
    double *swpspace,
    FFTDoublePlan pfR,
    FFTDoublePlan prR,
    double *omega,
    long iiters,
    long *p_iitert,
    double *p_res,
    void (*res_reporter)(long iiters, long iitert, double res, void*),
    void *res_reporter_closure,
    double laml2
);

/*
 * swpspace must point to space for 7 * G->LengthF doubles.  X, TF, R, and D
 * must point to G->LengthF doubles.
 */
void Res2D_old (
    const double *X,
    const double *TF,
    double Diagd,
    double *R,
    double *D,
    double lamf,
    const double *lamratio,
    double lampos,
    const double *beta,
    const double *alpha,
    const Geometry2D *G,
    double *swpspace,
    FFTDoublePlan pfR,
    FFTDoublePlan prR,
    double laml2
);
