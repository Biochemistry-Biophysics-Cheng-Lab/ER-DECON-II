/*
 *  regc.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 6/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


/*
 * swpspace must point to space for 5 * G->LengthF doubles. X, R, and D must
 * point to G->LengthF doubles.
 */
void Regsum(
    const double *X,
    double *R,
    double *D,
    const double *lamratio,
    double laml2,
    double lampos,
    double delfrac,
    const double *beta,
    const double *alpha,
    const Geometry *G, 
    double *swpspace
);

/*
 * swpspace must point to space for 3 * G->LengthF doubles.  X, R, and D must
 * point to G->LengthF doubles.
 */
void Reg(
    const double *X,
    double *R,
    double *D,
    const char *dir,
    long NF,
    const double *lamc,
    double laml2,
    double beta,
    double alpha,
    const Geometry *G, 
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X, R, and D must
 * point to G->LengthF doubles.
 */
void RegNL(
    const double *X,
    const double *D,
    double *R,
    const Geometry *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X and R must point
 * to G->LengthF doubles.
 */
void Regr(
    const double *X,
    double *R,
    const Geometry *G,
    const char *dir,
    int	NF,
    const double *lamc,
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X and R must point
 * to G->LengthF doubles.
 */
void Regs (
    const double *X,
    double *R,
    const Geometry *G,
    const char *dir,
    int	NF,
    const double *lamc,
    double *swpspace
);
