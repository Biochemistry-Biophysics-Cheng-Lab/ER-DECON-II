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
void Regsum2D(
    const double *X,
    double *R,
    double *D,
    const double *lamratio,
    double laml2,
    double lampos,
    const double *beta,
    const double *alpha,
    const Geometry2D *G,
    double *swpspace
);

/*
 * swpspace must point to space for 3 * G->LengthF doubles.  X, R, and D must
 * point to G->LengthF doubles.
 */
void Reg2D(
    const double *X,
    double *R,
    double *D,
    const char *dir,
    long NF,
    const double *lamc,
    double laml2,
    double beta,
    double alpha,
    const Geometry2D *G,
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X, R, and D must
 * point to G->LengthF doubles.
 */
void RegNL2D(
    const double *X,
    const double *D,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X and R must point
 * to G->LengthF doubles.
 */
void Regr2D(
    const double *X,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
);

/*
 * swpspace must point to space for 2 * G->LengthF doubles.  X and R must point
 * to G->LengthF doubles.
 */
void Regs2D(
    const double *X,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
);
