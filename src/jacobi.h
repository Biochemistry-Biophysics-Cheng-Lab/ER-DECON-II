/*
 *  jacobi.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 7/1/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include     "fft_double.h"

void         Jacobi(   double    *X,
					   double    *B,
					   double    *Hd,
					   double    *PCF,
					   double    *PCF2,
					   double    *Wm,
					   double    *alpham,
					   double     *Wd,
					   double     lamf,
					   double     *lam,
					   double     *alpha,
					   double     *swpspace,
					   double      *D,
					   double      tol,
					   long       NIter,
					   long       NIterT,
					   double     *omega,
					   Geometry   *G,
			                   void       (*res_reporter)(long iiters, long iitert, double res, void* res_reporter_closure),
			                   void*      res_reporter_closure );



void	Update (	double       *X,
				double       *R,
				double       *PCF,
				double       *Dr,
				FFTDoublePlan    pfR,
				FFTDoublePlan    prR,
				double       omega,
				long         Nx,
				long         NxSkip,
				long	     Ny,
				long         Nz        );


