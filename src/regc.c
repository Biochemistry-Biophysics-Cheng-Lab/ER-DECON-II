/*******************************************************************************************
 *  regc.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 11/23/09.
 *  Copyright 2009 Muthuvel Arigovindan. All rights reserved.
 * 
 *
 ******************************************************************************************/

#include  <stddef.h>
#include  <string.h>
#include  <math.h>
#include  <stdio.h>

/******************************************************************************************/

#include  "volumestruct.h"
#include  "derivatives.h"
#include "regc.h"

/******************************************************************************************/



void Regsum (
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
)
{
	
	
	long      J;
    double    *Rt, *Dt;
	double    lam1, lam2[7L], lamr1, lamr2;
	char      dir1[2L] = "  ";
	char      dir2[14L] = "  xxyyzzxyyzxz";
	
	long       Nx, NxPad, Ny, Nz;
	
	
	
	Nx = G->Nx;
	NxPad = G->Nx + G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	lamr1 =  lamratio[0L];
	lamr2 =  1.0 - lamratio[0L];
	
	lam1     = 1.0;
	
	lam2[0L] = lamratio[1L];
	lam2[1L] = 1.0;
	lam2[2L] = 1.0;
	lam2[3L] = (delfrac*delfrac*delfrac*delfrac);
	lam2[4L] = 2.0;
	lam2[5L] = 2.0*(delfrac*delfrac);
	lam2[6L] = 2.0*(delfrac*delfrac);
	
	
	Rt = swpspace;
	Dt = swpspace + (ptrdiff_t)(G->LengthF);
		
	
	Reg ( X,  Rt, Dt, dir1, 1L,  &lam1, laml2,
		 beta[0L], alpha[0L], G, swpspace+(ptrdiff_t)(2L)*(G->LengthF) );
	
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(J) shared(R, Rt, lamr1, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
	for ( J = 0; J < Ny * Nz; ++J )		{
	        const double* Rtc = Rt + J * (ptrdiff_t) NxPad;
		double* Rc = R + J * (ptrdiff_t) NxPad;
	        long K;

		for ( K = Nx-1L;  K > -1L; K-- )		{
		        *Rc++ = lamr1 * (*Rtc++);
		}
	}
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(J) shared(D, Dt, lamr1, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
	for ( J = 0; J < Ny * Nz; ++J ) {
	        const double* Dtc = Dt + J * (ptrdiff_t) NxPad;
		double* Dc = D + J * (ptrdiff_t) NxPad;
	        long K;

		for ( K = Nx-1L;  K > -1L; K-- )		{
			*Dc++ =  lamr1 * (*Dtc++);
		}
	}

	Reg ( X,  Rt, Dt, dir2, 7L,  lam2, laml2,
	  	 beta[1L], alpha[1L], G, swpspace+(ptrdiff_t)(2L)*(G->LengthF) );

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(J) shared(R, Rt, X, lamr2, lampos, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
	for ( J = 0; J < Ny * Nz; ++J ) {
	        const double* Rtc = Rt + J * (ptrdiff_t) NxPad;
		const double* Xc = X + J * (ptrdiff_t) NxPad;
		double* Rc = R + J * (ptrdiff_t) NxPad;
	        long K;

		for ( K = Nx-1L;  K > -1L; K-- )		{
		        *Rc++ +=  lamr2 * (*Rtc++) + ((*Xc < 0.0) ? lampos * (*Xc) : 0.0);
			++Xc;
		}
	}
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(J) shared(D, Dt, X, lamr2, lampos, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
	for ( J = 0; J < Ny * Nz; ++J ) {
	        const double* Dtc = Dt + J * (ptrdiff_t) NxPad;
		const double* Xc = X + J * (ptrdiff_t) NxPad;
		double* Dc = D + J * (ptrdiff_t) NxPad;
	        long K;

		for ( K = Nx-1L;  K > -1L; K-- )		{
		    *Dc++ +=  lamr2 * (*Dtc++) + ((*Xc++ < 0.0) ? lampos : 0.0);
		}
	}
}
		
		
	
/*
 */



void Reg (
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
)
{


	
	long Nx, NxPad, Ny, Nz;
	long  I;
	
	Nx = G->Nx;
	NxPad = G->Nx + G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	Regr(  X,  R, G, dir, NF, lamc,  swpspace );	

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(I) shared(D, R, alpha, laml2, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
	for ( I = 0; I < Ny * Nz; ++I ) {
	        const double* Rc = R + I * (ptrdiff_t) NxPad;
		double* Dc = D + I * (ptrdiff_t) NxPad;
	        long J;

		for ( J = Nx-1L; J > -1L; J-- )	{
			*Dc++ = 1.0 / ( 1 + alpha * (*Rc++) ) + laml2;
		}
	}
	
	RegNL ( X, D, R, G, dir, NF, lamc,  swpspace );	
	
	
	Regs(  D,  swpspace+(ptrdiff_t)2L*(G->LengthF), G, dir, NF, lamc,  swpspace );
	memcpy(D,  swpspace+(ptrdiff_t)2L*(G->LengthF),  G->LengthF*sizeof(double));
	
	
	
}


	


/*
 */


void RegNL (
    const double *X,
    const double *D,
    double *R,
    const Geometry *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
)
{
	
	double   *Rs, lam;
	long Nx, NxPad, Ny, Nz;
	long J, K;
	
	Nx = G->Nx;
	NxPad = G->Nx + G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	
	Rs = swpspace;
	
	memset( R, 0, NxPad*Ny*Nz*sizeof(double) );
    
	
    for ( J = 0L; J < NF; J++ )		{
        
		lam = lamc[J];
		
		if ( lam > 0.0 )   {
			
				
			
			Derivative( X, Rs, G, dir+(ptrdiff_t)(2*J), 'f', swpspace+(ptrdiff_t)(G->LengthF));

#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(Rs, D, lam, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
			for ( K = 0; K < Ny * Nz; ++K ) {
			        const double* Dc = D + K * (ptrdiff_t) NxPad;
				double* Rsc = Rs + K * (ptrdiff_t) NxPad;
			        long I;

				for ( I = Nx - 1L; I > -1L; I-- ) {
					*Rsc = lam * (*Rsc)* (*Dc++);  Rsc++;
				}
			}
			
			Derivative( Rs, R, G, dir+(ptrdiff_t)(2*J), 'r', swpspace+(ptrdiff_t)(G->LengthF));
			
			
				
			
			
		}
		
		
	}
	
	
}





/*
 */


    
	

void Regr (
    const double *X,
    double *R,
    const Geometry*G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
)
{
        int Rset = 0;
	long Nx, NxPad, Ny, Nz, LengthF;
	long J, K;
	char dirc[2L];
	
	Nx = G->Nx;
	NxPad = G->Nx + G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	LengthF = G->LengthF;
	
	for ( J = 0L; J < NF; J++ ) {
		double lam = lamc[J];
		
		if ( lam > 0.0 )  {
			dirc[0L] = *(dir+(ptrdiff_t)(J*2));
			dirc[1L] = *(dir+(ptrdiff_t)(J*2L+1L));
			Derivative( X, swpspace, G, dirc, 'f', swpspace+(ptrdiff_t)(G->LengthF));
			
			if ( (dirc[0L] == ' ') && (dirc[1L] == ' ')  )   {
			    if (Rset) {
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(swpspace, R, Nx, NxPad, Ny, Nz, lam) schedule(dynamic,1)
#endif
			        for ( K = 0; K < Ny * Nz; ++K )	{
				        const double* Rs = swpspace + K * (ptrdiff_t) NxPad;
					double* Rc = R + K * (ptrdiff_t) NxPad;
					long I;

					for ( I = Nx - 1L; I > -1L; I-- )		{
					    double tmp = (*Rs > 0.0) ? *Rs : 0.0;

					    ++Rs;
					    *Rc++ += lam * tmp * tmp;
					}
				}
			    } else {
				Rset = 1;
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(swpspace, R, Nx, NxPad, Ny, Nz, lam) schedule(dynamic,1)
#endif
			        for ( K = 0; K < Ny * Nz; ++K )	{
				        const double* Rs = swpspace + K * (ptrdiff_t) NxPad;
					double* Rc = R + K * (ptrdiff_t) NxPad;
					long I;

					for ( I = Nx - 1L; I > -1L; I-- )		{
					    double tmp = (*Rs > 0.0) ? *Rs : 0.0;

					    ++Rs;
					    *Rc++ = lam * tmp * tmp;
					}
				}
			    }
			} else {
			    if (Rset) {
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(swpspace, R, Nx, NxPad, Ny, Nz, lam) schedule(dynamic,1)
#endif
			        for ( K = 0; K < Ny * Nz; ++K )	{
				        const double* Rs = swpspace + K * (ptrdiff_t) NxPad;
					double* Rc = R + K * (ptrdiff_t) NxPad;
					long I;

					for ( I = Nx - 1L; I > -1L; I-- )		{
					    *Rc++ += lam * (*Rs) * (*Rs);
					    ++Rs;
					}
				}
			    } else {
				Rset = 1;
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(swpspace, R, Nx, NxPad, Ny, Nz, lam) schedule(dynamic,1)
#endif
			        for ( K = 0; K < Ny * Nz; ++K )	{
				        const double* Rs = swpspace + K * (ptrdiff_t) NxPad;
					double* Rc = R + K * (ptrdiff_t) NxPad;
					long I;

					for ( I = Nx - 1L; I > -1L; I-- )		{
					    *Rc++ = lam * (*Rs) * (*Rs);
					    ++Rs;
					}
				}
			    }
			}
		}
	}

	if (! Rset) {
	    (void) memset( R, 0, LengthF*sizeof(double) );
	}
}



/*
*/




void Regs (
    const double *X,
    double *R,
    const Geometry *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
)
{
        int Rset = 0;
	double   *Rs;
	long Nx, Ny, NxPad, Nz, LengthF;
	long J, K;
	
	Nx = G->Nx;
	NxPad = G->Nx + G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	LengthF = G->LengthF;
	
	Rs = swpspace;
   
	
        for ( J = 0L; J < NF; J++ )		{
	        double lam = lamc[J];
		
		
		if ( lam > 0.0 )  {
			
			Derivative( X, Rs, G, dir+(ptrdiff_t)(J*2), 's', swpspace+(ptrdiff_t)(G->LengthF));

			if (Rset) {
#if defined(_OPENMP) && defined(USE_OPENMP)
#pragma omp parallel for private(K) shared(R, Rs, lam, Nx, NxPad, Ny, Nz) schedule(dynamic,1)
#endif
			    for ( K = 0; K < Ny * Nz; ++K ) {
				const double* Rsc = Rs + K * (ptrdiff_t) NxPad;
				double* Rc = R + K * (ptrdiff_t) NxPad;
				long I;

				for ( I = Nx - 1L; I > -1L; I-- ) {
					*Rc++ += lam * (*Rsc++);
				}
			    }
			} else {
			    Rset = 1;
			    for ( K = 0; K < Ny * Nz; ++K ) {
				const double* Rsc = Rs + K * (ptrdiff_t) NxPad;
				double* Rc = R + K * (ptrdiff_t) NxPad;
				long I;

				for ( I = Nx - 1L; I > -1L; I-- ) {
					*Rc++ = lam * (*Rsc++);
				}
			    }
			}
		}
	}

	if (! Rset) {
	    (void) memset( R, 0, LengthF*sizeof(double) );
	}
}




