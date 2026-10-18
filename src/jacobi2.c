/*
 *  jacobi2.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 6/9/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include        <stdlib.h>
#include        <string.h>
#include        <stddef.h>
#include        <math.h>

#include        "fft_double.h"
#include		"volumestruct.h"
#include        "derivatives.h"
#include		"regc.h"
#include        "misc.h"
#include		"jacobi2.h"





void JacobiNL3D (
    double *X,
    const double *TF,
    double Diagd,
    const double *B,
    const double *PCF,
    double lamf,
    const double *lamratio,
    double lampos,
    double delfrac,
    const double *beta,
    const double *alpha,
    const Geometry *G,
    double *swpspace,
    double *omega,
    long NIter,
    double tol,
    void (*res_reporter)(
	long iiters, long ittert, double res, void* res_reporter_closure
    ),
    void *res_reporter_closure,
    double laml2
)
{
	
	
	
	double         resx, res, nrmb = 0.0;
	long            J, K,  Nx, NxSkip, Ny, Nz, Itind,  NIteracc = 0L;
	double         *R, *D, *U, arrsz;
	FFTDoublePlan  pfR, prR, pfU,  prU;
	
	
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	arrsz = 1.0/(Nx*Ny*Nz);
	
	R = swpspace;
	D = swpspace + (ptrdiff_t)(G->LengthF);
	U = swpspace + (ptrdiff_t)2L*(G->LengthF);
	
	
	pfR = plan_real_complex_3d_fft_double(Nx, Ny, Nz, R, R);
	prR = plan_complex_real_3d_fft_double(Nx, Ny, Nz, R, R);
	pfU = plan_real_complex_3d_fft_double(Nx, Ny, Nz, U, U);
	prU = plan_complex_real_3d_fft_double(Nx, Ny, Nz, U, U);
	
	
	memset(R, 0, sizeof(double)*(G->LengthF));
	memset(D, 0, sizeof(double)*(G->LengthF));
	memset(U, 0, sizeof(double)*(G->LengthF));

	
	
    for ( J = Nz*Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			nrmb +=  (*B)*(*B);
			B++;
		}
		B += NxSkip;
	}
	B    -= (ptrdiff_t)G->LengthF;	
	
	
	Res3D_old(  X, TF, Diagd, R, D, lamf, lamratio,
		    lampos,  delfrac, beta, alpha, G, swpspace  + (ptrdiff_t)3L*(G->LengthF),
			pfR, prR, laml2   );
	
	/*
	for ( J = 0L;  J < 10L; J++ )	{
		printf("%3.10f   ", R[J]);
	}
	printf("\n");  */

	
	
	
	res = 0.0;
	for ( J = Nz*Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			resx  = *B++ - *R;
			*R++  = resx;
			res  += resx*resx;				
		}
		R += NxSkip;
		B += NxSkip;
	}
	R    -= (ptrdiff_t)G->LengthF;
	B    -= (ptrdiff_t)G->LengthF;
	res  /=  nrmb;
	res = sqrt(res);
	if (res_reporter != 0) {
	    (*res_reporter)(0, 0, res, res_reporter_closure);
	}
	
	
	
	for ( Itind = 0L;  Itind < NIter; Itind++  )	{
	       long NIteracco = NIteracc;

		/*
		 printf("%d\n", (int)Itind);
		 printf("%3.10f\n", res);

		 for ( J = 0L;  J < 10L; J++ )	{
		 printf("%3.10f   ", R[J]);
		 }
		 printf("\n");  */
		
		if ( res <= tol  )  {
			break;
		}
		
		memcpy(U, R,  (G->LengthF)*sizeof(double));
		compute_real_complex_3d_fft_double(pfU, U, U);
		ApplyComplexTFp (  U, PCF,  G->LengthF  );
		compute_complex_real_3d_fft_double(prU, U, U);
		for ( J = Nz*Ny-1L; J > -1L; J-- )		{
			for ( K = Nx-1L;  K > -1L; K-- )		{
				*U++ *=  arrsz*arrsz/(*D++);
			}
			U += NxSkip;
			D += NxSkip;
		}
		U   -= (ptrdiff_t)G->LengthF;
		D   -= (ptrdiff_t)G->LengthF;
		compute_real_complex_3d_fft_double(pfU, U, U);
		ApplyComplexTFpr (  U, PCF,  G->LengthF  );
		compute_complex_real_3d_fft_double(prU, U, U);

		/*
		for ( J = 0L;  J < 10L; J++ )	{
			printf("%3.10f   ", U[J]);
		}
		printf("\n\n\n");  */
		

		
		
		JacobiNL3Dstep  (  X, U, TF, Diagd, R, B, D, lamf,
				   lamratio, lampos, delfrac, beta,
				   alpha, G, swpspace +  (ptrdiff_t)3L*(G->LengthF),
				   pfR, prR, omega, Itind, &NIteracc, &res, res_reporter, res_reporter_closure, laml2 );
		
		/*
		for ( J = 0L;  J < 10L; J++ )	{
			printf("%3.10f   ", R[J]);
		}
		printf("\n\n\n");	*/	
		

		if ( NIteracc - NIteracco  >  50  )	{
			break;
		}
		
		
		
	}
	
	destroy_fft_double_plan(prU);
	destroy_fft_double_plan(pfU);
	destroy_fft_double_plan(prR);
	destroy_fft_double_plan(pfR);
}




/*
*/
	
	
	

void JacobiNL3Dstep (
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
    double delfrac,
    const double *beta,
    const double *alpha,
    const Geometry *G, 
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
)
{
  
	
        double  res2, nrmb = 0.0, resx, omegan;
	double  *Xt;
	long     Nx, NxSkip, Ny, Nz,   J, K;
	int     count = 0;
	
	
	omegan = *omega;
	

	
	Xt  = swpspace;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	for ( J = Nz*Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			nrmb +=  (*B)*(*B);
			B++;
		}
		B += NxSkip;
	}
	B    -= (ptrdiff_t)G->LengthF;	
		
		
	
	while ( 1 )  {
		memcpy(Xt, X, sizeof(double)*(G->LengthF));
		
		
		for ( J = Nz*Ny-1L; J > -1L; J-- )		{
			for ( K = Nx-1L;  K > -1L; K-- )		{
				*Xt++  +=  omegan*(*U++);
			}
			Xt += NxSkip;
			U += NxSkip;
		}
		Xt   -= (ptrdiff_t)G->LengthF;
		U    -= (ptrdiff_t)G->LengthF;  
		
		
		/*
		for ( J = 0L;  J < 10L;  J++ )	{
			printf("%2.10f    ",Xt[J]);
		}
		printf("\n");  */
		
			
		
		Res3D_old(  Xt, TF, Diagd, R, D, lamf, lamratio,
				lampos,  delfrac, beta, alpha, G, swpspace+(ptrdiff_t)G->LengthF,
				pfR, prR, laml2   );
		
		
		/*
		for ( J = 0L;  J < 10L;  J++ )	{
			printf("%2.10f    ",R[J]);
		}
		printf("\n"); */
		
 		
		
		res2 = 0.0;
		for ( J = Nz*Ny-1L; J > -1L; J-- )		{
			for ( K = Nx-1L;  K > -1L; K-- )		{
				resx  = *B++ - *R;
				*R++  = resx;
				res2 += resx*resx;				
			}
			R += NxSkip;
			B += NxSkip;
		}
		R    -= (ptrdiff_t)G->LengthF;
		B    -= (ptrdiff_t)G->LengthF;
		res2  /=  nrmb;
		res2   = sqrt( res2 );
		
		
		/* printf("%f   %f\n\n\n\n\n", *p_res, res2); */
		
		/*
		 for ( I = 0L; I < 5L;  I++ )	{
		 printf("%3.10f   ", R[I] );
		 }
		 printf("\n\n"); */
		
		++count;
		if (res_reporter != 0) {
		    (*res_reporter)(
			iiters + ((res2 < *p_res) ? 1 : 0),
			*p_iitert + count,
			res2,
			res_reporter_closure
		    );
		}
		
		if ( res2 >= *p_res )		{
			omegan *= 0.7;
			if (count  > 50 )   {
				break;
			}
		}
		else  {
			break;
		}  
					
		
				
		
		
	}
	
	memcpy(X, Xt, sizeof(double)*(G->LengthF));
	*p_iitert += count;
	*p_res  = res2;
	*omega  = omegan;
}




/*
*/



void Res3D_old (
    const double *X,
    const double *TF,
    double Diagd,
    double *R,
    double *D,
    double lamf, 
    const double *lamratio,
    double lampos,
    double delfrac,
    const double *beta,
    const double *alpha,
    const Geometry *G,
    double *swpspace,
    FFTDoublePlan pfR,
    FFTDoublePlan prR,
    double laml2
)
{

	long Nx, NxSkip, Ny, Nz;
	double  arrsz;
	long     J, K;
		
	double  *Rt, *Dt;
	
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	arrsz = 1.0/(Nx*Ny*Nz);
	
 
	
	
	memcpy(R, X,  (G->LengthF)*sizeof(double));
	compute_real_complex_3d_fft_double(pfR, R, R);
	for ( J = G->LengthF/2L - 1; J > -1L; J-- ) {
	    double scl = ((*TF) * (*TF) + (*(TF + 1)) * (*(TF + 1))) * arrsz;

	    *R++ *= scl;
	    *R++ *= scl;
	    TF += 2;
	}
	TF -= G->LengthF;
	R -= G->LengthF;
        compute_complex_real_3d_fft_double(prR, R, R);
	
    /*
	for ( J = 0L;  J < 10L;  J++ )	{
		printf("%2.10f    ",R[J]);
	}
	printf("\n"); */

	
	Rt = swpspace;
	Dt = swpspace + (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	
	Regsum ( X, Rt, Dt,  lamratio, laml2, lampos, delfrac, beta, alpha, G, swpspace + (ptrdiff_t)2L*(Nx+NxSkip)*Ny*Nz );
	
    /*
	for ( J = 0L;  J < 10L;  J++ )	{
		printf("%2.10f    ",Rt[J]);
	}
	printf("\n"); */
	
    /*
	void  Regsum	(	 double     *X,
				     double     *R,
				     double     *D,
				     double     *lamratio,
					 double     laml2,
				     double     lampos,
				     double     delfrac,
				     double     *beta,
				     double     *alpha,
				     Geometry	*G, 
				     double     *swpspace
					 )  */
	
	
	
	for ( J = Nz*Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			*R++ +=  lamf*(*Rt++);
			*D++ =  lamf*(*Dt++) + Diagd;
		}
		R += NxSkip;
		D += NxSkip;
		Rt += NxSkip;
		Dt += NxSkip;
	}
	R     -= (ptrdiff_t)G->LengthF;
	D     -= (ptrdiff_t)G->LengthF;
	Rt    -= (ptrdiff_t)G->LengthF;
	Dt    -= (ptrdiff_t)G->LengthF;
}
