/*
 *  applyreg.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include     <stddef.h>
#include     <stdlib.h>
#include     <string.h>
#include     <math.h>

#include     "fft_double.h"

#include     "volumestruct.h"
#include     "applyreg.h"
#include     "jacobi.h"


void            Jacobi(   double    *Xini,
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
			                      void*      res_reporter_closure )
{
	
	
	FFTDoublePlan    pf_X, pr_X, pr_R,  pf_R, pr_Dr;
	FFTDoublePlan    pf_F, pf_Fsq, pf_Rtmp, pf_W, pr_Rtmp;
	double       *R, *R2, *Dr, *F, *Fsq, *Rtmp, *W, *X;
	long         ItInd, ItIndt, I, J;
	long         Nx, NxSkip, Ny, Nz;
	double       nrmb = 0.0f, res, resp;
	double       arrsz;
	double       omegan;
	
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	omegan = *omega;
	
	arrsz = 1.0f/((double)(Nx*Ny*Nz));
	
	Dr     = swpspace;
	R      = swpspace + (ptrdiff_t)(G->LengthF);
	R2     = swpspace + (ptrdiff_t)2L*(G->LengthF);
	
	Rtmp   = swpspace + (ptrdiff_t)3L*(G->LengthF);
	F      = swpspace + (ptrdiff_t)4L*(G->LengthF);
    Fsq    = swpspace + (ptrdiff_t)5L*(G->LengthF);
	W      = swpspace + (ptrdiff_t)6L*(G->LengthF);
	
	X     = swpspace + (ptrdiff_t)7L*(G->LengthF);
	
	
	
	pf_X = plan_real_complex_3d_fft_double(Nx, Ny, Nz, X, X);
	pr_X = plan_complex_real_3d_fft_double(Nx, Ny, Nz, X, X);
	pf_R = plan_real_complex_3d_fft_double(Nx, Ny, Nz, R, R);
	pr_R = plan_complex_real_3d_fft_double(Nx, Ny, Nz, R, R);

	pf_Rtmp = plan_real_complex_3d_fft_double(Nx, Ny, Nz, Rtmp, Rtmp);
	pr_Rtmp = plan_complex_real_3d_fft_double(Nx, Ny, Nz, Rtmp, Rtmp);
	pf_F = plan_real_complex_3d_fft_double(Nx, Ny, Nz, F, F);
	pf_Fsq = plan_real_complex_3d_fft_double(Nx, Ny, Nz, Fsq, Fsq);
	pf_W = plan_real_complex_3d_fft_double(Nx, Ny, Nz, W, W);
	pr_Dr = plan_complex_real_3d_fft_double(Nx, Ny, Nz, Dr, Dr);
	
	memcpy(  X, Xini, (G->LengthF)*sizeof(double));
	resp = 1.0e16;
	
	
	for ( J = Ny*Nz-1L; J > -1L;  J-- )	{
		for ( I = Nx - 1L;  I > -1L; I-- )		{
			nrmb += (*B)*(*B);
			B++;
		}
		B += NxSkip;
	}
	B   -= (ptrdiff_t)(G->LengthF);
	
	ItIndt = 0;
	ItInd = 0;
	
	while (1)  {
		
		
	        compute_real_complex_3d_fft_double(pf_X, X, X);
		
		memcpy(R, X, (G->LengthF)*sizeof(double));
		ApplyComplexTFp(R, Hd, G->LengthF);
		compute_complex_real_3d_fft_double(pr_R, R, R);
		for (J = G->LengthF-1L;  J > -1L;  J--)		{
			*R++ *= arrsz;
		}
		R -= (ptrdiff_t)G->LengthF;
		
		for ( J = (Nx+NxSkip)*Ny*Nz-1L; J > -1L;  J-- )  {
			*R++ *= *Wd++;
		}
		R   -= (ptrdiff_t)(G->LengthF);
		Wd  -= (ptrdiff_t)(G->LengthF);
		compute_real_complex_3d_fft_double(pf_R, R, R);
		ApplyComplexTFpr(R, Hd, G->LengthF);
		
		
		memset(R2, 0, G->LengthF*sizeof(double));
		memset(Dr, 0, G->LengthF*sizeof(double));
		
		ApplyReg(  X,  R2, Dr, PCF, Wm,  alpham, lam, alpha, lamf, G, F, Fsq, Rtmp, W,
				 pf_F, pf_Fsq, pf_Rtmp, pf_W, pr_Rtmp );
		
		
		compute_complex_real_3d_fft_double(pr_X, X, X);
		compute_complex_real_3d_fft_double(pr_Dr, Dr, Dr);
		for (J = G->LengthF-1L;  J > -1L;  J--)		{
			*X++  *= arrsz;
			*Dr++ *= arrsz;
		}
		X  -= (ptrdiff_t)G->LengthF;
		Dr -= (ptrdiff_t)G->LengthF; 
		
		
		for ( J = G->LengthF-1L; J > -1L;  J-- )		{
			*Dr++ +=  *D++;
			*R++  +=  *R2++;
		}
		D  -= (ptrdiff_t)G->LengthF;
		Dr -= (ptrdiff_t)G->LengthF;
		R  -= (ptrdiff_t)G->LengthF;
		R2 -= (ptrdiff_t)G->LengthF; 
		
		compute_complex_real_3d_fft_double(pr_R, R, R);
		for (J = G->LengthF-1L;  J > -1L;  J--)		{
			*R++ *= arrsz;
		}
		R -= (ptrdiff_t)G->LengthF;		
		
		
		
		res = 0.0f;
		
		for ( J = Ny*Nz-1L; J > -1L; J-- )		{
			for ( I = Nx-1L; I > -1L;  I-- )		{
				*R = *B - *R;
				res +=  (*R)*(*R);
				R++;B++;
			}
			R += NxSkip;
			B += NxSkip;
		}
		R  -= (ptrdiff_t)G->LengthF;
		B  -= (ptrdiff_t)G->LengthF;
		
		res /= nrmb;
		res = sqrt(res);

		if (res_reporter != 0) {
		    (*res_reporter)(
			ItInd + ((ItIndt != 0 && res - resp <= 1e-8) ? 1 : 0),
			ItIndt,
			res,
			res_reporter_closure
		    );
		}

		if ( res <= tol )   {
			ItIndt++;
			ItInd++;
			break;
		}
		else if  ( (res - resp) > 1e-8  )		{
			memcpy(X, Xini, G->LengthF*sizeof(double));
			omegan *= 0.7;
		}
		else {
		    
			if ( (res - resp) <= 1e-6 )  {
			    ItInd++;
		    }
			if ( (ItIndt >= NIterT) || (ItInd >= NIter) )  {
				break;
			}
			
			memcpy(Xini, X, G->LengthF*sizeof(double));
			resp = res;
			/*
			for ( J = Ny*Nz-1L; J > -1L; J-- )		{
				for ( I = Nx-1L; I > -1L;  I-- )		{
					*X += omegan*((*R)/(*Dr) );
					X++;R++;Dr++;
				}
				X += NxSkip;
				R += NxSkip;
				Dr += NxSkip;
			}
			X   -= (ptrdiff_t)G->LengthF;
			R   -= (ptrdiff_t)G->LengthF;
			Dr  -= (ptrdiff_t)G->LengthF;  */
			
			Update( X, R, PCF2, Dr, pf_R, pr_R, omegan, Nx, NxSkip, Ny, Nz ); 
			
						
			
			
		}  
		
		ItIndt++;
		
	}
	
	*omega = omegan;
	
	memcpy(Xini, X, G->LengthF*sizeof(double));
	
	destroy_fft_double_plan(pf_X);
	destroy_fft_double_plan(pr_X);
	destroy_fft_double_plan(pf_R);
	destroy_fft_double_plan(pr_R);
	destroy_fft_double_plan(pf_Rtmp);
	destroy_fft_double_plan(pr_Rtmp);
	destroy_fft_double_plan(pf_F);
	destroy_fft_double_plan(pf_Fsq);
	destroy_fft_double_plan(pf_W);	
	
	
	
}





/*
*/





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
					long         Nz        )


{
	
	
	long I, J;
	double   arrsz;
	
	arrsz = 1.0/(Nx*Ny*Nz);
	
	compute_real_complex_3d_fft_double(pfR, R, R);
	for ( I = (Nx+NxSkip)*Ny*Nz/2L-1L;  I > -1L;  I-- )	{
	    *R++  *= *PCF * arrsz; *R++  *= (*PCF++) * arrsz;PCF++;	
    }
	R   -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	PCF -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	compute_complex_real_3d_fft_double(prR, R, R);
	
	
	for ( I = Ny*Nz-1L;  I > -1L;  I-- )	{
		for ( J = Nx-1L;  J > -1L;  J-- )      {
			*R++  /=  *Dr++;
		}
		R += NxSkip;
		Dr += NxSkip;
	}
	R    -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	Dr   -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	
	compute_real_complex_3d_fft_double(pfR, R, R);
	for ( I = (Nx+NxSkip)*Ny*Nz/2L-1L;  I > -1L;  I-- )	{
	    *R++  *= *PCF * arrsz;*R++  *= (*PCF++) * arrsz ;PCF++;
    }
	R   -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	PCF -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	compute_complex_real_3d_fft_double(prR, R, R);
	
	
	for ( I = Ny*Nz-1L;  I > -1L;  I-- )	{
		for ( J = Nx-1L;  J > -1L;  J-- )      {
			*X++  +=  omega*(*R++);
		}
		X += NxSkip;
		R += NxSkip;
	}
	
	
	
}
	
	
	
			
		     	
		
		
		
		
		
		
	
	
	
				    
				
				
				
				
				
				
				
				
				
	
