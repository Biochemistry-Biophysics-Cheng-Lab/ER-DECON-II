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
#include     <math.h>
#include     <string.h>

#include     "fft_double.h"
#include     "volumestruct.h"
#include     "derivatives.h"
#include     "applyreg.h"





void		ApplyReg(  double       *X,
					   double       *R,
					   double       *D,
					   double       *PCF,
					   double       *Wm,
					   double       *alpham,
					   double       *lam,
					   double       *alpha,
					   double       lamf,
					   Geometry     *G,
					   double       *F,
					   double       *Fsq,  
					   double       *Rtmp,
					   double       *W,
					   FFTDoublePlan    pf_F,
					   FFTDoublePlan    pf_Fsq,  
					   FFTDoublePlan    pf_Rtmp,
					   FFTDoublePlan    pf_W,
					   FFTDoublePlan    pr_Rtmp  )
{
	
	char  dir[16] = "    xxyyzzxyyzxz";
	int   type[8L] = {0, 1, 2, 2, 2, 2, 2, 2};
	long   Nx, NxSkip, Ny, Nz;
	double tmpv, arrsz;
	
	long  I, J, K;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	arrsz = 1.0f/((double)(Nx*Ny*Nz));
	
	for ( I = 0L; I < 8L; I++ )			{
		
		if ( lam[I] != 0.0f )		{
		
			Derivative(PCF,  F, G, dir + (ptrdiff_t)(I*2L), 'f', Fsq );
			
			for ( J = Nz*Ny-1L;  J > -1L;  J--  )	{
				for ( K = Nx-1L;  K > -1L;  K--  )     {
					*Fsq++ = (*F)*(*F);
					F++;
				}
				Fsq += NxSkip;
				F += NxSkip;
			}
			F   -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			Fsq -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;

			compute_real_complex_3d_fft_double(pf_Fsq, Fsq, Fsq);
			compute_real_complex_3d_fft_double(pf_F, F, F);
			
			
			memcpy(Rtmp, X,  (Nx+NxSkip)*Ny*Nz*sizeof(double));
			ApplyComplexTFp(  Rtmp, F,  (Nx+NxSkip)*Ny*Nz  );
			compute_complex_real_3d_fft_double(pr_Rtmp, Rtmp, Rtmp);
			
			for ( J = (Nx+NxSkip)*Ny*Nz-1L;  J > -1L;  J-- ) {
				*Rtmp++ *= arrsz;
			}
			Rtmp -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			
			
			if ( type[I] == 0  )		{
			
				for ( J = Nz*Ny-1L;  J > -1L;  J--  )	{
					for ( K = Nx-1L;  K > -1L;  K--  )     {
						*W        =  (*Rtmp < 0.0f) ?  1.0f :  0.0f;
						*Rtmp++  *= *W++;
					}
					W += NxSkip;
					Rtmp += NxSkip;
				}
				W    -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				Rtmp -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				
			}
			else if ( type[I] == 1  )		{
				for ( J = Nz*Ny-1L;  J > -1L;  J--  )	{
					for ( K = Nx-1L;  K > -1L;  K--  )     {
						tmpv = *Rtmp;
						*W        =  1.0f/(1.0f+alpha[I]*(tmpv*tmpv));
						*Rtmp++  *= *W++;
					}
					W += NxSkip;
					Rtmp += NxSkip;
				}
				W    -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				Rtmp -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			}
			else		{
				for ( J = Nz*Ny-1L;  J > -1L;  J--  )	{
					for ( K = Nx-1L;  K > -1L;  K--  )     {
						tmpv = *Rtmp;
						*W        =  (*Wm++)/(1.0f+alpha[I]*(*alpham++)*(tmpv*tmpv));
						*Rtmp++  *= *W++;
					}
					W += NxSkip;
					Rtmp += NxSkip;
					Wm += NxSkip;
					alpham += NxSkip;
				}
				W       -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				Rtmp    -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				Wm      -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				alpham  -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
				
			}
			
			compute_real_complex_3d_fft_double(pf_W, W, W);
			compute_real_complex_3d_fft_double(pf_Rtmp, Rtmp, Rtmp);
			ApplyComplexTFpr(  Rtmp, F,  (Nx+NxSkip)*Ny*Nz  );
			ApplyComplexTFpr(  W, Fsq,  (Nx+NxSkip)*Ny*Nz  );
			
			for ( J = (Nx+NxSkip)*Ny*Nz-1L;  J > -1L;  J-- ) {
				*R++ +=  lam[I]*(*Rtmp++);
				*D++ +=  lam[I]*(*W++);
			}
			Rtmp   -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			R      -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			D      -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			W      -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
			
			
		}
		
		
	}
		
    
	for ( J = (Nx+NxSkip)*Ny*Nz-1L;  J > -1L;  J-- )	{
		*R++ *=  lamf;
		*D++ *=  lamf;
	}
	R      -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	D      -=  (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;		
    
	
}
		


/*
 */
			
			
			
void      ApplyComplexTFp (    double   *X,
							   double   *TF,
							   long     Length  )
	
	
{
	
	long J;
	double re, im;
	
	for ( J = Length/2L - 1L; J > -1L; J-- )		{
		re    =  *X++;               
		im    =  *X--;                
		*X       =  re*(*TF++);
		*X++    -=  im*(*TF);
		*X       =  re*(*TF--);
		*X++    +=  im*(*TF++);
		TF++;
	}
	
	
}
	
	
/*
 */
	
	
	
void      ApplyComplexTFpr       (  double   *X,
									double   *TF,
									long     Length  )
	
	
{
	
	long   J;
	double re, im;
	
	for ( J = Length/2L - 1L; J > -1L; J-- )		{
		re    =  *X++;               
		im    =  *X--;                
		*X       =  re*(*TF++);
		*X++    +=  im*(*TF);
		*X       =  -re*(*TF--);
		*X++    +=  im*(*TF++);
		TF++;
	}
	
	
	
}
			
			
			


