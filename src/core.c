/*
 *  applyreg.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#if defined(_OPENMP) && defined(USE_OPENMP)
#include     <omp.h>         /* omp_set_num_threads */
#endif
#include     <stddef.h>
#include     <stdlib.h>
#include     <math.h>
#include     <string.h>

#include     "fft_double.h"
#include     "volumestruct.h"
#include     "applyreg.h"
#include     "misc.h"
#include     "jacobi.h"

#include     "core.h"



int         core       (  double   *Imn,
						  double   *TF,
						  double   *TFc,
						  double   *X,
						  long     *Sz,
						  double   lamratio,
						  double   delfac,
						  double   lamf,
						  double   lampos,
						  double   lampc,
						  double   lamsmooth,
						  double   alpha,
						  long     NIter,
						  long     NIterT,
						  double   tol,
						  double   omega,
	   		                          int      nthread,
			                          void     (*res_reporter)(long iiters, long iitert, double res, void*),
			                          void*    res_reporter_closure )

/*
 
 
 Imn:   input data
 TF:     OTF
 TFc:    missing cone filter in Fourier that is derived from TF.  Function getmcf.c does the job.  Not used if lampc is less than or equal to zero.
 X:      output
 Sz:     [Nx Ny Nz];
 lamratio:   optional parameter. default value = 0.5;
 lamf:   smoothing parameter.  User has to choose it.   Will be calibrated soon.
 lampos:     optional parameter enforces postivity. Default: 100.0
 lampc:      weight for the missing cone filter. Default:  100.0
 lamsmooth:  aux. smoohting paramter. Default:  100.0
 alpha:  nonlinearity. Default: 1000
 NIter:  Number of iterations. Default: 50.
 tol.:   Termination tolerance.  Default:  0.001.
 omega:  relaxation parameter.  Default: 0.8.
 
 
 */

{
	
	Geometry     *G;
	
	long        *mempool;
	double      *mempoold;	
	double      *swpspace;
	
	double      *B, *PCF,  *PCF2, *Wm, *alpham, *Wd,  *D, *TFsq;
	
	char     dir[16L]  = "    xxyyzzxyyzxz";
	
	double   lam[8L];
	
	long	 Nx, Ny, Nz;
	long     I, J;
	
	double   arrsz;
	double   tmpval, maxval, factor = 1.0;
	double   re1, re2, re3, im1, tmpvalre, tmpvalim;
	double   red, imd;
	
	double   alphas[8L], Rmax[8L] = {1,1,1,1,1,1,1,1};
	
	
	FFTDoublePlan  pf_PCF, pr_PCF, pf_Wd, pr_Wd,  pf_Imn, pr_Imn,  pf_X, pr_X,  pf_B, pr_B, pf_TFsq,  pr_TFsq, pf_D, pr_D;
	
	if (initialize_fft_double() == 0) {
	    return -1;
	}
	set_fft_double_thread_count(nthread);
#if defined(_OPENMP) && defined(USE_OPENMP)
	omp_set_num_threads(nthread);
#endif

	Nx = Sz[0L];
	Ny = Sz[1L];
	Nz = Sz[2L];
	arrsz = 1.0/((double)(Nx*Ny*Nz));
	
	lam[0L] = lampos;
	lam[1L] = lamratio;
	for ( I = 2L; I < 8L;  I++ )	{
	     lam[I] = (1.0f - lamratio);
	}
	lam[4L] *= delfac*delfac;
	lam[6L] *= delfac;
	lam[7L] *= delfac;
	
	
	
	G = (Geometry  *)malloc(sizeof(Geometry));
	if (G == (Geometry*) 0) {
		return -1;
	}
	mempool = (long *)malloc( 7L*(2L*Nx*Ny + (Nz-2L)*(2L*Nx+2L*Ny-4L))*sizeof(long) );
	if ( mempool == (long *)NULL )	{
		free(G);
		return -1;
	}
	ConstructVolumeStruct( G, Nx, Ny, Nz, mempool );
	mempoold = (double *)malloc(15L*(Nx+G->NxFSkip)*Ny*Nz*sizeof(double));
	if ( mempoold == (double *)NULL )	{
		free(mempool);
		free(G);
		return -1;
	}
	B      = mempoold;
	PCF    = mempoold + (ptrdiff_t)1L*(Nx+G->NxFSkip)*Ny*Nz;
	Wm     = mempoold + (ptrdiff_t)2L*(Nx+G->NxFSkip)*Ny*Nz;
	alpham = mempoold + (ptrdiff_t)3L*(Nx+G->NxFSkip)*Ny*Nz;
	Wd     = mempoold + (ptrdiff_t)4L*(Nx+G->NxFSkip)*Ny*Nz;
	D      = mempoold + (ptrdiff_t)5L*(Nx+G->NxFSkip)*Ny*Nz;
	PCF2   = mempoold + (ptrdiff_t)6L*(Nx+G->NxFSkip)*Ny*Nz;
	
	
	swpspace =  mempoold + (ptrdiff_t)7L*(Nx+G->NxFSkip)*Ny*Nz;
	
	TFsq     =  swpspace;
	
	
	
	/* FFTs */
	pf_PCF = plan_real_complex_3d_fft_double(Nx, Ny, Nz, PCF, PCF);
	pr_PCF = plan_complex_real_3d_fft_double(Nx, Ny, Nz, PCF, PCF);
	pf_Wd = plan_real_complex_3d_fft_double(Nx, Ny, Nz, Wd, Wd);
	pr_Wd = plan_complex_real_3d_fft_double(Nx, Ny, Nz, Wd, Wd);
	pf_Imn =
	    plan_no_overwrite_real_complex_3d_fft_double(Nx, Ny, Nz, Imn, Imn);
	pr_Imn =
	    plan_no_overwrite_complex_real_3d_fft_double(Nx, Ny, Nz, Imn, Imn);
	pf_X = plan_real_complex_3d_fft_double(Nx, Ny, Nz, X, X);
	pr_X = plan_complex_real_3d_fft_double(Nx, Ny, Nz, X, X);
	pf_B = plan_real_complex_3d_fft_double(Nx, Ny, Nz, B, B);
	pr_B = plan_complex_real_3d_fft_double(Nx, Ny, Nz, B, B);
	pf_TFsq = plan_real_complex_3d_fft_double(Nx, Ny, Nz, TFsq, TFsq);
	pr_TFsq = plan_complex_real_3d_fft_double(Nx, Ny, Nz, TFsq, TFsq);
	pf_D = plan_real_complex_3d_fft_double(Nx, Ny, Nz, D, D);
	pr_D = plan_complex_real_3d_fft_double(Nx, Ny, Nz, D, D);
	if (pf_PCF == 0 || pr_PCF == 0 || pf_Wd == 0 || pr_Wd == 0 ||
	    pf_Imn == 0 || pr_Imn == 0 || pf_X == 0 || pr_X == 0 ||
	    pf_B == 0 || pr_B == 0 || pf_TFsq == 0 || pr_TFsq == 0 ||
	    pf_D == 0 || pr_D == 0) {
	    destroy_fft_double_plan(pf_PCF);
	    destroy_fft_double_plan(pr_PCF);
	    destroy_fft_double_plan(pf_Wd);
	    destroy_fft_double_plan(pr_Wd);
	    destroy_fft_double_plan(pf_Imn);
	    destroy_fft_double_plan(pr_Imn);
	    destroy_fft_double_plan(pf_X);
	    destroy_fft_double_plan(pr_X);
	    destroy_fft_double_plan(pf_B);
	    destroy_fft_double_plan(pr_B);
	    destroy_fft_double_plan(pf_TFsq);
	    destroy_fft_double_plan(pr_TFsq);
	    destroy_fft_double_plan(pf_D);
	    destroy_fft_double_plan(pr_D);
	    free(mempoold);
	    free(mempool);
	    free(G);
	    return -1;
	}
	
	
	
	
	/* Prepare Poisson weights */
	memset(PCF, 0, (G->LengthF)*sizeof(double));
	AddDuchon3D(PCF, Nx, Ny, Nz, delfac);
	compute_real_complex_3d_fft_double(pf_PCF, PCF, PCF);
	
	memcpy(Wd, Imn,  (G->LengthF)*sizeof(double));
	compute_real_complex_3d_fft_double(pf_Wd, Wd, Wd);
	for ( I = G->LengthF/2L-1L;  I > -1L;  I-- )		{
		*Wd++ /=  (1.0f+lamf*lamsmooth*(1.0f-lamratio)*(*PCF));
		*Wd++ /=  (1.0f+lamf*lamsmooth*(1.0f-lamratio)*(*PCF++));
		PCF++;
	}
	Wd    -= (ptrdiff_t)G->LengthF;
	PCF   -= (ptrdiff_t)G->LengthF;
	compute_complex_real_3d_fft_double(pr_Wd, Wd, Wd);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*Wd++ *=  arrsz;
	}
	Wd    -= (ptrdiff_t)G->LengthF;
	
	


	
	/* Square and normalize Poissons weights */
	maxval = 0.0f;
	for ( I = G->Nz*Ny-1L;  I > -1L;  I-- )	{
		for ( J = (Nx)-1L;  J > -1L;  J-- )	{
			tmpval = *Wd;
		    tmpval *= tmpval;
		    *Wd++ = tmpval;
		    if ( tmpval > maxval  )	{  maxval = tmpval; }
		}
		Wd++;Wd++;
	}
	Wd -= (ptrdiff_t)G->LengthF;
	for ( I = Nz*Ny-1L;  I > -1L;  I-- )	{
		for ( J = (Nx)-1L;  J > -1L;  J-- )	{
		    *Wd   *= 10.0f/maxval;
			*Wd++ += 1.0f;
		}
	    Wd++;Wd++;
	}
	Wd -= (ptrdiff_t)G->LengthF;
	
	
	
	/* Compute Tikhonov solution */
	compute_real_complex_3d_fft_double(pf_Imn, Imn, Imn);
	
	if (lampc > 0.0) {
	    for ( I = G->LengthF/2L - 1L; I > -1L;  I-- )		{
		re1 = *TF++;
		im1 = *TF++;
		red = *Imn++;
		imd = *Imn++;
		re2 = *TFc++;  TFc++;
		re3 = *PCF;
		tmpval = 1.0f/(re1*re1 + im1*im1 + lamf*(lampc*re2*re2 + (1.0f-lamratio)*re3  + lamratio));
		*PCF++ = tmpval;  PCF++;
		tmpvalre = re1*red  + im1*imd;
		tmpvalim = re1*imd  - im1*red;
		*X++ = tmpvalre*tmpval;
		*X++ = tmpvalim*tmpval;
	    }
	    TFc  -= (ptrdiff_t)G->LengthF;
	} else {
	    for ( I = G->LengthF/2L - 1L; I > -1L;  I-- )		{
		re1 = *TF++;
		im1 = *TF++;
		red = *Imn++;
		imd = *Imn++;
		re3 = *PCF;
		tmpval = 1.0f/(re1*re1 + im1*im1 + lamf*((1.0f-lamratio)*re3 + lamratio));
		*PCF++ = tmpval;  PCF++;
		tmpvalre = re1*red  + im1*imd;
		tmpvalim = re1*imd  - im1*red;
		*X++ = tmpvalre*tmpval;
		*X++ = tmpvalim*tmpval;
	    }
	}
	TF   -= (ptrdiff_t)G->LengthF;
	Imn  -= (ptrdiff_t)G->LengthF;
	PCF  -= (ptrdiff_t)G->LengthF;
	X    -= (ptrdiff_t)G->LengthF;
	
	compute_complex_real_3d_fft_double(pr_Imn, Imn, Imn);
	compute_complex_real_3d_fft_double(pr_X, X, X);
	for ( I = G->LengthF-1L;  I > -1L; I-- )	{
		*X++   *= arrsz;
		*Imn++ *= arrsz;
	}
	X   -=  (ptrdiff_t)(G->LengthF);
	Imn -=  (ptrdiff_t)(G->LengthF);

	
	for ( I = G->LengthF-1L;  I > -1L;  I-- )    {
		tmpval = *X;
		*X++ =  (tmpval >= 0) ?  tmpval :  0;
	}
	X -=  (ptrdiff_t)(G->LengthF);	

	
	
	/* compute alpham and Wm */
	for (    I = G->LengthF-1L;  I > -1L;  I-- )		{
		tmpval =  *X++;
		tmpval *= tmpval;
		*alpham++  = tmpval;
		*Wm++ =  (alpha*tmpval > 0) ? tanh(alpha*tmpval)/(alpha*tmpval)  :  1.0f;
	}
	X         -=  (ptrdiff_t)(G->LengthF);
	alpham    -=  (ptrdiff_t)(G->LengthF);
	Wm        -=  (ptrdiff_t)(G->LengthF);

	for (    I = G->LengthF-1L;  I > -1L;  I-- )		{
		*Wm    *=  100.0f;
		*Wm++  +=  1.0f;
	}
	Wm -= (ptrdiff_t)(G->LengthF); 	
	
	
	
	/* Compute alphas */
	RegRMax	(  X,  (double *)Rmax,  G, dir,  8L, swpspace   );
	for ( I = 0L;  I < 8L;  I++ )   {
		Rmax[I] =  ( Rmax[I] >= 1.0f ) ?  Rmax[I] : 1.0f;
	}
	for ( I =0L;  I < 8L;  I++  )  {
		alphas[I] =  alpha/(Rmax[I]);
	}
	
	if ( factor != 1.0 )  alphas[1L] *= 100.0;
	
	
	/*****************************************************************/
	
	
	memcpy(PCF2,PCF,G->LengthF*sizeof(double));
	
	if ( lamf < 0.02 )   {
		
		for ( I = G->LengthF/2L-1L; I > -1L;  I--  )	{
			
			*PCF++ = 1.0;  PCF++;
			
	    }
		PCF -= (ptrdiff_t)G->LengthF;
		
	}
	else {
		
	   
		for ( I = G->LengthF/2L-1L; I > -1L;  I--  )	{
			
			*PCF2++ = 1.0;  PCF2++;
			
	    }
		PCF2 -= (ptrdiff_t)G->LengthF;
		
		
	}
	
	
	/* Normalize */
	compute_real_complex_3d_fft_double(pf_X, X, X);
	
	for ( I = G->LengthF/2L-1L;  I > -1L; I--  )	{
		tmpval  = sqrt(*PCF); 
		*PCF++  = tmpval;   PCF++;
		*TF++  *= tmpval;
		*TF++  *= tmpval;
		*X++   /= tmpval;
		*X++   /= tmpval;
	}
	TF   -= (ptrdiff_t)G->LengthF;
	PCF  -= (ptrdiff_t)G->LengthF;
	X    -= (ptrdiff_t)G->LengthF;
	
	
	for ( I = G->LengthF/2L-1L; I > -1L;  I--  )	{
		
		*PCF2 = sqrt(*PCF2);  PCF2++;  PCF2++;
		
	}
	PCF2 -= (ptrdiff_t)G->LengthF;


	compute_complex_real_3d_fft_double(pr_X, X, X);
	
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*X++ *=  arrsz;
	}
	
	X    -= (ptrdiff_t)G->LengthF;
	
	
	
	/* Compute RHS */
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*B++ = (*Imn++)*(*Wd++);
	}
	B    -= (ptrdiff_t)G->LengthF;
	Imn  -= (ptrdiff_t)G->LengthF;
	Wd   -= (ptrdiff_t)G->LengthF;
	compute_real_complex_3d_fft_double(pf_B, B, B);
    ApplyComplexTFpr(B, TF, G->LengthF);
	compute_complex_real_3d_fft_double(pr_B, B, B);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*B++ *=  arrsz;
	}
	B    -= (ptrdiff_t)G->LengthF;
	
	

	
	
	/* Compute diagonal of data part */
	memcpy(TFsq, TF, sizeof(double)*G->LengthF);
	memcpy(D, Wd,  sizeof(double)*G->LengthF);
	compute_complex_real_3d_fft_double(pr_TFsq, TFsq, TFsq);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*TFsq *=  arrsz*arrsz*(*TFsq);
		++TFsq;
	}
	TFsq    -= (ptrdiff_t)G->LengthF;
	compute_real_complex_3d_fft_double(pf_TFsq, TFsq, TFsq);
	compute_real_complex_3d_fft_double(pf_D, D, D);
	ApplyComplexTFpr(D, TFsq, G->LengthF);
	compute_complex_real_3d_fft_double(pr_D, D, D);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*D++ *=  arrsz;
	}
	D    -= (ptrdiff_t)G->LengthF;
	
	
	
	/* call Jacobi */
	compute_complex_real_3d_fft_double(pr_PCF, PCF, PCF);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*PCF++ *=  arrsz;
	}
	PCF  -= (ptrdiff_t)G->LengthF;


	
	
	Jacobi(X, B, TF, PCF, PCF2, Wm, alpham, Wd, lamf, lam, alphas, swpspace, D, tol, NIter, NIterT, &omega, G, res_reporter, res_reporter_closure );
	
	/* Undo normalization */
	compute_real_complex_3d_fft_double(pf_PCF, PCF, PCF);
	compute_real_complex_3d_fft_double(pf_X, X, X);

	for ( I = G->LengthF/2L-1L;  I > -1L; I--  )	{
		tmpval = *PCF++;   PCF++;
		*X++   *= tmpval;
		*X++   *= tmpval;
	}
	
	PCF  -= (ptrdiff_t)G->LengthF;
	X    -= (ptrdiff_t)G->LengthF;
	
	
	compute_complex_real_3d_fft_double(pr_X, X, X);
	for ( I = G->LengthF-1L;  I > -1L;  I-- )	{
		*X++ *=  arrsz;
	}
	X  -= (ptrdiff_t)G->LengthF;
	
	
	
	/* Exit */
	destroy_fft_double_plan(pf_PCF);
	destroy_fft_double_plan(pr_PCF);
	destroy_fft_double_plan(pf_Wd);
	destroy_fft_double_plan(pr_Wd);
	destroy_fft_double_plan(pf_Imn);
	destroy_fft_double_plan(pr_Imn);
	destroy_fft_double_plan(pf_X);
	destroy_fft_double_plan(pr_X);
	destroy_fft_double_plan(pf_B);
	destroy_fft_double_plan(pr_B);
	destroy_fft_double_plan(pf_TFsq);
	destroy_fft_double_plan(pr_TFsq);
	destroy_fft_double_plan(pf_D);
	destroy_fft_double_plan(pr_D);

	
	free(mempoold);
	free(mempool);
	free(G);
	
	
	return 0;
	
	
}

