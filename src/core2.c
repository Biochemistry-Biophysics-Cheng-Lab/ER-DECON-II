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
#include     <limits.h>      /* INT_MAX */
#include     <math.h>
#include     <string.h>
#include     <stdio.h>

#include     "fft_double.h"
#include     "volumestruct.h"
#include     "regc.h"
#include     "misc.h"

#include     "jacobi2.h"
#include     "JacobiNL3Ddu.h"
#include     "core2.h"



int         core2      (  double   *Imn,         /* Input image in 2*(Nx/2+1)*Ny*Nz  format */
						  double   *TF,          /* OTF   */
						  double   *TFc,         /* Missing cone filter; not used if lampc is less than or equal to zero */
						  double   *X,           /* Allocation for solution 2*(Nx/2+1)*Ny*Nz   */
			                          double   *mask,        /* Spatial mask dimensioned 2*(Nx/2+1)*Ny*Nz; may be null if lammask is not positive */
						  long     *Sz,          /* [Nx Ny Nz]  */ 
						  double   delfrac,       /* dxy/dz  */  
			                          const Core2Params* p_params,
			                          ResReporter* p_resrep,
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
			                          void* mask_updater_closure
)
{
	Geometry G;
	long *mempool = 0;
	double *mempoold = 0;
	double *swpspace;

	double      *B, *PCF, Diagd, regt;

	long	 Nx, Ny, Nz;
	long     I, J;
	
	double   lampc;
	double   arrsz;
	double   tmpval;
	double   re1, re2, im1, im2, tmpvalre, tmpvalim;
	double   red, imd;
	
	FFTDoublePlan pf_PCF = 0, pf_Imn = 0, pr_Imn = 0, pr_X = 0, pf_B = 0,
	    pr_B = 0;
	int status = -1;

	if (initialize_fft_double() == 0) {
	    goto cleanup;
	}
	set_fft_double_thread_count(p_params->nthread);
#if defined(_OPENMP) && defined(USE_OPENMP)
	omp_set_num_threads(p_params->nthread);
#endif

	Nx = Sz[0L];
	Ny = Sz[1L];
	Nz = Sz[2L];
	arrsz = 1.0/((double)(Nx*Ny*Nz));

	mempool = (long *)malloc( 7L*(2L*Nx*Ny + (Nz-2L)*(2L*Nx+2L*Ny-4L))*sizeof(long) );
	if ( mempool == (long *)NULL )	{
	    goto cleanup;
	}
	ConstructVolumeStruct( &G, Nx, Ny, Nz, mempool );
	mempoold = (double *)malloc(13L*(Nx+G.NxFSkip)*Ny*Nz*sizeof(double));  
	if ( mempoold == (double *)NULL )	{
	    goto cleanup;
	}
	B      = mempoold;
	PCF    = mempoold + (ptrdiff_t)1L*(Nx+G.NxFSkip)*Ny*Nz;
	
	swpspace =  mempoold + (ptrdiff_t)2L*(Nx+G.NxFSkip)*Ny*Nz;

	lampc = p_params->lampc / p_params->lamf;
	
	/* FFTs */
	pf_PCF    =  plan_real_complex_3d_fft_double(Nx, Ny, Nz, PCF, PCF);
	pf_Imn    =
	    plan_no_overwrite_real_complex_3d_fft_double(Nx, Ny, Nz, Imn, Imn);
	pr_Imn    =
	    plan_no_overwrite_complex_real_3d_fft_double(Nx, Ny, Nz, Imn, Imn);
	if (p_params->skips_initial_guess) {
	    pr_X = plan_no_overwrite_complex_real_3d_fft_double(
		Nx, Ny, Nz, X, X
	    );
	} else {
	    pr_X = plan_complex_real_3d_fft_double(Nx, Ny, Nz, X, X);
	}
	pf_B      =  plan_real_complex_3d_fft_double(Nx, Ny, Nz, B, B);
	pr_B      =  plan_complex_real_3d_fft_double(Nx, Ny, Nz, B, B);
	if (pf_PCF == 0 || pf_Imn == 0 || pr_Imn == 0 || pr_X == 0 ||
	    pf_B == 0 || pr_B == 0) {
	    goto cleanup;
	}
	memset(PCF, 0, (G.LengthF)*sizeof(double));
	AddDuchon3D(PCF, Nx, Ny, Nz, delfrac);
	compute_real_complex_3d_fft_double(pf_PCF, PCF, PCF);	

	/* Compute Tikhonov solution */
	{
	    double lr0, lr1;

	    if (p_params->uses_old_line_search) {
		lr0 = p_params->lamratio[0];
		lr1 = p_params->lamratio[1];
	    } else {
		lr0 = 0.0;
		lr1 = (p_params->regularization_type == 1) ?
		    p_params->lamratio[1] : 0.0;
	    }
	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- ) {
		*PCF = (1.0 + p_params->laml2) *
		    (lr0 + (1.0 - lr0) * (lr1 + (*PCF)));
		PCF += 2;
	    }
	    PCF -= (ptrdiff_t)G.LengthF;
	}
	
	/*
	for ( I = 0L;  I < 5L; I++ )	{
	 printf("%3.10f   ", PCF[(Nz/2-I)*(Nx+G.NxFSkip)*Ny]);
	}
	printf("\n"); */
	
	if (lampc > 0.0) {
	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- )		{
		
		re1  = *TF++;
		im1  = *TF++;
		
		re2  = *TFc++;  
		im2  = *TFc++;
		
		regt  = *PCF;
		
		tmpval = 1.0 /
		    (re1*re1 + im1*im1 + p_params->lamf *
		     (lampc * (re2*re2 + im2*im2) + regt));
		*PCF++ = tmpval;  PCF++;
	    }
	    TFc  -= (ptrdiff_t)G.LengthF;
	} else {
	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- )		{
		
		re1  = *TF++;
		im1  = *TF++;
		
		regt  = *PCF;
		
		tmpval = 1.0/(re1*re1 + im1*im1 + p_params->lamf * regt);
		
		*PCF++ = tmpval;  PCF++;
	
	    }
	}
	TF   -= (ptrdiff_t)G.LengthF;
	PCF  -= (ptrdiff_t)G.LengthF;
	
	/*
	for ( I = 0L;  I < 5L; I++ )	{
		printf("%3.10e   ", PCF[(Nz/2-I)*(Nx+G.NxFSkip)*Ny]);
	}
	printf("\n"); */
	
	
	if (p_params->skips_initial_guess) {
	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- )		{
		*PCF = sqrt(*PCF);
		PCF += 2;
	    }
	    PCF  -= (ptrdiff_t)G.LengthF;
	} else {
	    compute_real_complex_3d_fft_double(pf_Imn, Imn, Imn);

	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- )		{
		re1  = *TF++;
		im1  = *TF++;
		*Imn *= arrsz;
		red  = *Imn++;
		*Imn *= arrsz;
		imd  = *Imn++;
		tmpvalre = re1*red  + im1*imd;
		tmpvalim = re1*imd  - im1*red;
		*X++ = tmpvalre;
		*X++ = tmpvalim;
	    }
	    TF   -= (ptrdiff_t)G.LengthF;
	    Imn  -= (ptrdiff_t)G.LengthF;
	    X    -= (ptrdiff_t)G.LengthF;

	    for ( I = G.LengthF/2L - 1L; I > -1L;  I-- )		{
		tmpval = *PCF;
		*PCF++ = sqrt(tmpval);  PCF++;
		*X++ *= tmpval;
		*X++ *= tmpval;
	    }
	    PCF  -= (ptrdiff_t)G.LengthF;
	    X    -= (ptrdiff_t)G.LengthF;
	    compute_complex_real_3d_fft_double(pr_Imn, Imn, Imn);
	    compute_complex_real_3d_fft_double(pr_X, X, X);
	}

	/* Compute alphas */
	/*
    Rmax[0L] = 1.0;  Rmax[1L] = 1.0;	
	for ( I = 0L;  I < 2L;  I++ )   {
		Rmax[I] =  ( Rmax[I] >= 1.0f ) ?  Rmax[I] : 1.0f;
	}
	for ( I =0L;  I < 2L;  I++  )  {
		alpha[I] =  alpha[I]/(Rmax[I]);
	}
	*/	
		
	
	/* Compute RHS */
	for ( I = G.LengthF-1L;  I > -1L;  I-- )	{
		*B++ = (*Imn++)*arrsz;
	}
	B    -= (ptrdiff_t)G.LengthF;
	Imn  -= (ptrdiff_t)G.LengthF;
	compute_real_complex_3d_fft_double(pf_B, B, B);
    ApplyComplexTFpr(B, TF, G.LengthF);
        compute_complex_real_3d_fft_double(pr_B, B, B);
	
	/*
	for ( I = 0L;  I < 5L;  I++  )		{
		printf("%2.10f    ", B[I]);
	}
	printf("\n");  */
	

	
	
	/* Compute diagonal of data part */
	Diagd = 0.0;
	for ( I = G.Nz*Ny-1L; I > -1L; I-- ) {
	    double v0, vnyq, sum_other = 0.0;

	    v0 = (*TF) * (*TF) + (*(TF + 1)) * (*(TF + 1));
	    TF += 2;
	    for (J = (Nx-G.NxFSkip-2L)/2; J > -1L; --J) {
		sum_other += (*TF) * (*TF) + (*(TF + 1)) * (*(TF + 1));
		TF += 2;
	    }
	    if (G.NxFSkip == 2L) {
		vnyq = (*TF) * (*TF) + (*(TF + 1)) * (*(TF + 1));
		TF += 2;
	    } else {
		vnyq = 0.0;
	    }
	    Diagd += v0 + vnyq + 2.0 * sum_other;
	}
	TF -= G.LengthF;
	Diagd *= arrsz;

	/* call Jacobi */
	if (p_params->uses_old_line_search) {
	    double omega = p_params->omega;

	    JacobiNL3D  ( X, TF, Diagd, B, PCF, p_params->lamf,
			  p_params->lamratio, p_params->lampos,
			  delfrac, p_params->beta, p_params->alpha, &G,
			  swpspace, &omega, p_params->niter,
			  p_params->tol, p_resrep->old_res_reporter,
			  p_resrep->old_res_closure, p_params->laml2 );
	} else {
	    double cost;
	    int result = JacobiNL3Ddu(
		Nx,
		Ny,
		Nz,
		p_params->lamf,
		(p_params->regularization_type != 1) ?
		    0.0 : p_params->lamratio[1],
		p_params->lampos,
		p_params->alpha[0],
		p_params->laml2,
		p_params->lammask,
		delfrac,
		p_params->regularization_type,
		p_resrep->res_reporter,
		p_resrep->res_closure,
		mask_updater,
		mask_updater_closure,
		p_params->outputs_rfactor,
		(p_params->niter > INT_MAX) ? INT_MAX : (int) p_params->niter,
		B,
		TF,
		Imn,
		Diagd,
		PCF,
		X,
		mask,
		&cost,
		swpspace
	    );

	    if (result != 0) {
		goto cleanup;
	    }
	}

	/* Exit */
	status = 0;
 cleanup:
	if (pr_B != 0) {
	    destroy_fft_double_plan(pr_B);
	}
	if (pf_B != 0) {
	    destroy_fft_double_plan(pf_B);
	}
	if (pr_X != 0) {
	    destroy_fft_double_plan(pr_X);
	}
	if (pr_Imn != 0) {
	    destroy_fft_double_plan(pr_Imn);
	}
	if (pf_Imn != 0) {
	    destroy_fft_double_plan(pf_Imn);
	}
	if (pf_PCF != 0) {
	    destroy_fft_double_plan(pf_PCF);
	}

	if (mempoold != 0) {
	    free(mempoold);
	}
	if (mempool != 0) {
	    free(mempool);
	}

	return status;
}
