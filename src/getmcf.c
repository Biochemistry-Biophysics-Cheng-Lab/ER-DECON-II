/*************************************************************************************************
 *
 *  getmcf.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 4/13/10.
 *  
 *
 ************************************************************************************************/



#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include  "fft_double.h"


#include "getmcf.h"


int		getmcf(	 double       *OTF,  
			     double       *OTFc,
			     long	      Nx,
			     long         Ny,
			     long         Nz,
			     double       NA,
			     double       nimm,
			     double       delxy,             
			     double       delz,
                             int          nthread      )
{
	
	long    I, J, K;
	double  re, im, r;
	double  mx;
	double  ang,  angc;
	double  delfx, delfy, delfz;
	double  arrsz;
	FFTDoublePlan planf, planr;		
	
	delfx = 1.0 / (delxy * Nx);
	delfy = 1.0 / (delxy * Ny);
	delfz = 1.0 / (delz * Nz);
	arrsz = 1.0/((double)(Nx*Ny*Nz));
	
	if (initialize_fft_double() == 0) {
	    return -1;
	}
	set_fft_double_thread_count(nthread);

	planf = plan_real_complex_3d_fft_double(Nx, Ny, Nz, OTFc, OTFc);
	planr = plan_complex_real_3d_fft_double(Nx, Ny, Nz, OTFc, OTFc);
	if (planf == 0 || planr == 0) {
	    if (planf != 0) {
		destroy_fft_double_plan(planf);
	    }
	    if (planr != 0) {
		destroy_fft_double_plan(planr);
	    }
	    return -1;
	}

	ang = M_PI / 2.0 - asin(NA / nimm);
	
	mx = OTF[0L]*OTF[0L];
	
	for ( K = 0L;  K < Nz;  K++ )  {
	        double z;

		if ( K < Nz/2L+1L )  {
		    z = delfz * K;
		}
		else {
		    z = delfz * (Nz - K);
		}
		for ( J = 0L;  J < Ny;  J++ )   {
		        double y;

		        if ( J < Ny/2L+1L )  {
			    y = delfy * J;
			}
			else {
			    y = delfy * (Ny - J);
			}
			for ( I = 0L;  I < Nx/2+1L;  I++ )		{
			        double x = delfx * I;
				
				r = sqrt(x*x + y*y);
				angc = atan(r/z);
				
				re = *OTF++;
				im = *OTF++;
				
				*OTFc++  = angc < ang ? sqrt(mx - re*re - im*im) :   0.0;
				*OTFc++  = 0;
				
			}
		}
	}
	OTFc -=  (ptrdiff_t)((Nz)*Ny*2L*(Nx/2L+1L));
	OTF  -=  (ptrdiff_t)((Nz)*Ny*2L*(Nx/2L+1L));
	
   compute_complex_real_3d_fft_double(planr, OTFc, OTFc);
							   
							
   for ( K = 0L;  K < Nz;  K++ )  {
           double z;

	   if ( K < Nz/2L+1L )  {
	       z  = K;
	   }
	   else {
	       z = Nz - K;
	   }
	   for ( J = 0L;  J < Ny;  J++ )   {
	           double y;

		   if ( J < Ny/2L+1L )  {
		       y = J;
		   }
		   else {
		       y = Ny - J;
		   }
		   for ( I = 0L;  I < Nx;  I++ )		{
		           double x;   

			   if ( I < Nx/2L+1L )  {
				   x  = I;
			   }
			   else {
				   x = Nx - I;
			   }
			   
			   r = x*x + y*y + z*z;
			   
			   *OTFc++ *= arrsz*exp(-r/100.0);
			   
		   }
		   OTFc++;OTFc++;
	   }
   }
   OTFc -=  (ptrdiff_t)(Nz*Ny*2L*(Nx/2L+1L));  
	
	
							   
   compute_real_complex_3d_fft_double(planf, OTFc, OTFc);
	
   destroy_fft_double_plan(planf);
   destroy_fft_double_plan(planr);

   return 0;
	

}

