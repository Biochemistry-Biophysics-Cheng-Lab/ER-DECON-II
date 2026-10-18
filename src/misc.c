/*
 *  misc.c
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include  <stdlib.h>
#include  <stddef.h>
#include  <string.h>

/*#include "fft_double.h"*/
#include  "volumestruct.h"
#include  "derivatives.h"
#include "misc.h"


/*
int BackGroundSubtract(     double    *Imn,
							double    *Wienr,
							long      Nx,
							long      Ny,
							long      Nz,
							double    lamd       )
{
	
	
	double    *Ims;
	fftw_plan  pf_Ims,  pr_Ims;
	long      NxSkip, I, J, K;
	double  bg = 1.0e16, tmpv;
	

	NxSkip = (Nx % 2L == 0L) ? 2L : 1L;
	Ims = (double *)malloc((Nx+NxSkip)*Ny*Nz*sizeof(double));
	if ( Ims = (double *)NULL )		{
		return -1;
	}
	
	pf_Ims = plan_real_complex_3d_fft_double(Nx, Ny, Nz, Ims, Ims);
	pr_Ims = plan_complex_real_3d_fft_double(Nx, Ny, Nz, Ims, Ims);
	
	
	memcpy(Ims, Imn,  (Nx+NxSkip)*Ny*Nz*sizeof(double));
	
	compute_real_complex_3d_fft_double(pf_Ims, Ims, Ims);
	
	for ( I = (Nx/2+1L)*Ny*Nz-1L;  I > -1L;  I-- )		{
		*Ims++ /=  (1.0f+lamd*(*Wienr));
		*Ims++ /=  (1.0f+lamd*(*Wienr++));
		Wienr++;
	}
	Ims   -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	Wienr -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	
	compute_complex_real_3d_fft_double(pr_Ims, Ims, Ims);
	
	for ( K = Nz-1L;  K > -1L;  K-- )	{
		for ( J = Ny-1L;  J > -1L;  J--  )   {
			for ( I = Nx-1L;  I > -1L;  I-- )	{
				tmpv = *Ims++;
				if ( tmpv < 1.0e16 )		{
					bg = tmpv;
				}
			}
			Ims += NxSkip;
		}
	}
	Ims -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	
	
	for ( K = (Nx+NxSkip)*Ny*Nz-1L;  K > -1L; K-- )		{
		+Imn++ -=  bg;
	}
	Imn -= (ptrdiff_t)(Nx+NxSkip)*Ny*Nz;
	
	free(Ims);
	
	return 0;
	
	destroy_fft_double_plan(pf_Ims);
	destroy_fft_double_plan(pr_Ims);
}



void		ThreshOTF(	double  *TF,
					    long     Nx,
					    long     Ny,
						long     Nz,
					    double   dfr,
						double   phi,
					    double   thr       )   
{
	
	
	long  I, J, K;
	double  x, y, z;
	double  r,  ang;
	double  re, im;
	double  dc;
	
	dc = *TF;
	
	for ( I = Nz-1L; I > -1L;  I-- )	{
		for ( J = Ny-1L;  J > -1L;  J-- )	{
			for ( K = Nx/2L;  K > -1L;  K-- )	{
				
				
				if (  I <  (Nz/2L+1L)  )  {
					z = I;
				}
				else {  z = Nz - I;  }
				
				if (  J <  (Ny/2L+1L)  )  {
					y = J;
				}
				else {  y = Ny - J;  }
				x = K;
				
				
				r = sqrt(x*x + y*y);
				ang = atan(r/(dfr*z));
				re = *TF++;
				im = *TF--;
				r = sqrt(re*re+im*im)/dc;
				
				if (  (ang < phi) & (r < thr)  )   {
					*TF++ = 0.0f;
					*TF++ = 0.0f;
				}
				
			}
		}
	}
	TF -= (ptrdiff_t)(Nx+((Nx % 2L == 0L) ? 2L : 1L))*Ny*Nz;
	
	
}
	
				
				

				
				


void		RegSq	(	double      *X,
					    double      *R,
					    Geometry	*G, 
					    char	    *dir,
					    int		    NF, 
					    double      *lam,
					    double      *swpspace   )


{
	
	double   *Rs, lamtmp;
	long Nx, NxSkip, Ny, Nz;
	long I, J, K;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	Rs = swpspace;
	
	memset( R, 0, (Nx+NxSkip)*Ny*Nz*sizeof(double) );
    
	
    for ( J = NF-1L; J > -1L; J-- )		{
        
		lamtmp = lam[J];
		
	    Derivative( X, Rs, G, dir+(ptrdiff_t)(2*J), 'f', swpspace+(ptrdiff_t)(G->LengthF));
		
		
		for ( K = Ny*Nz-1L; K > -1L; K-- )	{
			for ( I = Nx - 1L; I > -1L; I-- )		{
				*Rs++ = lamtmp*(*Rs);
			}
			Rs += NxSkip;
		}
		Rs -= (ptrdiff_t)((Nx+NxSkip)*Ny*Nz);
		
		
		Derivative( Rs, R, G, dir+(ptrdiff_t)(2*J), 'r', swpspace+(ptrdiff_t)(G->LengthF));
		
		
	}
	
	
}




 */




void	RegRMax	(	double      *X,
					 double     *Rmax,
					 Geometry	*G, 
					 char	    *dir,
					 int		NF, 
					 double      *swpspace   )


{
	
	double   *Rs;
	long Nx, NxSkip, Ny, Nz;
	long I, J, K;
	double  tmpval,  maxval;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	Nz = G->Nz;
	
	Rs = swpspace;    
	
    for ( J = NF-1L; J > -1L; J-- )		{
		
	    Derivative( X, Rs, G, dir+(ptrdiff_t)(2*J), 'f', swpspace+(ptrdiff_t)(G->LengthF));
		
		maxval = 0.0f;
		for ( K = Ny*Nz-1L; K > -1L; K-- )	{
			for ( I = Nx - 1L; I > -1L; I-- )		{
				tmpval = *Rs++;
				tmpval *= tmpval;
				if ( maxval < tmpval )  { maxval = tmpval; }
			}
			Rs += NxSkip;
		}
		Rs -= (ptrdiff_t)(G->LengthF);
		Rmax[J] = maxval;
		
		
	}
	
	
}



/*
 */




void		AddDuchon3D(  double  *V,
						  long    Nx,
						  long    Ny,
						  long    Nz,
						  double  delfac   )

{
	
	
	double   V1[5L] = {1, -4, 6, -4, 1};
	double    V2[3L] = {1, -2, 1};
	
	
	long shx[5L],  shy[5L],  shz[5L];
	
	
	long I, J;
	
	long  BSY, BSZ;
	
	BSY = Nx+((Nx % 2L == 0L) ? 2L : 1L);
	BSZ = Ny*BSY;
	
	
	
	shx[0L] = Nx-2L; shx[1L] = Nx-1L; shx[2L] = 0L;  shx[3L] = 1L;  shx[4L] = 2L;
	shy[0L] = BSY*(Ny-2L); shy[1L] = BSY*(Ny-1L); shy[2L] = 0L;  shy[3L] = BSY;  shy[4L] = 2L*BSY;
	shz[0L] = BSZ*(Nz-2L); shz[1L] = BSZ*(Nz-1L); shz[2L] = 0L;  shz[3L] = BSZ;  shz[4L] = 2L*BSZ;


	for ( I = 0L; I < 5L; I++ )		{
		*(V + (ptrdiff_t)(shx[I])) += V1[I];
		*(V + (ptrdiff_t)(shy[I])) += V1[I];
		*(V + (ptrdiff_t)(shz[I])) += delfac*delfac*delfac*delfac*V1[I];
	}
	
	
	for ( I = 0L; I < 3L;  I++ )	{
		for ( J = 0L;  J < 3L; J++ )	{
			*(V + (ptrdiff_t)(shx[I+1L]+shy[J+1L])) += 2.0f*V2[I]*V2[J];
			*(V + (ptrdiff_t)(shy[I+1L]+shz[J+1L])) += delfac*delfac*2.0f*V2[I]*V2[J];
			*(V + (ptrdiff_t)(shx[I+1L]+shz[J+1L])) += delfac*delfac*2.0f*V2[I]*V2[J];
		}
	}
	
	
	
	
}

