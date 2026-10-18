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

#include  "volumestruct2d.h"
#include  "derivatives2d.h"
#include "misc2d.h"





void	RegRMax2D	(	double      *X,
					 double     *Rmax,
					 Geometry2D	*G, 
					 char	    *dir,
					 int		NF, 
					 double     *swpspace   )


{
	
	double   *Rs;
	long Nx, NxSkip, Ny;
	long I, J, K;
	double  tmpval,  maxval;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	Rs = swpspace;    
	
    for ( J = NF-1L; J > -1L; J-- )		{
		
	    Derivative2D( X, Rs, G, dir+(ptrdiff_t)(2*J), 'f', swpspace+(ptrdiff_t)(G->LengthF));
		
		maxval = 0.0f;
		for ( K = Ny-1L; K > -1L; K-- )	{
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




void		AddDuchon2D(  double  *V,
						  long    Nx,
						  long    Ny        )

{
	
	
	double   V1[5L] = {1, -4, 6, -4, 1};
	double    V2[3L] = {1, -2, 1};
	
	
	long shx[5L],  shy[5L];
	
	
	long I, J;
	
	long  BSY;
	
	BSY = Nx+((Nx % 2L == 0L) ? 2L : 1L);
	
	
	
	shx[0L] = Nx-2L; shx[1L] = Nx-1L; shx[2L] = 0L;  shx[3L] = 1L;  shx[4L] = 2L;
	shy[0L] = BSY*(Ny-2L); shy[1L] = BSY*(Ny-1L); shy[2L] = 0L;  shy[3L] = BSY;  shy[4L] = 2L*BSY;


	for ( I = 0L; I < 5L; I++ )		{
		*(V + (ptrdiff_t)(shx[I])) += V1[I];
		*(V + (ptrdiff_t)(shy[I])) += V1[I];
	}
	
	
	for ( I = 0L; I < 3L;  I++ )	{
		for ( J = 0L;  J < 3L; J++ )	{
			*(V + (ptrdiff_t)(shx[I+1L]+shy[J+1L])) += 2.0f*V2[I]*V2[J];		}
	}
	
	
	
	
}

