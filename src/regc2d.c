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

#include  "volumestruct2d.h"
#include  "derivatives2d.h"
#include "regc2d.h"

/******************************************************************************************/



void Regsum2D (
    const double *X,
    double *R,
    double *D,
    const double *lamratio,
    double laml2,
    double lampos,
    const double *beta,
    const double *alpha,
    const Geometry2D *G,
    double *swpspace
)
{
	
	
	long      J, K;
    double    *Rt, *Dt;
	double    lam1, lam2[4L], lamr1, lamr2;
	char      dir1[2L] = "  ";
	char      dir2[8L] = "  xxyyxy";
	
	long       Nx, NxSkip, Ny;
	
	
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	lamr1 =  lamratio[0L];
	lamr2 =  1.0 - lamratio[0L];
	
	lam1     = 1.0;
	
	lam2[0L] = lamratio[1L];
	lam2[1L] = 1.0;
	lam2[2L] = 1.0;
	lam2[3L] = 2.0;
	
	
	Rt = swpspace;
	Dt = swpspace + (ptrdiff_t)(G->LengthF);
		
	
	 Reg2D ( X,  Rt, Dt, dir1, 1L,  &lam1, laml2,
		 beta[0L], alpha[0L], G, swpspace+(ptrdiff_t)(2L)*(G->LengthF) );
	
	
	
	
	
	for ( J = Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			*R++ =  lamr1*(*Rt++);
			*D++ =  lamr1*(*Dt++);
		}
		R += NxSkip;
		D += NxSkip;
		Rt += NxSkip;
		Dt += NxSkip;
	}
	R  -= (ptrdiff_t)G->LengthF;
	Rt -= (ptrdiff_t)G->LengthF;
	D  -= (ptrdiff_t)G->LengthF;
	Dt -= (ptrdiff_t)G->LengthF;
	
	

	
	 Reg2D ( X,  Rt, Dt, dir2, 4L,  lam2, laml2,
	  	 beta[1L], alpha[1L], G, swpspace+(ptrdiff_t)(2L)*(G->LengthF) );
	
	
    /*
	for ( J = 0L;  J < 5L; J++ )	{
		printf("%2.10f   ", Rt[J]);
	}
	printf("\n");
	for ( J = 0L;  J < 5L; J++ )	{
		printf("%2.10f   ", Dt[J]);
	}
	printf("\n");  */
	
	
	for ( J = Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			*R++ +=  lamr2*(*Rt++);
			*D++ +=  lamr2*(*Dt++);
		}
		R += NxSkip;
		D += NxSkip;
		Rt += NxSkip;
		Dt += NxSkip;
	}
	R  -= (ptrdiff_t)G->LengthF;
	Rt -= (ptrdiff_t)G->LengthF;
	D  -= (ptrdiff_t)G->LengthF;
	Dt -= (ptrdiff_t)G->LengthF;
	
	
	/*
	for ( J = 0L;  J < 5L; J++ )	{
		printf("%2.10f   ", R[J]);
	}
	printf("\n");
	for ( J = 0L;  J < 5L; J++ )	{
		printf("%2.10f   ", D[J]);
	}
	printf("\n");  */
	
	
	for ( J = Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			*Dt =  (*X < 0) ? 1.0 : 0.0;
			*Rt =   (*Dt)*(*X);
			Dt++;Rt++;X++;
		}
		Rt += NxSkip;
		Dt += NxSkip;
		X += NxSkip;
	}
	Rt -= (ptrdiff_t)G->LengthF;
	X  -= (ptrdiff_t)G->LengthF;
	Dt -= (ptrdiff_t)G->LengthF;
	
	
		
	
	for ( J = Ny-1L; J > -1L; J-- )		{
		for ( K = Nx-1L;  K > -1L; K-- )		{
			*R++ +=  lampos*(*Rt++);
			*D++ +=  lampos*(*Dt++);
		}
		R += NxSkip;
		D += NxSkip;
		Rt += NxSkip;
		Dt += NxSkip;
	}
	R  -= (ptrdiff_t)G->LengthF;
	Rt -= (ptrdiff_t)G->LengthF;
	D  -= (ptrdiff_t)G->LengthF;
	Dt -= (ptrdiff_t)G->LengthF;
	
	
	
}
		
		
	
/*
 */



void Reg2D (
    const double *X,
    double *R,
    double *D,
    const char *dir,
    long NF,
    const double *lamc,
    double laml2,
    double beta,
    double alpha,
    const Geometry2D *G, 
    double *swpspace
)
{


	
	long Nx, NxSkip, Ny, LengthF;
	long  I, J;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	LengthF = G->LengthF;
	
	Regr2D(  X,  R, G, dir, NF, lamc,  swpspace );	
	
	
	
	
	for ( I = Ny-1L;  I > -1L; I-- )  {
		for ( J = Nx-1L; J > -1L; J-- )	{
			*D++ = 1.0/( 1+alpha*(*R++) ) + laml2;
		}
		D += NxSkip;
		R += NxSkip;
	}
	D  -= (ptrdiff_t)LengthF;
	R  -= (ptrdiff_t)LengthF;
	
	
			
	
	RegNL2D ( X, D, R, G, dir, NF, lamc,  swpspace );	
	
	
	
	
	Regs2D(  D,  swpspace+(ptrdiff_t)2L*G->LengthF, G, dir, NF, lamc,  swpspace );
	memcpy(D,  swpspace+(ptrdiff_t)2L*G->LengthF,  G->LengthF*sizeof(double));
	
	
	
	
	
	
}


	
/*
 */



void RegNL2D (
    const double *X,
    const double *D,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF, 
    const double *lamc,
    double *swpspace
)
{
	
	double   *Rs, lam;
	long Nx, NxSkip, Ny;
	long I, J, K;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	
	Rs = swpspace;
	
	memset( R, 0, (Nx+NxSkip)*Ny*sizeof(double) );
    
	
    for ( J = 0L; J < NF; J++ )		{
        
		lam = lamc[J];
		
		if ( lam > 0.0 )   {
			
				
			
			 Derivative2D( X, Rs, G, dir+(ptrdiff_t)(2*J), 'f', swpspace+(ptrdiff_t)(G->LengthF));
			
			


			
			for ( K = Ny-1L; K > -1L; K-- )	{
				for ( I = Nx - 1L; I > -1L; I-- )		{
					*Rs = lam*(*Rs)*(*D++);  Rs++;
				}
				Rs += NxSkip;
				D += NxSkip;
			}
			Rs -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			D  -= (ptrdiff_t)((Nx+NxSkip)*Ny); 
			
			
				
			
			Derivative2D( Rs, R, G, dir+(ptrdiff_t)(2*J), 'r', swpspace+(ptrdiff_t)(G->LengthF));
			
			
			
				
			
			
		}
		
		
	}
	
	
}





/*
 */


    
	

void Regr2D (
    const double *X,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
)
{
	
	double   *Rs, lam, tmp;
	long Nx, NxSkip, Ny, LengthF;
	long I, J, K;
	char dirc[2L];
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	LengthF = G->LengthF;
	
	Rs = swpspace;	    
	memset( R, 0, LengthF*sizeof(double) );
	
	
    for ( J = 0L; J < NF; J++ )		{
        
		lam = lamc[J];
		
		if ( lam > 0.0 )  {
			
			dirc[0L] = *(dir+(ptrdiff_t)(J*2));
			dirc[1L] = *(dir+(ptrdiff_t)(J*2L+1L));

			Derivative2D( X, Rs, G, dirc, 'f', swpspace+(ptrdiff_t)(G->LengthF));
			
			if ( (dirc[0L] == ' ') && (dirc[1L] == ' ')  )   {
				for ( K = Ny-1L; K > -1L; K-- )	{
					for ( I = Nx - 1L; I > -1L; I-- )		{
						tmp    = *Rs;
						*Rs++  = tmp > 0 ?  tmp : 0; 
					}
					Rs += NxSkip;
				}
				Rs -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			}
			
			for ( K = Ny-1L; K > -1L; K-- )	{
				for ( I = Nx - 1L; I > -1L; I-- )		{
					*R++ += lam*(*Rs)*(*Rs);
					Rs++;
				}
				R += NxSkip;
				Rs += NxSkip;
			}
			Rs -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			R  -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			
		}
		
		
		
	}
	
	
}



/*
*/




void Regs2D (
    const double *X,
    double *R,
    const Geometry2D *G,
    const char *dir,
    int NF,
    const double *lamc,
    double *swpspace
)
{
	
	double   *Rs, lam;
	long Nx, NxSkip, Ny,  LengthF;
	long I, J, K;
	
	Nx = G->Nx;
	NxSkip = G->NxFSkip;
	Ny = G->Ny;
	
	LengthF = G->LengthF;
	
	Rs = swpspace;
	memset( R, 0, LengthF*sizeof(double) );
   
	
    for ( J = 0L; J < NF; J++ )		{
        
		lam = lamc[J];
		
		
		if ( lam > 0.0 )  {
			
			Derivative2D( X, Rs, G, dir+(ptrdiff_t)(J*2), 's', swpspace+(ptrdiff_t)(G->LengthF));
			
			
			
			for ( K = Ny-1L; K > -1L; K-- )	{
				for ( I = Nx - 1L; I > -1L; I-- )		{
					*R++ += lam*(*Rs);
					Rs++;
				}
				R += NxSkip;
				Rs += NxSkip;
			}
			Rs -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			R  -= (ptrdiff_t)((Nx+NxSkip)*Ny);
			
			
			
		}
		
		
		
	}
	
	
}
