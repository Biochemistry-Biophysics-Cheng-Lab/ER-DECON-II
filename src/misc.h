/*
 *  misc.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
/*
void		RegSq	(	double      *X,
					    double      *R,
					    Geometry	*G, 
					    char	    *dir,
					    int		    NF, 
					    double      *lam,
					    double      *swpspace   );  */





void	RegRMax	(	double      *X,
				    double     *Rmax,
				    Geometry	*G, 
				    char	    *dir,
				    int  		NF, 
				    double      *swpspace   );






void      AddDuchon3D(  double  *V,
						long    Nx,
						long    Ny,
						long    Nz,
						double  delfac    );


