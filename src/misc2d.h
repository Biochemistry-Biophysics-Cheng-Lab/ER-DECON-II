/*
 *  misc.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
/*
void		RegSq2D	(	double      *X,
					    double      *R,
					    Geometry2D	*G, 
					    char	    *dir,
					    int		    NF, 
					    double      *lam,
					    double      *swpspace   );  */





void	RegRMax2D	(	double      *X,
				    double      *Rmax,
				    Geometry2D	*G, 
				    char	    *dir,
				    int  		NF, 
				    double      *swpspace   );






void      AddDuchon2D(  double  *V,
						long    Nx,
						long    Ny        );


