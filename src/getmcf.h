/*
 *  getmcf.h
 *  clgetmcf
 *
 *  Created by Muthuvel Arigovindan on 4/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

int		getmcf(	 double       *OTF,  
			   double       *OTFc,
			   long	      Nx,
			   long         Ny,
			   long         Nz,
			   double       NA,
			   double       nimm,
			   double       delxy,             
			   double       delz,
			   int          nthread );
