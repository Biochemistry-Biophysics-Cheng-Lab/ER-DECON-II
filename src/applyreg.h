/*
 *  applyreg.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 5/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

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
					 FFTDoublePlan    pr_Rtmp  );



void      ApplyComplexTFp (    double   *X,
						   double   *TF,
						   long     Length  );





void      ApplyComplexTFpr       (  double   *X,
								  double   *TF,
								  long     Length  );
