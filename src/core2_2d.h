/*
 *  core2.h
 *  
 *
 *  Created by Muthuvel Arigovindan on 6/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "core2_params.h"
#include "res_reporter.h"

int         core2_2d   (  double   *Imn,
						double   *TF,
			                        double   *TFc,
						double   *X,
			                        double   *mask,
						long     *Sz,
			                        const Core2Params* p_params,
			                        ResReporter* p_resrep,
			                        int (*mask_updater)(
						    int iter,
						    int nx,
						    int ny,
						    const double* guess,
						    double* mask,
						    double* work,
						    void* closure
						),
			                        void* mask_updater_closure);
