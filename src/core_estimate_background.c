#include "core_estimate_background.h"
#include "volumestruct.h"
#include "misc.h"
#include "volumestruct2d.h"
#include "misc2d.h"
#include "fft_double.h"
#include <stdlib.h>  /* free, malloc */
#include <string.h>  /* memset */


int
core_estimate_background(
    double* measurement,
    const long size[3],
    double delfac,
    double lamf,
    double lamsmooth,
    int nthread,
    double* p_background
)
{
    FFTDoublePlan pf_PCF, pf_m, pr_m;
    double* PCF;
    const double* in;
    double* out;
    long nxh;
    size_t i;
    double background;

    nxh = size[0] / 2 + 1;

    PCF = (double*)
	malloc(2 * nxh * (size_t) size[1] * (size_t) size[2] * sizeof(double));
    if (PCF == 0) {
	return 1;
    }

    if (initialize_fft_double() == 0) {
	free(PCF);
	return 1;
    }
    set_fft_double_thread_count(nthread);

    pf_PCF =
	plan_real_complex_3d_fft_double(size[0], size[1], size[2], PCF, PCF);
    pf_m = plan_no_overwrite_real_complex_3d_fft_double(
	size[0], size[1], size[2], measurement, measurement
    );
    pr_m = plan_no_overwrite_complex_real_3d_fft_double(
	size[0], size[1], size[2], measurement, measurement
    );
    if (pf_PCF == 0 || pf_m == 0 || pr_m == 0) {
	if (pr_m != 0) {
  	    destroy_fft_double_plan(pr_m);
	}
	if (pf_m != 0) {
  	    destroy_fft_double_plan(pf_m);
	}
	if (pf_PCF != 0) {
  	    destroy_fft_double_plan(pf_PCF);
	}
	free(PCF);
	return 1;
    }

    (void) memset(
	PCF, 0, 2 * nxh * (size_t) size[1] * (size_t) size[2] * sizeof(double)
    );
    if (size[2] == 1) {
	AddDuchon2D(PCF, size[0], size[1]);
    } else {
	AddDuchon3D(PCF, size[0], size[1], size[2], delfac);
    }
    compute_real_complex_3d_fft_double(pf_PCF, PCF, PCF);

    compute_real_complex_3d_fft_double(pf_m, measurement, measurement);
    out = measurement;
    in = PCF;
    for (i = 0; i < nxh * (size_t) size[1] * (size_t) size[2]; ++i) {
	double filt = 1.0 + lamf * lamsmooth * in[i + i];

	out[i + i] /= filt;
	out[i + i + 1] /= filt;
    }

    compute_complex_real_3d_fft_double(pr_m, measurement, measurement);
    background = measurement[0];
    in = measurement;
    for (i = 0; i < (size_t) size[1] * (size_t) size[2]; ++i) {
	int j;

	for (j = 0; j < size[0]; ++j) {
	    if (background > in[j]) {
		background = in[j];
	    }
	}
	in += nxh + nxh;
    }
    *p_background = background /
        (size[0] * (double) size[1] * (double) size[2]);

    destroy_fft_double_plan(pf_PCF);
    destroy_fft_double_plan(pf_m);
    destroy_fft_double_plan(pr_m);
    free(PCF);

    return 0;
}
