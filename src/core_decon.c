/*
 * Implements a command-line front-end to Muthuvel Arigovindan's ER Type I
 * deconvolution algorithm.  Uses Priism's libraries for input and output of
 * MRC format files.
 *
 * Usage:
 *     core_decon measurement result CTF [options]
 *
 * where options can include any of the following in any order:
 *
 * -x=start[:end]
 *     Limits the processed data in the x dimension to be the points with
 *     x indices, i, such that start <= i <= end.  start and end must be
 *     non-negative integers.  If end is omitted, it is equivalent to setting
 *     end to the last possible index.  When no -x option is given,
 *     the entire range in x is processed.
 *
 * -y=start[:end]
 *     Limits the processed data in the y dimension.  Behaves the same as the
 *     -x option but for the y dimension.
 *
 * -z=start[:end]
 *     Limits the processed data in the z dimension.  Behaves the same as the
 *     -x option but for the z dimension.
 *
 * -t=start[:end[:step]]
 *     Limits the range of time points to be processed.  Only time points
 *     equal to start + i * step where i is a non-negative integer and
 *     start + i * step is less than or equal to end are processed.  start,
 *     end, and step must be non-negative integers.  If step is omitted, it
 *     is equivalent to setting step to be one.  If end is omitted, it is
 *     equivalent to setting end to be the last possible index.  When no -t
 *     option is given, all time points are processed.
 *
 * -w=w1[:w2...]
 *     Limits the wavelengths to process.  A wavelength may be specified by
 *     its zero-based index in the input file or by its wavelength.  By
 *     default, all wavelengths are processed.
 *
 * -abs_background
 *     If you supply this option, the program then interprets the values
 *     specified with the -sub option as absolute values which will be
 *     subtracted from the input data.  If you do not supply this option,
 *     the values specified with -sub option act as multipliers for the
 *     background level that is estimated internally.
 *
 * -alpha=a
 *     Sets the nonlinearity.  If this option is omitted, the default value is
 *     1000.
 *
 * -cuth=t
 *     Internally, the intensities in each volume are rescaled to reduce
 *     the sensitivity of the optimal parameter values to the range of values
 *     in the data.  After deconvolution, the intensities in the result are
 *     rescaled by the inverse of the factor applied to original value.  The
 *     rescaled intensity of one prior to deconvolution matches the intensity
 *     in the original volume at which the fraction, t, of the values have
 *     greater intensities.  t must be greater than or equal to zero and less
 *     than one.  If this option is omitted, the default value of t is 0.001.
 *
 * -lamratio=r
 *     If omitted the default value is 0.5.
 *
 * -lamf=f
 *     Sets the smoothing parameter.  The user will have to tune this for
 *     the data set.
 *
 * -lampc=p
 *     Sets the weight for the missing cone filter.  If this option is omitted,
 *     the weighting for the missing cone filter is one hundred.
 *
 * -lampos=p
 *     Controls the positivity constraint.  If this option is omitted, the
 *     default value is one hundred.
 *
 * -lamsmooth=s
 *     Sets an auxiliary smoothing parameter.  If this option is omitted,
 *     the auxiliary smoothing parameter is one hundred.
 *
 * -na=v
 *     Sets the numerical aperture (used in the computation of the missing
 *     cone filter) to be v.  If this option is omitted, the numerical aperture
 *     will be assumed to be 1.4.
 *
 * -nimm=n
 *     Sets the immersion media refractive index (used in the computation of
 *     the missing cone filter) to be n.  If this option is omitted, the
 *     refractive index of the immersion media is assumed to be 1.512.  The
 *     refractive index of the immersion media must be greater than or equal
 *     to the numerical aperture.
 *
 * -ncycl=n
 *     Set ths number of iterations.  n must be a positive integer.  If this
 *     option is omitted, fifty iteration will be used.
 *
 * -nonorm_otf
 *     Do not rescale the tranfer function to have a zero frequency of response
 *     of one.  Note that radially averaged transfer functions are never
 *     rescaled.
 *
 * -np=n
 *     Sets the number of threads to be n.  This option is ignored if the
 *     environment variable, NSLOTS, is set.  The default value for the
 *     number of threads is the larger of one or the value of the environment
 *     variable, NSLOTS.  If the environment variable, NSLOTS, is not set or is
 *     not set to an integer value, the default number of threads is the
 *     number of processors on the system.
 *
 * -nzpad=n
 *     Specifies the padded size (in samples) to use in the z dimension
 *     during deconvolution.  If this option is omitted, the padded size will
 *     be the smallest power of two greater than or equal to 1.3 times the
 *     size in z of the processed region.
 *
 * -omega=o
 *     Sets the relaxation parameter.  If this option is omitted, the
 *     relaxation parameter will be 0.8.
 *
 * -oplotfile=pltname
 *     Specifies that histories of the termination criteria as a function of
 *     the iteration number be written to the file named pltname.
 *
 * -sub=b1[:b2...]
 *     The behavior of this option depends on whether the -abs_background
 *     has also been specified on the command line.  When -abs_background
 *     does not appear on the command line, the program computes an internal
 *     estimate of the background for each wavelength, multiplies that
 *     estimate by the value (b1 for wavelength 1, b2 for wavelength 2, ...)
 *     supplied by the -sub option, and then subtracts the result from the
 *     incoming data.  If the -sub option does not appear on the command line
 *     or does not set the value for a wavelength, the program acts as if the
 *     value for that wavelength is one.  When the -abs_background option
 *     appears on the command line, the program directly subtracts the value
 *     supplied by the -sub option (b1 for wavelength 1, b2 for wavelength2,
 *     ...) from the incoming data.  In that case, if the -sub option does not
 *     appear on the command line or does not set the value for a wavelength,
 *     the program acts as if the value for that wavelength is zero.
 *
 * -tol=t
 *     Specifies the termination tolerance.  If this option is omitted,
 *     the termination tolerance will be 0.001.
 */

#include "ive_standards.h"  /* for ive_thread.h */
#include "ive_thread.h"
#include "core.h"
#include "core_estimate_background.h"
#include "getmcf.h"
#include "interp_util.h"
#include "binary_heap.h"
#include "IMInclude.h"
#include <errno.h>
#include <float.h>       /* DBL_MAX, FLT_MAX */
#include <limits.h>      /* CHAR_BIT, INT_MAX, SIZE_MAX */
#include <math.h>        /* fabs, floor */
#include <stdarg.h>      /* va_end, va_start, vfprintf */
#include <stdio.h>       /* fprintf, printf, sscanf */
#include <stdlib.h>      /* free, getenv, malloc, strtol */
#include <string.h>      /* memcpy, strerror, strncmp */
#ifndef SIZE_MAX
#define SIZE_MAX ((((((size_t) 1) << (sizeof(size_t) * CHAR_BIT - 1)) - 1) << 1) + 1)
#endif


typedef struct COREParameters {
    double lamratio;
    double lamf;
    double lampos;
    double lampc;
    double lamsmooth;
    double alpha;
    double tol;
    double omega;
    long iteration_count;
} COREParameters;

typedef struct FrontEndParameters {
    char* measurement_file_name;
    char* result_file_name;
    char* ctf_file_name;
    char* plot_file_name;
    double background[IW_MAX_WAVE];
    double cuth;
    double na;
    double nimm;
    int xbounds[2];
    int ybounds[2];
    int zbounds[2];
    int tbounds[3];
    int wavelength_indices[IW_MAX_WAVE + 1];
    int xpadded, ypadded, zpadded;
    int absolute_background;
    int normalize_otf;
    int nthread;
} FrontEndParameters;

typedef struct ResReporterPt {
    double res;
    /*
     * The absolute value of this is the total number of iterations at the
     * time res was computed.  If iteration is negative, the iteration was
     * successful.
     */
    long iteration;
} ResReporterPt;

typedef struct ResReporterData {
    char* filename;
    FILE* fo;
    ResReporterPt* values;
    double last_successful_res;
    long last_successful_iteration;
    long ncount, nalloc;
    int status;
} ResReporterData;


static int complete_plot(
    int iw, int it, int iz, int in_im, ResReporterData* prrd
);
static int complete_plot_file(ResReporterData* prrd);
static int extract_doubles(const char* s, int max_count, double* values);
static int extract_ints(const char* s, int max_count, int* values);
static double* fetch_ctf(
    const FrontEndParameters* pfp, int in_im, int ctf_im, int iw
);
static double* fetch_measurement(
    const FrontEndParameters* pfp, int in_im, int iw, int it
);
static int find_free_stream(void);
static int find_matching_wavelength_mrc(int mrc_stream, float wavelength_nm);
static void force_conj_sym_c(int nxh, int ny, int nz, int even_x, double* io);
static int get_default_thread_count(void);
static int interpolate_3d_cc(
    const int nxyzi[3],
    const float dkxyzi[3],
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    float scl,
    double op[]
);
static int interpolate_3d_rc(
    const int nxyzi[3],
    const float dkxyzi[3],
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    float scl,
    double op[]
);
static int interpolate_cyl_to_rect_cc(
    int nri,
    float dkri,
    int nzi,
    float dkzi,
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    double op[]
);
static int interpolate_cyl_to_rect_rc(
    int nri,
    float dkri,
    int nzi,
    float dkzi,
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    double op[]
);
static void pad_ramp_3d(
    int n0,
    int n1,
    int n2,
    int n0p,
    int n1p,
    int n2p,
    int w0,
    int w1,
    double* inout
);
static int parse_command_line(
    int argc, char* argv[], COREParameters* pcp, FrontEndParameters* pfp
);
static void print_parameters(
    const COREParameters* pcp, FrontEndParameters* pfp
);
static void print_usage(void);
static void report_error(const char* fmt, ...);
static int rescale(
    const long size[3], double cuth, double* measurement, double* p_scale
);
static int setup_ctf(const FrontEndParameters* pfp, int* p_im_stream);
static int setup_input(FrontEndParameters* pfp, int* p_im_stream);
static int setup_output(
    const FrontEndParameters* pfp, int in_im, int *p_im_stream
);
static int setup_plot(
    int iw, int it, int iz, int in_im, ResReporterData* prrd
);
static int setup_plot_file(
    const FrontEndParameters* pfp, ResReporterData* prrd
);
static void update_plot_file(
    long iiters, long iitert, double res, void* closure
);
static int write_result_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    const double* volume,
    double scale,
    float mmm[3]
);


int
main(int argc, char* argv[])
{
    COREParameters cp;
    FrontEndParameters fp;
    ResReporterData rrd = { 0, 0, 0, 0.0, 0, 0, 0, 0 };
    int in_im = -1, ctf_im = -1, out_im = -1;
    double* deconvoluted = 0;
    double* ctf_copy = 0;
    int status = 1;
    int iw;
    float delta[3];

    if (parse_command_line(argc, argv, &cp, &fp) != 0) {
	goto cleanup;
    }
    if (setup_input(&fp, &in_im) != 0 || setup_ctf(&fp, &ctf_im) != 0) {
	goto cleanup;
    }
    if (setup_output(&fp, in_im, &out_im) != 0 ||
	setup_plot_file(&fp, &rrd) != 0) {
	goto cleanup;
    }
    (void) printf("\n");
    print_parameters(&cp, &fp);
    (void) printf("\n");

    deconvoluted = (double*) malloc(
	2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded * (size_t) fp.zpadded *
	sizeof(double)
    );
    ctf_copy = (double*) malloc(
	2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded * (size_t) fp.zpadded *
	sizeof(double)
    );
    if (deconvoluted == 0 || ctf_copy == 0) {
	report_error(
	    "memory allocation for result failed: %s", strerror(errno)
	);
	goto cleanup;
    }

    IMRtDel(in_im, delta);
    if (delta[0] <= 0.0f || delta[1] <= 0.0f || delta[2] <= 0.0f) {
	report_error(
	    "one or more of the pixel spacings for the measurement are not positive"
	);
	goto cleanup;
    }
    if (fabs(delta[0] - delta[1]) > 1e-5 * delta[0]) {
	report_error("expected data with equal x and y sampling");
	goto cleanup;
    }

    iw = 0;
    while (fp.wavelength_indices[iw] != -1) {
	double* ctf = fetch_ctf(&fp, in_im, ctf_im, iw);
	double* cone_filter = 0;
	float mmm[3] = { FLT_MAX, -FLT_MAX, 0.0f };
	float first_wave_mmm[3];
	long size[3];
	int had_error = 0, it;
	int result;

	if (ctf == 0) {
	    goto cleanup;
	}

	if (cp.lampc > 0.0) {
	    cone_filter = (double*) malloc(
		2 * (size_t) (fp.xpadded / 2 + 1) * (size_t) fp.ypadded *
		(size_t) fp.zpadded * sizeof(double)
	    );
	    if (cone_filter == 0) {
		report_error(
		    "memory allocation for cone filter failed: %s",
		    strerror(errno)
		);
		free(ctf);
		goto cleanup;
	    }
	    result = getmcf(
		ctf,
		cone_filter,
		fp.xpadded,
		fp.ypadded,
		fp.zpadded,
		fp.na,
		fp.nimm,
		delta[0],
		delta[2],
		fp.nthread
	    );
	    if (result != 0) {
		report_error("cone filter computation failed");
		free(cone_filter);
		free(ctf);
		goto cleanup;
	    }
	}

	size[0] = fp.xpadded;
	size[1] = fp.ypadded;
	size[2] = fp.zpadded;

	if (! fp.absolute_background) {
	    double background_scale = fp.background[iw];
	    double background = DBL_MAX;

	    /*
	     * Set this to zero so the subraction in fetch_measurement() has
	     * no effect.
	     */
	    fp.background[iw] = 0.0;
	    it = 0;
	    while (1) {
		double* measurement;
		double vback;

		if (it >=
		    (fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1) {
		    break;
		}

		measurement = fetch_measurement(&fp, in_im, iw, it);
		if (measurement == 0) {
		    had_error = 1;
		    break;
		}
		result = core_estimate_background(
		    measurement,
		    size,
		    delta[0] / (double) delta[2],
		    cp.lamf,
		    cp.lamsmooth,
		    fp.nthread,
		    &vback
		);

		free(measurement);
		if (result != 0) {
		    had_error = 1;
		    break;
		}
		if (background > vback) {
		    background = vback;
		}
		++it;
	    }

	    if (had_error) {
		free(ctf);
		free(cone_filter);
		goto cleanup;
	    }

	    (void) printf("unscaled background iw %d : %g\n", iw, background);
	    fp.background[iw] = background * background_scale;
	}

	it = 0;
	while (1) {
	    double* measurement;
	    double scale;

	    if (it >= (fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1) {
		if (iw == 0) {
		    first_wave_mmm[0] = mmm[0];
		    first_wave_mmm[1] = mmm[1];
		    first_wave_mmm[2] = mmm[2];
		} else {
		    IMAlWavMM(out_im, iw, mmm[0], mmm[1]);
		}
		if (fp.wavelength_indices[iw + 1] == -1) {
		  char title[80];

		  IMWrHdr(
		      out_im,
		      title,
		      -1,
		      first_wave_mmm[0],
		      first_wave_mmm[1],
		      first_wave_mmm[2] /
		      ((fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1)
		    );
		}
		break;
	    }

	    measurement = fetch_measurement(&fp, in_im, iw, it);
	    if (measurement == 0) {
		had_error = 1;
		break;
	    }
	    if (rescale(size, fp.cuth, measurement, &scale) != 0) {
		report_error(
		    "unable to allocate temporary space for rescaling: %s",
		    strerror(errno)
		);
		had_error = 1;
		break;
	    }
	    (void) printf("scale iw %d it %d : %g\n", iw, it, 1.0 / scale);

	    (void) memcpy(
		ctf_copy,
		ctf,
		2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded *
		    (size_t) fp.zpadded * sizeof(double)
	    );
	    result = setup_plot(iw, it, -1, in_im, &rrd);
	    if (result != 0) {
		had_error = 1;
		goto cleanup;
	    }
	    result = core(
		measurement,
		ctf_copy,
		cone_filter,
		deconvoluted,
		size,
		cp.lamratio,
		delta[0] / (double) delta[2],
		cp.lamf,
		cp.lampos,
		cp.lampc,
		cp.lamsmooth,
		cp.alpha,
		cp.iteration_count,
		20 * cp.iteration_count,
		cp.tol,
		cp.omega,
		fp.nthread,
		&update_plot_file,
		&rrd
	    );
	    free(measurement);
	    if (result != 0) {
		report_error("deconvolution computation failed");
		had_error = 1;
		break;
	    }
	    result = complete_plot(iw, it, -1, in_im, &rrd);
	    if (result != 0) {
		had_error = 1;
		goto cleanup;
	    }

	    result = write_result_and_collect_stats(
		&fp, out_im, iw, it, deconvoluted, scale, mmm
	    );
	    if (result != 0) {
		had_error = 1;
		break;
	    }

	    ++it;
	}

	free(ctf);
	free(cone_filter);
	if (had_error) {
	    goto cleanup;
	}
	++iw;
    }

    status = 0;

 cleanup:
    if (ctf_copy != 0) {
	free(ctf_copy);
    }
    if (deconvoluted != 0) {
	free(deconvoluted);
    }
    if (complete_plot_file(&rrd) != 0) {
	if (status == 0) {
	    status = 1;
	}
    }
    if (out_im != -1) {
	IMClose(out_im);
    }
    if (ctf_im != -1) {
	IMClose(ctf_im);
    }
    if (in_im != -1) {
	IMClose(in_im);
    }
    return status;
}


static
int
complete_plot(int iw, int it, int iz, int in_im, ResReporterData* prrd)
{
    char name[80];
    long i, j;

    if (prrd->fo == 0) {
	return 0;
    }
    if (prrd->status != 0) {
	report_error(
	    "no data for termination criteria plot: memory allocation error\n"
	);
	return 1;
    }

    if (iw >= 0) {
	int offset = 0;
	float wavelengths[IW_MAX_WAVE];
	int nw;

	offset += sprintf(name + offset, "Wave %d", iw);
	IMRtWav(in_im, &nw, wavelengths);
	if (iw < nw &&
	    wavelengths[iw] >= 200.0f && wavelengths[iw] < 2000.0f) {
	    offset += sprintf(name + offset, " (%.0f nm)", wavelengths[iw]);
	}
	if (it >= 0) {
	    offset += sprintf(name + offset, " Time %d", it);
	    if (iz >= 0) {
		offset += sprintf(name + offset, " Z %d", iz);
	    }
	}
    } else {
	(void) strcpy(name, "Full Dataset");
    }

    (void) fprintf(
	prrd->fo, "TTermination Criteria For %s: Successful Iterations\n", name
    );
    (void) fprintf(
	prrd->fo, "Xiteration\nYtermination criteria\nLsuccessful iterations\n"
    );
    for (i = 0, j = 0; i < prrd->ncount; ++i) {
	if (prrd->values[i].iteration > 0) {
	    continue;
	}
	(void) fprintf(prrd->fo, "%ld %.6f\n", j, prrd->values[i].res);
	++j;
    }
    (void) fprintf(
	prrd->fo, "TTermination Criteria For %s: All Iterations\n", name
    );
    (void) fprintf(
	prrd->fo, "Xiteration\nYtermination criteria\nLall iterations\n"
    );
    for (i = 0; i < prrd->ncount; ++i) {
	(void) fprintf(prrd->fo, "%ld %.6f\n", i, prrd->values[i].res);
    }

    if (ferror(prrd->fo) != 0) {
	report_error(
	    "an error occurred while writing the termination criteria plot, %s",
	    prrd->filename
	);
	return 1;
    }
    return 0;
}


static
int
complete_plot_file(ResReporterData* prrd)
{
    int status = 0;

    if (prrd->values != 0) {
	free(prrd->values);
    }
    if (prrd->fo != 0) {
	if (fclose(prrd->fo) != 0) {
	    status = 1;
	    report_error(
		"an error occurred while closing the plot file, %s: %s",
		prrd->filename,
		strerror(errno)
	    );
	}
    }
    return status;
}


static
int
extract_doubles(const char* s, int max_count, double* values)
{
    int count = 0;

    while (s[0] != '\0' && count >= 0) {
	double value;
	int offset;

	if (sscanf(s, "%lf%n", &value, &offset) == 1) {
	    switch (s[offset]) {
	    case ':':
		++s;
		/* Fall through. */

	    case '\0':
		s += offset;
		if (count < max_count) {
		    values[count] = value;
		}
		++count;
		break;

	    default:
		count = -count - 1;
		break;
	    }
	} else {
	    count = -count - 1;
	}
    }

    return count;
}


static
int
extract_ints(const char* s, int max_count, int* values)
{
    int count = 0;

    while (s[0] != '\0' && count >= 0) {
	int value, offset;

	if (sscanf(s, "%d%n", &value, &offset) == 1) {
	    switch (s[offset]) {
	    case ':':
		++s;
		/* Fall through. */

	    case '\0':
		s += offset;
		if (count < max_count) {
		    values[count] = value;
		}
		++count;
		break;

	    default:
		count = -count - 1;
		break;
	    }
	} else {
	    count = -count - 1;
	}
    }

    return count;
}


static
double*
fetch_ctf(
    const FrontEndParameters* pfp, int in_im, int ctf_im, int iw
)
{
    int measurement_nw;
    float measurement_wavelengths[IW_MAX_WAVE], measurement_delta[3];
    float measurement_nm;
    int nxhlf = pfp->xpadded / 2 + 1;
    int nxyz[3], mxyz[3], nxyzst[3], nw, pixel_format;
    float delta[3], wavelengths[IW_MAX_WAVE], wave_ratio;
    int wavelength_slot, is_real;
    double* out;

    IMRtWav(in_im, &measurement_nw, measurement_wavelengths);
    IMRtDel(in_im, measurement_delta);
    measurement_nm = measurement_wavelengths[pfp->wavelength_indices[iw]];

    wavelength_slot = find_matching_wavelength_mrc(ctf_im, measurement_nm);
    IMRtSiz(ctf_im, nxyz, mxyz, nxyzst);
    IMRtWav(ctf_im, &nw, wavelengths);
    IMRtMode(ctf_im, &pixel_format);
    IMRtDel(ctf_im, delta);

    if (nxyz[0] < 1 || nxyz[1] < 1 || nxyz[2] < 1) {
	report_error(
	    "one or more number of samples for the transfer function are not positive"
	);
	return 0;
    }
    if (measurement_nm <= 0.0f || wavelengths[wavelength_slot] <= 0.0f) {
	wave_ratio = 1.0f;
	if (measurement_nm <= 0.0f && wavelengths[wavelength_slot] <= 0.0f) {
	    report_error("warning: non-positive wavelengths for input data and OTF; assuming wavelengths are the same");
	} else {
	    report_error(
		(measurement_nm <= 0.0f) ?
		"invalid wavelength for input data" :
		"invalid wavelength for OTF"
	    );
	    return 0;
	}
    } else {
	wave_ratio = wavelengths[wavelength_slot] / measurement_nm;
    }

    is_real = pixel_format != IW_COMPLEX && pixel_format != IW_COMPLEX_SHORT;

    out = (double*) malloc(
	2 * (size_t) nxhlf * (size_t) pfp->ypadded * (size_t) pfp->zpadded *
	sizeof(double)
    );
    if (out == 0) {
	report_error("memory allocation for CTF failed: %s", strerror(errno));
	return 0;
    }

    if (nxyz[2] == nw) {
	float* section;
	int nxyzo[3];
	float dkxyzo[3];
	int status;

	if (delta[0] <= 0.0f || delta[1] <= 0.0f) {
	    report_error(
		"the x or y pixel spacing for the CTF is not positive"
	    );
	    free(out);
	    return 0;
	}
	if (nxyz[0] > SIZE_MAX / (((is_real) ? 1 : 2) * sizeof(float)) ||
	    nxyz[1] > SIZE_MAX / (((is_real) ? 1 : 2) * (size_t) nxyz[0] *
				  sizeof(float))) {
	    report_error(
		"the size of on-file CTF exceeds the machine limits"
	    );
	    free(out);
	    return 0;
	}
	section = (float*) malloc(
	    ((is_real) ? 1 : 2) * nxyz[0] * (size_t) nxyz[1] * sizeof(float)
	);
	if (section == 0) {
	    report_error(
		"memory allocation for CTF buffer failed: %s", strerror(errno)
	    );
	    free(out);
	    return 0;
	}
	if (IMPosnZWT(ctf_im, 0, wavelength_slot, 0) != 0) {
	    report_error(
		"unable to position CTF file at the wavelength index, %d",
		wavelength_slot
	    );
	    free(section);
	    free(out);
	    return 0;
	}
	IMRdSec(ctf_im, section);

	nxyzo[0] = nxhlf;
	nxyzo[1] = pfp->ypadded;
	nxyzo[2] = pfp->zpadded;
	dkxyzo[0] = 1.0f / (pfp->xpadded * measurement_delta[0]);
	dkxyzo[1] = 1.0f / (pfp->ypadded * measurement_delta[1]);
	dkxyzo[2] = 1.0f / (pfp->zpadded * measurement_delta[2]);
	if (is_real) {
	    status = interpolate_cyl_to_rect_rc(
		nxyz[1],
		delta[1],
		nxyz[0],
		delta[0],
		section,
		nxyzo,
		dkxyzo,
		pfp->xpadded % 2 == 0,
		wave_ratio,
		out
	    );
	} else {
	    status = interpolate_cyl_to_rect_cc(
		nxyz[1],
		delta[1],
		nxyz[0],
		delta[0],
		section,
		nxyzo,
		dkxyzo,
		pfp->xpadded % 2 == 0,
		wave_ratio,
		out
	    );
	}
	free(section);
	if (status != 0) {
	    report_error(
		"memory allocation failed in CTF interpolation: %s",
		strerror(errno)
	    );
	    free(out);
	    return 0;
	}
    } else {
	int ntin, intlv;
	float* otfin;
	int iz;
	int nxyzo[3];
	float dkxyzo[3];
	int status;

	if (delta[0] <= 0.0f || delta[1] <= 0.0f || delta[2] <= 0.0f) {
	    report_error(
		"the x, y, or z pixel spacing for the CTF is not positive"
	    );
	    free(out);
	    return 0;
	}
	IMRtZWT(ctf_im, &nxyz[2], &nw, &ntin, &intlv);
	if (nxyz[0] > SIZE_MAX / (((is_real) ? 1 : 2) * sizeof(float)) ||
	    nxyz[1] > SIZE_MAX / (((is_real) ? 1 : 2) * (size_t) nxyz[0] *
				  sizeof(float)) ||
	    nxyz[2] > SIZE_MAX / (((is_real) ? 1 : 2) * (size_t) nxyz[0] *
				  (size_t) nxyz[1] * sizeof(float))) {
	    report_error(
		"the size of the on-file CTF exceeds the machine limits"
	    );
	    free(out);
	    return 0;
	}
	otfin = (float*) malloc(
	    ((is_real) ? 1 : 2) * nxyz[0] * (size_t) nxyz[1] *
	    (size_t) nxyz[2] * sizeof(float)
	);
	if (otfin == 0) {
	    report_error(
		"memory allocation for CTF buffer failed: %s", strerror(errno)
	    );
	    free(out);
	    return 0;
	}

	for (iz = 0; iz < nxyz[2]; ++iz) {
	    if (IMPosnZWT(ctf_im, iz, wavelength_slot, 0) != 0) {
		report_error(
		    "unable to position CTF file at wavelength index, %d, and z index, %d",
		    wavelength_slot,
		    iz
		);
		free(otfin);
		free(out);
		return 0;
	    }
	    IMRdSec(
		ctf_im,
		otfin + iz * (size_t) nxyz[0] * (size_t) nxyz[1] *
		    (size_t) ((is_real) ? 1 : 2)
	    );
	}

	delta[0] = 1.0f / (2 * (nxyz[0] - 1) * delta[0]);
	delta[1] = 1.0f / (nxyz[1] * delta[1]);
	delta[2] = 1.0f / (nxyz[2] * delta[2]);
	nxyzo[0] = nxhlf;
	nxyzo[1] = pfp->ypadded;
	nxyzo[2] = pfp->zpadded;
	dkxyzo[0] = 1.0f / (pfp->xpadded * measurement_delta[0]);
	dkxyzo[1] = 1.0f / (pfp->ypadded * measurement_delta[1]);
	dkxyzo[2] = 1.0f / (pfp->zpadded * measurement_delta[2]);
	if (is_real) {
	    status = interpolate_3d_rc(
		nxyz,
		delta,
		otfin,
		nxyzo,
		dkxyzo,
	        pfp->xpadded % 2 == 0,
		wave_ratio,
		1.0f / ((pfp->normalize_otf) ? otfin[0] : 1.0f),
		out
	    );
	} else {
	    status = interpolate_3d_cc(
		nxyz,
		delta,
		otfin,
		nxyzo,
		dkxyzo,
		pfp->xpadded % 2 == 0,
		wave_ratio,
		1.0f / ((pfp->normalize_otf) ? otfin[0] : 1.0f),
		out
	    );
	}
	free(otfin);
	if (status != 0) {
	    report_error("memory allocation failed in CTF interpolation");
	    free(out);
	    return 0;
	}
    }

    return out;
}


static
double*
fetch_measurement(
    const FrontEndParameters* pfp, int in_im, int iw, int it
)
{
    int nxp = 2 * (pfp->xpadded / 2 + 1);
    double* volume = (double*) malloc(
	nxp * (size_t) pfp->ypadded * (size_t) pfp->zpadded * sizeof(double)
    );
    float* section = (float*) malloc(
	(size_t) (pfp->xbounds[1] - pfp->xbounds[0] + 1) *
	(size_t) (pfp->ybounds[1] - pfp->ybounds[0] + 1) * sizeof(float)
    );
    double* out;
    int iz;

    if (volume == 0 || section == 0) {
	report_error(
	    "memory allocation failed for measurement or measurement buffer: %s\n",
	    strerror(errno)
	);
	if (section != 0) {
	    free(section);
	}
	if (volume != 0) {
	    free(volume);
	}
	return 0;
    }

    out = volume;
    for (iz = 0; iz < pfp->zbounds[1] - pfp->zbounds[0] + 1; ++iz) {
	const float* cursor = section;
	int iy;

	if (IMPosnZWT(in_im, pfp->zbounds[0] + iz, pfp->wavelength_indices[iw], pfp->tbounds[0] + it * pfp->tbounds[2]) != 0) {
	    report_error(
		"unable to position measurement file at z, wavelength, and time indices, %d, %d, %d",
		pfp->zbounds[0] + iz,
		pfp->wavelength_indices[iw],
		pfp->tbounds[0] + it * pfp->tbounds[2]
	    );
	    free(section);
	    free(volume);
	    return 0;
	}
	IMRdPas(
	    in_im,
	    section,
	    pfp->xbounds[1] - pfp->xbounds[0] + 1,
	    pfp->ybounds[1] - pfp->ybounds[0] + 1,
	    pfp->xbounds[0],
	    pfp->xbounds[1],
	    pfp->ybounds[0],
	    pfp->ybounds[1]
	);

	for (iy = 0; iy < pfp->ybounds[1] - pfp->ybounds[0] + 1; ++iy) {
	    int ix;

	    for (ix = 0; ix < pfp->xbounds[1] - pfp->xbounds[0] + 1; ++ix) {
		out[ix] = cursor[ix] - pfp->background[iw];
	    }
	    for (ix = pfp->xbounds[1] - pfp->xbounds[0] + 1; ix < nxp; ++ix) {
		out[ix] = 0.0f;
	    }
	    cursor += pfp->xbounds[1] - pfp->xbounds[0] + 1;
	    out += nxp;
	}
	for (iy = pfp->ybounds[1] - pfp->ybounds[0] + 1;
	     iy < pfp->ypadded;
	     ++iy) {
	    int ix;

	    for (ix = 0; ix < nxp; ++ix) {
		out[ix] = 0.0f;
	    }
	    out += nxp;
	}
    }

    free(section);

    pad_ramp_3d(
	pfp->xpadded,
	pfp->ypadded,
	pfp->zbounds[1] - pfp->zbounds[0] + 1,
	0,
	0,
	pfp->zpadded - (pfp->zbounds[1] - pfp->zbounds[0] + 1),
	nxp,
	pfp->ypadded,
	volume
    );

    return volume;
}


static
int
find_free_stream(void)
{
    int i = 1;

    while (i < INT_MAX && (IMUnit(i) != -1 || IMRtWID(i) != 0)) {
	++i;
    }

    return (i != INT_MAX) ? i : -1;
}


static
int
find_matching_wavelength_mrc(int mrc_stream, float wavelength_nm)
{
    int nw;
    float wavelengths[IW_MAX_WAVE];
    float min_diff;
    int iw_closest;
    int iw;

    IMRtWav(mrc_stream, &nw, wavelengths);
    if (nw < 1) {
	return -1;
    }
    iw_closest = 0;
    min_diff = fabs(wavelength_nm - wavelengths[0]);
    for (iw = 1; iw < nw; ++iw) {
	float diff = fabs(wavelength_nm - wavelengths[iw]);

	if (min_diff > diff) {
	    iw_closest = iw;
	    min_diff = diff;
	}
    }

    return iw_closest;
}


static
void
force_conj_sym_c(int nxh, int ny, int nz, int even_x, double* io)
{
    int iyo, izo;

    io[1] = 0.0f;
    if (even_x) {
	io[2 * (size_t) nxh - 1] = 0.0f;
    }
    if (ny % 2 == 0) {
	io[ny * (size_t) nxh + 1] = 0.0f;
	if (even_x) {
	    io[(ny + 2) * (size_t) nxh - 1] = 0.0f;
	}
    }
    if (even_x) {
	for (iyo = ny / 2 + 1; iyo < ny; ++iyo) {
	    io[iyo * (size_t) (nxh + nxh)] =
		io[(ny - iyo) * (size_t) (nxh + nxh)];
	    io[iyo * (size_t) (nxh + nxh) + 1] =
		-io[(ny - iyo) * (size_t) (nxh + nxh) + 1];
	    io[(iyo + 1) * (size_t) (nxh + nxh) - 2] =
		io[(ny - iyo + 1) * (size_t) (nxh + nxh) - 2];
	    io[(iyo + 1) * (size_t) (nxh + nxh) - 1] =
		-io[(ny - iyo + 1) * (size_t) (nxh + nxh) - 1];
	}
    } else {
	for (iyo = ny / 2 + 1; iyo < ny; ++iyo) {
	    io[iyo * (size_t) (nxh + nxh)] =
		io[(ny - iyo) * (size_t) (nxh + nxh)];
	    io[iyo * (size_t) (nxh + nxh) + 1] =
		-io[(ny - iyo) * (size_t) (nxh + nxh) + 1];
	}
    }
    if (nz % 2 == 0) {
	io[nz * (size_t) nxh * (size_t) ny + 1] = 0.0f;
	if (even_x) {
	    io[(nz * (size_t) ny + 2) * (size_t) nxh - 1] = 0.0f;
	}
	if (ny % 2 == 0) {
	    io[(nz + 1) * (size_t) nxh * (size_t) ny + 1] = 0.0f;
	    if (even_x) {
		io[((nz + 1) * (size_t) ny + 2) * (size_t) nxh - 1] = 0.0f;
	    }
	}
    }
    if (even_x) {
	for (izo = nz / 2 + 1; izo < nz; ++izo) {
	    io[izo * (size_t) (nxh + nxh) * (size_t) ny] = 
		io[(nz - izo) * (size_t) (nxh + nxh) * (size_t) ny];
	    io[izo * (size_t) (nxh + nxh) * (size_t) ny + 1] = 
		-io[(nz - izo) * (size_t) (nxh + nxh) * (size_t) ny + 1];
	    io[(izo * (size_t) ny + 1) * (size_t) (nxh + nxh) - 2] = io[
		((nz - izo) * (size_t) ny + 1) * (size_t) (nxh + nxh) - 2
	    ];
	    io[(izo * (size_t) ny + 1) * (size_t) (nxh + nxh) - 1] = -io[
		((nz - izo) * (size_t) ny + 1) * (size_t) (nxh + nxh) - 1
	    ];
	    for (iyo = 1; iyo < ny; ++iyo) {
		io[(izo * (size_t) ny + iyo) * (size_t) (nxh + nxh)] = io[
		    ((nz - izo + 1) * (size_t) ny - iyo) * (size_t) (nxh + nxh)
		];
		io[(izo * (size_t) ny + iyo) * (size_t) (nxh + nxh) + 1] = -io[
		    ((nz - izo + 1) * (size_t) ny - iyo) *
		    (size_t) (nxh + nxh) + 1
		];
		io[
		    (izo * (size_t) ny + iyo + 1) * (size_t) (nxh + nxh) - 2
		] = io[
		    ((nz - izo + 1) * (size_t) ny - iyo + 1) *
		    (size_t) (nxh + nxh) - 2
		];
		io[
		    (izo * (size_t) ny + iyo + 1) * (size_t) (nxh + nxh) - 1
		] = -io[
		    ((nz - izo + 1) * (size_t) ny - iyo + 1) *
		    (size_t) (nxh + nxh) - 1
		];
	    }
	}
    } else {
	for (izo = nz / 2 + 1; izo < nz; ++izo) {
	    io[izo * (size_t) (nxh + nxh) * (size_t) ny] = 
		io[(nz - izo) * (size_t) (nxh + nxh) * (size_t) ny];
	    io[izo * (size_t) (nxh + nxh) * (size_t) ny + 1] = 
		-io[(nz - izo) * (size_t) (nxh + nxh) * (size_t) ny + 1];
	    for (iyo = 1; iyo < ny; ++iyo) {
		io[(izo * (size_t) ny + iyo) * (size_t) (nxh + nxh)] = io[
		    ((nz - izo + 1) * (size_t) ny - iyo) * (size_t) (nxh + nxh)
		];
		io[(izo * (size_t) ny + iyo) * (size_t) (nxh + nxh) + 1] = -io[
		    ((nz - izo + 1) * (size_t) ny - iyo) *
		    (size_t) (nxh + nxh) + 1
		];
	    }
	}
    }
}


static
int
get_default_thread_count(void)
{
    /*
     * For integration with the Sun's Grid Engine software, allow the
     * value of the NSLOTS environment variable to control the maximum number
     * of threads to use; otherwise use a number of threads equal to the
     * number of processors reported by the operating system.
     */
    char* envvar = getenv("NSLOTS");
    int default_value;

    if (envvar != 0) {
	char* endp;
	long slots;

	slots = strtol(envvar, &endp, 10);

	if (*endp == '\0' && *envvar != '\0') {
	    if (slots >= INT_MIN) {
		if (slots <= INT_MAX) {
		    default_value = slots;
		} else {
		    default_value = INT_MAX;
		}
	    } else {
		default_value = INT_MIN;
	    }
	} else {
	    /*
	     * NSLOTS did not contain a decimal integer value.  Instead, use
	     * the number of processors reported by the operating system.
	     */
	    default_value = IVEThreadCountProcessors();
	}
    } else {
	default_value = IVEThreadCountProcessors();
    }

    return (default_value < 1) ? 1 : default_value;
}


static
int
interpolate_3d_cc(
    const int nxyzi[3],
    const float dkxyzi[3],
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    float scl,
    double op[]
)
{
    Interp1DSourceWeights* xweights =
	setup_1d_interp(nxyzi[0], dkxyzi[0] * mag, nxyzo[0], dkxyzo[0]);
    Interp1DSourceWeights* yweights = setup_1d_interp_split(
	nxyzi[1], dkxyzi[1] * mag, nxyzo[1], dkxyzo[1]
    );
    Interp1DSourceWeights* zweights = setup_1d_interp_split(
	nxyzi[2], dkxyzi[2] * mag, nxyzo[2], dkxyzo[2]
    );
    int i;
    int izo;

    if (xweights == 0 || yweights == 0 || zweights == 0) {
	if (xweights != 0) {
	    free(xweights);
	}
	if (yweights != 0) {
	    free(yweights);
	}
	if (zweights != 0) {
	    free(zweights);
	}
	return 1;
    }
    for (i = 0; i < nxyzo[0]; ++i) {
	xweights[i].i1 *= 2;
	xweights[i].i2 *= 2;
    }

    for (izo = 0; izo < nxyzo[2]; ++izo) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) izo;
	int iyo;

	for (iyo = 0; iyo < nxyzo[1]; ++iyo) {
	    const float* ip11 = ip + 2 * (size_t) nxyzi[0] *
		(yweights[iyo].i1 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i1);
	    const float* ip21 = ip + 2 * (size_t) nxyzi[0] *
		(yweights[iyo].i2 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i1);
	    const float* ip12 = ip + 2 * (size_t) nxyzi[0] *
		(yweights[iyo].i1 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i2);
	    const float* ip22 = ip + 2 * (size_t) nxyzi[0] *
		(yweights[iyo].i2 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i2);
	    float w11 = yweights[iyo].w1 * zweights[izo].w1;
	    float w21 = yweights[iyo].w2 * zweights[izo].w1;
	    float w12 = yweights[iyo].w1 * zweights[izo].w2;
	    float w22 = yweights[iyo].w2 * zweights[izo].w2;
	    int ixo;

	    for (ixo = 0; ixo < nxyzo[0]; ++ixo) {
		oc[ixo + ixo] = scl * (
		    xweights[ixo].w1 * (
			w11 * ip11[xweights[ixo].i1] +
			w21 * ip21[xweights[ixo].i1] +
			w12 * ip12[xweights[ixo].i1] +
			w22 * ip22[xweights[ixo].i1]
		    ) + xweights[ixo].w2 * (
			w11 * ip11[xweights[ixo].i2] +
			w21 * ip21[xweights[ixo].i2] +
			w12 * ip12[xweights[ixo].i2] +
			w22 * ip22[xweights[ixo].i2]
		    )
		);
		oc[ixo + ixo + 1] = scl * (
		    xweights[ixo].w1 * (
			w11 * ip11[xweights[ixo].i1 + 1] +
			w21 * ip21[xweights[ixo].i1 + 1] +
			w12 * ip12[xweights[ixo].i1 + 1] +
			w22 * ip22[xweights[ixo].i1 + 1]
		    ) + xweights[ixo].w2 * (
			w11 * ip11[xweights[ixo].i2 + 1] +
			w21 * ip21[xweights[ixo].i2 + 1] +
			w12 * ip12[xweights[ixo].i2 + 1] +
			w22 * ip22[xweights[ixo].i2 + 1]
		    )
		);
	    }
	    oc += nxyzo[0] + nxyzo[0];
	}
    }

    force_conj_sym_c(nxyzo[0], nxyzo[1], nxyzo[2], even_x, op);
    free(xweights);
    free(yweights);
    free(zweights);

    return 0;
}


static
int
interpolate_3d_rc(
    const int nxyzi[3],
    const float dkxyzi[3],
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    float scl,
    double op[]
)
{
    Interp1DSourceWeights* xweights =
	setup_1d_interp(nxyzi[0], dkxyzi[0] * mag, nxyzo[0], dkxyzo[0]);
    Interp1DSourceWeights* yweights = setup_1d_interp_split(
	nxyzi[1], dkxyzi[1] * mag, nxyzo[1], dkxyzo[1]
    );
    Interp1DSourceWeights* zweights = setup_1d_interp_split(
	nxyzi[2], dkxyzi[2] * mag, nxyzo[2], dkxyzo[2]
    );
    int izo;

    if (xweights == 0 || yweights == 0 || zweights == 0) {
	if (xweights != 0) {
	    free(xweights);
	}
	if (yweights != 0) {
	    free(yweights);
	}
	if (zweights != 0) {
	    free(zweights);
	}
	return 1;
    }

    for (izo = 0; izo < nxyzo[2]; ++izo) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) izo;
	int iyo;

	for (iyo = 0; iyo < nxyzo[1]; ++iyo) {
	    const float* ip11 = ip + (size_t) nxyzi[0] *
		(yweights[iyo].i1 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i1);
	    const float* ip21 = ip + (size_t) nxyzi[0] *
		(yweights[iyo].i2 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i1);
	    const float* ip12 = ip + (size_t) nxyzi[0] *
		(yweights[iyo].i1 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i2);
	    const float* ip22 = ip + (size_t) nxyzi[0] *
		(yweights[iyo].i2 + (size_t) nxyzi[1] *
		 (size_t) zweights[izo].i2);
	    float w11 = yweights[iyo].w1 * zweights[izo].w1;
	    float w21 = yweights[iyo].w2 * zweights[izo].w1;
	    float w12 = yweights[iyo].w1 * zweights[izo].w2;
	    float w22 = yweights[iyo].w2 * zweights[izo].w2;
	    int ixo;

	    for (ixo = 0; ixo < nxyzo[0]; ++ixo) {
		oc[ixo + ixo] = scl * (
		    xweights[ixo].w1 * (
			w11 * ip11[xweights[ixo].i1] +
			w21 * ip21[xweights[ixo].i1] +
			w12 * ip12[xweights[ixo].i1] +
			w22 * ip22[xweights[ixo].i1]
		    ) + xweights[ixo].w2 * (
			w11 * ip11[xweights[ixo].i2] +
			w21 * ip21[xweights[ixo].i2] +
			w12 * ip12[xweights[ixo].i2] +
			w22 * ip22[xweights[ixo].i2]
		    )
		);
		oc[ixo + ixo + 1] = 0.0;
	    }
	    oc += nxyzo[0] + nxyzo[0];
	}
    }

    force_conj_sym_c(nxyzo[0], nxyzo[1], nxyzo[2], even_x, op);
    free(xweights);
    free(yweights);
    free(zweights);

    return 0;
}


static
int
interpolate_cyl_to_rect_cc(
    int nri,
    float dkri,
    int nzi,
    float dkzi,
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    double op[]
)
{
    int nzhlf = nxyzo[2] / 2 + 1;
    float krmax = nri * dkri * mag;
    float krmax2 = krmax * krmax;
    Interp1DSourceWeights* zweights =
	setup_1d_interp_split(nzi, dkzi * mag, nxyzo[2], dkxyzo[2]);
    int iz;

    if (zweights == 0) {
	return 1;
    }

    for (iz = 0; iz < nzhlf; ++iz) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) iz;
	const float* ip1 = ip + 2 * (size_t) zweights[iz].i1;
	const float* ip2 = ip + 2 * (size_t) zweights[iz].i2;
	int iy;

	for (iy = 0; iy < nxyzo[1]; ++iy) {
	    float ky = ((iy <= nxyzo[1] / 2) ? iy : iy - nxyzo[1]) * dkxyzo[1];

	    if (fabs(ky) <= krmax) {
		float ky2 = ky * ky;
		float xhigh = sqrt(krmax2 - ky2) / dkxyzo[0];
		int ixhigh = (xhigh <= INT_MAX) ? (int) floor(xhigh) : INT_MAX;
		int ix, ix2;

		if (ixhigh > nxyzo[0] - 1) {
		    ixhigh = nxyzo[0] - 1;
		}
		for (ix = 0, ix2 = 0; ix <= ixhigh; ++ix, ix2 += 2) {
		    float r = sqrt(ky2 + ix * dkxyzo[0] * ix * dkxyzo[0]) /
			(dkri * mag);
		    int il = floor(r);
		    float wl = 1.0f - (r - il);

		    if (il < nri - 1) {
			float wh = r - il;

			il *= nzi + nzi;
			oc[ix2] =
			    wl * (zweights[iz].w1 * ip1[il] +
				  zweights[iz].w2 * ip2[il]) +
			    wh * (zweights[iz].w1 * ip1[il + nzi + nzi] +
				  zweights[iz].w2 * ip2[il + nzi + nzi]);
			oc[ix2 + 1] =
			    wl * (zweights[iz].w1 * ip1[il + 1] +
				  zweights[iz].w2 * ip2[il + 1]) +
			    wh * (zweights[iz].w1 * ip1[il + nzi + nzi + 1] +
				  zweights[iz].w2 * ip2[il + nzi + nzi + 1]);
		    } else {
			il *= nzi + nzi;
			oc[ix2] =
			    wl * (zweights[iz].w1 * ip1[il] +
				  zweights[iz].w2 * ip2[il]);
			oc[ix2 + 1] =
			    wl * (zweights[iz].w1 * ip1[il + 1] +
				  zweights[iz].w2 * ip2[il + 1]);
		    }
		}
		for (ix = ixhigh + ixhigh + 2;
		     ix < nxyzo[0] + nxyzo[0];
		     ++ix) {
		    oc[ix] = 0.0f;
		}
		oc += nxyzo[0] + nxyzo[0];
	    } else {
		int ix;

		for (ix = 0; ix < nxyzo[0] + nxyzo[0]; ++ix) {
		    oc[ix] = 0.0f;
		}
		oc += nxyzo[0] + nxyzo[0];
	    }
	}
    }
    /*
     * Enforce the consequences of conjugate symmetric and cylindrical
     * symmetry (f(x,y,z) = f*(-x,-y,-z) but f*(-x,-y,-z) = f*(-x,y,-z) =
     * f*(x,-y,-z) = f*(x,y,-z) by cylindrical symmetry so f(x,y,z) =
     * f*(x,y,-z).
     */
    for (iz = nzhlf; iz < nxyzo[2]; ++iz) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) iz;
	const double* ic = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) (nxyzo[2] - iz);
	int iy;

	for (iy = 0; iy < nxyzo[1]; ++iy) {
	    int ix;

	    for (ix = 0; ix < nxyzo[0]; ++ix) {
		oc[ix + ix] = ic[ix + ix];
		oc[ix + ix + 1] = -ic[ix + ix + 1];
	    }
	    oc += nxyzo[0] + nxyzo[0];
	    ic += nxyzo[0] + nxyzo[0];
	}
    }

    force_conj_sym_c(nxyzo[0], nxyzo[1], nxyzo[2], even_x, op);
    free(zweights);
    return 0;
}


static
int
interpolate_cyl_to_rect_rc(
    int nri,
    float dkri,
    int nzi,
    float dkzi,
    const float ip[],
    const int nxyzo[3],
    const float dkxyzo[3],
    int even_x,
    float mag,
    double op[]
)
{
    int nzhlf = nxyzo[2] / 2 + 1;
    float krmax = nri * dkri * mag;
    float krmax2 = krmax * krmax;
    Interp1DSourceWeights* zweights =
	setup_1d_interp_split(nzi, dkzi * mag, nxyzo[2], dkxyzo[2]);
    int iz;

    if (zweights == 0) {
	return 1;
    }

    for (iz = 0; iz < nzhlf; ++iz) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) iz;
	const float* ip1 = ip + zweights[iz].i1;
	const float* ip2 = ip + zweights[iz].i2;
	int iy;

	for (iy = 0; iy < nxyzo[1]; ++iy) {
	    float ky = ((iy <= nxyzo[1] / 2) ? iy : iy - nxyzo[1]) * dkxyzo[1];

	    if (fabs(ky) <= krmax) {
		float ky2 = ky * ky;
		float xhigh = sqrt(krmax2 - ky2) / dkxyzo[0];
		int ixhigh = (xhigh <= INT_MAX) ? (int) floor(xhigh) : INT_MAX;
		int ix, ix2;

		if (ixhigh > nxyzo[0] - 1) {
		    ixhigh = nxyzo[0] - 1;
		}
		for (ix = 0, ix2 = 0; ix <= ixhigh; ++ix, ix2 += 2) {
		    float r = sqrt(ky2 + ix * dkxyzo[0] * ix * dkxyzo[0]) /
			(dkri * mag);
		    int il = floor(r);
		    float wl = 1.0f - (r - il);

		    if (il < nri - 1) {
			float wh = r - il;

			il *= nzi;
			oc[ix2] =
			    wl * (zweights[iz].w1 * ip1[il] +
				  zweights[iz].w2 * ip2[il]) +
			    wh * (zweights[iz].w1 * ip1[il + nzi] +
				  zweights[iz].w2 * ip2[il + nzi]);
			oc[ix2 + 1] = 0.0;
		    } else {
			il *= nzi;
			oc[ix2] = 
			    wl * (zweights[iz].w1 * ip1[il] +
				  zweights[iz].w2 * ip2[il]);
			oc[ix2 + 1] = 0.0;
		    }
		}
		for (ix = ixhigh + ixhigh + 2;
		     ix < nxyzo[0] + nxyzo[0];
		     ++ix) {
		    oc[ix] = 0.0f;
		}
		oc += nxyzo[0] + nxyzo[0];
	    } else {
		int ix;

		for (ix = 0; ix < nxyzo[0] + nxyzo[0]; ++ix) {
		    oc[ix] = 0.0f;
		}
		oc += nxyzo[0] + nxyzo[0];
	    }
	}
    }
    /*
     * Enforce the consequences of conjugate symmetric and cylindrical
     * symmetry (f(x,y,z) = f*(-x,-y,-z) but f*(-x,-y,-z) = f*(-x,y,-z) =
     * f*(x,-y,-z) = f*(x,y,-z) by cylindrical symmetry so f(x,y,z) =
     * f*(x,y,-z).
     */
    for (iz = nzhlf; iz < nxyzo[2]; ++iz) {
	double* oc = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) iz;
	const double* ic = op + 2 * (size_t) nxyzo[0] * (size_t) nxyzo[1] *
	    (size_t) (nxyzo[2] - iz);
	int iy;

	for (iy = 0; iy < nxyzo[1]; ++iy) {
	    int ix;

	    for (ix = 0; ix < nxyzo[0]; ++ix) {
		oc[ix + ix] = ic[ix + ix];
		oc[ix + ix + 1] = 0.0f;
	    }
	    oc += nxyzo[0] + nxyzo[0];
	    ic += nxyzo[0] + nxyzo[0];
	}
    }

    force_conj_sym_c(nxyzo[0], nxyzo[1], nxyzo[2], even_x, op);
    free(zweights);
    return 0;
}


static
void
pad_ramp_3d(
    int n0,
    int n1,
    int n2,
    int n0p,
    int n1p,
    int n2p,
    int w0,
    int w1,
    double* inout
)
{
    double norm0 = 1.0 / (n0p + 1);
    double* out = inout;
    int i2;

    for (i2 = 0; i2 < n2; ++i2) {
	double* section = out;
	const double* line = out;
	int i1;

	for (i1 = 0; i1 < n1; ++i1) {
	    int i0;

	    out += n0;
	    for (i0 = n0; i0 < n0 + n0p; ++i0) {
		(*out++) = line[n0 - 1] +
		    (i0 - n0 + 1) * norm0 * (line[0] - line[n0 - 1]);
	    }
	    out += w0 - n0 - n0p;
	    line += w0;
	}
	for (i1 = n1; i1 < n1 + n1p; ++i1) {
	    const double* fl = section;
	    const double* ll = section + w0 * (size_t) (n1 - 1);
	    double weight1 = (double) (i1 - n1 + 1) / (n1p + 1);
	    int i0;

	    for (i0 = 0; i0 < n0; ++i0) {
		(*out++) = (*ll) + weight1 * ((*fl++) - (*ll));
		++ll;
	    }
	    for (i0 = n0; i0 < n0 + n0p; ++i0) {
		(*out++) = line[n0 - 1] +
		    (i0 - n0 + 1) * norm0 * (line[0] - line[n0 - 1]);
	    }
	    out += w0 - n0 - n0p;
	    line += w0;
	}
	out += w0 * (size_t) (w1 - n1 - n1p);
    }

    for (i2 = n2; i2 < n2 + n2p; ++i2) {
	const double* fs = inout;
	const double* ls = inout + w0 * (size_t) w1 * (size_t) (n2 - 1);
	double* section = out;
	const double* line = out;
	double weight3 = (double) (i2 - n2 + 1) / (n2p + 1);
	int i1;

	for (i1 = 0; i1 < n1; ++i1) {
	    int i0;

	    for (i0 = 0; i0 < n0; ++i0) {
		(*out++) = (*ls) + weight3 * ((*fs++) - (*ls));
		++ls;
	    }
	    for (i0 = n0; i0 < n0 + n0p; ++i0) {
		(*out++) = line[n0 - 1] +
		    (i0 - n0 + 1) * norm0 * (line[0] - line[n0 - 1]);
	    }
	    out += w0 - n0 - n0p;
	    fs += w0 - n0;
	    ls += w0 - n0;
	    line += w0;
	}

	for (i1 = n1; i1 < n1 + n1p; ++i1) {
	    const double* fl = section;
	    const double* ll = section + w0 * (size_t) (n1 - 1);
	    double weight1 = (double) (i1 - n1 + 1) / (n1p + 1);
	    int i0;

	    for (i0 = 0; i0 < n0; ++i0) {
		(*out++) = (*ll) + weight1 * ((*fl++) - (*ll));
		++ll;
	    }
	    for (i0 = n0; i0 < n0 + n0p; ++i0) {
		(*out++) = line[n0 - 1] +
		    (i0 - n0 + 1) * norm0 * (line[0] - line[n0 - 1]);
	    }
	    out += w0 - n0 - n0p;
	    line += w0;
	}
	out += w0 * (size_t) (w1 - n1 - n1p);
    }
}


static
int
parse_command_line(
    int argc, char* argv[], COREParameters* pcp, FrontEndParameters* pfp
)
{
    int iarg, had_error, background_set;

    pcp->lamratio = 0.5;
    pcp->lamf = 1.0;
    pcp->lampos = 100.0;
    pcp->lampc = 100.0;
    pcp->lamsmooth = 100.0;
    pcp->alpha = 1000.0;
    pcp->tol = 0.001;
    pcp->omega = 0.8;
    pcp->iteration_count = 50;

    pfp->measurement_file_name = 0;
    pfp->result_file_name = 0;
    pfp->ctf_file_name = 0;
    pfp->cuth = .001;
    pfp->na = 1.4;
    pfp->nimm = 1.512;
    pfp->xbounds[0] = 0;
    pfp->xbounds[1] = -1;
    pfp->ybounds[0] = 0;
    pfp->ybounds[1] = -1;
    pfp->zbounds[0] = 0;
    pfp->zbounds[1] = -1;
    pfp->tbounds[0] = 0;
    pfp->tbounds[1] = -1;
    pfp->tbounds[2] = 1;
    pfp->wavelength_indices[0] = -1;
    pfp->xpadded = -1;
    pfp->ypadded = -1;
    pfp->zpadded = -1;
    pfp->absolute_background = 0;
    pfp->normalize_otf = 1;
    pfp->nthread = get_default_thread_count();

    if (argc < 4) {
	print_usage();
	return 1;
    }

    pfp->measurement_file_name = argv[1];
    pfp->result_file_name = argv[2];
    pfp->ctf_file_name = argv[3];
    iarg = 4;
    had_error = 0;
    background_set = 0;
    while (iarg < argc) {
	int recognized = 1;

	if (argv[iarg][0] == '-') {
	    switch (argv[iarg][1]) {
	    case 'a':
		if (strcmp(argv[iarg] + 2, "bs_background") == 0) {
		    pfp->absolute_background = 1;
		} else if (strncmp(argv[iarg] + 2, "lpha=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->alpha) != 1) {
			report_error(
			    "expected one floating-point value for -alpha"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'c':
		if (strncmp(argv[iarg] + 2, "uth=", 4) == 0) {
		    if (extract_doubles(argv[iarg] + 6, 1, &pfp->cuth) != 1 ||
			pfp->cuth < 0.0 || pfp->cuth >= 1.0) {
			report_error(
			    "expected one floating-point value greater than or equal to zero and less than one for -cuth"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'l':
		if (strncmp(argv[iarg] + 2, "amratio=", 8) == 0) {
		    if (extract_doubles(argv[iarg] + 10, 1, &pcp->lamratio) != 1) {
			report_error(
			    "expected one floating-point value for -lamratio"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "amf=", 4) == 0) {
		    if (extract_doubles(argv[iarg] + 6, 1, &pcp->lamf) != 1) {
			report_error(
			    "expected one floating-point value for -lamf"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "ampc=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->lampc) != 1) {
			report_error(
			    "expected one floating-point value for -lampc"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "ampos=", 6) == 0) {
		    if (extract_doubles(argv[iarg] + 8, 1, &pcp->lampos) != 1) {
			report_error(
			    "expected one floating-point value for -lampos"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "amsmooth=", 9) == 0) {
		    if (extract_doubles(argv[iarg] + 11, 1, &pcp->lamsmooth) != 1) {
			report_error(
			    "expected one floating-point value for -lamsmooth"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'n':
		if (strncmp(argv[iarg] + 2, "a=", 2) == 0) {
		    if (extract_doubles(argv[iarg] + 4, 1, &pfp->na) != 1 ||
			pfp->na <= 0.0) {
			report_error(
			    "expected one positive floating-point value for -na"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "cycl=", 5) == 0) {
		    int ncycl;

		    if (extract_ints(argv[iarg] + 7, 1, &ncycl) == 1 &&
			ncycl > 0) {
			pcp->iteration_count = ncycl;
		    } else {
			report_error(
			    "expected one integer value greater than zero for -ncycl"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "imm=", 4) == 0) {
		    if (extract_doubles(argv[iarg] + 6, 1, &pfp->nimm) != 1 ||
			pfp->nimm <= 0.0) {
			report_error(
			    "expected one positive floating-point values for -nimm"
			);
			had_error = 1;
		    }
		} else if (strcmp(argv[iarg] + 2, "onorm_otf") == 0) {
		    pfp->normalize_otf = 0;
		} else if (strncmp(argv[iarg] + 2, "p=", 2) == 0) {
		    if (getenv("NSLOTS") == 0) {
			if (extract_ints(argv[iarg] + 4, 1, &pfp->nthread) != 1 ||
			    pfp->nthread < 1) {
			    report_error(
				"expected one positive integer value for -np"
			    );
			    had_error = 1;
			}
		    } else {
			(void) printf(
			    "NSLOTS environment variable set; ignoring -np command-line option\n"
			);
		    }
		} else if (strncmp(argv[iarg] + 2, "zpad=", 5) == 0) {
		    if (extract_ints(argv[iarg] + 7, 1, &pfp->zpadded) != 1 ||
			pfp->zpadded < 1) {
			report_error(
			    "expected one integer value greater than zero for -nzpad"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'o':
		if (strncmp(argv[iarg] + 2, "mega=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->omega) != 1) {
			report_error(
			    "expected one floating-point value for -omega"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "plotfile=", 9) == 0) {
		    pfp->plot_file_name = argv[iarg] + 11;
		} else {
		    recognized = 0;
		}
		break;

	    case 's':
	        if (strncmp(argv[iarg] + 2, "ub=", 3) == 0) {
		    int extracted = extract_doubles(
		        argv[iarg] + 5, IW_MAX_WAVE, pfp->background
		    );

		    if (extracted < 1 || extracted > IW_MAX_WAVE) {
		        report_error(
			    "expected one to %d floating-point values for -sub",
			    IW_MAX_WAVE
			);
		        had_error = 1;
		    }
		    background_set = 1;
	        } else {
		    recognized = 0;
		}
		break;

	    case 't':
		if (argv[iarg][2] == '=') {
		    int extracted =
			extract_ints(argv[iarg] + 3, 3, pfp->tbounds);

		    if (extracted >= 1 && extracted <= 3) {
			if (pfp->tbounds[0] < 0) {
			    report_error(
				"the starting time index must not be negative"
			    );
			    had_error = 1;
			}
			if (extracted > 1 &&
			    pfp->tbounds[1] < pfp->tbounds[0]) {
			    report_error(
				"the final time index must not be smaller than the first time index"
			    );
			    had_error = 1;
			    
			}
			if (extracted > 2 && pfp->tbounds[2] < 1) {
			    report_error(
				"the increment between time indices must be positive"
			    );
			    had_error = 1;
			}
		    } else {
			report_error("-t expects one to three integer values");
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "ol=", 3) == 0) {
		    if (extract_doubles(argv[iarg] + 5, 1, &pcp->tol) != 1) {
			report_error(
			    "expected one floating-point value for -tol"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'w':
		if (argv[iarg][2] == '=') {
		    int extracted = extract_ints(
			argv[iarg] + 3, IW_MAX_WAVE, pfp->wavelength_indices
		    );

		    if (extracted >= 1 && extracted <= IW_MAX_WAVE) {
			pfp->wavelength_indices[extracted] = -1;
		    } else {
			report_error(
			    "expected 1 to %d values for -w", IW_MAX_WAVE
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'x':
		if (argv[iarg][2] == '=') {
		    int extracted =
			extract_ints(argv[iarg] + 3, 2, pfp->xbounds);

		    if (extracted == 1 || extracted == 2) {
			if (pfp->xbounds[0] < 0) {
			    report_error(
				"the starting index of the x range must not be negative"
			    );
			    had_error = 1;
			}
			if (extracted > 1 &&
			    pfp->xbounds[1] < pfp->xbounds[0]) {
			    report_error(
				"the final index of the x range must not be smaller than the starting index"
			    );
			    had_error = 1;
			}
		    } else {
			report_error("expected two integers for -x");
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'y':
		if (argv[iarg][2] == '=') {
		    int extracted =
			extract_ints(argv[iarg] + 3, 2, pfp->ybounds);

		    if (extracted == 1 || extracted == 2) {
			if (pfp->xbounds[0] < 0) {
			    report_error(
				"the starting index of the y range must not be negative"
			    );
			    had_error = 1;
			}
			if (extracted > 1 &&
			    pfp->ybounds[1] < pfp->ybounds[0]) {
			    report_error(
				"the final index of the y range must not be smaller than the starting index"
			    );
			    had_error = 1;
			}
		    } else {
			report_error("expected two integers for -y");
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'z':
		if (argv[iarg][2] == '=') {
		    int extracted =
			extract_ints(argv[iarg] + 3, 2, pfp->zbounds);

		    if (extracted == 1 || extracted == 2) {
			if (pfp->zbounds[0] < 0) {
			    report_error(
				"the starting index of the y range must not be negative"
			    );
			    had_error = 1;
			}
			if (extracted > 1 &&
			    pfp->zbounds[1] < pfp->zbounds[0]) {
			    report_error(
				"the final index of the z range must not be smaller than the starting index"
			    );
			    had_error = 1;
			}
		    } else {
			report_error("expected two integers for -z");
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;
	    }
	}
	if (! recognized) {
	    report_error("unrecognized option: %s", argv[iarg]);
	    had_error = 1;
	}
	++iarg;
    }

    if (pcp->lampc > 0.0 && pfp->na > pfp->nimm) {
	report_error(
	    "the numerical aperture, %f, exceeds the immersion media refractive index, %f",
	    pfp->na,
	    pfp->nimm
	);
	had_error = 1;
    }
    if (! background_set) {
	int i;

	if (pfp->absolute_background) {
	    for (i = 0; i < IW_MAX_WAVE; ++i) {
		pfp->background[i] = 0.0;
	    }
	} else {
	    for (i = 0; i < IW_MAX_WAVE; ++i) {
		pfp->background[i] = 1.0;
	    }
	}
    }
    if (had_error) {
	print_usage();
	return 1;
    }

    return 0;
}


static
void
print_parameters(const COREParameters* pcp, FrontEndParameters* pfp)
{
    int i;

    (void) printf("Measurement file:     %s\n", pfp->measurement_file_name);
    (void) printf("Result file:          %s\n", pfp->result_file_name);
    (void) printf("CTF file:             %s\n", pfp->ctf_file_name);
    (void) printf(
	"Plot file:            %s\n",
	(pfp->plot_file_name == 0 || pfp->plot_file_name[0] == '\0') ?
	    "not written" : pfp->plot_file_name
    );
    (void) printf(
	"X range:              %d %d\n", pfp->xbounds[0], pfp->xbounds[1]
    );
    (void) printf(
	"Y range:              %d %d\n", pfp->ybounds[0], pfp->ybounds[1]
    );
    (void) printf(
	"Z range:              %d %d\n", pfp->zbounds[0], pfp->zbounds[1]
    );
    (void) printf(
	"Time range:           %d %d %d\n",
	pfp->tbounds[0],
	pfp->tbounds[1],
	pfp->tbounds[2]
    );
    (void) printf("Wavelengths:         ");
    i = 0;
    while (pfp->wavelength_indices[i] != -1) {
	(void) printf(" %d", pfp->wavelength_indices[i]);
	++i;
    }
    (void) printf("\nBackground:          ");
    i = 0;
    while (pfp->wavelength_indices[i] != -1) {
	(void) printf(" %.6g", pfp->background[i]);
	++i;
    }
    (void) printf(
	" %s\n", (pfp->absolute_background) ? "(absolute)" : "(relative)"
    );
    (void) printf(
	"Normalize CTF:        %s\n", (pfp->normalize_otf) ? "yes" : "no"
    );
    (void) printf("Nonlinearity:         %.6g\n", pcp->alpha);
    (void) printf("Ratio:                %.6g\n", pcp->lamratio);
    (void) printf("Smoothing:            %.6g\n", pcp->lamf);
    (void) printf("Missing cone weight:  %.6g\n", pcp->lampc);
    (void) printf("Positivity:           %.6g\n", pcp->lampos);
    (void) printf("Auxiliary smoothing:  %.6g\n", pcp->lamsmooth);
    (void) printf("Numerical aperture:   %.6g\n", pfp->na);
    (void) printf("Imm. ref. index:      %.6g\n", pfp->nimm);
    (void) printf("High cutoff:          %.6g\n", pfp->cuth);
    (void) printf("Iteration count:      %ld\n", pcp->iteration_count);
    (void) printf("Padded z size:        %d\n", pfp->zpadded);
    (void) printf("Relaxation:           %.6g\n", pcp->omega);
    (void) printf("Tolerance:            %.6g\n", pcp->tol);
    (void) printf("Thread count:         %d\n", pfp->nthread);
}


static
void
print_usage(void)
{
    (void) fprintf(
	stderr,
	"usage: core_decon measurement result CTF \\\n"
	"    [-x=start[:end]] [-y=start[:end]] [-z=start[:end]] \\\n"
	"    [-t=start[:end[:step]]] [-w=w1[:w2...]] [-abs_background] \\\n"
	"    [-alpha=a] [-cuth=t] [-lamratio=r] [-lamf=f] [-lampc=p] \\\n"
	"    [-lampos=p] [-lamsmooth=s] [-na=v] [-ncycl=n] [-nimm=n] \\\n"
	"    [-nonorm_ctf] [-np=n] [-nzpad=n] [-omega=o] \\\n"
	"    [-oplotfile=pltname] [-sub=b1[:b2...]] [-tol=t] \n"
    );
}


static
void
report_error(const char* fmt, ...)
{
    va_list ap;

    (void) fprintf(stderr, "core_decon: ");
    va_start(ap, fmt);
    (void) vfprintf(stderr, fmt, ap);
    va_end(ap);
    (void) fprintf(stderr, "\n");
}


/*
 * When tuning parameters, Muthuvel had rescaled the data to be in the
 * range of zero to one.  Here rescale by finding the highest intensity, h,
 * after removing the highest floor(cuth * size[0] * size[1] * size[2])
 * elements and then multiply the entire volume by 1/h.  Assumes size[0],
 * size[1], and size[2] are positive and there product fits in a size_t.
 * Assumes cuth falls in the range [0, 1).
 *
 * If successful, returns zero, the modified volume values are in measurement,
 * and *p_scale is set to h.  If unsuccessful, returns one, measurement is
 * not modified, and *p_scale is not modified.
 */
static
int
rescale(
    const long size[3], double cuth, double* measurement, double* p_scale
)
{
    double* work;
    const double* cursor;
    double inv_h;
    size_t nmax, ncut, nochild_lim, nskip, i, ix;
    int keep_high;

    nmax = size[0] * (size_t) size[1] * (size_t) size[2];
    ncut = (size_t) floor(cuth * nmax);
    if (ncut > nmax / 2) {
	ncut = nmax - ncut;
	keep_high = 0;
    } else {
	++ncut;
	keep_high = 1;
    }
    if (ncut > SIZE_MAX / sizeof(double)) {
	return 1;
    }
    work = (double*) malloc(ncut * sizeof(double));
    if (ncut == 0) {
	return 1;
    }

    nochild_lim = (ncut - 1) / 2;
    nskip = 2 * (size[0] / 2 + 1) - size[0];
    ix = 0;
    cursor = measurement;
    for (i = 0; i < nochild_lim; ++i) {
	work[ncut - 1 - i] = *cursor;
	++cursor;
	++ix;
	if (ix == size[0]) {
	    ix = 0;
	    cursor += nskip;
	}
    }
    if (keep_high) {
	for (i = nochild_lim; i < ncut; ++i) {
	    work[ncut - 1 - i] = *cursor;
	    adjust_min_double_heap(ncut - 1 - i, i + 1, work);
	    ++cursor;
	    ++ix;
	    if (ix == size[0]) {
		ix = 0;
		cursor += nskip;
	    }
	}
	for (i = ncut; i < nmax; ++i) {
	    if (work[0] < *cursor) {
		work[0] = *cursor;
		adjust_min_double_heap(0, ncut, work);
	    }
	    ++cursor;
	    ++ix;
	    if (ix == size[0]) {
		ix = 0;
		cursor += nskip;
	    }
	}
    } else {
	for (i = nochild_lim; i < ncut; ++i) {
	    work[ncut - 1 - i] = *cursor;
	    adjust_max_double_heap(ncut - 1 - i, i + 1, work);
	    ++cursor;
	    ++ix;
	    if (ix == size[0]) {
		ix = 0;
		cursor += nskip;
	    }
	}
	for (i = ncut; i < nmax; ++i) {
	    if (work[0] > *cursor) {
		work[0] = *cursor;
		adjust_max_double_heap(0, ncut, work);
	    }
	    ++cursor;
	    ++ix;
	    if (ix == size[0]) {
		ix = 0;
		cursor += nskip;
	    }
	}
    }

    *p_scale = work[0];
    inv_h = 1.0 / work[0];
    free(work);
    for (i = 0; i < size[1] * (size_t) size[2]; ++i) {
	for (ix = 0; ix < size[0]; ++ix) {
	    measurement[ix] *= inv_h;
	}
	measurement += size[0] + nskip;
    }

    return 0;
}


static
int
setup_ctf(const FrontEndParameters* pfp, int* p_im_stream)
{
    int nxyz[3], mxyz[3], fmt;
    float fmin, fmax, fmean;

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for CTF"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->ctf_file_name, "ro") != 0) {
	report_error(
	    "unable to open CTF file for reading", pfp->ctf_file_name
	);
	*p_im_stream = -1;
	return 1;
    }

    IMRdHdr(*p_im_stream, nxyz, mxyz, &fmt, &fmin, &fmax, &fmean);

    if (nxyz[0] < 1 || nxyz[1] < 1 || nxyz[2] < 1) {
	report_error(
	    "non-positive dimensions for CTF file, %s: %d %d %d",
	    pfp->ctf_file_name,
	    nxyz[0],
	    nxyz[1],
	    nxyz[2]
	);
	return 1;
    }

    return 0;
}


static
int
setup_input(FrontEndParameters* pfp, int* p_im_stream)
{
    int nxyz[3], mxyz[3], nw, nt, intlv, fmt, status, i;
    float wavelengths[IW_MAX_WAVE], fmin, fmax, fmean;

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for measurement"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->measurement_file_name, "ro") != 0) {
	report_error(
	    "unable to open measurement, %s, for reading",
	    pfp->measurement_file_name
	);
	*p_im_stream = -1;
	return 1;
    }

    IMRdHdr(*p_im_stream, nxyz, mxyz, &fmt, &fmin, &fmax, &fmean);

    if (nxyz[0] < 1 || nxyz[1] < 1 || nxyz[1] < 1) {
	report_error(
	    "non-positive dimensions for measurement file, %s: %d %d %d",
	    pfp->measurement_file_name,
	    nxyz[0],
	    nxyz[1],
	    nxyz[2]
	);
	return 1;
    }
    if (fmt == IW_COMPLEX || fmt == IW_COMPLEX_SHORT) {
	report_error(
	    "measurement file, %s, has complex-valued data",
	    pfp->measurement_file_name
	);
	return 1;
    }

    IMRtZWT(*p_im_stream, &nxyz[2], &nw, &nt, &intlv);
    IMRtWav(*p_im_stream, &nw, wavelengths);

    status = 0;

    if (pfp->xbounds[1] < 0) {
	pfp->xbounds[1] = nxyz[0] - 1;
    }
    if (pfp->ybounds[1] < 0) {
	pfp->ybounds[1] = nxyz[1] - 1;
    }
    if (pfp->zbounds[1] < 0) {
	pfp->zbounds[1] = nxyz[2] - 1;
    }
    if (pfp->tbounds[1] < 0) {
	pfp->tbounds[1] = nt - 1;
    }

    if (pfp->xbounds[0] > pfp->xbounds[1] || pfp->xbounds[1] > nxyz[0] ||
	pfp->ybounds[0] > pfp->ybounds[1] || pfp->ybounds[1] > nxyz[1] ||
	pfp->zbounds[0] > pfp->zbounds[1] || pfp->zbounds[1] > nxyz[2] ||
	pfp->tbounds[0] > pfp->tbounds[1] || pfp->tbounds[1] > nt) {
	report_error(
	    "bounds to process (%d:%d,%d:%d,%d:%d,%d:%d) have a minimum that exceeds a maximum or fall outside the measurement (size %d,%d,%d,%d)",
	    pfp->xbounds[0],
	    pfp->xbounds[1],
	    pfp->ybounds[0],
	    pfp->ybounds[1],
	    pfp->zbounds[0],
	    pfp->zbounds[1],
	    pfp->tbounds[0],
	    pfp->tbounds[1],
	    nxyz[0],
	    nxyz[1],
	    nxyz[2],
	    nt
	);
	status = 1;
    }

    if (pfp->wavelength_indices[0] == -1) {
	for (i = 0; i < nw; ++i) {
	    pfp->wavelength_indices[i] = i;
	}
	if (nw >= 1) {
	    pfp->wavelength_indices[nw] = -1;
	} else {
	    report_error(
		"measurement file, %s, has zero wavelengths",
		pfp->measurement_file_name
	    );
	    status = 1;
	}
    } else {
	i = 0;
	while (pfp->wavelength_indices[i] != -1) {
	    if (pfp->wavelength_indices[i] >= 200 &&
		pfp->wavelength_indices[i] < 2000) {
		int j = 0;

		while (j < nw &&
		       (int) wavelengths[j] != pfp->wavelength_indices[i]) {
		    ++j;
		}
		if (j < nw) {
		    pfp->wavelength_indices[i] = j;
		} else {
		    report_error(
			"selected wavelength, %d, not present in measurement",
			pfp->wavelength_indices[i]
		    );
		    status = 1;
		}
	    } else if (pfp->wavelength_indices[i] < 0 ||
		       pfp->wavelength_indices[i] >= nw) {
		report_error(
		    "selected wavelength, %d, not present in measurement",
		    pfp->wavelength_indices[i]
		);
		status = 1;
	    }
	    ++i;
	}
    }

    if (pfp->xpadded < 0) {
	pfp->xpadded = pfp->xbounds[1] - pfp->xbounds[0] + 1;
    } else if (pfp->xpadded < pfp->xbounds[1] - pfp->xbounds[0] + 1) {
	report_error(
	    "padded x size, %d, is smaller than size, %d, of processed region",
	    pfp->xpadded,
	    pfp->xbounds[1] - pfp->xbounds[0] + 1
	);
	status = 1;
    } else if (pfp->xpadded % 2 != 0) {
	report_error("padded x size, %d, is not even", pfp->xpadded);
	status = 1;
    }
    if (pfp->ypadded < 0) {
	pfp->ypadded = pfp->ybounds[1] - pfp->ybounds[0] + 1;
    } else if (pfp->ypadded < pfp->ybounds[1] - pfp->ybounds[0] + 1) {
	report_error(
	    "padded y size, %d, is smaller than size, %d, of processed region",
	    pfp->ypadded,
	    pfp->ybounds[1] - pfp->ybounds[0] + 1
	);
	status = 1;
    }
    if (pfp->zpadded < 0) {
	pfp->zpadded = 1;
	while (pfp->zpadded < 1.3f * (pfp->zbounds[1] - pfp->zbounds[0] + 1) &&
	       pfp->zpadded < INT_MAX / 2) {
	    pfp->zpadded += pfp->zpadded;
	}
    } else if (pfp->zpadded < pfp->zbounds[1] - pfp->zbounds[0] + 1) {
	report_error(
	    "padded z size, %d, is smaller than size, %d, of processed region",
	    pfp->zpadded,
	    pfp->zbounds[1] - pfp->zbounds[0] + 1
	);
	status = 1;
    }

    if (pfp->xpadded / 2 > INT_MAX / 2 - 1 ||
	2 * (pfp->xpadded / 2 + 1) > SIZE_MAX / sizeof(double) ||
	pfp->ypadded > SIZE_MAX / (2 * (pfp->xpadded / 2 + 1) *
				   sizeof(double)) ||
	pfp->zpadded > SIZE_MAX / (2 * (pfp->xpadded / 2 + 1) *
				   (size_t) pfp->ypadded * sizeof(double))) {
	report_error("size of padded measurement exceeds machine limits");
	status = 1;
    }

    return status;
}


static
int
setup_output(const FrontEndParameters* pfp, int in_im, int *p_im_stream)
{
    int nw;
    int nxyz[3], mxyz[3], nxyzst[3], itstart;
    float delta[3];
    float xyz[3];

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for output"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->result_file_name, "new") != 0) {
	report_error("unable to open %s for writing", pfp->result_file_name);
	*p_im_stream = -1;
	return 1;
    }

    nw = 0;
    while (pfp->wavelength_indices[nw] != -1) {
	++nw;
    }
    IMTrHdr(*p_im_stream, in_im);
    IMRtSiz(*p_im_stream, nxyz, mxyz, nxyzst);
    nxyz[0] = pfp->xbounds[1] - pfp->xbounds[0] + 1;
    nxyz[1] = pfp->ybounds[1] - pfp->ybounds[0] + 1;
    nxyzst[0] += pfp->xbounds[0];
    nxyzst[1] += pfp->ybounds[0];
    nxyzst[2] += pfp->zbounds[0];
    IMAlSiz(*p_im_stream, nxyz, nxyzst);
    IMFixExHdrZWT(
	*p_im_stream,
	pfp->zbounds[0],
	pfp->zbounds[1],
	nw,
	pfp->wavelength_indices,
	pfp->tbounds[0],
	pfp->tbounds[1],
	pfp->tbounds[2],
	1
    );

    IMAlMode(*p_im_stream, IW_FLOAT);

    /*
     * Nothing is done to set up lower resolution versions of the results
     * so set the number of resolutions to one.
     */
    IMAlRes(*p_im_stream, 1, 2);

    IMRtTSt(*p_im_stream, &itstart);
    itstart += pfp->tbounds[0];
    IMAlTSt(*p_im_stream, itstart);

    IMRtDel(*p_im_stream, delta);
    IMRtOrig(*p_im_stream, &xyz[0], &xyz[1], &xyz[2]);
    xyz[0] += delta[0] * pfp->xbounds[0];
    xyz[1] += delta[1] * pfp->ybounds[0];
    xyz[2] += delta[2] * pfp->zbounds[0];
    IMAlOrig(*p_im_stream, xyz[0], xyz[1], xyz[2]);

    return 0;
}


static
int
setup_plot(int iw, int it, int iz, int in_im, ResReporterData* prrd)
{
    prrd->last_successful_res = DBL_MAX;
    prrd->last_successful_iteration = 0;
    prrd->ncount = 0;
    prrd->status = 0;
    return 0;
}


static
int
setup_plot_file(const FrontEndParameters* pfp, ResReporterData* prrd)
{
    if (pfp->plot_file_name == 0 || pfp->plot_file_name[0] == '\0') {
	return 0;
    }
    prrd->filename = pfp->plot_file_name;
    prrd->fo = fopen(prrd->filename, "w");
    if (prrd->fo == 0) {
	report_error(
	    "unable to open plot file, %s, for writing: %s",
	    prrd->filename,
	    strerror(errno)
	);
	return 1;
    }
    return 0;
}


static
void
update_plot_file(long iiters, long iitert, double res, void* closure)
{
    ResReporterData* prrd = (ResReporterData*) closure;
    long last_successful_iteration = prrd->last_successful_iteration;

    if (iitert == 0 || iiters != prrd->last_successful_iteration) {
	prrd->last_successful_res = res;
	prrd->last_successful_iteration = iiters;
    }
    if (prrd->fo != 0) {
	if (prrd->status != 0) {
	    /*
	     * An earlier call to update_plot_file() for the same plot
	     * encountered a memory allocation error so don't do anything in
	     * this call.
	     */
	    return;
	}
	if (prrd->ncount == prrd->nalloc) {
	    if (prrd->ncount < LONG_MAX / 2 &&
		prrd->ncount < SIZE_MAX / (2 * sizeof(ResReporterPt))) {
		long new_count = (prrd->ncount > 0) ?
		    prrd->ncount + prrd->ncount : 16;
		ResReporterPt* new_values;

		new_values =
		    (ResReporterPt*) malloc(new_count * sizeof(ResReporterPt));
		if (new_values != 0) {
		    if (prrd->ncount > 0) {
			(void) memcpy(
			    new_values,
			    prrd->values,
			    prrd->ncount * sizeof(ResReporterPt)
			);
		    }
		    if (prrd->values != 0) {
			free(prrd->values);
		    }
		    prrd->values = new_values;
		    prrd->nalloc = new_count;
		} else {
		    prrd->status = 1;
		    return;
		}
	    } else {
		prrd->status = 1;
		return;
	    }
	}
	prrd->values[prrd->ncount].res = res;
	prrd->values[prrd->ncount].iteration =
	    (iitert == 0 || iiters != last_successful_iteration) ?
	    -iitert : iitert;
	++prrd->ncount;
    } else {
	prrd->ncount = iitert;
    }
}


static
int
write_result_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    const double* volume,
    double scale,
    float mmm[3]
)
{
    int nxp = 2 * (pfp->xpadded / 2 + 1);
    float vol_sum = 0.0f;
    const double* in = volume;
    float* section = (float*) malloc(
	(size_t) (pfp->xbounds[1] - pfp->xbounds[0] + 1) *
	(size_t) (pfp->ybounds[1] - pfp->ybounds[0] + 1) * sizeof(float)
    );
    int iz;

    if (section == 0) {
	report_error(
	    "memory allocation for output buffer failed: %s", strerror(errno)
	);
	return 1;
    }

    for (iz = 0; iz < pfp->zbounds[1] - pfp->zbounds[0] + 1; ++iz) {
	float sec_sum = 0.0f;
	float* out = section;
	int iy;

	for (iy = 0; iy < pfp->ybounds[1] - pfp->ybounds[0] + 1; ++iy) {
	    float line_sum = 0.0f;
	    int ix;

	    for (ix = 0; ix < pfp->xbounds[1] - pfp->xbounds[0] + 1; ++ix) {
		if (scale * in[ix] >= -FLT_MAX) {
		    if (scale * in[ix] <= FLT_MAX) {
			out[ix] = scale * in[ix];
			if (mmm[0] > out[ix]) {
			    mmm[0] = out[ix];
			}
			if (mmm[1] < out[ix]) {
			    mmm[1] = out[ix];
			}
		    } else {
			out[ix] = FLT_MAX;
			mmm[1] = FLT_MAX;
		    }
		} else {
		    out[ix] = -FLT_MAX;
		    mmm[0] = -FLT_MAX;
		}
		line_sum += out[ix];
	    }
	    in += nxp;
	    out += pfp->xbounds[1] - pfp->xbounds[0] + 1;
	    sec_sum += line_sum;
	}
	in += nxp *
	    (size_t) (pfp->ypadded - (pfp->ybounds[1] - pfp->ybounds[0] + 1));
	vol_sum += sec_sum;

	if (IMPosnZWT(out_im, iz, iw, it) != 0) {
	    report_error(
		"unable to position output file at z, wavelength, and time indices, %d, %d, %d",
		iz,
		iw,
		it
	    );
	    free(section);
	    return 1;
	}
	IMWrSec(out_im, section);
    }

    mmm[2] += vol_sum /
        ((pfp->xbounds[1] - pfp->xbounds[0] + 1) *
	 (size_t) (pfp->ybounds[1] - pfp->ybounds[0] + 1) *
	 (size_t) (pfp->zbounds[1] - pfp->zbounds[0] + 1)
        );
    free(section);

    return 0;
}
