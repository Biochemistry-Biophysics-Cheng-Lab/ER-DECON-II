/*
 * Implements a command-line front-end to Muthuvel Arigovindan's ER Type II
 * deconvolution algorithm.  Uses Priism's libraries for input and output of
 * MRC format files.
 *
 * Usage:
 *     core2_decon measurement result CTF [options]
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
 *     10000.
 *
 * -cuth=t
 *     By default, the intensities in each volume are rescaled to reduce
 *     the sensitivity of the optimal parameter values to the range of values
 *     in the data.  After deconvolution, the intensities in the result are
 *     rescaled by the inverse of the factor applied to original value.  The
 *     rescaled intensity of one prior to deconvolution matches the intensity
 *     in the original volume at which the fraction, t, of the values have
 *     greater intensities.  If t is not greater than or equal to zero and less
 *     than one, the rescaling is not done.
 *
 * -guess=filename
 *     Causes the initial guess for the deconvolution result to be read from
 *     MRC-formatted file specified rather than use an internally calculated
 *     initial guess.  The dimension of the file must match the region
 *     selected with the -x, -y, -z, -t, and -w options except for the z
 *     dimension for a 3D deconvolution which must either be the size selected
 *     with -z or the padded size set with -nzpad.  When the initial guess
 *     does not match the padded size set with -nzpad, the additional sections
 *     will be filled zeros.
 *
 * -imask=filename
 *     Causes the initial spatial constraint mask to be read from the
 *     MRC-formatted file specified.  The data in the file must be real-valued,
 *     and it must have x and y dimensions that match the padded dimensions of
 *     the processed data.  If the deconvolution is done in 2D, the z dimension
 *     of the mask file must be one or match the number of z sections
 *     selected from the input with -z.  If the z dimension is one, the same
 *     initial mask is used for all z sections.  If the deconvolution is done
 *     in 3D, the z dimension of the mask file must match the padded z
 *     dimension of the processed data.  The number of wavelengths in the mask
 *     file must be one or match the number of wavelengths selected from the
 *     input with -w.  If the number of wavelengths is one, the same initial
 *     mask will be used for each wavelength.  The number of time points in
 *     the mask file must be one or match the number of time points selected
 *     from the input with -t.  If the number of time points is one, the same
 *     initial mask will be used for each time point.  This option will be
 *     parsed but not otherwise used if the weight set with -lammask is zero.
 *     If the weight is positive and no -imask option has been set, the
 *     initial mask will zero everywhere indicating that all positions are
 *     completely in the the region of support.
 *
 * -imcf=filename
 *     The CORE algorithm uses a filter to compensate for the limited region
 *     of support for a 3D optical transfer function.  Normally, that filter
 *     is generated automatically from the supplied optical transfer function.
 *     Instead, you can use the -imcf option to have the filter read from
 *     the MRC-formatted file whose name you supply with the -imcf option.  The
 *     data in the filter file must be complex-valued and it must have x, y,
 *     and z dimensions of nx / 2 + 1 x ny x nz where nx, ny, and nz are the
 *     padded dimensions of the processed data.  The ordering of the frequency
 *     components in the file is the same as for the file with the optical
 *     transfer function.  When the filter is read from the file, the
 *     -na and -nimm options no longer affect the deconvolution result.
 *     
 * -include_padding
 *     If specified and the deconvolution is done in 3D, the additional z
 *     sections for padding will not be stripped out when the result is
 *     written to the output file.  By default, those sections will be
 *     stripped.  -include_padding is useful if you want to later use the
 *     result as the initial guess for additional deconvolution iterations
 *     and have the result of that better match what would have happened
 *     if all the iterations had been done with one run from core2_decon.
 *
 * -lamratio=r1:r2
 *     If omitted, the default values are zero for r1 and 1 for r2.  Both
 *     r1 must be between zero and one.  r2 must be non-negative.  For the 2014
 *     line search algorithm, a valid value for r1 will have no effect as that
 *     algorithm always takes r1 to be zero.  The value for r2 will also be
 *     ignored if you use the 2014 line search with -regtype=tv2.
 *
 * -lamf=f
 *     Sets the smoothing parameter.  The user will have to tune this for
 *     the data set.
 *
 * -laml2=f
 *     During development of the ER Type II deconvolution Muthuvel experimented
 *     with a term in the entropy that was, in the notation used in the 2013
 *     PNAS paper, a scalar times the sum over r of Eg(r).  He abandoned the
 *     use of that term in the PNAS paper and in the Matlab implementation of
 *     his algorithm, but it remained in this C implementation where the
 *     scalar had a fixed value of 0.0001.  The -laml2 option allows you to
 *     set the scalar for the additional term in the entropy.  The value of f
 *     must be a non-negative floating-point value.  If you do not specify
 *     this option, zero will be used as the value of f.  That gives the
 *     entropy described in the 2013 PNAS paper.  If you want to reproduce or
 *     mimic the results by the original C version of core2_decon, use a value
 *     of 0.0001 for f.
 *
 * -lammask=f
 *     Sets the weight for the spatial constraint mask.  The weight must be a
 *     a non-negative values.  If this option is  omitted, the weight for the
 *     spatial constraint mask is zero.  When the weight is zero, the other
 *     options controlling the mask (-imask, -omask, -umask) are parsed but
 *     are not used.
 *
 * -lampc=p
 *     Sets the weight for the missing cone filter.  If this option is omitted,
 *     the weighting for the missing cone filter is one hundred.
 *
 * -lampos=p
 *     Controls the positivity constraint.  If this option is omitted, the
 *     default value is one.
 *
 * -lamsmooth=s
 *     Sets an auxiliary smoothing parameter.  If this option is omitted,
 *     the auxiliary smoothing parameter is one hundred.
 *
 * -linesearch=m
 *     Selects the line search algorithm.  m may be 2013 which selects the
 *     algorithm used in the 2013 PNAS paper or it may be 2014 which selects
 *     the line search algorithm from Muthuvel Arigovindan's 4D deconvolution.
 *     The 2014 algorithm is the default as it is expected to be more robust.
 *     The -omega and -tol options control parameters that are used with the
 *     2013 line search algorithm.  The -regtype option controls a parameter
 *     used with the 2014 line search algorithm.  The choice of the line search
 *     algorithm has a side effect of altering the contents of the file
 *     generated when you use the -oplotfile option.  With the 2014 line
 *     search algorithm, using the -rfactor option will include R-factor
 *     results in the file generated by -oplotfile.
 *
 * -logging=level
 *     Selects which messages will be generated by the applicaton as it runs.
 *     level may be error, warning, info, progress, or debug.  If you do not
 *     set a logging level, the level will be info.  As currently implemented,
 *     the only change from altering the log level will be if the level is set
 *     to progress or debug.  Then messages will be written at each iteration
 *     to indicate how far the algorithm has progressed.
 *
 * -na=v
 *     Sets the numerical aperture (used in the computation of the missing
 *     cone filter) to be v.  If this option is omitted, the numerical aperture
 *     will be assumed to be 1.4.  This option only has an effect when you
 *     do not use -imcf option or give that option an empty file name.
 *
 * -nimm=n
 *     Sets the immersion media refractive index (used in the computation of
 *     the missing cone filter) to be n.  If this option is omitted, the
 *     refractive index of the immersion media is assumed to be 1.512.  The
 *     refractive index of the immersion media must be greater than or equal
 *     to the numerical aperture.  This option only has an effect when you
 *     do not use -imcf option or give that option an empty file name.
 *
 * -ncycl=n
 *     Set ths number of iterations.  n must be a non-negative integer.  If
 *     this option is omitted, fifty iteration will be used.
 *
 * -nonorm_otf
 *     Do not rescale the transfer function to have a zero frequency response
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
 * -omask=filename
 *     Causes the final value for the spatial constraint mask to be written to
 *     the MRC-formatted file specified.  This option will be parsed but not
 *     otherwise used if the weight set with -lammask is zero.
 *
 * -omcf=filename
 *     Causes the filter used to compensate for the limited region of support
 *     to be written in MRC format to the file with the given name.
 *
 * -omega=o
 *     Sets the relaxation parameter used with the 2013 line search algorithm
 *     (see -linesearch).  If this option is omitted, the relaxation parameter
 *     will be 0.8.
 *
 * -oplotfile=pltname
 *     With the 2013 line search algorithm (see -linesearch), causes a history
 *     of the termination criteria as a function of the iteration number to be
 *     written to the file named pltname.  With the 2014 line search algorithm,
 *     causes a history of the cost and, optionally, the R-factor to be
 *     written to the file named pltname.
 *
 * -regtype=t
 *     Selects the form of the regularization to use with the 2014 line
 *     search (see -linesearch).  t may be ma or tv2.  The default is ma.
 *
 * -rfactor
 *     Causes R-factor results to be included in the file specified by
 *     -oplotfile if you are using the 2014 line search algorithm.
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
 *     Specifies the termination tolerance used with the 2013 line search
 *     algorithm (see -linesearch).  If this option is omitted, the
 *     termination tolerance will be 0.0001.
 *
 * -umask=plugin_name[,[interface_name],...]
 *     Specifies the plugin to use for updating the spatial constraint mask and
 *     the options to use with that plugin.  The plugin's name may be either
 *     a path name for the module to load, or one of the short names listed
 *     in core2_masku_plugins.def.  core2_masku_plugins.def may either be
 *     in the CONFIG directory of the Priism distribution or in the .iveprefs
 *     subdirectory within your home directory.  After the plugin's name and
 *     separated from it by a comma, you may specify the name of the
 *     initialization routine to call for the plugin.  That is typically only
 *     necessary if the plugin supports multiple update algorithms each with a
 *     different initialization routine or the plugin does not use the default
 *     initialization routine, initialize_plugin.  If you use a short name for
 *     a plugin, the name of the initialization routine will be ignored and
 *     what was specified in core2_masku_plugins.def will be used.  After the
 *     name of the initialization routine and separated from it by a comma,
 *     specify any of the options for the plugin.  If you are specifying
 *     options but want to use the default plugin initialization routine,
 *     use -umask=plugin_name,,opts where opts are the options you want to
 *     specify.  The format of the options is specific to each plugin.  The
 *     -umask option will be parsed but have no other effect if the weight set
 *     with -lammask is zero.
 */

#include "ive_standards.h"  /* for ive_thread.h */
#include "ive_thread.h"
#include "core2.h"
#include "core2_2d.h"
#include "core_estimate_background.h"
#include "getmcf.h"
#include "masku_plugin.h"
#include "interp_util.h"
#include "binary_heap.h"
#include "IMInclude.h"
#include <assert.h>
#include <errno.h>
#include <float.h>       /* DBL_MAX, FLT_MAX */
#include <limits.h>      /* CHAR_BIT, INT_MAX, LONG_MAX, SIZE_MAX */
#include <math.h>        /* fabs, floor */
#include <stdarg.h>      /* va_end, va_start, vfprintf */
#include <stdio.h>       /* fclose, ferror, fflush, fopen, fprintf, printf,
			    sscanf, stderr, stdout, vfprintf */
#include <stdlib.h>      /* free, getenv, malloc, strtol */
#include <string.h>      /* memcpy, strchr, strcpy, strerror, strncmp */
#ifndef SIZE_MAX
#define SIZE_MAX ((((((size_t) 1) << (sizeof(size_t) * CHAR_BIT - 1)) - 1) << 1) + 1)
#endif


typedef enum LogLevel {
    LOG_ERROR, LOG_WARN, LOG_INFO, LOG_PROGRESS, LOG_DEBUG
} LogLevel;

typedef struct FrontEndParameters {
    char* measurement_file_name;
    char* result_file_name;
    char* ctf_file_name;
    char* in_filter_name;
    char* out_filter_name;
    char* plot_file_name;
    char* in_mask_name;
    char* out_mask_name;
    char* masku_plugin_name;
    char* masku_intfc_name;
    char* guess_name;
    MaskuPlugin masku_plug;
    MaskuUpdaterClosure masku_closure;
    double background[IW_MAX_WAVE];
    double lamsmooth;
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
    int as_2d;
    /*
     * After the initial setup and error checking, this will only be true if
     * as_2d is also true.
     */
    int as_2d_ctf;
    int include_padding;
    LogLevel log_level;
} FrontEndParameters;

typedef struct OldStylePt {
    double res;
    /*
     * The absolute value of this is the total number of iterations at the
     * time res was computed.  If iteration is negative, the iteration was
     * successful.
     */
    long iteration;
} OldStylePt;

typedef struct NewStylePt {
    double cost;
    double cost_fidel;
    double cost_reg;
    double cost_pos;
    double cost_laml2;
    double cost_spa;
    double rfactor;
    int iteration;
} NewStylePt;

typedef struct ResReporterData {
    char* filename;
    FILE* fo;
    OldStylePt* old_values;
    NewStylePt* new_values;
    double last_successful_res;
    long last_successful_iteration;
    long ncount, nalloc;
    int uses_old_style_results;
    int outputs_rfactor;
    int updates_stdout;
    int status;
} ResReporterData;

static int complete_plot(
    int iw,
    int it,
    int iz,
    int in_im,
    int includes_laml2,
    int includes_mask,
    ResReporterData* prrd
);
static int complete_plot_file(ResReporterData* prrd);
static void dzero(size_t n, double* v);
static int extract_doubles(const char* s, int max_count, double* values);
static int extract_ints(const char* s, int max_count, int* values);
static double* fetch_ctf(
    const FrontEndParameters* pfp, int in_im, int ctf_im, int iw
);
static int fetch_guess(
    const FrontEndParameters* pfp,
    int in_im,
    int iw,
    int it,
    int iz,
    double scale,
    double* values
);
static int fetch_mask(
    const FrontEndParameters* pfp,
    int in_im,
    int iw,
    int it,
    int iz,
    double* values
);
static int fetch_mcf(
    const FrontEndParameters* pfp, int in_im, int iw, double* values
);
static double* fetch_measurement(
    const FrontEndParameters* pfp, int in_im, int iw, int it, int iz
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
    int argc, char* argv[], Core2Params* pcp, FrontEndParameters* pfp
);
static int parse_masku_command_line(
    FrontEndParameters* pfp,
    char* arg,
    AvailableMaskuPlugin* plugin_list,
    int list_count
);
static void print_for_masku(void* completion, const char* fmt, ...);
static void print_parameters(
    const Core2Params* pcp, FrontEndParameters* pfp
);
static void print_usage(void);
static void report_error(const char* fmt, ...);
static int rescale(
    const long size[3], double cuth, double* measurement, double* p_scale
);
static int setup_ctf(const FrontEndParameters* pfp, int* p_im_stream);
static int setup_guess(const FrontEndParameters* pfp, int* p_im_stream);
static int setup_imcf(FrontEndParameters* pfp, int* p_im_stream);
static int setup_input(FrontEndParameters* pfp, int* p_im_stream);
static int setup_mask_input(FrontEndParameters* pfp, int* p_im_stream);
static int setup_mask_output(
    FrontEndParameters* pfp, int template_im, int* p_im_stream
);
static int setup_omcf(
    FrontEndParameters* pfp, int template_im, int* p_im_stream
);
static int setup_output(
    const FrontEndParameters* pfp, int in_im, int *p_im_stream
);
static int setup_plot(
    int iw, int it, int iz, int in_im, ResReporterData* prrd
);
static int setup_plot_file(
    const FrontEndParameters* pfp, ResReporterData* prrd
);
static void update_new_plot(
    int iter, const ResReporterCosts* p_costs, double rfactor, void* closure
);
static void update_old_plot(
    long iiters, long iitert, double res, void* closure
);
static int write_mcf(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    const double* values,
    float mmm[3]
);
static int write_mask_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    int iz,
    const double* mask,
    float mmm[3]
);
static int write_result_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    int iz,
    const double* volume,
    double scale,
    float mmm[3]
);


int
main(int argc, char* argv[])
{
    Core2Params cp;
    FrontEndParameters fp;
    ResReporterData rrd = { 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0 };
    ResReporter rr = { &update_new_plot, &rrd, &update_old_plot, &rrd };
    const char old_col_labels[] =
	"time   scale          resf           NiterS     NIterAcc";
    const char old_prg_labels[] = "nsucc  ntot   resf";
    const char new_col_labels[] =
	"time   scale          init. cost    final cost";
    const char new_prg_labels[] =
	"iter   cost          cost_fid      cost_reg      cost_pos";
    const char new_rf_col_labels[] =
	"time   scale          init. cost    final cost    init. rfactor  final rfactor";
    const char new_rf_prg_labels[] =
	"iter   cost          cost_fid      cost_reg      cost_pos      rfactor";
    int in_im = -1, ctf_im = -1, in_filt_im = -1, out_im = -1, out_filt_im = -1;
    int in_mask_im = -1, out_mask_im = -1, guess_im = -1;
    double* deconvoluted = 0;
    double* ctf_copy = 0;
    double* mask = 0;
    int status = 1;
    int iw;
    float delta[3];
    const char* col_labels;
    const char* prg_labels;

    if (parse_command_line(argc, argv, &cp, &fp) != 0) {
	goto cleanup;
    }
    if (setup_input(&fp, &in_im) != 0 || setup_ctf(&fp, &ctf_im) != 0) {
	goto cleanup;
    }
    /*
     * If processing as 2D, do not internally calculate the cone filter so
     * disable the cone filter unless one is provided by the user.
     */
    if (fp.as_2d && (fp.in_filter_name == 0 || fp.in_filter_name[0] == '\0')) {
	cp.lampc = 0.0;
    }
    /*
     * If the cone filter has been disabled, suppress input and output of the
     * cone filter.
     */
    if (cp.lampc <= 0.0) {
	fp.in_filter_name = 0;
	fp.out_filter_name = 0;
    }
    if ((fp.in_filter_name == 0 || fp.in_filter_name[0] == '\0') &&
	cp.lampc > 0.0 && fp.na > fp.nimm) {
	report_error(
	    "the numerical aperture, %f, exceeds the immersion media refractive index, %f",
	    fp.na,
	    fp.nimm
	);
	goto cleanup;
    }
    if (setup_imcf(&fp, &in_filt_im) != 0) {
	goto cleanup;
    }
    if (setup_guess(&fp, &guess_im) != 0) {
	goto cleanup;
    }
    if (setup_output(&fp, in_im, &out_im) != 0 ||
	setup_omcf(&fp, in_im, &out_filt_im) != 0 ||
	setup_plot_file(&fp, &rrd) != 0) {
	goto cleanup;
    }
    if (cp.lammask > 0.0 && ! cp.uses_old_line_search &&
	(setup_mask_input(&fp, &in_mask_im) != 0 ||
	 setup_mask_output(&fp, in_im, &out_mask_im) != 0)) {
	goto cleanup;
    }
    rrd.uses_old_style_results = cp.uses_old_line_search;
    rrd.outputs_rfactor = cp.outputs_rfactor;
    rrd.updates_stdout = (fp.log_level == LOG_PROGRESS ||
			  fp.log_level == LOG_DEBUG);
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
    if (cp.lammask > 0.0 && ! cp.uses_old_line_search) {
	mask = (double*) malloc(
	    2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded *
	    (size_t) fp.zpadded * sizeof(double)
	);
	if (mask == 0) {
	    report_error(
		"memory allocation for mask failed: %s", strerror(errno)
	    );
	    goto cleanup;
	}
    }

    IMRtDel(in_im, delta);
    if (delta[0] <= 0.0f || delta[1] <= 0.0f ||
	(delta[2] <= 0.0f && ! fp.as_2d)) {
	report_error(
	    "one or more of the pixel spacings for the measurement are not positive"
	);
	goto cleanup;
    }
    if (fabs(delta[0] - delta[1]) > 1e-5 * delta[0]) {
	report_error("expected data with equal x and y sampling");
	goto cleanup;
    }
    if (fp.as_2d && delta[2] <= 0.0f) {
	/*
	 * Put in a reasonable value since some parts below
	 * (core2_estimate_background() and core2_2d()) expect the ratio
	 * of delta[0] to delta[2] but then do not use the value.
	 */
	delta[2] = delta[0];
    }

    iw = 0;
    while (fp.wavelength_indices[iw] != -1) {
	double* ctf = fetch_ctf(&fp, in_im, ctf_im, iw);
	double* cone_filter = 0;
	float mmm[3] = { FLT_MAX, -FLT_MAX, 0.0f };
	float mask_mmm[3] = { FLT_MAX, -FLT_MAX, 0.0f };
	float first_wave_mmm[3];
	float first_wave_mask_mmm[3];
	float first_wave_mcf_mmm[3];
	long size[3];
	int had_error = 0, it, iz, zbounds[2];
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
	    if (fp.in_filter_name == 0 || fp.in_filter_name[0] == '\0') {
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
		    cp.nthread
		);
		if (result != 0) {
		    report_error("cone filter computation failed");
		    free(cone_filter);
		    free(ctf);
		    goto cleanup;
		}
	    } else {
		if (fetch_mcf(&fp, in_filt_im, iw, cone_filter) != 0) {
		    free(cone_filter);
		    free(ctf);
		    goto cleanup;
		}
	    }
	    if (fp.out_filter_name != 0 && fp.out_filter_name[0] != '\0') {
		float mcf_mmm[3];

		if (write_mcf(&fp, out_filt_im, iw, cone_filter, (iw == 0) ? first_wave_mcf_mmm : mcf_mmm) != 0) {
		    free(cone_filter);
		    free(ctf);
		    goto cleanup;
		}
		if (fp.wavelength_indices[iw + 1] == -1) {
		    char title[80];

		    IMWrHdr(
			out_filt_im,
			title,
			-1,
			first_wave_mcf_mmm[0],
			first_wave_mcf_mmm[1],
			first_wave_mcf_mmm[2]
		    );
		}
	    }
	}

	size[0] = fp.xpadded;
	size[1] = fp.ypadded;
	size[2] = fp.zpadded;
	if (fp.as_2d) {
	    zbounds[0] = 0;
	    zbounds[1] = fp.zbounds[1] - fp.zbounds[0];
	} else {
	    zbounds[0] = -1;
	    zbounds[1] = -1;
	}

	if (! fp.absolute_background) {
	    double background_scale = fp.background[iw];
	    double background = DBL_MAX;

	    /*
	     * Set this to zero so the subtraction in fetch_measurement() has
	     * no effect.
	     */
	    fp.background[iw] = 0.0;
	    it = 0;
	    iz = zbounds[0];
	    while (1) {
		double* measurement;
		double vback;

		if (it >=
		    (fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1) {
		    break;
		}

		measurement = fetch_measurement(&fp, in_im, iw, it, iz);
		if (measurement == 0) {
		    had_error = 1;
		    break;
		}
		result = core_estimate_background(
		    measurement,
		    size,
		    delta[0] / (double) delta[2],
		    cp.lamf,
		    fp.lamsmooth,
		    cp.nthread,
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
		++iz;
		if (iz > zbounds[1]) {
		    iz = zbounds[0];
		    ++it;
		}
	    }

	    if (had_error) {
		free(ctf);
		free(cone_filter);
		goto cleanup;
	    }

	    (void) printf("iw %d: unscaled background %g\n", iw, background);
	    fp.background[iw] = background * background_scale;
	} else {
	    (void) printf(
		"iw %d: absolute background %g\n", iw, fp.background[iw]
	    );
	}

	if (rrd.uses_old_style_results) {
	    col_labels = old_col_labels;
	    prg_labels = old_prg_labels;
	} else {
	    if (rrd.outputs_rfactor) {
		col_labels = new_rf_col_labels;
		prg_labels = new_rf_prg_labels;
	    } else {
		col_labels = new_col_labels;
		prg_labels = new_prg_labels;
	    }
	}
	if (! rrd.updates_stdout) {
	    (void) printf("%s\n", col_labels);
	}

	it = 0;
	iz = zbounds[0];
	while (1) {
	    double* measurement;
	    double scale;

	    if (it >= (fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1) {
		if (iw == 0) {
		    first_wave_mmm[0] = mmm[0];
		    first_wave_mmm[1] = mmm[1];
		    first_wave_mmm[2] = mmm[2] / (zbounds[1] - zbounds[0] + 1);
		    if (out_mask_im != -1) {
			first_wave_mask_mmm[0] = mask_mmm[0];
			first_wave_mask_mmm[1] = mask_mmm[1];
			first_wave_mask_mmm[2] = mask_mmm[2] /
			    (zbounds[1] - zbounds[0] + 1);
		    }
		} else {
		    IMAlWavMM(out_im, iw, mmm[0], mmm[1]);
		    if (out_mask_im != -1) {
			IMAlWavMM(out_mask_im, iw, mask_mmm[0], mask_mmm[1]);
		    }
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
		  if (out_mask_im != -1) {
		      IMWrHdr(
			  out_mask_im,
			  title,
			  -1,
			  first_wave_mask_mmm[0],
			  first_wave_mask_mmm[1],
			  first_wave_mask_mmm[2] /
			  ((fp.tbounds[1] - fp.tbounds[0]) / fp.tbounds[2] + 1)
		      );
		  }
		}
		break;
	    }

	    measurement = fetch_measurement(&fp, in_im, iw, it, iz);
	    if (measurement == 0) {
		had_error = 1;
		break;
	    }
	    if (fp.cuth >= 0.0 && fp.cuth < 1.0) {
		if (rescale(size, fp.cuth, measurement, &scale) != 0) {
		    report_error(
			"unable to allocate temporary space for rescaling: %s",
			strerror(errno)
		    );
		    had_error = 1;
		    break;
		}
	    } else {
		scale = 1.0;
	    }

	    if (guess_im != -1) {
		cp.skips_initial_guess = 1;
		if (fetch_guess(
			&fp, guess_im, iw, it, iz, scale, deconvoluted) != 0) {
		    had_error = 1;
		    break;
		}
	    }

	    if (mask != 0) {
		if (in_mask_im != -1) {
		    if (fetch_mask(&fp, in_mask_im, iw, it, iz, mask) != 0) {
			had_error = 1;
			break;
		    }
		} else {
		    dzero(
			2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded *
			    (size_t) fp.zpadded,
			mask
		    );
		}
	    }

	    (void) memcpy(
		ctf_copy,
		ctf,
		2 * (fp.xpadded / 2 + 1) * (size_t) fp.ypadded *
		    (size_t) fp.zpadded * sizeof(double)
	    );
	    result = setup_plot(iw, it, iz, in_im, &rrd);
	    if (result != 0) {
		had_error = 1;
		goto cleanup;
	    }
	    if (rrd.updates_stdout) {
		(void) printf("%s\n", prg_labels);
	    }
	    if (! cp.uses_old_line_search && fp.masku_closure.attr != 0 &&
		fp.masku_closure.attr->pass_initializer != 0) {
		if ((*fp.masku_closure.attr->pass_initializer)(
			fp.masku_closure.attr,
			fp.masku_closure.options) != 0) {
		    had_error = 1;
		    goto cleanup;
		}
	    }
	    if (fp.zpadded == 1) {
		result = core2_2d(
		    measurement,
		    ctf_copy,
		    cone_filter,
		    deconvoluted,
		    mask,
		    size,
		    &cp,
		    &rr,
		    (fp.masku_closure.attr != 0) ?
		        fp.masku_closure.attr->mask_updater_2d : 0,
		    &fp.masku_closure
		);
	    } else {
		result = core2(
		    measurement,
		    ctf_copy,
		    cone_filter,
		    deconvoluted,
		    mask,
		    size,
		    delta[0] / (double) delta[2],
		    &cp,
		    &rr,
		    (fp.masku_closure.attr != 0) ?
		        fp.masku_closure.attr->mask_updater_3d : 0,
		    &fp.masku_closure
		);
	    }
	    free(measurement);
	    if (! cp.uses_old_line_search && fp.masku_closure.attr != 0 &&
		fp.masku_closure.attr->pass_finalizer != 0) {
		(*fp.masku_closure.attr->pass_finalizer)(
		    fp.masku_closure.attr, fp.masku_closure.options
		);
	    }
	    if (result < 0) {
		report_error("deconvolution computation failed");
		had_error = 1;
		break;
	    }
	    result = complete_plot(
		iw, it, iz, in_im, (cp.laml2 != 0.0), (cp.lammask > 0.0), &rrd
	    );
	    if (result != 0) {
		had_error = 1;
		goto cleanup;
	    }
	    if (rrd.updates_stdout) {
		(void) printf("%s\n", col_labels);
	    }
	    if (rrd.uses_old_style_results) {
		(void) printf(
		    "%5d  %13.6e  %13.6e  %8ld  %8ld\n",
		    it,
		    1.0 / scale,
		    rrd.last_successful_res,
		    rrd.last_successful_iteration,
		    rrd.ncount
		);
	    } else {
		if (rrd.ncount > 0) {
		    if (rrd.outputs_rfactor) {
			(void) printf(
			    "%5d  %13.6e  %12.6e  %12.6e  %13.6e  %13.6e\n",
			    it,
			    1.0 / scale,
			    rrd.new_values[0].cost,
			    rrd.new_values[rrd.ncount - 1].cost,
			    rrd.new_values[0].rfactor,
			    rrd.new_values[rrd.ncount - 1].rfactor
			);
		    } else {
			(void) printf(
			    "%5d  %13.6e  %12.6e  %12.6e\n",
			    it,
			    1.0 / scale,
			    rrd.new_values[0].cost,
			    rrd.new_values[rrd.ncount - 1].cost
			);
		    }
		} else {
		    if (rrd.outputs_rfactor) {
			(void) printf(
			    "%5d  %13.6e  %12.6e  %12.6e  %13.6e  %13.6e\n",
			    it,
			    1.0 / scale,
			    0.0,
			    0.0,
			    0.0,
			    0.0
			);
		    } else {
			(void) printf(
			    "%5d  %13.6e  %12.6e  %12.6e\n",
			    it,
			    1.0 / scale,
			    0.0,
			    0.0
			);
		    }
		}
	    }
	    result = write_result_and_collect_stats(
		&fp, out_im, iw, it, iz, deconvoluted, scale, mmm
	    );
	    if (result != 0) {
		had_error = 1;
		break;
	    }

	    if (out_mask_im != -1) {
		result = write_mask_and_collect_stats(
		    &fp, out_mask_im, iw, it, iz, mask, mask_mmm
		);
	    }

	    ++iz;
	    if (iz > zbounds[1]) {
		iz = zbounds[0];
		++it;
	    }
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
    if (mask != 0) {
	free(mask);
    }
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
    if (out_mask_im != -1) {
	IMClose(out_mask_im);
    }
    if (out_filt_im != -1) {
	IMClose(out_filt_im);
    }
    if (guess_im != -1) {
	IMClose(guess_im);
    }
    if (out_im != -1) {
	IMClose(out_im);
    }
    if (in_mask_im != -1) {
	IMClose(in_mask_im);
    }
    if (in_filt_im != -1) {
	IMClose(in_filt_im);
    }
    if (ctf_im != -1) {
	IMClose(ctf_im);
    }
    if (in_im != -1) {
	IMClose(in_im);
    }
    if (fp.masku_closure.attr != 0 &&
	fp.masku_closure.attr->options_destroyer != 0 &&
	fp.masku_closure.options != 0) {
	(*fp.masku_closure.attr->options_destroyer)(
	    fp.masku_closure.attr,
	    fp.masku_closure.options
	);
    }
    if (fp.masku_plug != 0) {
	unload_masku_plugin(fp.masku_plug);
    }
    return status;
}


static
int
complete_plot(
    int iw,
    int it,
    int iz,
    int in_im,
    int includes_laml2,
    int includes_mask,
    ResReporterData* prrd
)
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

    if (prrd->uses_old_style_results) {
	(void) fprintf(
	    prrd->fo,
	    "TTermination Criteria For %s: Successful Iterations\n",
	    name
	);
	(void) fprintf(
	    prrd->fo,
	    "Xiteration\nYtermination criteria\nLsuccessful iterations\n"
	);
	for (i = 0, j = 0; i < prrd->ncount; ++i) {
	    if (prrd->old_values[i].iteration > 0) {
		continue;
	    }
	    (void) fprintf(prrd->fo, "%ld %.6f\n", j, prrd->old_values[i].res);
	    ++j;
	}
	(void) fprintf(
	    prrd->fo, "TTermination Criteria For %s: All Iterations\n", name
	);
	(void) fprintf(
	    prrd->fo, "Xiteration\nYtermination criteria\nLall iterations\n"
	);
	for (i = 0; i < prrd->ncount; ++i) {
	    (void) fprintf(prrd->fo, "%ld %.6f\n", i, prrd->old_values[i].res);
	}
    } else {
	(void) fprintf(prrd->fo, "TCost For %s\n", name);
	(void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	for (i = 0; i < prrd->ncount; ++i) {
	    (void) fprintf(
		prrd->fo,
		"%d %.6f\n",
		prrd->new_values[i].iteration,
		prrd->new_values[i].cost
	    );
	}
	if (prrd->outputs_rfactor) {
	    (void) fprintf(prrd->fo, "TR-factor For %s\n", name);
	    (void) fprintf(prrd->fo, "Xiteration\nYr-factor\nLr-factor\n");
	    for (i = 0; i < prrd->ncount; ++i) {
		(void) fprintf(
		    prrd->fo,
		    "%d %.6f\n",
		    prrd->new_values[i].iteration,
		    prrd->new_values[i].rfactor
		);
	    }
	}
	(void) fprintf(
	    prrd->fo, "TData Fidelity Cost Component For %s\n", name
	);
	(void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	for (i = 0; i < prrd->ncount; ++i) {
	    (void) fprintf(
		prrd->fo,
		"%d %.6f\n",
		prrd->new_values[i].iteration,
		prrd->new_values[i].cost_fidel
	    );
	}
	(void) fprintf(
	    prrd->fo, "TRegularization Cost Component For %s\n", name
	);
	(void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	for (i = 0; i < prrd->ncount; ++i) {
	    (void) fprintf(
		prrd->fo,
		"%d %.6f\n",
		prrd->new_values[i].iteration,
		prrd->new_values[i].cost_reg
	    );
	}
	(void) fprintf(
	    prrd->fo, "TPositivity Cost Component For %s\n", name
	);
	(void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	for (i = 0; i < prrd->ncount; ++i) {
	    (void) fprintf(
		prrd->fo,
		"%d %.6f\n",
		prrd->new_values[i].iteration,
		prrd->new_values[i].cost_pos
	    );
	}
	if (includes_laml2) {
	    (void) fprintf(
		prrd->fo, "Tlaml2 Cost Component For %s\n", name
	    );
	    (void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	    for (i = 0; i < prrd->ncount; ++i) {
		(void) fprintf(
		    prrd->fo,
		    "%d %.6f\n",
		    prrd->new_values[i].iteration,
		    prrd->new_values[i].cost_laml2
		);
	    }
	}
	if (includes_mask) {
	    (void) fprintf(
		prrd->fo, "TSpatial Constraint Cost Component For %s\n", name
	    );
	    (void) fprintf(prrd->fo, "Xiteration\nYcost\nLcost\n");
	    for (i = 0; i < prrd->ncount; ++i) {
		(void) fprintf(
		    prrd->fo,
		    "%d %.6f\n",
		    prrd->new_values[i].iteration,
		    prrd->new_values[i].cost_spa
		);
	    }
	}
    }

    if (ferror(prrd->fo) != 0) {
	report_error(
	    "an error occurred while writing the summary plot, %s",
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

    if (prrd->old_values != 0) {
	free(prrd->old_values);
	prrd->old_values = 0;
    }
    if (prrd->new_values != 0) {
	free(prrd->new_values);
	prrd->new_values = 0;
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
	prrd->fo = 0;
    }
    return status;
}


static
void
dzero(size_t n, double* v)
{
    size_t i;

    for (i = 0; i < n; ++i) {
	v[i] = 0.0;
    }
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

    if (nxyz[2] == nw && ! pfp->as_2d_ctf) {
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
	if (pfp->zpadded > 1) {
	    dkxyzo[2] = 1.0f / (pfp->zpadded * measurement_delta[2]);
	} else {
	    dkxyzo[2] = 1.0f;
	}
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
	int zbounds[2];

	if (pfp->zpadded > 1) {
	    zbounds[0] = 0;
	    zbounds[1] = nxyz[2] - 1;
	} else {
	    zbounds[0] = 0;
	    zbounds[1] = 0;
	}
	if (delta[0] <= 0.0f || delta[1] <= 0.0f ||
	    (delta[2] <= 0.0f && pfp->zpadded > 1)) {
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
	    zbounds[1] - zbounds[0] + 1 > SIZE_MAX /
	        (((is_real) ? 1 : 2) * (size_t) nxyz[0] * (size_t) nxyz[1] *
		 sizeof(float))) {
	    report_error(
		"the size of the on-file CTF exceeds the machine limits"
	    );
	    free(out);
	    return 0;
	}
	otfin = (float*) malloc(
	    ((is_real) ? 1 : 2) * nxyz[0] * (size_t) nxyz[1] *
	    (size_t) (zbounds[1] - zbounds[0] + 1) * sizeof(float)
	);
	if (otfin == 0) {
	    report_error(
		"memory allocation for CTF buffer failed: %s", strerror(errno)
	    );
	    free(out);
	    return 0;
	}

	for (iz = zbounds[0]; iz <= zbounds[1]; ++iz) {
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
	nxyzo[0] = nxhlf;
	nxyzo[1] = pfp->ypadded;
	nxyzo[2] = pfp->zpadded;
	dkxyzo[0] = 1.0f / (pfp->xpadded * measurement_delta[0]);
	dkxyzo[1] = 1.0f / (pfp->ypadded * measurement_delta[1]);
	if (pfp->zpadded > 1) {
	    delta[2] = 1.0f / (nxyz[2] * delta[2]);
	    dkxyzo[2] = 1.0f / (pfp->zpadded * measurement_delta[2]);
	} else {
	    delta[2] = 1.0f;
	    dkxyzo[2] = 1.0f;
	    nxyz[2] = 1;
	}
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
int
fetch_guess(
    const FrontEndParameters* pfp,
    int in_im,
    int iw,
    int it,
    int iz,
    double scale,
    double* values
)
{
    float* section = (float*) malloc(
	pfp->xpadded * (size_t) pfp->ypadded * sizeof(float)
    );
    int nxv = 2 * (pfp->xpadded / 2 + 1);
    int needs_z_fill = 0;
    int nz, nw, nt, intlv;
    int zbounds[2];

    if (section == 0) {
	report_error(
	    "memory allocation failed for initial guess buffer: %s",
	    strerror(errno)
	);
	return 1;
    }

    IMRtZWT(in_im, &nz, &nw, &nt, &intlv);
    if (iz == -1) {
	zbounds[0] = 0;
	if (nz >= pfp->zpadded) {
	    zbounds[1] = pfp->zpadded - 1;
	} else {
	    zbounds[1] = nz - 1;
	    needs_z_fill = 1;
	}
    } else {
	zbounds[0] = iz;
	zbounds[1] = iz;
    }
    for (iz = zbounds[0]; iz <= zbounds[1]; ++iz) {
	const float* cursor = section;
	int iy;

	if (IMPosnZWT(in_im, iz, iw, it) != 0) {
	    report_error(
		"unable to position initial guess file at z, wavelength, and time indices, %d, %d, %d",
		iz,
		iw,
		it
	    );
	    free(section);
	    return 1;
	}
	IMRdSec(in_im, section);
	for (iy = 0; iy < pfp->ypadded; ++iy) {
	    int ix;

	    for (ix = 0; ix < pfp->xpadded; ++ix) {
		values[ix] = cursor[ix] / scale;
	    }
	    for (ix = pfp->xpadded; ix < nxv; ++ix) {
		values[ix] = 0.0;
	    }
	    values += nxv;
	    cursor += pfp->xpadded;
	}
    }
    if (needs_z_fill) {
	dzero(
	    nxv * (size_t) pfp->ypadded *
	        (size_t) (pfp->zpadded - 1 - zbounds[1]),
	    values
	);
    }

    free(section);
    return 0;
}


static
int
fetch_mask(
    const FrontEndParameters* pfp,
    int in_im,
    int iw,
    int it,
    int iz,
    double* values
)
{
    float* section = (float*) malloc(
	pfp->xpadded * (size_t) pfp->ypadded * sizeof(float)
    );
    int nxv = 2 * (pfp->xpadded / 2 + 1);
    int nz, nw, nt, intlv;
    int zbounds[2];

    if (section == 0) {
	report_error(
	    "memory allocation failed for mask buffer: %s", strerror(errno)
	);
	return 1;
    }

    IMRtZWT(in_im, &nz, &nw, &nt, &intlv);
    if (iz == -1) {
	zbounds[0] = 0;
	zbounds[1] = pfp->zpadded - 1;
    } else {
	if (nz == 1) {
	    zbounds[0] = 0;
	    zbounds[1] = 0;
	} else {
	    zbounds[0] = iz;
	    zbounds[1] = iz;
	}
    }
    if (nw == 1) {
	iw = 0;
    }
    if (nt == 1) {
	it = 0;
    }
    for (iz = zbounds[0]; iz <= zbounds[1]; ++iz) {
	const float* cursor;
	int iy;

	if (IMPosnZWT(in_im, iz, iw, it) != 0) {
	    report_error(
		"unable to position input spatial mask file at z, wavelength, and time indices, %d, %d, %d",
		iz,
		iw,
		it
	    );
	    free(section);
	    return 1;
	}
	IMRdSec(in_im, section);
	cursor = section;
	for (iy = 0; iy < pfp->ypadded; ++iy) {
	    int ix;

	    for (ix = 0; ix < pfp->xpadded; ++ix) {
		values[ix] = cursor[ix];
	    }
	    values += nxv;
	    cursor += pfp->xpadded;
	}
    }

    free(section);
    return 0;
}


static
int
fetch_mcf(const FrontEndParameters* pfp, int in_im, int iw, double* values)
{
    int nxp = 2 * (pfp->xpadded / 2 + 1);
    float* section =
	(float*) malloc(nxp * (size_t) pfp->ypadded * sizeof(float));
    int iz;

    if (section == 0) {
	report_error(
	    "memory allocation failed for cone filter buffer: %s",
	    strerror(errno)
	);
	return 1;
    }
    for (iz = 0; iz < pfp->zpadded; ++iz) {
	const float* cursor;
	int iy;

	if (IMPosnZWT(in_im, iz, iw, 0) != 0) {
	    report_error(
		"unable to postion cone filter file at z, wavelength, and time indices, %d, %d, 0",
		iz,
		iw
	    );
	    free(section);
	    return 1;
	}
	IMRdSec(in_im, section);
	cursor = section;
	for (iy = 0; iy < pfp->ypadded; ++iy) {
	    int ix;

	    for (ix = 0; ix < nxp; ++ix) {
		values[ix] = cursor[ix];
	    }
	    values += nxp;
	    cursor += nxp;
	}
    }

    free(section);
    return 0;
}


static
double*
fetch_measurement(
    const FrontEndParameters* pfp, int in_im, int iw, int it, int iz
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
    int zbounds[2];

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
    if (iz == -1) {
	zbounds[0] = pfp->zbounds[0];
	zbounds[1] = pfp->zbounds[1];
    } else {
	zbounds[0] = iz;
	zbounds[1] = iz;
    }
    for (iz = 0; iz < zbounds[1] - zbounds[0] + 1; ++iz) {
	const float* cursor = section;
	int iy;

	if (IMPosnZWT(in_im, zbounds[0] + iz, pfp->wavelength_indices[iw], pfp->tbounds[0] + it * pfp->tbounds[2]) != 0) {
	    report_error(
		"unable to position measurement file at z, wavelength, and time indices, %d, %d, %d",
		zbounds[0] + iz,
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
	zbounds[1] - zbounds[0] + 1,
	0,
	0,
	pfp->zpadded - (zbounds[1] - zbounds[0] + 1),
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
    int argc, char* argv[], Core2Params* pcp, FrontEndParameters* pfp
)
{
    AvailableMaskuPlugin* masku_plugins_list = 0;
    int masku_plugins_count = 0;
    int iarg, had_error, background_set;

    pcp->lamratio[0] = 0.0;
    pcp->lamratio[1] = 1.0;
    pcp->lamf = 1.0;
    pcp->lampos = 1.0;
    pcp->lampc = 100.0;
    pcp->laml2 = 0.0;
    pcp->lammask = 0.0;
    pcp->alpha[0] = 10000.0;
    pcp->alpha[1] = 10000.0;
    /*
     * core2 currently takes values in beta but never uses them in
     * computations so do not expose it via the command-line interface.
     */
    pcp->beta[0] = 1.0;
    pcp->beta[1] = 1.0;
    pcp->tol = 0.0001;
    pcp->omega = 0.8;
    pcp->niter = 50;
    pcp->nthread = get_default_thread_count();
    pcp->uses_old_line_search = 0;
    pcp->regularization_type = 1;
    pcp->outputs_rfactor = 0;
    pcp->skips_initial_guess = 0;

    pfp->measurement_file_name = 0;
    pfp->result_file_name = 0;
    pfp->ctf_file_name = 0;
    pfp->in_filter_name = 0;
    pfp->out_filter_name = 0;
    pfp->plot_file_name = 0;
    pfp->in_mask_name = 0;
    pfp->out_mask_name = 0;
    pfp->masku_plugin_name = 0;
    pfp->masku_intfc_name = 0;
    pfp->guess_name = 0;
    pfp->masku_plug = 0;
    pfp->masku_closure.attr = 0;
    pfp->masku_closure.options = 0;
    pfp->cuth = .001;
    pfp->lamsmooth = 100.0;
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
    pfp->as_2d = -1;
    pfp->as_2d_ctf = 0;
    pfp->include_padding = 0;
    pfp->log_level = LOG_INFO;

    if (argc < 4) {
	print_usage();
	return 1;
    }

    get_masku_plugin_list(
	"core2_masku_plugins.def",
	(const char*) 0,
	&masku_plugins_list,
	&masku_plugins_count
    );
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
	    case '2':
		if (argv[iarg][2] == 'D' && argv[iarg][3] == '\0') {
		    pfp->as_2d = 1;
		} else {
		    recognized = 0;
		}
		break;

	    case 'a':
		if (strcmp(argv[iarg] + 2, "bs_background") == 0) {
		    pfp->absolute_background = 1;
		} else if (strncmp(argv[iarg] + 2, "lpha=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->alpha[0]) == 1) {
			pcp->alpha[1] = pcp->alpha[0];
		    } else {
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
		if (strcmp(argv[iarg] + 2, "tf2d") == 0) {
		    pfp->as_2d_ctf = 1;
		} else if (strncmp(argv[iarg] + 2, "uth=", 4) == 0) {
		    if (extract_doubles(argv[iarg] + 6, 1, &pfp->cuth) != 1) {
			report_error(
			    "expected one floating-point value for -cuth"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'g':
		if (strncmp(argv[iarg] + 2, "uess=", 5) == 0) {
		    pfp->guess_name = argv[iarg] + 7;
		} else {
		    recognized = 0;
		}
		break;

	    case 'i':
		if (strncmp(argv[iarg] + 2, "mask=", 5) == 0) {
		    pfp->in_mask_name = argv[iarg] + 7;
		} else if (strncmp(argv[iarg] + 2, "mcf=", 4) == 0) {
		    pfp->in_filter_name = argv[iarg] + 6;
		} else if (strcmp(argv[iarg] + 2, "nclude_padding") == 0) {
		    pfp->include_padding = 1;
		} else {
		    recognized = 0;
		}
		break;

	    case 'l':
		if (strncmp(argv[iarg] + 2, "amratio=", 8) == 0) {
		    if (extract_doubles(argv[iarg] + 10, 2, pcp->lamratio) != 2 ||
			pcp->lamratio[0] < 0.0 || pcp->lamratio[0] > 1.0 ||
			pcp->lamratio[1] < 0.0) {
			report_error(
			    "expected two non-negative floating-point values for -lamratio; first must be less than or equal to one"
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
		} else if (strncmp(argv[iarg] + 2, "aml2=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->laml2) != 1 ||
			pcp->laml2 < 0.0) {
			report_error(
			    "expected one non-negative floating-point value for -laml2"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "ammask=", 7) == 0) {
		    if (extract_doubles(argv[iarg] + 9, 1, &pcp->lammask) != 1 ||
			pcp->lammask < 0.0) {
			report_error(
			    "expected on non-negative floating-point value for -lammask"
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
		    if (extract_doubles(argv[iarg] + 11, 1, &pfp->lamsmooth) != 1) {
			report_error(
			    "expected one floating-point value for -lamsmooth"
			);
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "inesearch=", 10) == 0) {
		    if (strcmp(argv[iarg] + 12, "2013") == 0) {
			pcp->uses_old_line_search = 1;
		    } else if (strcmp(argv[iarg] + 12, "2014") == 0) {
			pcp->uses_old_line_search = 0;
		    } else {
			report_error("unrecognized method for -linesearch");
			had_error = 1;
		    }
		} else if (strncmp(argv[iarg] + 2, "ogging=", 7) == 0) {
		    if (strcmp(argv[iarg] + 9, "error") == 0) {
			pfp->log_level = LOG_ERROR;
		    } else if (strcmp(argv[iarg] + 9, "warning") == 0) {
			pfp->log_level = LOG_WARN;
		    } else if (strcmp(argv[iarg] + 9, "info") == 0) {
			pfp->log_level = LOG_INFO;
		    } else if (strcmp(argv[iarg] + 9, "progress") == 0) {
			pfp->log_level = LOG_PROGRESS;
		    } else if (strcmp(argv[iarg] + 9, "debug") == 0) {
			pfp->log_level = LOG_DEBUG;
		    } else {
			report_error("unrecognized value for -logging");
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
			ncycl >= 0) {
			pcp->niter = ncycl;
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
			if (extract_ints(argv[iarg] + 4, 1, &pcp->nthread) != 1 ||
			    pcp->nthread < 1) {
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
		if (strncmp(argv[iarg] + 2, "mask=", 5) == 0) {
		    pfp->out_mask_name = argv[iarg] + 7;
		} else if (strncmp(argv[iarg] + 2, "mcf=", 4) == 0) {
		    pfp->out_filter_name = argv[iarg] + 6;
		} else if (strncmp(argv[iarg] + 2, "plotfile=", 9) == 0) {
		    pfp->plot_file_name = argv[iarg] + 11;
		} else if (strncmp(argv[iarg] + 2, "mega=", 5) == 0) {
		    if (extract_doubles(argv[iarg] + 7, 1, &pcp->omega) != 1) {
			report_error(
			    "expected one floating-point value for -omega"
			);
			had_error = 1;
		    }
		} else {
		    recognized = 0;
		}
		break;

	    case 'r':
		if (strncmp(argv[iarg] + 2, "egtype=", 7) == 0) {
		    if (strcmp(argv[iarg] + 9, "ma") == 0) {
			pcp->regularization_type = 1;
		    } else if (strcmp(argv[iarg] + 9, "tv2") == 0) {
			pcp->regularization_type = 2;
		    } else {
			report_error("-regtype expects either ma or tv2");
			had_error = 1;
		    }
		} else if (strcmp(argv[iarg] + 2, "factor") == 0) {
		    pcp->outputs_rfactor = 1;
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

	    case 'u':
		if (strncmp(argv[iarg] + 2, "mask=", 5) == 0) {
		    if (parse_masku_command_line(
			    pfp,
			    argv[iarg] + 7,
			    masku_plugins_list,
			    masku_plugins_count) != 0) {
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

    if (pcp->lammask > 0.0 && pcp->uses_old_line_search) {
	report_error(
	    "warning: positive value set with -lammask but using 2013 line search; the spatial constraint will be ignored"
	);
    }

    if (masku_plugins_list != 0) {
	destroy_masku_plugin_list(masku_plugins_list);
    }

    if (had_error) {
	print_usage();
	return 1;
    }

    return 0;
}


static
int
parse_masku_command_line(
    FrontEndParameters* pfp,
    char* arg,
    AvailableMaskuPlugin* plugin_list,
    int list_count
)
{
    char* opts = 0;
    char* next_comma;
    AvailableMaskuPlugin* p_plugin;

    pfp->masku_plugin_name = arg;
    next_comma = strchr(pfp->masku_plugin_name, ',');
    if (next_comma != 0) {
	*next_comma = '\0';
	pfp->masku_intfc_name = next_comma + 1;
	next_comma = strchr(pfp->masku_intfc_name, ',');
	if (next_comma != 0) {
	    *next_comma = '\0';
	    if (*(next_comma + 1) != '\0') {
		opts = next_comma + 1;
	    }
	}
	if (pfp->masku_intfc_name[0] == '\0') {
	    pfp->masku_intfc_name = 0;
	}
    }
    if (pfp->masku_plugin_name[0] == '\0') {
	report_error("no name specified for mask update plugin");
	return 1;
    }

    p_plugin = search_masku_plugin_list_by_short_name(
	plugin_list, list_count, pfp->masku_plugin_name
    );
    if (p_plugin != 0) {
	if (pfp->masku_intfc_name != 0) {
	    report_error(
		"warning: short name specified for mask update plugin; ignoring the interface name specified: %s",
		pfp->masku_intfc_name
	    );
	    pfp->masku_intfc_name = 0;
	}
	pfp->masku_plug = load_masku_plugin_from_list(p_plugin);
    } else {
	pfp->masku_plug = load_masku_plugin_directly(
	    pfp->masku_plugin_name, pfp->masku_intfc_name
	);
    }
    if (pfp->masku_plug == 0) {
	report_error(
	    "failed to load mask update plugin %s", pfp->masku_plugin_name
	);
	return 1;
    }
    pfp->masku_closure.attr = get_masku_plugin_attr(pfp->masku_plug);
    if (pfp->masku_closure.attr == 0) {
	report_error(
	    "mask update plugin, %s, did not return pointer to attributes",
	    pfp->masku_plugin_name
	);
	return 1;
    }

    pfp->masku_closure.attr->error_handler = &print_for_masku;
    pfp->masku_closure.attr->error_handler_arg = (void*) stderr;
    if (pfp->masku_closure.attr->options_creator != 0) {
	pfp->masku_closure.options =
	    (*pfp->masku_closure.attr->options_creator)(
		pfp->masku_closure.attr
	    );
	if (pfp->masku_closure.options == 0) {
	    report_error(
		"mask update plugin, %s, failed to allocate options",
		pfp->masku_plugin_name
	    );
	    return 1;
	}
	if (opts != 0) {
	    if (pfp->masku_closure.attr->arg_parser != 0) {
		if ((*pfp->masku_closure.attr->arg_parser)(
			pfp->masku_closure.attr,
			opts,
			pfp->masku_closure.options) != 0) {
		    if (pfp->masku_closure.attr->usage_printer != 0) {
			(*pfp->masku_closure.attr->usage_printer)(
			    pfp->masku_closure.attr,
			    &print_for_masku,
			    (void*) stderr
			);
		    }
		    return 1;
		}
	    } else {
		report_error(
		    "ignored options, %s, for mask update; plugin, %s, does not supply argument parser",
		    opts,
		    pfp->masku_plugin_name
		);
	    }
	}
    } else {
	if (opts != 0) {
	    report_error(
		"ignored options, %s, for mask update; plugin, %s, does not provide an options allocator",
		opts,
		pfp->masku_plugin_name
	    );
	}
    }

    return 0;
}


static
void
print_for_masku(void* completion, const char* fmt, ...)
{
    FILE* fo = (FILE*) completion;
    va_list ap;

    (void) fprintf(fo, "    ");
    va_start(ap, fmt);
    (void) vfprintf(fo, fmt, ap);
    va_end(ap);
    (void) fprintf(fo, "\n");
}


static
void
print_parameters(const Core2Params* pcp, FrontEndParameters* pfp)
{
    int i;

    (void) printf("Measurement file:     %s\n", pfp->measurement_file_name);
    (void) printf("Result file:          %s\n", pfp->result_file_name);
    (void) printf("CTF file:             %s\n", pfp->ctf_file_name);
    (void) printf(
	"Initial guess:        %s\n",
	(pfp->guess_name == 0 || pfp->guess_name[0] == '\0') ?
	"internally computed" : pfp->guess_name
    );
    (void) printf(
	"Cone filter:          %s\n",
	(pfp->in_filter_name == 0 || pfp->in_filter_name[0] == '\0') ?
	"internally computed" : pfp->in_filter_name
    );
    (void) printf(
	"Cone filter output:   %s\n",
	(pfp->out_filter_name == 0 || pfp->out_filter_name[0] == '\0') ?
	"not written" : pfp->out_filter_name
    );
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
    (void) printf("Dimensionality:       %dD\n", (pfp->as_2d) ? 2 : 3);
    (void) printf(
	"Assume 2D CTF:        %s\n", (pfp->as_2d_ctf) ? "yes" : "no"
    );
    (void) printf(
	"Normalize CTF:        %s\n", (pfp->normalize_otf) ? "yes" : "no"
    );
    (void) printf("Nonlinearity:         %.6g\n", pcp->alpha[0]);
    (void) printf("Smoothing:            %.6g\n", pcp->lamf);
    (void) printf("Missing cone weight:  %.6g\n", pcp->lampc);
    (void) printf("Positivity:           %.6g\n", pcp->lampos);
    (void) printf("Auxiliary smoothing:  %.6g\n", pfp->lamsmooth);
    (void) printf("Extra entropy scalar: %.6g\n", pcp->laml2);
    (void) printf("Numerical aperture:   %.6g\n", pfp->na);
    (void) printf("Imm. ref. index:      %.6g\n", pfp->nimm);
    (void) printf("High cutoff:          %.6g\n", pfp->cuth);
    (void) printf("Iteration count:      %ld\n", pcp->niter);
    (void) printf("Padded z size:        %d\n", pfp->zpadded);
    if (pcp->uses_old_line_search) {
	(void) printf("Line search:          2013 algorithm\n");
	(void) printf(
	    "Ratio:                %.6g %.6g\n",
	    pcp->lamratio[0],
	    pcp->lamratio[1]
	);
	(void) printf("Relaxation:           %.6g\n", pcp->omega);
	(void) printf("Tolerance:            %.6g\n", pcp->tol);
    } else {
	(void) printf("Line search:          2014 algorithm\n");
	if (pcp->regularization_type == 1) {
	    (void) printf("Ratio:                %.6g\n", pcp->lamratio[1]);
	} else {
	    (void) printf("Ratio:                %.6g\n", 0.0);
	}
	(void) printf(
	    "Regularization type:  %s\n",
	    (pcp->regularization_type == 1) ? "Muthuvel Arigovindan" : "Other"
	);
	(void) printf(
	    "Outputs R-factor:     %s\n", (pcp->outputs_rfactor) ? "yes" : "no"
	);
    }
    if (! pfp->as_2d) {
	(void) printf(
	    "Output padded:        %s\n", (pfp->include_padding) ? "yes" : "no"
	);
    }
    (void) printf("Thread count:         %d\n", pcp->nthread);
    (void) printf("Log level:            ");
    switch (pfp->log_level) {
    case LOG_ERROR:
	(void) printf("error");
	break;

    case LOG_WARN:
	(void) printf("warning");
	break;

    case LOG_INFO:
	(void) printf("info");
	break;

    case LOG_PROGRESS:
	(void) printf("progress");
	break;

    case LOG_DEBUG:
	(void) printf("debug");
	break;
    }
    (void) printf("\n");
    (void) printf("Spatial mask weight:  %.6g\n", pcp->lammask);
    if (pcp->lammask > 0.0) {
	(void) printf(
	    "Input mask:           %s\n",
	    (pfp->in_mask_name != 0) ? pfp->in_mask_name : "all zeros"
	);
	(void) printf(
	    "Output mask:          %s\n",
	    (pfp->out_mask_name != 0) ? pfp->out_mask_name : "not written"
	);
	(void) printf(
	    "Mask update plugin:   %s\n",
	    (pfp->masku_plugin_name != 0) ? pfp->masku_plugin_name : "none"
	);
	(void) printf(
	    "Mask upd. interface:  %s\n",
	    (pfp->masku_intfc_name != 0) ? pfp->masku_intfc_name : "default"
	);
	if (pfp->masku_closure.attr != 0 &&
	    pfp->masku_closure.attr->options_summarizer != 0 &&
	    pfp->masku_closure.options != 0) {
	    (void) printf("Mask update options:\n");
	    (*pfp->masku_closure.attr->options_summarizer)(
		pfp->masku_closure.attr,
		pfp->masku_closure.options,
		&print_for_masku,
		(void*) stdout
	    );
	}
    }
}


static
void
print_usage(void)
{
    (void) fprintf(
	stderr,
	"usage: core2_decon measurement result CTF \\\n"
	"    [-x=start[:end]] [-y=start[:end]] [-z=start[:end]] \\\n"
	"    [-t=start[:end[:step]]] [-w=w1[:w2...]] \\\n"
	"    [-2D] [-abs_background] [-alpha=a] [-ctf2d] [-cuth=t] \\\n"
	"    [-guess=filename] [-imcf=filename] [-imask=filename] \\\n"
	"    [-include_padding] [-lamratio=r1:r2] [-lamf=f] [-laml2=f] \\\n"
	"    [-lammask=f] [-lampc=p] [-lampos=p] [-lamsmooth=s] \\\n"
	"    [-linesearch=2013|2014] \\\n"
	"    [-logging=error|warning|info|progress|debug] [-na=v] [-ncycl=n] \\\n"
	"    [-nimm=n] [-nonorm_otf] [-np=n] [-nzpad=n] [-omask=filename] \\\n"
	"    [-omcf=filename] [-omega=o] [-oplotfile=pltname] \\\n"
	"    [-sub=b1[:b2...]] [-regtype=t] [-rfactor] [-tol=t] \\\n"
	"    [-umask=plugin_name[,[interface_name],...]]\n"
    );
}


static
void
report_error(const char* fmt, ...)
{
    va_list ap;

    (void) fprintf(stderr, "core2_decon: ");
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
setup_guess(const FrontEndParameters* pfp, int* p_im_stream)
{
    int nxyz[3], mxyz[3], fmt, nz, nw, nt, intlv;
    float fmin, fmax, fmean;
    int status;

    if (pfp->guess_name == 0 || pfp->guess_name[0] == '\0') {
	/* Will use an internally computed initial guess. */
	return 0;
    }

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for initial guess"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->guess_name, "ro") != 0) {
	report_error(
	    "unable to open initial guess, %s, for reading", pfp->guess_name
	);
	*p_im_stream = -1;
	return 1;
    }
    IMRdHdr(*p_im_stream, nxyz, mxyz, &fmt, &fmin, &fmax, &fmean);
    IMRtZWT(*p_im_stream, &nz, &nw, &nt, &intlv);
    status = 0;
    if (nxyz[0] != pfp->xpadded || nxyz[1] != pfp->ypadded) {
	report_error(
	    "initial guess must have same x and y dimensions as the padded data to process"
	);
	status = 1;
    }
    if (pfp->as_2d) {
	if (nz != pfp->zbounds[1] - pfp->zbounds[0] + 1) {
	    report_error(
		"z dimension of initial guess must have the same size as the data to process"
	    );
	    status = 1;
	}
    } else {
	if (nz != pfp->zbounds[1] - pfp->zbounds[0] + 1 &&
	    nz != pfp->zpadded) {
	    report_error(
		"z dimension of initial guess must match the z size of the data to process or the padded z size"
	    );
	    status = 1;
	}
    }
    if (fmt == IW_COMPLEX || fmt == IW_COMPLEX_SHORT) {
	report_error("initial guess must be real-valued");
	status = 1;
    }
    return status;
}


static
int
setup_imcf(FrontEndParameters* pfp, int* p_im_stream)
{
    int nxyz[3], mxyz[3], fmt, nz, nw, nt, intlv;
    float fmin, fmax, fmean;

    if (pfp->in_filter_name == 0 || pfp->in_filter_name[0] == '\0') {
	/*
	 * Will use the filter computed from the transfer function if the cone
	 * filter is needed.
	 */
	return 0;
    }

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error("internal error: unable to find free IM stream for cone filter"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->in_filter_name, "ro") != 0) {
	report_error(
	    "unable to open cone filter, %s, for reading", pfp->in_filter_name
	);
	*p_im_stream = -1;
	return 1;
    }
    IMRdHdr(*p_im_stream, nxyz, mxyz, &fmt, &fmin, &fmax, &fmean);
    IMRtZWT(*p_im_stream, &nz, &nw, &nt, &intlv);
    if (nxyz[0] != pfp->xpadded / 2 + 1 || nxyz[1] != pfp->ypadded ||
	nz != pfp->zpadded || (fmt != IW_COMPLEX && fmt != IW_COMPLEX_SHORT)) {
	report_error(
	    "cone filter must be complex-valued with the same dimensions as the padded data to process."
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

    if (pfp->as_2d == -1) {
	pfp->as_2d = (pfp->zbounds[1] > pfp->zbounds[0]) ? 0 : 1;
    }
    if (! pfp->as_2d && pfp->as_2d_ctf) {
	report_error("forcing (-ctf2d) the CTF to be interpreted as 2D while processing is done as 3D");
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
	if (! pfp->as_2d) {
	    while (pfp->zpadded <
		       1.3f * (pfp->zbounds[1] - pfp->zbounds[0] + 1) &&
		   pfp->zpadded < INT_MAX / 2) {
		pfp->zpadded += pfp->zpadded;
	    }
	}
    } else if (pfp->zpadded > 1 && pfp->as_2d) {
	pfp->zpadded = 1;
    } else if (! pfp->as_2d &&
	       pfp->zpadded < pfp->zbounds[1] - pfp->zbounds[0] + 1) {
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
setup_mask_input(FrontEndParameters* pfp, int* p_im_stream)
{
    int nxyz[3], mxyz[3], fmt, nz, nw, nt, intlv;
    float fmin, fmax, fmean;
    int status, nw_in;

    if (pfp->in_mask_name == 0 || pfp->in_mask_name[0] == '\0') {
	return 0;
    }

    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error("internal errror: unable to find free IM stream for spatial mask input");
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->in_mask_name, "ro") != 0) {
	report_error(
	    "unable to open spatial mask, %s, for reading", pfp->in_mask_name
	);
	*p_im_stream = -1;
	return 1;
    }
    IMRdHdr(*p_im_stream, nxyz, mxyz, &fmt, &fmin, &fmax, &fmean);
    IMRtZWT(*p_im_stream, &nz, &nw, &nt, &intlv);
    status = 0;
    if (fmt == IW_COMPLEX || fmt == IW_COMPLEX_SHORT) {
	report_error("spatial mask must be real valued");
	status = 1;
    }
    if (nxyz[0] != pfp->xpadded || nxyz[1] != pfp->ypadded) {
	report_error(
	    "x and y dimensions of spatial mask do not match x and y dimensions of padded data"
	);
	status = 1;
    }
    if (pfp->zpadded != nz &&
	(pfp->zpadded == 1 && nz != pfp->zbounds[1] - pfp->zbounds[0] + 1)) {
	report_error(
	    "z dimension of spatial mask is not consistent with padded data"
	);
	status = 1;
    }
    if (nt != 1 &&
	nt != (pfp->tbounds[1] - pfp->tbounds[0]) / pfp->tbounds[2] + 1) {
	report_error(
	    "time dimension of spatial mask is not consistent with padded data"
	);
	status = 1;
    }
    nw_in = 0;
    while (pfp->wavelength_indices[nw_in] != -1) {
	++nw_in;
    }
    if (nw != 1 && nw != nw_in) {
	report_error(
	    "wavelength dimension of spatial mask in not consistent with padded data"
	);
	status = 1;
    }

    return status;
}


static
int
setup_mask_output(FrontEndParameters* pfp, int template_im, int* p_im_stream)
{
    int nz, nt, nw;
    int nxyz[3], mxyz[3], nxyzst[3], itstart;
    float wavelengths_in[IW_MAX_WAVE], wavelengths_out[IW_MAX_WAVE];
    float delta[3], xyz[3];

    if (pfp->out_mask_name == 0 || pfp->out_mask_name[0] == '\0') {
	return 0;
    }
    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for spatial mask output"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->out_mask_name, "new") != 0) {
	report_error("unable to open %s for writing", pfp->out_mask_name);
	*p_im_stream = -1;
	return 1;
    }

    IMTrHdr(*p_im_stream, template_im);
    IMRtSiz(*p_im_stream, nxyz, mxyz, nxyzst);
    IMRtWav(*p_im_stream, &nw, wavelengths_in);
    nz = (pfp->zpadded == 1) ?
	pfp->zbounds[1] - pfp->zbounds[0] + 1 : pfp->zpadded;
    nw = 0;
    while (pfp->wavelength_indices[nw] != -1) {
	wavelengths_out[nw] = wavelengths_in[pfp->wavelength_indices[nw]];
	++nw;
    }
    nt = (pfp->tbounds[1] - pfp->tbounds[0]) / pfp->tbounds[2] + 1;
    nxyz[0] = pfp->xpadded;
    nxyz[1] = pfp->ypadded;
    nxyz[2] = nz * nw * nt;
    nxyzst[0] += pfp->xbounds[0];
    nxyzst[1] += pfp->ybounds[0];
    nxyzst[2] += pfp->zbounds[0];
    IMAlSiz(*p_im_stream, nxyz, nxyzst);
    IMAlZWT(*p_im_stream, nz, nw, nt, ZTW_SEQUENCE);
    IMAlWav(*p_im_stream, nw, wavelengths_out);
    IMAlExHdrSize(*p_im_stream, 0, 0, nxyz[2]);

    IMAlMode(*p_im_stream, IW_FLOAT);    

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
setup_omcf(FrontEndParameters* pfp, int template_im, int* p_im_stream)
{
    int nw;
    int nxyz[3], nxyzst[3];
    float wavelengths_in[IW_MAX_WAVE], wavelengths_out[IW_MAX_WAVE];

    if (pfp->out_filter_name == 0 || pfp->out_filter_name[0] == '\0') {
	return 0;
    }
    *p_im_stream = find_free_stream();
    if (*p_im_stream == -1) {
	report_error(
	    "internal error: unable to find free IM stream for cone filter output"
	);
	return 1;
    }

    if (IMOpen(*p_im_stream, pfp->out_filter_name, "new") != 0) {
	report_error("unable to open %s for writing", pfp->out_filter_name);
	*p_im_stream = -1;
	return 1;
    }

    IMRtWav(template_im, &nw, wavelengths_in);
    IMTrHdr(*p_im_stream, template_im);
    nw = 0;
    while (pfp->wavelength_indices[nw] != -1) {
	wavelengths_out[nw] = wavelengths_in[pfp->wavelength_indices[nw]];
	++nw;
    }
    nxyz[0] = pfp->xpadded / 2 + 1;
    nxyz[1] = pfp->ypadded;
    nxyz[2] = pfp->zpadded * nw;
    nxyzst[0] = 0;
    nxyzst[1] = 0;
    nxyzst[2] = 0;
    IMAlSiz(*p_im_stream, nxyz, nxyzst);
    IMAlZWT(*p_im_stream, pfp->zpadded, nw, 1, ZTW_SEQUENCE);
    IMAlWav(*p_im_stream, nw, wavelengths_out);
    IMAlMode(*p_im_stream, IW_COMPLEX);
    IMAlExHdrSize(*p_im_stream, 0, 0, nxyz[2]);
    IMAlOrig(*p_im_stream, 0.0f, 0.0f, 0.0f);
    IMAlTSt(*p_im_stream, 0);
    IMAlRes(*p_im_stream, 1, 2);

    return 0;
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
	(! pfp->as_2d && pfp->include_padding) ?
	    pfp->zbounds[0] + pfp->zpadded - 1 : pfp->zbounds[1],
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
update_new_plot(
    int iter, const ResReporterCosts* p_costs, double rfactor, void* closure
)
{
    ResReporterData* prrd = (ResReporterData*) closure;

    assert(! prrd->uses_old_style_results);
    if (prrd->updates_stdout) {
	(void) printf(
	    "%6d %13.7e %13.7e %13.7e %13.7e",
	    iter,
	    p_costs->total,
	    p_costs->data_fidelity,
	    p_costs->reg,
	    p_costs->pos
	);
	if (prrd->outputs_rfactor) {
	    (void) printf(" %13.7e\n", rfactor);
	} else {
	    (void) printf("\n");
	}
	(void) fflush(stdout);
    }
    if (prrd->status != 0) {
	/*
	 * An earlier call to update_new_plot() for the same plot
	 * encountered a memory allocation error so don't do anything in
	 * this call.
	 */
	return;
    }
    if (prrd->ncount == prrd->nalloc) {
	if (prrd->ncount < LONG_MAX / 2 &&
	    prrd->ncount < SIZE_MAX / (2 * sizeof(NewStylePt))) {
	    long new_count = (prrd->ncount > 0) ?
		prrd->ncount + prrd->ncount : 16;
	    NewStylePt* values;

	    values = (NewStylePt*) malloc(new_count * sizeof(NewStylePt));
	    if (values != 0) {
		if (prrd->ncount > 0) {
		    (void) memcpy(
			values,
			prrd->new_values,
			prrd->ncount * sizeof(NewStylePt)
		    );
		}
		if (prrd->new_values != 0) {
		    free(prrd->new_values);
		}
		prrd->new_values = values;
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
    prrd->new_values[prrd->ncount].cost = p_costs->total;
    prrd->new_values[prrd->ncount].cost_fidel = p_costs->data_fidelity;
    prrd->new_values[prrd->ncount].cost_reg = p_costs->reg;
    prrd->new_values[prrd->ncount].cost_pos = p_costs->pos;
    prrd->new_values[prrd->ncount].cost_laml2 = p_costs->laml2;
    prrd->new_values[prrd->ncount].cost_spa = p_costs->spa;
    prrd->new_values[prrd->ncount].rfactor = rfactor;
    prrd->new_values[prrd->ncount].iteration = iter;
    ++prrd->ncount;
}


static
void
update_old_plot(long iiters, long iitert, double res, void* closure)
{
    ResReporterData* prrd = (ResReporterData*) closure;
    long last_successful_iteration = prrd->last_successful_iteration;

    assert(prrd->uses_old_style_results);
    if (iitert == 0 || iiters != prrd->last_successful_iteration) {
	prrd->last_successful_res = res;
	prrd->last_successful_iteration = iiters;
    }
    if (prrd->updates_stdout) {
	(void) printf("%6ld %6ld %13.7e\n", iiters, iitert, res);
	(void) fflush(stdout);
    }
    if (prrd->fo != 0) {
	if (prrd->status != 0) {
	    /*
	     * An earlier call to update_old_plot() for the same plot
	     * encountered a memory allocation error so don't do anything in
	     * this call.
	     */
	    return;
	}
	if (prrd->ncount == prrd->nalloc) {
	    if (prrd->ncount < LONG_MAX / 2 &&
		prrd->ncount < SIZE_MAX / (2 * sizeof(OldStylePt))) {
		long new_count = (prrd->ncount > 0) ?
		    prrd->ncount + prrd->ncount : 16;
		OldStylePt* values;

		values = (OldStylePt*) malloc(new_count * sizeof(OldStylePt));
		if (values != 0) {
		    if (prrd->ncount > 0) {
			(void) memcpy(
			    values,
			    prrd->old_values,
			    prrd->ncount * sizeof(OldStylePt)
			);
		    }
		    if (prrd->old_values != 0) {
			free(prrd->old_values);
		    }
		    prrd->old_values = values;
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
	prrd->old_values[prrd->ncount].res = res;
	prrd->old_values[prrd->ncount].iteration =
	    (iitert == 0 || iiters != last_successful_iteration) ?
	    -iitert : iitert;
	++prrd->ncount;
    } else {
	prrd->ncount = iitert;
    }
}


static
int
write_mcf(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    const double* values,
    float mmm[3]
)
{
    int nxp = 2 * (pfp->xpadded / 2 + 1);
    float* section =
	(float*) malloc(nxp * (size_t) pfp->ypadded * sizeof(float));
    int iz;

    if (section == 0) {
	report_error(
	    "memory allocation for cone filter output buffer failed: %s",
	    strerror(errno)
	);
	return 1;
    }

    mmm[0] = FLT_MAX;
    mmm[1] = -FLT_MAX;
    mmm[2] = 0.0f;
    for (iz = 0; iz < pfp->zpadded; ++iz) {
	float sec_mms[3] = { FLT_MAX, -FLT_MAX, 0.0f };
	float* out = section;
	int iy;

	for (iy = 0; iy < pfp->ypadded; ++iy) {
	    float line_sum = 0.0f;
	    int ix;

	    for (ix = 0; ix < nxp; ix += 2) {
		out[ix] = values[ix];
		out[ix + 1] = values[ix + 1];
		if (sec_mms[0] > out[ix]) {
		    sec_mms[0] = out[ix];
		}
		if (sec_mms[1] < out[ix]) {
		    sec_mms[1] = out[ix];
		}
		line_sum += out[ix];
	    }
	    out += nxp;
	    values += nxp;
	    sec_mms[2] += line_sum;

	    if (IMPosnZWT(out_im, iz, iw, 0) != 0) {
		report_error(
		    "unable to position output cone filter file at z, wavelength, and time index, %d, %d, 0",
		    iz,
		    iw
		);
		free(section);
		return 1;
	    }
	    IMWrSec(out_im, section);
	}
	if (mmm[0] > sec_mms[0]) {
	    mmm[0] = sec_mms[0];
	}
	if (mmm[1] < sec_mms[1]) {
	    mmm[1] = sec_mms[1];
	}
	mmm[2] += sec_mms[2];
    }
    mmm[2] /= pfp->xpadded * (float) pfp->ypadded * (float) pfp->zpadded;

    free(section);
    return 0;
}


static
int
write_mask_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    int iz,
    const double* mask,
    float mmm[3]
)
{
    float vol_sum = 0.0f;
    float* section = (float*) malloc(
	pfp->xpadded * (size_t) pfp->ypadded * sizeof(float)
    );
    int nxv = 2 * (pfp->xpadded / 2 + 1);
    int zbounds[2];

    if (section == 0) {
	report_error(
	    "memory allocation for output buffer failed: %s", strerror(errno)
	);
	return 1;
    }

    if (iz == -1) {
	zbounds[0] = 0;
	zbounds[1] = pfp->zpadded - 1;
    } else {
	zbounds[0] = iz;
	zbounds[1] = iz;
    }
    for (iz = zbounds[0]; iz <= zbounds[1]; ++iz) {
	float sec_sum = 0.0f;
	float* out = section;
	int iy;

	for (iy = 0; iy < pfp->ypadded; ++iy) {
	    float line_sum = 0.0;
	    int ix;

	    for (ix = 0; ix < pfp->xpadded; ++ix) {
		out[ix] = mask[ix];
		if (mmm[0] > out[ix]) {
		    mmm[0] = out[ix];
		}
		if (mmm[1] < out[ix]) {
		    mmm[1] = out[ix];
		}
		line_sum += out[ix];
	    }
	    mask += nxv;
	    out += pfp->xpadded;
	    sec_sum += line_sum;
	}
	vol_sum += sec_sum;

	if (IMPosnZWT(out_im, iz, iw, it) != 0) {
	    free(section);
	    return 1;
	}
	IMWrSec(out_im, section);
    }

    mmm[2] += vol_sum /
	(pfp->xpadded * (size_t) pfp->ypadded *
	 (size_t) (zbounds[1] - zbounds[0] + 1));
    free(section);

    return 0;
}


static
int
write_result_and_collect_stats(
    const FrontEndParameters* pfp,
    int out_im,
    int iw,
    int it,
    int iz,
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
    int zbounds[2];

    if (section == 0) {
	report_error(
	    "memory allocation for output buffer failed: %s", strerror(errno)
	);
	return 1;
    }

    if (iz == -1) {
	zbounds[0] = 0;
	if (pfp->include_padding) {
	    zbounds[1] = pfp->zpadded - 1;
	} else {
	    zbounds[1] = pfp->zbounds[1] - pfp->zbounds[0];
	}
    } else {
	zbounds[0] = iz;
	zbounds[1] = iz;
    }
    for (iz = zbounds[0]; iz <= zbounds[1]; ++iz) {
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
	 (size_t) (zbounds[1] - zbounds[0] + 1)
        );
    free(section);

    return 0;
}
