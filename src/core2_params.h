#ifndef CORE2_PARAMS_H
#define CORE2_PARAMS_H

typedef struct Core2Params {
    /*
     * First value is only used for the old line search.
     * Second is used for both line search methods (but will be ignored for
     * the new line search if the regularization_type member is one) and sets
     * the scale factor for the signal squared regularization term relative to
     * the other regularization terms.  The second value should be
     * non-negative.
     */
    double lamratio[2];
    /*
     * Scales the contribution of the regularization terms to the cost\
     * function.  Should be non-negative.
     */
    double lamf;
    /*
     * Scales the positivity term relative to the other regularization terms.
     * Should be non-negative.
     */
    double lampos;
    /*
     * Scales the cone filter contribution when computing the approximate
     * inverse.  Should be non-negative.
     */
    double lampc;
    /*
     * In the C version of the algorithm, Muthuvel experimented with a term in
     * the cost function that was sum over r of Eg^2(r).  That term remained
     * in the C version (though after October 2014 it became a parameter and
     * defaulted to zero) but was abandoned in the Matlab code and PNAS paper.
     * This value is the scale factor for that term.  Use .0001 to match what
     * was used in old versions of the C code.  The value should be
     * non-negative.
     */
    double laml2;
    /*
     * Scales the contribution from the spatial constraint term.  The value
     * should be non-negative.  A value of zero disables the spatial
     * constraint and the pointer to the spatial mask data will not be used.
     */
    double lammask;
    /*
     * Holds one over epsilon from the PNAS paper.  For the old line search,
     * the first element is used with the signal squared term in the
     * regularization, and the second element is used with the the second
     * derivative terms.  For the new line search, only the first value is
     * used.  It is customary to set both elements to the same value.  The
     * values should be positive and the customary default is 10000.
     */
    double alpha[2];
    /*
     * The old line search algorithm accepts these two values but they're not
     * applied in any calculations.  The recommended default is to use one for
     * both.
     */
    double beta[2];
    /*
     * Sets the termination tolerance for the old line search.  The value
     * should be positive.  A recommended default is .0001.
     */
    double tol;
    /*
     * Sets the relaxation parameter for the old line search.  The value
     * should be positive.  A recommended default is 0.8.
     */
    double omega;
    /*
     * Sets the number of iterations to perform.  For the new line search,
     * values larger than INT_MAX are coerced to INT_MAX.
     */
    long niter;
    /* Sets the number of threads of execution to use.  Must be positive. */
    int nthread;
    /*
     * If nonzero, uses the line search described in the PNAS paper.  If zero,
     * use the line search Muthuvel implemented in 2014.
     */
    int uses_old_line_search;
    /*
     * This is only used for the new line search.  If this is one with that
     * line search, uses log(epsilon + sum over r Eg^2(r)) in the cost
     * function.  Otherwise, uses sqrt(epsilon + sum over r Eg^2(r)) to
     * give a total variation-like regularization.
     */
    int regularization_type;
    /*
     * If nonzero, the progress reports through the ResReporter callback
     * for the new line search with the new line search will include the
     * R-factor between the measurement and the blurred guess.
     */
    int outputs_rfactor;
    /*
     * If nonzero, the space pointed to be X is assumed to contain the values
     * for the initial guess.  If zero, the initial contents of X are not
     * used and are replaced with an internally calculated initial guess.
     */
    int skips_initial_guess;
} Core2Params;

#endif /* include guard */
