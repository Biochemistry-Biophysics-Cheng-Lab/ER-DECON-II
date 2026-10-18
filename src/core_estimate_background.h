#ifndef UCSF_MSG_CORE_ESTIMATE_BACKGROUND_H
#define UCSF_MSG_CORE_ESTIMATE_BACKGROUND_H

/*
 * Returns zero if successful; returns a nonzero value if not successful.
 * Sets *p_background to the estimated background level on a successful return.
 * Destroys the contents of measurement.
 * lamf and lamsmooth are only used as the product lamf * lamsmooth; they are
 * kept separate in the funtion interface only to follow the interface of
 * core().
 */
int core_estimate_background(
    double* measurement,
    const long size[3],
    double delfac,
    double lamf,
    double lamsmooth,
    int nthread,
    double* p_background
);

#endif /* include guard */
