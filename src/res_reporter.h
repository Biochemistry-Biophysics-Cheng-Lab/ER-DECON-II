#ifndef RESREPORTER_H
#define RESREPORTER_H

/*
 * The old line search algorithm and the new one adapted from the 4D
 * deconvolution report progress in different ways.  Pass instances of this
 * structure to define how to report for both cases.  The old line search
 * has been retained primarily to facilitate comparison with previous results.
 */
typedef struct ResReporterCosts {
    double total;
    double data_fidelity;
    double reg; /* guess^2 and derivative contribution */
    double pos; /* positivty contribution */
    double laml2; /* contribution from term scaled by laml2 */
    double spa; /* spatial constraint contribution */
} ResReporterCosts;

typedef struct ResReporter {
    void (*res_reporter)(
	int iter,
	const ResReporterCosts* p_costs,
	double rfactor,
	void* closure
    );
    void* res_closure;
    void (*old_res_reporter)(
	long iiters, long iitert, double res, void* res_reporter_closure
    );
    void* old_res_closure;
} ResReporter;

#endif /* include guard */
