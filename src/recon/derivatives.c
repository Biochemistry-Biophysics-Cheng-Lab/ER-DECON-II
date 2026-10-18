/*
 * derivatives.c -- RECONSTRUCTED.  See derivatives.h and deriv_kernel.h.
 *
 * Periodic second differences of the PNAS 2013 SI (S1-S6), backward-aligned
 * mixed terms; 'f' forward, 'r' adjoint (accumulating), 's' squared stencil
 * in the adjoint direction.  The Geometry boundary lists (G->boundary[])
 * are not needed: wrap-around indexing is done with per-axis neighbour
 * tables, which is what those lists encoded in the original.
 *
 * ERDecon reference build, 2026-09-21.
 */
#include "derivatives.h"
#include "deriv_kernel.h"

void
Derivative(
    const double* X,
    double* R,
    const Geometry* G,
    const char* dir,
    char mode,
    double* swpspace
)
{
    long nxpad = G->Nx + G->NxFSkip;

    if (X == R) {
	if (swpspace == 0) {
	    (void) fprintf(stderr,
		"Derivative: in-place call without scratch space\n");
	    abort();
	}
	(void) memcpy(swpspace, X, (size_t) G->LengthF * sizeof(double));
	X = swpspace;
    }
    deriv_apply(X, R, G->Nx, G->Ny, G->Nz, nxpad,
		deriv_parse_slot(dir), mode, "Derivative");
}
