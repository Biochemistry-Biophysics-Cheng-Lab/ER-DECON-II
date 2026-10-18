/*
 * derivatives2d.c -- RECONSTRUCTED.  2D counterpart of derivatives.c; the
 * kernel is shared (deriv_kernel.h) and run with Nz = 1.  Slots "  ", "xx",
 * "yy", "xy"; a z-slot is rejected.
 *
 * ERDecon reference build, 2026-09-21.
 */
#include "derivatives2d.h"
#include "deriv_kernel.h"

void
Derivative2D(
    const double* X,
    double* R,
    const Geometry2D* G,
    const char* dir,
    char mode,
    double* swpspace
)
{
    long nxpad = G->Nx + G->NxFSkip;
    enum deriv_slot slot = deriv_parse_slot(dir);

    if (slot == DERIV_ZZ || slot == DERIV_YZ || slot == DERIV_XZ) {
	(void) fprintf(stderr, "Derivative2D: z direction in 2D\n");
	abort();
    }
    if (X == R) {
	if (swpspace == 0) {
	    (void) fprintf(stderr,
		"Derivative2D: in-place call without scratch space\n");
	    abort();
	}
	(void) memcpy(swpspace, X, (size_t) G->LengthF * sizeof(double));
	X = swpspace;
    }
    deriv_apply(X, R, G->Nx, G->Ny, 1L, nxpad, slot, mode, "Derivative2D");
}
