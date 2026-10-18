/*
 * volumestruct2d.c -- RECONSTRUCTED (not part of the delivered source).
 * See volumestruct2d.h.  ERDecon reference build, 2026-09-18.
 */
#include "volumestruct2d.h"

long
VolumeStructPoolLength2D(long Nx, long Ny)
{
    long perimeter;

    if (Nx >= 2 && Ny >= 2) {
	perimeter = 2L * Nx + 2L * Ny - 4L;
    } else {
	perimeter = Nx * Ny;
    }
    return 7L * perimeter;
}

void
ConstructVolumeStruct2D(Geometry2D* G, long Nx, long Ny, long* mempool)
{
    int k;

    G->Nx = Nx;
    G->Ny = Ny;
    G->NxFSkip = (Nx % 2L == 0L) ? 2L : 1L;
    G->LengthF = (Nx + G->NxFSkip) * Ny;
    G->NSurface = VolumeStructPoolLength2D(Nx, Ny) / 7L;
    G->pool = mempool;
    for (k = 0; k < 7; ++k) {
	G->boundary[k] = (mempool != 0) ? mempool + k * G->NSurface : 0;
    }
}
