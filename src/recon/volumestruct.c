/*
 * volumestruct.c -- RECONSTRUCTED (not part of the delivered source).
 * See volumestruct.h.  ERDecon reference build, 2026-09-18.
 */
#include "volumestruct.h"

long
VolumeStructPoolLength(long Nx, long Ny, long Nz)
{
    long surface;

    if (Nz >= 2) {
	surface = 2L * Nx * Ny + (Nz - 2L) * (2L * Nx + 2L * Ny - 4L);
    } else {
	/* single slice: every voxel is on the surface */
	surface = Nx * Ny;
    }
    return 7L * surface;
}

void
ConstructVolumeStruct(Geometry* G, long Nx, long Ny, long Nz, long* mempool)
{
    int k;

    G->Nx = Nx;
    G->Ny = Ny;
    G->Nz = Nz;
    G->NxFSkip = (Nx % 2L == 0L) ? 2L : 1L;
    G->LengthF = (Nx + G->NxFSkip) * Ny * Nz;
    G->NSurface = VolumeStructPoolLength(Nx, Ny, Nz) / 7L;
    G->pool = mempool;
    for (k = 0; k < 7; ++k) {
	G->boundary[k] = (mempool != 0) ? mempool + k * G->NSurface : 0;
    }
    /*
     * Boundary index lists intentionally left unfilled: their semantics
     * belong to Derivative() and are fixed together with it.
     */
}
