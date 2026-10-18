/*
 * volumestruct2d.h -- RECONSTRUCTED (not part of the delivered source).
 * 2D counterpart of volumestruct.h; see the notes there.  The caller's pool
 * is 7 * (2*Nx + 2*Ny - 4) longs = 7 x (perimeter of the rectangle).
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_VOLUMESTRUCT2D_H
#define ERDECON_RECON_VOLUMESTRUCT2D_H

typedef struct Geometry2D {
    long Nx;
    long Ny;
    long NxFSkip;    /* 2 if Nx even, 1 if Nx odd */
    long LengthF;    /* (Nx + NxFSkip) * Ny */

    long NSurface;   /* 2*Nx + 2*Ny - 4: pixels on the rectangle perimeter */
    long* pool;
    long* boundary[7];
} Geometry2D;

void ConstructVolumeStruct2D(Geometry2D* G, long Nx, long Ny, long* mempool);
long VolumeStructPoolLength2D(long Nx, long Ny);

#endif /* include guard */
