/*
 * volumestruct.h -- RECONSTRUCTED (not part of the delivered source).
 *
 * The original Priism file defined the Geometry struct and
 * ConstructVolumeStruct().  What is known for certain from the delivered
 * code:
 *
 *   - Only five members are read anywhere: Nx, Ny, Nz, NxFSkip, LengthF.
 *   - NxFSkip is the padding that turns an Nx-wide real row into the
 *     2*(Nx/2+1)-wide half-complex row FFTW uses: 2 if Nx is even, 1 if odd.
 *   - LengthF = (Nx + NxFSkip) * Ny * Nz, the padded array length.
 *   - Callers hand in a long[] pool of exactly
 *         7 * (2*Nx*Ny + (Nz-2)*(2*Nx + 2*Ny - 4))
 *     elements = 7 x (number of surface voxels of the box).  The original
 *     filled it with seven boundary index lists whose only consumer was
 *     Derivative(); nothing we hold dereferences them.
 *
 * This reconstruction keeps the pool pointer and reserves the seven list
 * slots so that derivatives.c can be designed against them later, but it
 * does not (yet) fill them.  The boundary rule is decided in derivatives.c.
 *
 * ERDecon reference build, 2026-09-18.
 */
#ifndef ERDECON_RECON_VOLUMESTRUCT_H
#define ERDECON_RECON_VOLUMESTRUCT_H

typedef struct Geometry {
    long Nx;         /* real-space extent along x */
    long Ny;
    long Nz;
    long NxFSkip;    /* 2 if Nx even, 1 if Nx odd */
    long LengthF;    /* (Nx + NxFSkip) * Ny * Nz */

    long NSurface;   /* 2*Nx*Ny + (Nz-2)*(2*Nx+2*Ny-4): voxels on the box surface */
    long* pool;      /* caller-owned long[7*NSurface]; may be NULL */
    long* boundary[7]; /* pool + k*NSurface, k = 0..6; reserved for Derivative() */
} Geometry;

/*
 * Fill *G for an Nx x Ny x Nz volume.  mempool may be NULL if the caller
 * does not need the boundary lists (they are not used by any delivered
 * code).  If non-NULL it must hold at least 7 * NSurface longs.
 */
void ConstructVolumeStruct(Geometry* G, long Nx, long Ny, long Nz, long* mempool);

/* Number of longs ConstructVolumeStruct expects mempool to provide. */
long VolumeStructPoolLength(long Nx, long Ny, long Nz);

#endif /* include guard */
