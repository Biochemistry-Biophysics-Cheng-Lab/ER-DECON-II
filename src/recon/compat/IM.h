/* IM.h -- ERDecon compatibility header: the IM calls used by the front ends. */
#ifndef ERDECON_COMPAT_IM_H
#define ERDECON_COMPAT_IM_H

#ifdef __cplusplus
extern "C" {
#endif

int IMOpen(int istream, const char* name, const char* attrib);
void IMClose(int istream);
int IMUnit(int istream);
int IMRtWID(int istream);
void IMRdHdr(int istream, int ixyz[3], int mxyz[3], int* imode,
             float* min, float* max, float* mean);
void IMWrHdr(int istream, const char title[80], int ntflag,
             float dmin, float dmax, float dmean);
void IMTrHdr(int istream, int jstream);
void IMAlDel(int istream, const float delta[3]);
void IMAlExHdrSize(int istream, int nints, int nreals, int nsecs);
void IMAlMode(int istream, int imode);
void IMAlOrig(int istream, float xorig, float yorig, float zorig);
void IMAlRes(int istream, int mres, int mzfact);
void IMAlSiz(int istream, const int inxyz[3], const int nxyzst[3]);
void IMAlTSt(int istream, int itst);
void IMAlWav(int istream, int nwave, const float wavelen[]);
void IMAlWavMM(int istream, int nwave, float dmin, float dmax);
void IMAlZWT(int istream, int nz, int nw, int nt, int wflag);
void IMRtDat(int istream, int* itype, int* lensnum, int* n1, int* n2,
             float* v1, float* v2);
void IMRtDel(int istream, float delta[3]);
void IMRtMode(int istream, int* mode);
void IMRtOrig(int istream, float* xorig, float* yorig, float* zorig);
void IMRtSiz(int istream, int inxyz[3], int mxyz[3], int nxyzst[3]);
void IMRtTSt(int istream, int* itst);
void IMRtWav(int istream, int* nwave, float wavelen[]);
void IMRtZWT(int istream, int* nz, int* nw, int* nt, int* wflag);
int IMPosnZWT(int istream, int iz, int iw, int it);
void IMRdPas(int istream, void* array, int mx, int my,
             int nx1, int nx2, int ny1, int ny2);
void IMRdSec(int istream, void* array);
void IMWrSec(int istream, const void* array);
void IMFixExHdrZWT(int istream, int izst, int izend, int nwt,
                   const int iwtable[], int itst, int itend, int itinc,
                   int iflag);
void IMCalcDen(const float* array, int mx, int my, int nx1, int nx2,
               int ny1, int ny2, float* dmin, float* dmax, float* dmean);

#ifdef __cplusplus
}
#endif

#endif
