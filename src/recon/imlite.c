/*
 * imlite.c -- ERDecon's replacement for the part of Priism's IM image-file
 * library (libimlib/libive) used by the delivered front ends core_decon.c
 * and core2_decon.c.
 *
 * Goal: the Priism-free build writes the same bytes and prints the same
 * header summaries as the build linked against Priism 4.7.2's libimlib.a.
 * The behaviour below was established by differential tests against that
 * library; where the library does something
 * surprising the comment says "as libimlib".
 *
 * Files: MRC with either Priism's header layout (dvid 0xc0a0 at byte 96:
 * wavelengths, time points, interleaving, origin z,x,y at bytes 208-219) or
 * the standard/EM layout (MRC2000/2014: origin x,y,z at 196-207, "MAP " at
 * 208, machine stamp, rms).  Both byte orders are read; files are written
 * in the order they were read from (IMTrHdr copies it) or native order.
 * Pixel modes 0 (int8, as libimlib), 1 (int16), 2 (float32), 3 (complex
 * int16), 4 (complex float32), 6 (uint16); sections are always delivered to
 * and accepted from the caller as float / complex float.
 *
 * ERDecon, 2026-09-23.
 */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L

#include "IMInclude.h"
#include "ive_thread.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>

#define IML_MAX_STREAM 64
#define PRIISM_DVID ((int16_t) -16224)   /* 0xc0a0 */

/* Priism header layout; every member is naturally aligned (1024 bytes). */
typedef struct {
    int32_t num[3];        /*   0 nx ny nsec */
    int32_t mode;          /*  12 */
    int32_t nxyzst[3];     /*  16 */
    int32_t mxyz[3];       /*  28 */
    float cell[3];         /*  40 */
    float angle[3];        /*  52 */
    int32_t axis[3];       /*  64 */
    float mmm1[3];         /*  76 min max mean */
    int32_t nspg;          /*  88 */
    int32_t next;          /*  92 extended header bytes */
    int16_t dvid;          /*  96 */
    int16_t nblank;        /*  98 */
    int32_t ntst;          /* 100 */
    char blank[24];        /* 104 */
    int16_t nint;          /* 128 */
    int16_t nreal;         /* 130 */
    int16_t sub;           /* 132 */
    int16_t zfac;          /* 134 */
    float mm2[2];          /* 136 */
    float mm3[2];          /* 144 */
    float mm4[2];          /* 152 */
    int16_t imtype;        /* 160 */
    int16_t lensnum;       /* 162 */
    int16_t n1, n2, v1, v2;/* 164 */
    float mm5[2];          /* 172 */
    int16_t ntime;         /* 180 */
    int16_t interleave;    /* 182 */
    float tilt[3];         /* 184 */
    int16_t nwave;         /* 196 */
    int16_t wave[5];       /* 198 */
    float zxy0[3];         /* 208 */
    int32_t ntitle;        /* 220 */
    char title[10][80];    /* 224 */
} PriismHeader;

typedef char iml_size_check[(sizeof(PriismHeader) == 1024) ? 1 : -1];

typedef struct {
    int open;
    int writable;
    FILE* fp;
    char* name;
    int swapped;            /* file byte order is not native */
    int em;                 /* standard (non-Priism) header layout */
    PriismHeader h;         /* native byte order */
    unsigned char* ext;     /* extended header, native order (4-byte words) */
    size_t next_alloc;
    long long cursec;       /* current section */
} Stream;

static Stream streams[IML_MAX_STREAM];


/* ---------------------------------------------------------------- utils */

static Stream*
get(int istream)
{
    if (istream < 0 || istream >= IML_MAX_STREAM || ! streams[istream].open) {
        return NULL;
    }
    return &streams[istream];
}

static void
swab4(void* p, size_t n)
{
    unsigned char* c = (unsigned char*) p;
    size_t i;
    for (i = 0; i < n; ++i, c += 4) {
        unsigned char t0 = c[0], t1 = c[1];
        c[0] = c[3]; c[1] = c[2]; c[2] = t1; c[3] = t0;
    }
}

static void
swab2(void* p, size_t n)
{
    unsigned char* c = (unsigned char*) p;
    size_t i;
    for (i = 0; i < n; ++i, c += 2) {
        unsigned char t = c[0];
        c[0] = c[1]; c[1] = t;
    }
}

/* Swap a header between file and native order (involution). */
static void
swab_header(PriismHeader* h, int em)
{
    unsigned char* b = (unsigned char*) h;
    swab4(b, 24);                    /* 0..95 */
    if (em) {
        /* standard layout: 96..195 "extra" words, origin, rms, nlabl;
           "MAP " and the machine stamp are bytes. */
        swab4(b + 96, 25);           /* 96..195 */
        swab4(b + 196, 3);           /* origin */
        swab4(b + 216, 2);           /* rms, nlabl */
    } else {
        swab2(b + 96, 2);            /* dvid nblank */
        swab4(b + 100, 1);           /* ntst */
        swab2(b + 128, 4);
        swab4(b + 136, 6);
        swab2(b + 160, 6);
        swab4(b + 172, 2);
        swab2(b + 180, 2);
        swab4(b + 184, 3);
        swab2(b + 196, 6);
        swab4(b + 208, 4);           /* zxy0, ntitle */
    }
}

static void
get_zwt(const Stream* s, int* nz, int* nw, int* nt, int* seq)
{
    int w, t;

    if (s->em) {
        /* as libimlib: an MRC2014 image stack of volumes (space group 401)
           has mz sections per volume, one volume per time point */
        if (s->h.nspg == 401 && s->h.mxyz[2] > 0) {
            *nz = s->h.mxyz[2];
            *nt = s->h.num[2] / s->h.mxyz[2];
        } else {
            *nz = s->h.num[2];
            *nt = 1;
        }
        *nw = 1;
        *seq = ZTW_SEQUENCE;
        return;
    }
    w = (s->h.nwave > 0) ? s->h.nwave : 1;
    t = (s->h.ntime > 0) ? s->h.ntime : 1;
    *nw = w;
    *nt = t;
    *nz = s->h.num[2] / (w * t);
    *seq = s->h.interleave;
}

static int
bytes_per_pixel(int mode)
{
    switch (mode) {
    case 0: return 1;
    case 1: case 6: return 2;
    case 2: return 4;
    case 3: return 4;
    case 4: return 8;
    default: return 0;
    }
}

static int
is_complex(int mode)
{
    return mode == 3 || mode == 4;
}

static off_t
data_offset(const Stream* s)
{
    return (off_t) 1024 + (off_t) s->h.next;
}

static off_t
section_bytes(const Stream* s)
{
    return (off_t) s->h.num[0] * (off_t) s->h.num[1] *
        (off_t) bytes_per_pixel(s->h.mode);
}


/* --------------------------------------------------------- header print */

static void
pr_g(int width, int prec, double v)
{
    printf("%c%-#*.*g", (v < 0.0) ? '-' : ' ', width, prec, fabs(v));
}

static const char*
pixel_type_name(int mode)
{
    switch (mode) {
    case 0: return "1 Byte Integer";
    case 1: return "2 Byte Integer";
    case 2: return "Floating-point";
    case 3: return "Complex, 2 Byte Integer";
    case 4: return "Complex, Floating-point";
    case 6: return "2 Byte Integer, Unsigned";
    case 7: return "4 Byte Integer";
    default: return "Unknown";
    }
}

static const char*
sequence_name(int seq)
{
    switch (seq) {
    case ZTW_SEQUENCE: return "Z -> T -> W";
    case WZT_SEQUENCE: return "W -> Z -> T";
    case ZWT_SEQUENCE: return "Z -> W -> T";
    default: return "Unknown";
    }
}

static void
print_titles(const PriismHeader* h)
{
    int n = h->ntitle, i;

    if (n < 0) n = 0;
    if (n > 10) n = 10;
    printf(" Titles (%2d).............................\n", n);
    for (i = 0; i < n; ++i) {
        printf(" %.*s\n", 80, h->title[i]);
    }
}

static void
print_header(const Stream* s)
{
    const PriismHeader* h = &s->h;
    int nz, nw, nt, seq, i;
    float d[3];

    get_zwt(s, &nz, &nw, &nt, &seq);
    for (i = 0; i < 3; ++i) {
        d[i] = (h->mxyz[i] != 0) ? h->cell[i] / h->mxyz[i] : 0.0f;
    }
    printf(" Image Size (X,Y,Sections)...............   %4d %4d %4d\n",
           h->num[0], h->num[1], h->num[2]);
    /* as libimlib: the stored wave and time counts are printed as they are
       (an old file with ntime 0 shows NT 0; IMRtZWT still reports 1) */
    printf(" Data Organization (NZ,NW,NT)............   %4d %4d %4d\n",
           nz, s->em ? nw : h->nwave, s->em ? nt : h->ntime);
    printf(" Image Sequence..........................  %s\n",
           sequence_name(s->em ? ZTW_SEQUENCE : seq));
    printf(" Multi-Resolution Settings...............   %2d %2d\n",
           h->sub, h->zfac);
    printf(" Data Type...............................  %s\n",
           pixel_type_name(h->mode));
    printf(" Start Pixel / Grid Spacing (X,Y,Z)......   %3d %3d %3d   %2d %2d %2d\n",
           h->nxyzst[0], h->nxyzst[1], h->nxyzst[2],
           h->mxyz[0], h->mxyz[1], h->mxyz[2]);
    printf(" Pixel Size..............................  ");
    for (i = 0; i < 3; ++i) {
        printf(" ");
        pr_g(9, 4, d[i]);
    }
    printf("\n");
    printf(" Cell Angles.............................  %7.2f%7.2f%7.2f\n",
           h->angle[0], h->angle[1], h->angle[2]);
    {
        float o[3];
        if (s->em) {
            memcpy(o, (const unsigned char*) h + 196, sizeof(o));
        } else {
            o[0] = h->zxy0[1]; o[1] = h->zxy0[2]; o[2] = h->zxy0[0];
        }
        printf(" Origin (X,Y,Z)..........................   ");
        pr_g(9, 4, o[0]); printf(" ");
        pr_g(9, 4, o[1]); printf(" ");
        pr_g(9, 4, o[2]); printf("\n");
    }
    printf(" Minimum Intensity.......................   "); pr_g(11, 6, h->mmm1[0]); printf("\n");
    printf(" Maximum Intensity.......................   "); pr_g(11, 6, h->mmm1[1]); printf("\n");
    printf(" Mean    Intensity.......................   "); pr_g(11, 6, h->mmm1[2]); printf("\n");
    if (s->em) {
        float rms;
        memcpy(&rms, (const unsigned char*) h + 216, sizeof(rms));
        printf(" RMS              .......................   "); pr_g(11, 6, rms); printf("\n");
    }
    printf(" Rotation Angles (X,Y,Z).................  %7.2f%7.2f%7.2f\n",
           h->tilt[0], h->tilt[1], h->tilt[2]);
    printf(" Crystallography Space Group Number......  %6d\n", h->nspg);
    printf(" Extended Header Size (Bytes)............  %6d\n", h->next);
    printf(" Image Type..............................  %6d\n", h->imtype);
    printf(" Lens ID Number..........................  %6d\n", h->lensnum);
    if (! s->em) {
        printf(" Number of Wavelengths...................  %3d\n", h->nwave);
        printf(" Wavelengths (nm)........................  ");
        for (i = 0; i < h->nwave && i < 5; ++i) {
            printf("%4d", h->wave[i]);
        }
        printf("\n");
    }
    print_titles(h);
    fflush(stdout);
}


/* ---------------------------------------------------------- open/close */

static int
read_header(Stream* s)
{
    PriismHeader h;
    unsigned char* b = (unsigned char*) &h;

    if (fread(&h, 1, sizeof(h), s->fp) != sizeof(h)) {
        return -1;
    }
    /* as libimlib (IMcheck_header_swabbed): with "MAP " at byte 208 the
       machine stamp decides; otherwise nx and ny both above 65535 mean the
       file was written with the other byte order.  The standard (EM)
       layout is recognised by "MAP " and the absence of Priism's dvid. */
    {
        uint32_t nx = (uint32_t) h.num[0], ny = (uint32_t) h.num[1];
        int has_map = memcmp(b + 208, "MAP ", 4) == 0;
        s->swapped = 0;
        if (has_map && b[212] == 0x11 && b[213] == 0x11) {
            s->swapped = 1;
        } else if (has_map && b[212] == 0x44 && (unsigned char) (b[213] - 0x40) <= 0xf) {
            s->swapped = 0;
        } else {
            s->swapped = nx > 0xffff && ny > 0xffff;
        }
        {
            int16_t dv = h.dvid;
            if (s->swapped) swab2(&dv, 1);
            s->em = has_map && dv != PRIISM_DVID;
        }
    }
    if (s->swapped) {
        swab_header(&h, s->em);
    }
    (void) b;
    s->h = h;
    if (h.next < 0) {
        return -1;
    }
    if (h.next > 0) {
        s->ext = (unsigned char*) malloc((size_t) h.next);
        if (s->ext == NULL) return -1;
        s->next_alloc = (size_t) h.next;
        if (fread(s->ext, 1, (size_t) h.next, s->fp) != (size_t) h.next) {
            return -1;
        }
        if (s->swapped) swab4(s->ext, (size_t) h.next / 4);
    }
    return 0;
}

int
IMOpen(int istream, const char* name, const char* attrib)
{
    Stream* s;
    const char* how;

    if (istream < 0 || istream >= IML_MAX_STREAM || streams[istream].open) {
        return 1;
    }
    s = &streams[istream];
    memset(s, 0, sizeof(*s));
    if (strcmp(attrib, "ro") == 0) {
        s->fp = fopen(name, "rb");
        how = "RO";
    } else if (strcmp(attrib, "new") == 0) {
        s->fp = fopen(name, "w+b");
        s->writable = 1;
        how = "NEW";
    } else if (strcmp(attrib, "old") == 0) {
        s->fp = fopen(name, "r+b");
        s->writable = 1;
        how = "OLD";
    } else {
        return 1;
    }
    if (s->fp == NULL) {
        return 1;
    }
    s->name = (char*) malloc(strlen(name) + 1);
    strcpy(s->name, name);
    s->open = 1;
    if (! (attrib[0] == 'n')) {
        if (read_header(s) != 0) {
            fclose(s->fp);
            free(s->ext);
            free(s->name);
            memset(s, 0, sizeof(*s));
            return 1;
        }
    } else {
        /* as libimlib: a new file starts with a blank Priism header */
        s->h.dvid = PRIISM_DVID;
        s->h.mode = IW_FLOAT;
        s->h.mxyz[0] = s->h.mxyz[1] = s->h.mxyz[2] = 1;
        s->h.axis[0] = 1; s->h.axis[1] = 2; s->h.axis[2] = 3;
        s->h.nwave = 1;
        s->h.ntime = 1;
    }
    printf("\n %s image file on unit: %2d %s\n\n", how, istream, name);
    fflush(stdout);
    return 0;
}

void
IMClose(int istream)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    fclose(s->fp);
    free(s->ext);
    free(s->name);
    memset(s, 0, sizeof(*s));
}

int
IMUnit(int istream)
{
    return (get(istream) != NULL) ? istream : -1;
}

int
IMRtWID(int istream)
{
    (void) istream;
    return 0;
}


/* --------------------------------------------------------- header I/O */

void
IMRdHdr(int istream, int ixyz[3], int mxyz[3], int* imode,
        float* min, float* max, float* mean)
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    for (i = 0; i < 3; ++i) {
        ixyz[i] = s->h.num[i];
        mxyz[i] = s->h.mxyz[i];
    }
    *imode = s->h.mode;
    *min = s->h.mmm1[0];
    *max = s->h.mmm1[1];
    *mean = s->h.mmm1[2];
    print_header(s);
}

void
IMWrHdr(int istream, const char title[80], int ntflag,
        float dmin, float dmax, float dmean)
{
    Stream* s = get(istream);
    PriismHeader h;
    int i;

    if (s == NULL || ! s->writable) return;
    s->h.mmm1[0] = dmin;
    s->h.mmm1[1] = dmax;
    s->h.mmm1[2] = dmean;
    if (ntflag >= 0 && title != NULL) {
        if (s->h.ntitle < 0) s->h.ntitle = 0;
        if (ntflag == 0) {
            memcpy(s->h.title[0], title, 80);
            s->h.ntitle = 1;
        } else if (ntflag == 1) {
            if (s->h.ntitle < 10) {
                memcpy(s->h.title[s->h.ntitle], title, 80);
                ++s->h.ntitle;
            } else {
                for (i = 0; i < 9; ++i) memcpy(s->h.title[i], s->h.title[i + 1], 80);
                memcpy(s->h.title[9], title, 80);
            }
        } else if (ntflag == 2) {
            int n = (s->h.ntitle < 10) ? s->h.ntitle : 9;
            for (i = n; i > 0; --i) memcpy(s->h.title[i], s->h.title[i - 1], 80);
            memcpy(s->h.title[0], title, 80);
            s->h.ntitle = n + 1;
        }
    }
    if (s->em) {
        /* as libimlib: the machine stamp is set to the writer's order */
        unsigned char* b = (unsigned char*) &s->h;
        b[212] = 0x44; b[213] = 0x41; b[214] = 0; b[215] = 0;
    }
    h = s->h;
    if (s->swapped) swab_header(&h, s->em);
    if (fseeko(s->fp, 0, SEEK_SET) != 0 ||
        fwrite(&h, 1, sizeof(h), s->fp) != sizeof(h)) {
        fprintf(stderr, "IMWrHdr: write failed for %s: %s\n", s->name, strerror(errno));
        return;
    }
    if (s->h.next > 0) {
        if (s->ext == NULL || s->next_alloc < (size_t) s->h.next) {
            unsigned char* e = (unsigned char*) realloc(s->ext, (size_t) s->h.next);
            if (e == NULL) return;
            memset(e + s->next_alloc, 0, (size_t) s->h.next - s->next_alloc);
            s->ext = e;
            s->next_alloc = (size_t) s->h.next;
        }
        if (fwrite(s->ext, 1, (size_t) s->h.next, s->fp) != (size_t) s->h.next) {
            fprintf(stderr, "IMWrHdr: write failed for %s: %s\n", s->name, strerror(errno));
        }
    }
    fflush(s->fp);
}

void
IMTrHdr(int istream, int jstream)
{
    Stream* s = get(istream);
    Stream* t = get(jstream);
    if (s == NULL || t == NULL) return;
    s->h = t->h;
    s->em = t->em;
    s->swapped = 0;          /* as libimlib: output is written in native order */
    free(s->ext);
    s->ext = NULL;
    s->next_alloc = 0;
    if (t->h.next > 0 && t->ext != NULL) {
        s->ext = (unsigned char*) malloc((size_t) t->h.next);
        if (s->ext != NULL) {
            memcpy(s->ext, t->ext, (size_t) t->h.next);
            s->next_alloc = (size_t) t->h.next;
        }
    }
}


/* ------------------------------------------------- header field access */

void
IMAlDel(int istream, const float delta[3])
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    for (i = 0; i < 3; ++i) s->h.cell[i] = delta[i] * s->h.mxyz[i];
}

void
IMAlExHdrSize(int istream, int nints, int nreals, int nsecs)
{
    Stream* s = get(istream);
    size_t n;
    if (s == NULL) return;
    s->h.nint = (int16_t) nints;
    s->h.nreal = (int16_t) nreals;
    n = (size_t) (nints + nreals) * 4u * (size_t) nsecs;
    s->h.next = (int32_t) n;
    free(s->ext);
    s->ext = (n > 0) ? (unsigned char*) calloc(n, 1) : NULL;
    s->next_alloc = (s->ext != NULL) ? n : 0;
}

void
IMAlMode(int istream, int imode)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    s->h.mode = imode;
}

void
IMAlOrig(int istream, float xorig, float yorig, float zorig)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    if (s->em) {
        float o[3];
        o[0] = xorig; o[1] = yorig; o[2] = zorig;
        memcpy((unsigned char*) &s->h + 196, o, sizeof(o));
    } else {
        s->h.zxy0[0] = zorig;
        s->h.zxy0[1] = xorig;
        s->h.zxy0[2] = yorig;
    }
}

void
IMAlRes(int istream, int mres, int mzfact)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    s->h.sub = (int16_t) mres;
    s->h.zfac = (int16_t) mzfact;
}

void
IMAlSiz(int istream, const int inxyz[3], const int nxyzst[3])
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    for (i = 0; i < 3; ++i) {
        s->h.num[i] = inxyz[i];
        s->h.nxyzst[i] = nxyzst[i];
    }
}

void
IMAlTSt(int istream, int itst)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    s->h.ntst = itst;
}

void
IMAlWav(int istream, int nwave, const float wavelen[])
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    if (s->em) return;
    s->h.nwave = (int16_t) nwave;
    for (i = 0; i < 5; ++i) {
        s->h.wave[i] = (i < nwave) ? (int16_t) (wavelen[i] + 0.5f) : 0;
    }
}

void
IMAlWavMM(int istream, int iwave, float dmin, float dmax)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    switch (iwave) {
    case 0: s->h.mmm1[0] = dmin; s->h.mmm1[1] = dmax; break;
    case 1: s->h.mm2[0] = dmin; s->h.mm2[1] = dmax; break;
    case 2: s->h.mm3[0] = dmin; s->h.mm3[1] = dmax; break;
    case 3: s->h.mm4[0] = dmin; s->h.mm4[1] = dmax; break;
    case 4: s->h.mm5[0] = dmin; s->h.mm5[1] = dmax; break;
    default: break;
    }
}

void
IMAlZWT(int istream, int nz, int nw, int nt, int wflag)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    s->h.num[2] = nz * nw * nt;
    if (s->em) return;
    s->h.nwave = (int16_t) nw;
    s->h.ntime = (int16_t) nt;
    s->h.interleave = (int16_t) wflag;
}

void
IMRtDat(int istream, int* itype, int* lensnum, int* n1, int* n2,
        float* v1, float* v2)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    *itype = s->h.imtype;
    *lensnum = s->h.lensnum;
    *n1 = s->h.n1;
    *n2 = s->h.n2;
    *v1 = s->h.v1;
    *v2 = s->h.v2;
}

void
IMRtDel(int istream, float delta[3])
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    for (i = 0; i < 3; ++i) {
        delta[i] = (s->h.mxyz[i] != 0) ? s->h.cell[i] / s->h.mxyz[i] : 0.0f;
    }
}

void
IMRtMode(int istream, int* mode)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    *mode = s->h.mode;
}

void
IMRtOrig(int istream, float* xorig, float* yorig, float* zorig)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    if (s->em) {
        float o[3];
        memcpy(o, (unsigned char*) &s->h + 196, sizeof(o));
        *xorig = o[0]; *yorig = o[1]; *zorig = o[2];
    } else {
        *xorig = s->h.zxy0[1];
        *yorig = s->h.zxy0[2];
        *zorig = s->h.zxy0[0];
    }
}

void
IMRtSiz(int istream, int inxyz[3], int mxyz[3], int nxyzst[3])
{
    Stream* s = get(istream);
    int i;
    if (s == NULL) return;
    for (i = 0; i < 3; ++i) {
        inxyz[i] = s->h.num[i];
        mxyz[i] = s->h.mxyz[i];
        nxyzst[i] = s->h.nxyzst[i];
    }
}

void
IMRtTSt(int istream, int* itst)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    *itst = s->h.ntst;
}

void
IMRtWav(int istream, int* nwave, float wavelen[])
{
    Stream* s = get(istream);
    int i, n;
    if (s == NULL) return;
    if (s->em) {
        *nwave = 1;
        wavelen[0] = 0.0f;
        return;
    }
    n = (s->h.nwave > 0) ? s->h.nwave : 1;
    if (n > IW_MAX_WAVE) n = IW_MAX_WAVE;
    *nwave = n;
    for (i = 0; i < n; ++i) wavelen[i] = (float) s->h.wave[i];
}

void
IMRtZWT(int istream, int* nz, int* nw, int* nt, int* wflag)
{
    Stream* s = get(istream);
    if (s == NULL) return;
    get_zwt(s, nz, nw, nt, wflag);
}


/* ------------------------------------------------ extended header fix */

/*
 * IMFixExHdrZWT: the output keeps only the selected z range, waves and time
 * points.  As libimlib: the section organisation (interleave) is kept; the
 * wavelength slots 0..nwt-1 are replaced by the selected ones and the
 * remaining slots are left as they were; the extended header keeps the
 * records of the selected sections, in the output's own section order, and
 * its size is rounded up to a multiple of 1024 bytes.  An extended header
 * without records (nint = nreal = 0) is dropped.
 */
static long long
sec_index(int seq, int nz, int nw, int nt, int iz, int iw, int it)
{
    switch (seq) {
    case WZT_SEQUENCE:
        return iw + (long long) nw * (iz + (long long) nz * it);
    case ZWT_SEQUENCE:
        return iz + (long long) nz * (iw + (long long) nw * it);
    default:
        return iz + (long long) nz * (it + (long long) nt * iw);
    }
}

void
IMFixExHdrZWT(int istream, int izst, int izend, int nwt,
              const int iwtable[], int itst, int itend, int itinc, int iflag)
{
    Stream* s = get(istream);
    int nz_in, nw_in, nt_in, seq;
    int nz, nt, iz, iw, it;
    size_t rec;

    if (s == NULL) return;
    get_zwt(s, &nz_in, &nw_in, &nt_in, &seq);
    if (itinc < 1) itinc = 1;
    nz = izend - izst + 1;
    nt = (itend - itst) / itinc + 1;

    rec = (size_t) (s->h.nint + s->h.nreal) * 4u;
    if (s->h.next > 0 && rec > 0) {
        size_t used = rec * (size_t) nz * (size_t) nwt * (size_t) nt;
        size_t next_new = (used + 1023u) / 1024u * 1024u;
        unsigned char* ext_new = (unsigned char*) calloc(next_new > 0 ? next_new : 1, 1);
        if (ext_new != NULL) {
            for (it = 0; it < nt; ++it) {
                for (iw = 0; iw < nwt; ++iw) {
                    for (iz = 0; iz < nz; ++iz) {
                        int zi = izst + iz, wi = iwtable[iw], ti = itst + it * itinc;
                        long long src = sec_index(seq, nz_in, nw_in, nt_in, zi, wi, ti);
                        long long dst = sec_index(seq, nz, nwt, nt, iz, iw, it);
                        if (zi >= 0 && zi < nz_in && src >= 0 && s->ext != NULL &&
                            (size_t) (src + 1) * rec <= (size_t) s->h.next) {
                            memcpy(ext_new + (size_t) dst * rec, s->ext + (size_t) src * rec, rec);
                        }
                    }
                }
            }
            free(s->ext);
            s->ext = ext_new;
            s->next_alloc = next_new;
            s->h.next = (int32_t) next_new;
        }
    } else if (s->h.next > 0) {
        /* as libimlib: an extended header without integer/float records
           (e.g. an FEI/EM one) is dropped */
        free(s->ext);
        s->ext = NULL;
        s->next_alloc = 0;
        s->h.next = 0;
    }
    if (iflag) {
        s->h.num[2] = nz * nwt * nt;
        if (! s->em) {
            int16_t wav[5];
            memcpy(wav, s->h.wave, sizeof(wav));
            for (iw = 0; iw < nwt && iw < 5; ++iw) {
                if (iwtable[iw] >= 0 && iwtable[iw] < 5) wav[iw] = s->h.wave[iwtable[iw]];
            }
            memcpy(s->h.wave, wav, sizeof(wav));
            s->h.nwave = (int16_t) nwt;
            s->h.ntime = (int16_t) nt;
        }
    }
}


/* ---------------------------------------------------------- data I/O */

int
IMPosnZWT(int istream, int iz, int iw, int it)
{
    Stream* s = get(istream);
    int nz, nw, nt, seq;
    long long sec;

    if (s == NULL) return 1;
    get_zwt(s, &nz, &nw, &nt, &seq);
    if (iz < 0 || iz >= nz || iw < 0 || iw >= nw || it < 0 || it >= nt) {
        return 1;
    }
    sec = sec_index(seq, nz, nw, nt, iz, iw, it);
    s->cursec = sec;
    return 0;
}

/* Convert n pixels of file mode `mode` (file byte order) to float. */
static void
to_float(const Stream* s, const unsigned char* in, float* out, size_t npix)
{
    size_t i;
    int mode = s->h.mode;

    switch (mode) {
    case 0:
        /* as libimlib 4.7.2: mode 0 is signed (MRC2014 int8) */
        for (i = 0; i < npix; ++i) out[i] = (float) (int8_t) in[i];
        break;
    case 1: case 3: {
        size_t n = (mode == 3) ? 2 * npix : npix;
        for (i = 0; i < n; ++i) {
            int16_t v;
            memcpy(&v, in + 2 * i, 2);
            if (s->swapped) swab2(&v, 1);
            out[i] = (float) v;
        }
        break;
    }
    case 6:
        for (i = 0; i < npix; ++i) {
            uint16_t v;
            memcpy(&v, in + 2 * i, 2);
            if (s->swapped) swab2(&v, 1);
            out[i] = (float) v;
        }
        break;
    case 2: case 4: {
        size_t n = (mode == 4) ? 2 * npix : npix;
        memcpy(out, in, 4 * n);
        if (s->swapped) swab4(out, n);
        break;
    }
    default:
        memset(out, 0, sizeof(float) * npix * (is_complex(mode) ? 2 : 1));
        break;
    }
}

static int
read_rows(Stream* s, int y0, int ny, unsigned char** p_buf)
{
    off_t row = (off_t) s->h.num[0] * bytes_per_pixel(s->h.mode);
    off_t off = data_offset(s) + (off_t) s->cursec * section_bytes(s) + (off_t) y0 * row;
    size_t n = (size_t) row * (size_t) ny;
    unsigned char* buf = (unsigned char*) malloc(n > 0 ? n : 1);

    if (buf == NULL) return 1;
    if (fseeko(s->fp, off, SEEK_SET) != 0 || fread(buf, 1, n, s->fp) != n) {
        fprintf(stderr, "IMRdSec: read failed for %s section %lld\n", s->name, s->cursec);
        memset(buf, 0, n);
    }
    *p_buf = buf;
    return 0;
}

void
IMRdSec(int istream, void* array)
{
    Stream* s = get(istream);
    unsigned char* buf;

    if (s == NULL) return;
    if (read_rows(s, 0, s->h.num[1], &buf) != 0) return;
    to_float(s, buf, (float*) array, (size_t) s->h.num[0] * (size_t) s->h.num[1]);
    free(buf);
    ++s->cursec;
}

void
IMRdPas(int istream, void* array, int mx, int my,
        int nx1, int nx2, int ny1, int ny2)
{
    Stream* s = get(istream);
    unsigned char* buf;
    float* row;
    int ny = ny2 - ny1 + 1, nx = nx2 - nx1 + 1, iy;
    int cplx;
    size_t bpp;

    (void) my;
    if (s == NULL || nx < 1 || ny < 1) return;
    cplx = is_complex(s->h.mode) ? 2 : 1;
    bpp = (size_t) bytes_per_pixel(s->h.mode);
    if (read_rows(s, ny1, ny, &buf) != 0) return;
    row = (float*) malloc(sizeof(float) * (size_t) cplx * (size_t) nx);
    if (row != NULL) {
        for (iy = 0; iy < ny; ++iy) {
            to_float(s, buf + ((size_t) iy * (size_t) s->h.num[0] + (size_t) nx1) * bpp,
                     row, (size_t) nx);
            memcpy((float*) array + (size_t) iy * (size_t) mx * cplx, row,
                   sizeof(float) * (size_t) cplx * (size_t) nx);
        }
        free(row);
    }
    free(buf);
    ++s->cursec;
}

void
IMWrSec(int istream, const void* array)
{
    Stream* s = get(istream);
    size_t npix, n;
    off_t off;

    if (s == NULL || ! s->writable) return;
    if (s->h.mode != IW_FLOAT && s->h.mode != IW_COMPLEX) {
        fprintf(stderr, "IMWrSec: only float and complex output are supported\n");
        return;
    }
    npix = (size_t) s->h.num[0] * (size_t) s->h.num[1];
    n = npix * (s->h.mode == IW_COMPLEX ? 2 : 1);
    off = data_offset(s) + (off_t) s->cursec * section_bytes(s);
    if (fseeko(s->fp, off, SEEK_SET) != 0) {
        fprintf(stderr, "IMWrSec: seek failed for %s\n", s->name);
        return;
    }
    if (s->swapped) {
        float* tmp = (float*) malloc(sizeof(float) * n);
        if (tmp == NULL) return;
        memcpy(tmp, array, sizeof(float) * n);
        swab4(tmp, n);
        if (fwrite(tmp, sizeof(float), n, s->fp) != n) {
            fprintf(stderr, "IMWrSec: write failed for %s\n", s->name);
        }
        free(tmp);
    } else if (fwrite(array, sizeof(float), n, s->fp) != n) {
        fprintf(stderr, "IMWrSec: write failed for %s\n", s->name);
    }
    ++s->cursec;
}

void
IMCalcDen(const float* array, int mx, int my, int nx1, int nx2,
          int ny1, int ny2, float* dmin, float* dmax, float* dmean)
{
    int ix, iy;
    double sum = 0.0;
    float lo = array[(size_t) ny1 * mx + nx1], hi = lo;

    (void) my;
    for (iy = ny1; iy <= ny2; ++iy) {
        for (ix = nx1; ix <= nx2; ++ix) {
            float v = array[(size_t) iy * mx + ix];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
            sum += v;
        }
    }
    *dmin = lo;
    *dmax = hi;
    *dmean = (float) (sum / ((double) (nx2 - nx1 + 1) * (double) (ny2 - ny1 + 1)));
}


/* --------------------------------------------------------- libive part */

int
IVEThreadCountProcessors(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int) n : 1;
}
