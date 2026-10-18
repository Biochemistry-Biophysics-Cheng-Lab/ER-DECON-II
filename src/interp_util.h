#ifndef UCSF_MSG_INTERP_UTIL_H
#define UCSF_MSG_INTERP_UTIL_H

enum {
    INTERP_X_SPLIT = 1 << 0,
    INTERP_Y_SPLIT = 1 << 1,
    INTERP_CONJ_SYM = 1 << 2
};

typedef struct Interp1DSourceWeights {
    int i1, i2; float w1, w2;
} Interp1DSourceWeights;


#ifdef __cplusplus
extern "C" {
#endif

int
interpolate_2d(
    int nxi,
    int nyi,
    float dkxi,
    float dkyi,
    const float ip[],
    int nxo,
    int nyo,
    float dkxo,
    float dkyo,
    float mag,
    float scl,
    int flags,
    float op[]
);

int interpolate_2d_r(
    int nxi,
    int nyi,
    float dkxi,
    float dkyi,
    const float ip[],
    int nxo,
    int nyo,
    float dkxo,
    float dkyo,
    float mag,
    float scl,
    int flags,
    float op[]
);

Interp1DSourceWeights* setup_1d_interp(
    int n1, float d1, int n2, float d2
);

Interp1DSourceWeights* setup_1d_interp_split(
    int n1, float d1, int n2, float d2
);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
