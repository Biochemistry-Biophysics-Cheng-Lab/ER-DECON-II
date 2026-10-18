#include "interp_util.h"
#include <math.h>    /* ceil, floor */
#include <stdlib.h>  /* free, malloc */


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
)
{
    /* Uses a bilinear interpolation method. */
    Interp1DSourceWeights* xweights;
    Interp1DSourceWeights* yweights;
    int nxif, nxof, nyof;
    int iyo, ixo;

    if (nxo < 1 || nyo < 1) {
	return 0;
    }
    if ((flags & INTERP_X_SPLIT) != 0) {
	nxif = nxi;
	nxof = nxo;
	xweights = setup_1d_interp_split(nxif, dkxi * mag, nxof, dkxo);
    } else {
	nxif = nxi / 2 + 1;
	nxof = nxo / 2 + 1;
	xweights = setup_1d_interp(nxif, dkxi * mag, nxof, dkxo);
    }
    if ((flags & INTERP_Y_SPLIT) != 0) {
	nyof = nyo;
	yweights = setup_1d_interp_split(nyi, dkyi * mag, nyof, dkyo);
    } else {
	nyof = nyo / 2 + 1;
	yweights = setup_1d_interp(nyi / 2 + 1, dkyi * mag, nyof, dkyo);
    }
    if (xweights == 0 || yweights == 0) {
	if (xweights != 0) {
	    free(xweights);
	} else if (yweights != 0) {
	    free(yweights);
	}
	return 1;
    }

    for (ixo = 0; ixo < nxof; ++ixo) {
	xweights[ixo].i1 *= 2;
	xweights[ixo].i2 *= 2;
    }
    for (iyo = 0; iyo < nyof; ++iyo) {
	yweights[iyo].i1 *= nxif + nxif;
	yweights[iyo].i2 *= nxif + nxif;
    }

    if ((flags & INTERP_CONJ_SYM) != 0) {
	int nxo_even = (nxo % 2 == 0);
	int nyo_even = (nyo % 2 == 0);
	int ixh = nxo / 2;
	int iyh = nyo / 2;
	int ixo2;

	op[0] = xweights[0].w1 * ip[yweights[0].i1 + xweights[0].i1] +
	    xweights[0].w2 * ip[yweights[0].i1 + xweights[0].i2];
	op[1] = 0.0f;
	for (ixo = 1, ixo2 = 2; ixo < ixh; ++ixo, ixo2 += 2) {
	    op[ixo2] =
		xweights[ixo].w1 * ip[yweights[0].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i1 + xweights[ixo].i2];
	    op[ixo2] *= yweights[0].w1;
	    op[ixo2] += yweights[0].w2 * (
		xweights[ixo].w1 * ip[yweights[0].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i2 + xweights[ixo].i2]
	    );
	    op[ixo2] *= scl;
	    op[ixo2 + 1] =
		xweights[ixo].w1 * ip[yweights[0].i1 + xweights[ixo].i1 + 1] +
		xweights[ixo].w2 * ip[yweights[0].i1 + xweights[ixo].i2 + 1];
	    op[ixo2 + 1] *= yweights[0].w1;
	    op[ixo2 + 1] += yweights[0].w2 * (
		xweights[ixo].w1 * ip[yweights[0].i2 + xweights[ixo].i1 + 1] +
		xweights[ixo].w2 * ip[yweights[0].i2 + xweights[ixo].i2 + 1]
	    );
	    op[ixo2 + 1] *= scl;
	}
	op[ixh + ixh] =
	    xweights[ixh].w1 * ip[yweights[0].i1 + xweights[ixh].i1] +
	    xweights[ixh].w2 * ip[yweights[0].i1 + xweights[ixh].i2];
	op[ixh + ixh] *= yweights[0].w1;
	op[ixh + ixh] += yweights[0].w2 * (
	    xweights[ixh].w1 * ip[yweights[0].i2 + xweights[ixh].i1] +
	    xweights[ixh].w2 * ip[yweights[0].i2 + xweights[ixh].i2]
	);
	op[ixh + ixh] *= scl;
	if (nxo_even) {
	    op[ixh + ixh + 1] = 0.0f;
	} else {
	    op[ixh + ixh + 1] =
		xweights[ixh].w1 * ip[yweights[0].i1 + xweights[ixh].i1 + 1] +
		xweights[ixh].w2 * ip[yweights[0].i1 + xweights[ixh].i2 + 1];
	    op[ixh + ixh + 1] *= yweights[0].w1;
	    op[ixh + ixh + 1] += yweights[0].w2 * (
		xweights[ixh].w1 * ip[yweights[0].i2 + xweights[ixh].i1 + 1] +
		xweights[ixh].w2 * ip[yweights[0].i2 + xweights[ixh].i2 + 1]
	    );
	    op[ixh + ixh + 1] *= scl;
	}
	for (ixo2 = ixh + ixh + 2; ixo2 < nxof + nxof; ixo2 += 2) {
	    op[ixo2] = op[nxof + nxof - ixo2];
	    op[ixo2 + 1] = -op[nxof + nxof - ixo2 + 1];
	}
	op += nxof + nxof;
	for (iyo = 1; iyo < iyh; ++iyo) {
	    for (ixo = 0, ixo2 = 0; ixo < nxof; ++ixo, ixo2 += 2) {
		op[ixo2] =
		    xweights[ixo].w1 * ip[yweights[iyo].i1 + xweights[ixo].i1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i1 + xweights[ixo].i2];
		op[ixo2] *= yweights[iyo].w1;
		op[ixo2] += yweights[iyo].w2 * (
		    xweights[ixo].w1 * ip[yweights[iyo].i2 + xweights[ixo].i1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i2 + xweights[ixo].i2]
		);
		op[ixo2] *= scl;
		op[ixo2 + 1] =
		    xweights[ixo].w1 * ip[yweights[iyo].i1 + xweights[ixo].i1 + 1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i1 + xweights[ixo].i2 + 1];
		op[ixo2 + 1] *= yweights[iyo].w1;
		op[ixo2 + 1] += yweights[iyo].w2 * (
		    xweights[ixo].w1 * ip[yweights[iyo].i2 + xweights[ixo].i1 + 1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i2 + xweights[ixo].i2 + 1]
		);
		op[ixo2 + 1] *= scl;
	    }
	    op += nxof + nxof;
	}
	op[0] =
	    xweights[0].w1 * ip[yweights[iyh].i1 + xweights[0].i1] +
	    xweights[0].w2 * ip[yweights[iyh].i1 + xweights[0].i2];
	op[0] *= yweights[iyh].w1;
	op[0] += yweights[iyh].w2 * (
	    xweights[0].w1 * ip[yweights[iyh].i2 + xweights[0].i1] +
	    xweights[0].w2 * ip[yweights[iyh].i2 + xweights[0].i2]
	);
	op[0] *= scl;
	if (nyo_even) {
	    op[1] = 0.0f;
	} else {
	    op[1] =
		xweights[0].w1 * ip[yweights[iyh].i1 + xweights[0].i1 + 1] +
		xweights[0].w2 * ip[yweights[iyh].i1 + xweights[0].i2 + 1];
	    op[1] *= yweights[iyh].w1;
	    op[1] += yweights[iyh].w2 * (
		xweights[0].w1 * ip[yweights[iyh].i2 + xweights[0].i1 + 1] +
		xweights[0].w2 * ip[yweights[iyh].i2 + xweights[0].i2 + 1]
	    );
	    op[1] *= scl;
	}
	for (ixo = 1, ixo2 = 2; ixo < ixh; ++ixo, ixo2 += 2) {
	    op[ixo2] =
		xweights[ixo].w1 * ip[yweights[iyh].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[iyh].i1 + xweights[ixo].i2];
	    op[ixo2] *= yweights[iyh].w1;
	    op[ixo2] += yweights[iyh].w2 * (
		xweights[ixo].w1 * ip[yweights[iyh].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[iyh].i2 + xweights[ixo].i2]
	    );
	    op[ixo2] *= scl;
	    op[ixo2 + 1] =
		xweights[ixo].w1 * ip[yweights[iyh].i1 + xweights[ixo].i1 + 1] +
		xweights[ixo].w2 * ip[yweights[iyh].i1 + xweights[ixo].i2 + 1];
	    op[ixo2 + 1] *= yweights[iyh].w1;
	    op[ixo2 + 1] += yweights[iyh].w2 * (
		xweights[ixo].w1 * ip[yweights[iyh].i2 + xweights[ixo].i1 + 1] +
		xweights[ixo].w2 * ip[yweights[iyh].i2 + xweights[ixo].i2 + 1]
	    );
	    op[ixo2 + 1] *= scl;
	}
	op[ixh + ixh] =
	    xweights[ixh].w1 * ip[yweights[iyh].i1 + xweights[ixh].i1] +
	    xweights[ixh].w2 * ip[yweights[iyh].i1 + xweights[ixh].i2];
	op[ixh + ixh] *= yweights[iyh].w1;
	op[ixh + ixh] += yweights[iyh].w2 * (
	    xweights[ixh].w1 * ip[yweights[iyh].i2 + xweights[ixh].i1] +
	    xweights[ixh].w2 * ip[yweights[iyh].i2 + xweights[ixh].i2]
	);
	op[ixh + ixh] *= scl;
	if (nxo_even && nyo_even) {
	    op[ixh + ixh + 1] = 0.0f;
	} else {
	    op[ixh + ixh + 1] =
		xweights[ixh].w1 * ip[yweights[iyh].i1 + xweights[ixh].i1 + 1] +
		xweights[ixh].w2 * ip[yweights[iyh].i1 + xweights[ixh].i2 + 1];
	    op[ixh + ixh + 1] *= yweights[iyh].w1;
	    op[ixh + ixh + 1] += yweights[iyh].w2 * (
		xweights[ixh].w1 * ip[yweights[iyh].i2 + xweights[ixh].i1 + 1] +
		xweights[ixh].w2 * ip[yweights[iyh].i2 + xweights[ixh].i2 + 1]
	    );
	    op[ixh + ixh + 1] *= scl;
	}
	if (ixh != nxof - 1) {
	    if (nxo_even && nyo_even) {
		for (ixo2 = ixh + ixh + 2; ixo2 < nxof + nxof; ixo2 += 2) {
		    op[ixo2] = op[nxof + nxof - ixo2];
		    op[ixo2 + 1] = -op[nxof + nxof - ixo2 + 1];
		}
	    } else {
		for (ixo = ixh + 1, ixo2 = ixh + ixh + 2;
		     ixo2 < nxof + nxof;
		     ++ixo, ixo2 += 2) {
		    op[ixo2] = xweights[ixo].w1 *
			ip[yweights[iyh].i1 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyh].i1 + xweights[ixo].i2];
		    op[ixo2] *= yweights[iyh].w1;
		    op[ixo2] += yweights[iyh].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iyh].i2 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyh].i2 + xweights[ixo].i2]
		    );
		    op[ixo2] *= scl;
		    op[ixo2 + 1] = xweights[ixo].w1 *
			ip[yweights[iyh].i1 + xweights[ixo].i1 + 1] +
			xweights[ixo].w2 *
			ip[yweights[iyh].i1 + xweights[ixo].i2 + 1];
		    op[ixo2 + 1] *= yweights[iyh].w1;
		    op[ixo2 + 1] += yweights[iyh].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iyh].i2 + xweights[ixo].i1 + 1] +
			xweights[ixo].w2 *
			ip[yweights[iyh].i2 + xweights[ixo].i2 + 1]
		    );
		    op[ixo2 + 1] *= scl;
		}
	    }
	    op += nxof + nxof;
	    if (nxo_even) {
		for (iyo = iyh + 1; iyo < nyof; ++iyo) {
		    const float* src = op + (nyof - iyo - iyo) * (nxof + nxof);

		    op[0] = src[0];
		    op[1] = -src[1];
		    for (ixo2 = 2; ixo2 < ixh + ixh; ixo2 += 2) {
			op[ixo2] = src[nxof + nxof - ixo2];
			op[ixo2 + 1] = -src[nxof + nxof - ixo2 + 1];
		    }
		    op[ixh + ixh] = xweights[ixh].w1 *
			ip[yweights[iyo].i1 + xweights[ixh].i1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i1 + xweights[ixh].i2];
		    op[ixh + ixh] *= yweights[iyo].w1;
		    op[ixh + ixh] += yweights[iyo].w2 * (
			xweights[ixh].w1 *
			ip[yweights[iyo].i2 + xweights[ixh].i1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i2 + xweights[ixh].i2]
		    );
		    op[ixh + ixh] *= scl;
		    op[ixh + ixh + 1] = xweights[ixh].w1 *
			ip[yweights[iyo].i1 + xweights[ixh].i1 + 1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i1 + xweights[ixh].i2 + 1];
		    op[ixh + ixh + 1] *= yweights[iyo].w1;
		    op[ixh + ixh + 1] += yweights[iyo].w2 * (
			xweights[ixh].w1 *
			ip[yweights[iyo].i2 + xweights[ixh].i1 + 1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i2 + xweights[ixh].i2 + 1]
		    );
		    op[ixh + ixh + 1] *= scl;
		    for (ixo2 = ixh + ixh + 2; ixo2 < nxof + nxof; ixo2 += 2) {
			op[ixo2] = src[nxof + nxof - ixo2];
			op[ixo2 + 1] = -src[nxof + nxof - ixo2 + 1];
		    }
		    op += nxof + nxof;
		}
	    } else {
		for (iyo = iyh + 1; iyo < nyof; ++iyo) {
		    const float* src = op + (nyof - iyo - iyo) * (nxof + nxof);

		    op[0] = src[0];
		    op[1] = -src[1];
		    for (ixo2 = 2; ixo2 < nxof + nxof; ixo2 += 2) {
			op[ixo2] = src[nxof + nxof - ixo2];
			op[ixo2 + 1] = -src[nxof + nxof - ixo2 + 1];
		    }
		    op += nxof + nxof;
		}
	    }
	} else {
	    op += nxof + nxof;
	    for (iyo = iyh + 1; iyo < nyof; ++iyo) {
		const float* src = op + (nyof - iyo - iyo) * (nxof + nxof);

		op[0] = src[0];
		op[1] = -src[1];
		for (ixo = 1, ixo2 = 2; ixo < nxof; ++ixo, ixo2 += 2) {
		    op[ixo2] = xweights[ixo].w1 *
			ip[yweights[iyo].i1 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i1 + xweights[ixo].i2];
		    op[ixo2] *= yweights[iyo].w1;
		    op[ixo2] += yweights[iyo].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iyo].i2 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i2 + xweights[ixo].i2]
		    );
		    op[ixo2] *= scl;
		    op[ixo2 + 1] = xweights[ixo].w1 *
			ip[yweights[iyo].i1 + xweights[ixo].i1 + 1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i1 + xweights[ixo].i2 + 1];
		    op[ixo2 + 1] *= yweights[iyo].w1;
		    op[ixo2 + 1] += yweights[iyo].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iyo].i2 + xweights[ixo].i1 + 1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i2 + xweights[ixo].i2 + 1]
		    );
		    op[ixo2 + 1] *= scl;
		}
		op += nxof + nxof;
	    }
	}
    } else {
	for (iyo = 0; iyo < nyof; ++iyo) {
	    int ixo2;

	    for (ixo = 0, ixo2 = 0; ixo < nxof; ++ixo, ixo2 += 2) {
		op[ixo2] =
		    xweights[ixo].w1 * ip[yweights[iyo].i1 + xweights[ixo].i1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i1 + xweights[ixo].i2];
		op[ixo2] *= yweights[iyo].w1;
		op[ixo2] += yweights[iyo].w2 * (
		    xweights[ixo].w1 * ip[yweights[iyo].i2 + xweights[ixo].i1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i2 + xweights[ixo].i2]
		);
		op[ixo2] *= scl;
		op[ixo2 + 1] =
		    xweights[ixo].w1 * ip[yweights[iyo].i1 + xweights[ixo].i1 + 1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i1 + xweights[ixo].i2 + 1];
		op[ixo2 + 1] *= yweights[iyo].w1;
		op[ixo2 + 1] += yweights[iyo].w2 * (
		    xweights[ixo].w1 * ip[yweights[iyo].i2 + xweights[ixo].i1 + 1] +
		    xweights[ixo].w2 * ip[yweights[iyo].i2 + xweights[ixo].i2 + 1]
		);
		op[ixo2 + 1] *= scl;
	    }
	    op += nxof + nxof;
	}
    }

    free(xweights);
    free(yweights);

    return 0;
}


int
interpolate_2d_r(
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
)
{
    /* Uses a bilinear interpolation method. */
    Interp1DSourceWeights* xweights;
    Interp1DSourceWeights* yweights;
    const float* src = op;
    int ixh = nxo / 2;
    int nxif, nxof, nyof;
    int iyo, ixo;
    int iymax;

    if ((flags & INTERP_X_SPLIT) != 0) {
	nxif = nxi;
	nxof = nxo;
	xweights = setup_1d_interp_split(nxif, dkxi * mag, nxof, dkxo);
    } else {
	nxif = nxi / 2 + 1;
	nxof = nxo / 2 + 1;
	xweights = setup_1d_interp(nxif, dkxi * mag, nxof, dkxo);
    }
    if ((flags & INTERP_Y_SPLIT) != 0) {
	nyof = nyo;
	yweights = setup_1d_interp_split(nyi, dkyi * mag, nyof, dkyo);
    } else {
	nyof = nyo / 2 + 1;
	yweights = setup_1d_interp(nyi / 2 + 1, dkyi * mag, nyof, dkyo);
    }
    if (xweights == 0 || yweights == 0) {
	if (xweights != 0) {
	    free(xweights);
	} else if (yweights != 0) {
	    free(yweights);
	}
	return 1;
    }

    for (iyo = 0; iyo < nyof; ++iyo) {
	yweights[iyo].i1 *= nxif;
	yweights[iyo].i2 *= nxif;
    }

    if (nxof > 0) {
	for (ixo = 0; ixo < nxof / 2 + 1; ++ixo) {
	    op[ixo] =
		xweights[ixo].w1 * ip[yweights[0].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i1 + xweights[ixo].i2];
	    op[ixo] *= yweights[0].w1;
	    op[ixo] += yweights[0].w2 * (
		xweights[ixo].w1 * ip[yweights[0].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i2 + xweights[ixo].i2]
	    );
	    op[ixo] *= scl;
	}
    }
    if ((flags & (INTERP_X_SPLIT | INTERP_CONJ_SYM)) ==
	(INTERP_X_SPLIT | INTERP_CONJ_SYM)) {
	for (ixo = nxof / 2 + 1; ixo < nxof; ++ixo) {
	    op[ixo] = op[nxof - ixo];
	}
    } else {
	for (ixo = nxof / 2 + 1; ixo < nxof; ++ixo) {
	    op[ixo] =
		xweights[ixo].w1 * ip[yweights[0].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i1 + xweights[ixo].i2];
	    op[ixo] *= yweights[0].w1;
	    op[ixo] += yweights[0].w2 * (
		xweights[ixo].w1 * ip[yweights[0].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[0].i2 + xweights[ixo].i2]
	    );
	    op[ixo] *= scl;
	}
    }
    op += nxof;

    if ((flags & INTERP_CONJ_SYM) != 0) {
	iymax = nyo / 2;
    } else {
	iymax = nyof;
    }

    for (iyo = 1; iyo < iymax; ++iyo) {
	for (ixo = 0; ixo < nxof; ++ixo) {
	    op[ixo] =
		xweights[ixo].w1 * ip[yweights[iyo].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[iyo].i1 + xweights[ixo].i2];
	    op[ixo] *= yweights[iyo].w1;
	    op[ixo] += yweights[iyo].w2 * (
		xweights[ixo].w1 * ip[yweights[iyo].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[iyo].i2 + xweights[ixo].i2]
	    );
	    op[ixo] *= scl;
	}
	op += nxof;
    }

    if ((flags & INTERP_CONJ_SYM) != 0) {
	for (ixo = 0; ixo <= ixh; ++ixo) {
	    op[ixo] = xweights[ixo].w1 *
		ip[yweights[iymax].i1 + xweights[ixo].i1] +
		xweights[ixo].w2 * ip[yweights[iymax].i1 + xweights[ixo].i2];
	    op[ixo] *= yweights[iymax].w1;
	    op[ixo] += yweights[iymax].w2 * (
		xweights[ixo].w1 *
		ip[yweights[iymax].i2 + xweights[ixo].i1] +
		xweights[ixo].w2 *
		ip[yweights[iymax].i2 + xweights[ixo].i2]
	    );
	    op[ixo] *= scl;
	}
	if ((flags & INTERP_X_SPLIT) != 0) {
	    if (nxo % 2 == 0 && nyo % 2 == 0) {
		for (ixo = ixh + 1; ixo < nxof; ++ixo) {
		    op[ixo] = op[nxof - ixo];
		}
	    } else {
		for (ixo = ixh + 1; ixo < nxof; ++ixo) {
		    op[ixo] = xweights[ixo].w1 *
			ip[yweights[iymax].i1 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iymax].i1 + xweights[ixo].i2];
		    op[ixo] *= yweights[iymax].w1;
		    op[ixo] += yweights[iymax].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iymax].i2 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iymax].i2 + xweights[ixo].i2]
		    );
		    op[ixo] *= scl;
		}
	    }
	    op += nxof;
	    if (nxo % 2 == 0) {
		for (iyo = iymax + 1; iyo < nyof; ++iyo) {
		    const float* line = src + nxof * (nyof - iyo);

		    op[0] = line[0];
		    for (ixo = 1; ixo < ixh; ++ixo) {
			op[ixo] = line[nxof - ixo];
		    }
		    op[ixh] = xweights[ixh].w1 *
			ip[yweights[iyo].i1 + xweights[ixh].i1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i1 + xweights[ixh].i2];
		    op[ixh] *= yweights[iyo].w1;
		    op[ixh] += yweights[iyo].w2 * (
			xweights[ixh].w1 *
			ip[yweights[iyo].i2 + xweights[ixh].i1] +
			xweights[ixh].w2 *
			ip[yweights[iyo].i2 + xweights[ixh].i2]
		    );
		    op[ixh] *= scl;
		    for (ixo = ixh + 1; ixo < nxof; ++ixo) {
			op[ixo] = line[nxof - ixo];
		    }
		    op += nxof;
		}
	    } else {
		for (iyo = iymax + 1; iyo < nyof; ++iyo) {
		    const float* line = src + nxof * (nyof - iyo);

		    op[0] = line[0];
		    for (ixo = 1; ixo < nxof; ++ixo) {
			op[ixo] = line[nxof - ixo];
		    }
		    op += nxof;
		}
	    }
	} else {
	    op += nxof;
	    for (iyo = iymax + 1; iyo < nyof; ++iyo) {
		op[0] = src[nxof * (size_t) (nyof - iyo)];
		for (ixo = 1; ixo < nxof; ++ixo) {
		    op[ixo] = xweights[ixo].w1 *
			ip[yweights[iyo].i1 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i1 + xweights[ixo].i2];
		    op[ixo] *= yweights[iyo].w1;
		    op[ixo] += yweights[iyo].w2 * (
			xweights[ixo].w1 *
			ip[yweights[iyo].i2 + xweights[ixo].i1] +
			xweights[ixo].w2 *
			ip[yweights[iyo].i2 + xweights[ixo].i2]
		    );
		    op[ixo] *= scl;
		}
		op += nxof;
	    }
	}
    }

    free(xweights);
    free(yweights);

    return 0;
}


Interp1DSourceWeights*
setup_1d_interp(int n1, float d1, int n2, float d2)
{
    Interp1DSourceWeights* weights =
	(Interp1DSourceWeights*) malloc(n2 * sizeof(Interp1DSourceWeights));
    float clim;
    int i, ilim;

    if (weights == 0) {
	return 0;
    }

    clim = (n1 - 1) * d1 / d2;
    if (clim >= n2) {
	ilim = n2;
    } else {
	ilim = ceil(clim);
    }
    for (i = 0; i < ilim; ++i) {
	float c = i * d2 / d1;

	weights[i].i1 = floor(c);
	weights[i].i2 = weights[i].i1 + 1;
	weights[i].w2 = c - weights[i].i1;
	weights[i].w1 = 1.0f - weights[i].w2;
    }

    clim = n1 * d1 / d2;
    if (clim >= n2) {
	ilim = n2;
    } else {
	ilim = ceil(clim);
    }
    for (; i < ilim; ++i) {
	weights[i].i1 = n1 - 1;
	weights[i].i2 = weights[i].i1;
	weights[i].w1 = 1.0f - (i * d2 / d1 - weights[i].i1);
	weights[i].w2 = 0.0f;
    }

    for (; i < n2; ++i) {
	weights[i].i1 = n1 - 1;
	weights[i].i2 = weights[i].i1;
	weights[i].w1 = 0.0f;
	weights[i].w2 = 0.0f;
    }

    return weights;
}


Interp1DSourceWeights*
setup_1d_interp_split(int n1, float d1, int n2, float d2)
{
    Interp1DSourceWeights* weights =
	(Interp1DSourceWeights*) malloc(n2 * sizeof(Interp1DSourceWeights));
    float clim;
    int i, ilim;

    if (weights == 0) {
	return 0;
    }

    clim = (n1 / 2) * d1 / d2;
    if (clim > n2 / 2) {
	ilim = n2 / 2 + 1;
    } else {
	ilim = ceil(clim);
    }
    for (i = 0; i < ilim; ++i) {
	float c = i * d2 / d1;

	weights[i].i1 = floor(c);
	weights[i].i2 = weights[i].i1 + 1;
	weights[i].w2 = c - weights[i].i1;
	weights[i].w1 = 1.0f - weights[i].w2;
    }

    clim = (n1 / 2 + 1) * d1 / d2;
    if (clim > n2 / 2) {
	ilim = n2 / 2 + 1;
    } else {
	ilim = ceil(clim);
    }
    for (; i < ilim; ++i) {
	weights[i].i1 = n1 / 2;
	weights[i].i2 = weights[i].i1;
	weights[i].w1 = 1.0f - (i * d2 / d1 - weights[i].i1);
	weights[i].w2 = 0.0f;
    }

    clim = n2 + (n1 / 2 - n1) * d1 / d2;
    if (clim < n2 / 2 + 1) {
	ilim = n2 / 2 + 1;
    } else {
	ilim = ceil(clim);
    }
    for (; i < ilim; ++i) {
	weights[i].i1 = n1 / 2;
	weights[i].i2 = weights[i].i1;
	weights[i].w1 = 0.0f;
	weights[i].w2 = 0.0f;
    }

    if (clim < n2 / 2 + 1) {
	ilim = n2 / 2 + 1;
    } else {
	ilim = ceil(clim);
    }
    for (; i < ilim; ++i) {
	weights[i].i1 = n1 / 2 + 1;
	weights[i].i2 = weights[i].i1;
	weights[i].w1 = 1.0f - (weights[i].i1 - n1 - (i - n2) * d2 / d1);
	weights[i].w2 = 0.0f;
    }

    clim = n2 - d1 / d2;
    if (clim < n2 / 2 + 1) {
	ilim = n2 / 2 + 1;
    } else {
	ilim = ceil(clim);
    }
    for (; i < ilim; ++i) {
	float c = n1 + (i - n2) * d2 / d1;

	weights[i].i1 = floor(c);
	weights[i].i2 = weights[i].i1 + 1;
	weights[i].w2 = c - weights[i].i1;
	weights[i].w1 = 1.0f - weights[i].w2;
    }

    for (; i < n2; ++i) {
	weights[i].i1 = n1 - 1;
	weights[i].i2 = 0;
	weights[i].w1 = (n2 - i) * d2 / d1;
	weights[i].w2 = 1.0f - weights[i].w1;
    }

    return weights;
}
