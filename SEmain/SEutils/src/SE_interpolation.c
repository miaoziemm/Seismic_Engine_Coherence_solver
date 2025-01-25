#include "../include/SE_interpolation.h"


#define NCOEF_LAGRANGE 12
#define NCOEF_P9LAGRANGE (NCOEF_LAGRANGE>10?NCOEF_LAGRANGE:10)
#define SC_LAGRANGE NCOEF_LAGRANGE-2
#define SC_P9 5

typedef enum {
    zero,
    nn,
    lgr,
    p9
} intptype;

typedef struct intpop_s {
    int beg;
    int end;
    intptype op;
    union {
        double lgrcoef[NCOEF_LAGRANGE];
    } opdata;
} intpop_t;

typedef struct interp_s {
    float* zz;
    int nz;

    float dzout;
    int nout;
    
    intpop_t* op;

    int flag;
} interp_t;

SE_inline void intp_find_range_idx(int i, int nidx, int nmax,
                                        int* beg, int* end, int sc)
{
    int nh;

    nh = nidx/2;
    if(sc < 0) sc = 0;

    if(i-nh < 0) {
        *beg = 0;
        *end = i + nidx - nh - 1;
        if(*end > 2*i + sc) *end = 2*i;
        if(*end > nmax-1) *end = nmax-1;
    } else {
        *beg = i - nh;
        *end = i + nidx - nh - 1;
        if(*end > nmax-1) {
            *end = nmax-1;
            if(i-*beg > *end-i+sc) {
                *beg = 2*i - *end;
                if(*beg < 0) *beg = 0;
            }
        }
    }
}

SE_inline void intp_init_idx(interp_t* intp)
{
    int i, k2 = 0;
    double z0 = intp->zz[0];
    intp->op[0].beg = 0;
    intp->op[0].end = 0;
    intp->op[0].op  = nn;

    for(i = 1; i < intp->nz; ++i) {
        int k, k1;

        k1 = (int)((intp->zz[i-1] - z0) / intp->dzout);
        k2 = (int)((intp->zz[i]   - z0) / intp->dzout);
        if(k2 >= intp->nout) k2 = intp->nout - 1;

        for( k = k1; k < k2; k++ ) {
            double z = z0 + k*intp->dzout;
            if( FEQUAL(z, intp->zz[i-1]) ) {
                intp->op[k].beg = i-1;
                intp->op[k].end = i-1;
                intp->op[k].op  = nn;
            } else if( FEQUAL(z, intp->zz[i]) ) {
                intp->op[k].beg = i;
                intp->op[k].end = i;
                intp->op[k].op  = nn;
            } else {
                int beg, end;
                int do_p9;
                intp_find_range_idx(i, 10, intp->nz, &beg, &end, SC_P9);
                do_p9 = ( end-beg+1 == 10 );
                if( do_p9 ) {
                    int j;
                    double dz;
                    dz = intp->zz[beg+1] - intp->zz[beg];
                    for(j = beg+2; j <= end; ++j) {
                        if(!FEQUAL(dz, (intp->zz[j] - intp->zz[j-1]))) {
                            do_p9 = 0;
                            break;
                        }
                    }
                    if(do_p9) {
                        intp->op[k].beg = beg;
                        intp->op[k].end = end;
                        intp->op[k].op  = p9;
                    }
                }
                if(!do_p9) {
                    int j;
                    const float* zz;

                    intp_find_range_idx(i, NCOEF_LAGRANGE,
                                            intp->nz, &beg, &end, SC_LAGRANGE);
                    intp->op[k].beg = beg;
                    intp->op[k].end = end;
                    intp->op[k].op  = lgr;

                    zz = intp->zz + beg;
                    for(j = 0; j <= end - beg; ++j) {
                        int l;
                        double c = 1.0;
                        double zzj = zz[j];
                        for(l = 0; l <= end - beg; ++l) {
                            double zzl = zz[l];
                            if(l != j) c *= (z - zzl)/(zzj - zzl);
                        }
                        intp->op[k].opdata.lgrcoef[j] = c;
                    }
                }
            }
        }
    }

    if(k2 >= 0) {
        for(i = k2; i < intp->nout; ++i) {
            intp->op[i].beg = -1;
            intp->op[i].end = -1;
            intp->op[i].op  = zero;
        }
    }
}

SE_inline void intp_do_p9_and_lagrange_one_trace(interp_t* intp,
                                                      const float* inp,
                                                      float* out)
{
    int i;
    double z0 = intp->zz[0];

    for(i = 0; i < intp->nout; ++i) {
        const intpop_t* op = intp->op+i;
        switch(op->op) {
        case p9: {
            double z = z0 + i*intp->dzout;
            out[i] = (float)intp_poly9fit(z-intp->zz[op->beg+4],
                                              inp+op->beg, 
                                              intp->zz[op->beg+1]-
                                              intp->zz[op->beg]);
            break;
        }
        case lgr: {
            double v = 0.0;
            const float* in = inp + op->beg;
            int j;
            for(j = 0; j <= op->end - op->beg; ++j) {
                v += in[j]*op->opdata.lgrcoef[j];
            }
            out[i] = (float)v;
            break;
        }
        case nn:
            out[i] = inp[op->beg];
            break;
        case zero:
            out[i] = 0.0;
            break;
        default:
            SE_ERROR("Unknown op: %d\n", op->op);
        }
    }
}


trace_interpolator create_trace_interpolator(const float* zz, int nz,
                                                     float dzout, int nout,
                                                     int flag )
{
    interp_t* intp;
    int i;

    /* first verify the arguments */
    if(zz == NULL) {
        SE_WARNING("Trace interpolator called with NULL zz");
        return NULL;
    }
    if(nz <= 0) {
        SE_WARNING("Trace interpolator called with nz=%d <= 0", nz);
        return NULL;
    }
    if(nout <= 0) {
        SE_WARNING("Trace interpolator called with nout=%d <= 0", nout);
        return NULL;
    }
    if(FISZERO(dzout)) {
        SE_WARNING("Trace interpolator called with zero dzout");
        return NULL;
    }
    if(dzout < 0) {
        SE_WARNING("Trace interpolator called with negative dzout");
        return NULL;
    }

    for(i = 1 ; i < nz; ++i) {
        float z1 = zz[i-1];
        float z2 = zz[i];
        if( FEQUAL(z1, z2) ) {
            SE_WARNING("Trace interpolator called with equal zz values; i=%d, z1=%f, z2=%f",
                  i, z1, z2);
            return NULL;
        }
        if( z2 < z1 ) {
            SE_WARNING("Trace interpolator called with unsorted zz values; i=%d, z1=%f, z2=%f",
                  i, z1, z2);
            return NULL;
        }
    }

    intp = (interp_t*)malloc(sizeof(*intp));
    memset(intp, 0, sizeof(*intp));

    intp->dzout = dzout;
    intp->nout  = nout;

    intp->nz = nz;
    intp->zz = alloc1float(intp->nz);
    memcpy( intp->zz, zz, intp->nz*sizeof(float) );

    intp->flag = flag;

    intp->op = (intpop_t*)malloc(intp->nout*sizeof(intpop_t));
    memset(intp->op, 0, intp->nout*sizeof(intpop_t));
    intp_init_idx(intp);

    return intp;
}


void interpolate_trace( trace_interpolator _intp,
                            const float* input,
                            float* output,
                            int ntraces )
{
    interp_t* intp = (interp_t*)_intp;
    int it;

    for(it = 0; it < ntraces; ++it) {
        const float* inp = input  + it*intp->nz;
        float* out = output +  it*intp->nout;
        intp_do_p9_and_lagrange_one_trace(intp, inp, out);
    }
}

void destroy_trace_interpolator( trace_interpolator _intp)
{
    interp_t* intp = (interp_t*)_intp;

    free(intp->zz);
    free(intp->op);
    free(intp);
}

void interp1_regular(const axa_t ain, const real *in,
                         const axa_t aout, real *out,
                         real xmin, real xmax) {
    int64_t i, imin, imax;
    int64_t jmin, jmax, j1, j2;

    imin = axa_closest_indx(&aout, xmin);
    imax = axa_closest_indx(&aout, xmax);
    jmin = axa_closest_indx(&ain, xmin);
    jmax = axa_closest_indx(&ain, xmax);
    j1 = jmin;
    j2 = jmin+1;

    for(i = imin; i <= imax; ++i) {
        while (j1 < jmax-1 && axa_x_from_idx(&ain, j1+1) < axa_x_from_idx(&aout, i)) j1++;
        while (j2 < jmax && axa_x_from_idx(&ain, j2) < axa_x_from_idx(&aout, i)) j2++;

        out[i] = in[j1] + (in[j2] - in[j1]) *
            ( axa_x_from_idx(&aout, i) - axa_x_from_idx(&ain, j1) ) / ain.d;
    }
}

void interp1_regular_float(const axa_t ain, const reala *in,
                               const axa_t aout, reala *out,
                               real xmin, real xmax) {
    int64_t i, imin, imax;
    int64_t jmin, jmax, j1, j2;

    imin = axa_closest_indx(&aout, xmin);
    imax = axa_closest_indx(&aout, xmax);
    jmin = axa_closest_indx(&ain, xmin);
    jmax = axa_closest_indx(&ain, xmax);
    j1 = jmin;
    j2 = jmin+1;

    for(i = imin; i <= imax; ++i) {
        while (j1 < jmax-1 && axa_x_from_idx(&ain, j1+1) < axa_x_from_idx(&aout, i)) j1++;
        while (j2 < jmax && axa_x_from_idx(&ain, j2) < axa_x_from_idx(&aout, i)) j2++;

        out[i] = (reala)( in[j1] + (in[j2] - in[j1]) *
                        ( axa_x_from_idx(&aout, i) - axa_x_from_idx(&ain, j1) ) / ain.d);
    }
}

int interp1_irregular(const real *x, const real *v, const real *xq, real *vq, int64_t lx, int64_t lxq) {
    int64_t i, j1=0, j2=1, istrt=0;

    while (xq[istrt] < x[0]) istrt++;
    for(i=istrt;i<lxq;i++) {
        while (j1 < lx-2 && x[j1+1] < xq[i]) j1++;
        while (j2 < lx-1 && x[j2] < xq[i]) j2++;

        // x[j1] and x[j2] bracket xq[i]
        vq[i]=v[j1]+(xq[i]-x[j1])/(x[j2]-x[j1])*(v[j2]-v[j1]);
    }
    return istrt;
}

int interp1_irregular_float(const reala *x, const reala *v, const reala *xq, reala *vq, int64_t lx, int64_t lxq) {
    int64_t i, j1=0, j2=1, istrt=0;

    while (xq[istrt] < x[0]) istrt++;
    for(i=istrt;i<lxq;i++) {
        while (j1 < lx-2 && x[j1+1] < xq[i]) j1++;
        while (j2 < lx-1 && x[j2] < xq[i]) j2++;

        // x[j1] and x[j2] bracket xq[i]
        vq[i]=v[j1]+(xq[i]-x[j1])/(x[j2]-x[j1])*(v[j2]-v[j1]);
    }
    return istrt;
}

int check_coeff(const real x) {
	if((x <0.0)||(x>1.0)) {
		SE_WARNING("Wrong value for convex combination coefficient %f",x);
		return 0;
	} else return 1;
}
