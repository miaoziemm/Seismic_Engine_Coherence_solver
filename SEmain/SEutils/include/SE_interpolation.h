#ifndef __SE_INTERPOLATION_H__
#define __SE_INTERPOLATION_H__

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEminmax.h"
#include "SE_rsf.h"


/**
 * Returns the value from a 2D grid by bilinear interpolation.
 */
SE_inline real rsf2d_value_binterp( real x, real y,
                                         const rsf2d_t* rsf )
{
    real    rx,ry;
    int64_t jx,jy;
    int64_t jx1, jy1;
    real    fx;

    if( !isfinite(x) ) x = rsf->g.x.o;
    if( !isfinite(y) ) y = rsf->g.y.o;


    /* the eps terms should make sure that (d*n)/d is evaluated at n */

    x = maxd(x, rsf->g.x.o + EPSILON);
    x = mind(x, rsf->g.x.o + (rsf->g.x.n-2)*rsf->g.x.d - EPSILON);
    rx = (x - rsf->g.x.o)/rsf->g.x.d;
    jx = (int64_t)(rx + EPSILON);
    fx = rx-jx;

    ASSERT(jx >=0 && jx <= rsf->g.x.n - 2);
    jx1 = jx+1; if(jx1 >= rsf->g.x.n) jx1 = rsf->g.x.n - 1;


    if(rsf->g.y.n > 1) {
        real  fy;
        real  v0, vx, vy, vxy;
        y = maxd(y, rsf->g.y.o + EPSILON);
        y = mind(y, rsf->g.y.o + (rsf->g.y.n-2)*rsf->g.y.d - EPSILON);
        ry = (y - rsf->g.y.o)/rsf->g.y.d;
        jy = (int64_t)(ry + EPSILON);
        fy = ry-jy;
        ASSERT(jy >=0 && jy <= rsf->g.y.n - 2);
        jy1 = jy+1; if(jy1 >= rsf->g.y.n) jy1 = rsf->g.y.n - 1;

        v0   = *I2(rsf, jx  ,jy  );
        vx   = *I2(rsf, jx1, jy  );
        vy   = *I2(rsf, jx  ,jy1);
        vxy  = *I2(rsf, jx1, jy1);
        
        return 
            (1-fy) *  ( (1-fx)*v0 + fx*vx  ) +
            fy     *  ( (1-fx)*vy + fx*vxy ) ;

    } else {
        real  v0, vx;
        v0   = *I2(rsf, jx , 0  );
        vx   = *I2(rsf, jx1, 0  );
        return (1-fx)*v0 + fx*vx;
    }
}

/**
 * Returns the value from a 3D grid by bilinear interpolation.
 */
typedef struct rsf3d_binterp_pt_s {
    real fx, fy, fz;
    int64_t jx, jx1;
    int64_t jy, jy1;
    int64_t jz, jz1;
} rsf3d_binterp_pt_t;

SE_inline void rsf3d_binterp_get_pt( real x, real y, real z,
                                          const rsf3d_t* rsf,
                                          rsf3d_binterp_pt_t *pt )
{
    if( !isfinite(x) ) x = rsf->g.x.o;
    if( !isfinite(y) ) y = rsf->g.y.o;
    if( !isfinite(z) ) z = rsf->g.z.o;

    ASSERT(rsf->g.x.n > 0);
    ASSERT(rsf->g.y.n > 0);
    ASSERT(rsf->g.z.n > 0);

    /* the eps terms should make sure that (d*n)/d is evaluated at n */

    pt->jx = (int64_t)floor( (x - rsf->g.x.o)/rsf->g.x.d );
    if(pt->jx < 0) {
        pt->jx = 0;
        if(rsf->g.x.n == 1) pt->jx1 = 0;
        else pt->jx1 = 1;
        pt->fx = 0;
        x = rsf->g.x.o;
    } else if(pt->jx >= rsf->g.x.n - 1) {
        pt->jx = pt->jx1 = rsf->g.x.n - 1;
        pt->fx = 0;
        x = rsf->g.x.o + pt->jx*rsf->g.x.d;
    } else {
        pt->jx1 = pt->jx + 1;
        pt->fx = (x - rsf->g.x.o)/rsf->g.x.d - pt->jx;
    }

    pt->jy = (int64_t)floor( (y - rsf->g.y.o)/rsf->g.y.d );

    if(pt->jy < 0) {
        pt->jy = 0;
        if(rsf->g.y.n == 1) pt->jy1 = 0;
        else pt->jy1 = 1;
        pt->fy = 0;
        y = rsf->g.y.o;
    } else if(pt->jy >= rsf->g.y.n - 1) {
        pt->jy = pt->jy1 = rsf->g.y.n - 1;
        pt->fy = 0;
        y = rsf->g.y.o + pt->jy*rsf->g.y.d;
    } else {
        pt->jy1 = pt->jy + 1;
        pt->fy = (y - rsf->g.y.o)/rsf->g.y.d - pt->jy;
    }

    pt->jz = (int64_t)floor( (z - rsf->g.z.o)/rsf->g.z.d );
    if(pt->jz < 0) {
        pt->jz = 0;
        if(rsf->g.z.n == 1) pt->jz1 = 0;
        else pt->jz1 = 1;
        pt->fz = 0;
        z = rsf->g.z.o;
    } else if(pt->jz >= rsf->g.z.n - 1) {
        pt->jz = pt->jz1 = rsf->g.z.n - 1;
        pt->fz = 0;
        z = rsf->g.z.o + pt->jz*rsf->g.z.d;
    } else {
        pt->jz1 = pt->jz + 1;
        pt->fz = (z - rsf->g.z.o)/rsf->g.z.d - pt->jz;
    }
}

SE_inline real rsf3d_value_binterp( real x, real y, real z,
                                         const rsf3d_t* rsf )
{
    rsf3d_binterp_pt_t pt;
    real    v0, vx, vy, vxy, vz, vxz, vyz, vxyz;

    rsf3d_binterp_get_pt(x, y, z, rsf, &pt);

    v0   = *I3(rsf, pt.jx  , pt.jy  , pt.jz  );
    vx   = *I3(rsf, pt.jx1 , pt.jy  , pt.jz  );
    vy   = *I3(rsf, pt.jx  , pt.jy1 , pt.jz  );
    vxy  = *I3(rsf, pt.jx1 , pt.jy1 , pt.jz  );
    vz   = *I3(rsf, pt.jx  , pt.jy  , pt.jz1 );
    vxz  = *I3(rsf, pt.jx1 , pt.jy  , pt.jz1 );
    vyz  = *I3(rsf, pt.jx  , pt.jy1 , pt.jz1 );
    vxyz = *I3(rsf, pt.jx1 , pt.jy1 , pt.jz1 );

    return
        (1-pt.fz) * ( (1-pt.fy) * ( (1-pt.fx)*v0  + pt.fx*vx   ) +
                       pt.fy     * ( (1-pt.fx)*vy  + pt.fx*vxy  ) ) +
        pt.fz     * ( (1-pt.fy) * ( (1-pt.fx)*vz  + pt.fx*vxz  ) +
                       pt.fy     * ( (1-pt.fx)*vyz + pt.fx*vxyz ) );
}

SE_inline real rsf3d_value_nbinterp( real x, real y, real z,
                                          const rsf3d_t* rsf )
{
    int64_t ix = axa_nidx_trim(&rsf->g.x, x);
    int64_t iy = axa_nidx_trim(&rsf->g.y, y);
    int64_t iz = axa_nidx_trim(&rsf->g.z, z);
    return *I3(rsf, ix, iy, iz);
}

/**
 * The S array contains 10 equally spaced points with the origin of
 * the X coordinates located at the 5th element.
 *
 * DX is the spacing in X between the points.
 *
 * Xo is the X location of the output value.
 *
 * output = a + b x + c x^2 + d x^3 + e x^4 + f x^5 + g x^6 + h x^7 + p x^8 + q x^9
 *
 * Inspired by a method worked out by Joe Higginbotham
 */
SE_inline double intp_poly9fit(double Xo, const float* S, float DX)
{
    static const double O8  = 1.0/8.0;
    static const double O16 = 1.0/16.0;
    static const double O24 = 1.0/24.0;
    static const double O32 = 1.0/32.0;
    static const double O40 = 1.0/40.0;
    static const double O48 = 1.0/48.0;
    static const double O56 = 1.0/56.0;
    static const double O72 = 1.0/72.0;
    static const double O80 = 1.0/80.0;
    double X;
    double S1, S3, S5, S7, S9;
    double T3, T5, T7, T9;
    double U5, U7, U9;
    double V7, V9;
    double p,g,e,c,a;
    double q,h,f,d,b;

    X=(2.0*Xo/DX)-1.0;
    S1=0.5*(S[5+0]+S[5-1]);
    S3=0.5*(S[5+1]+S[5-2]);
    S5=0.5*(S[5+2]+S[5-3]);
    S7=0.5*(S[5+3]+S[5-4]);
    S9=0.5*(S[5+4]+S[5-5]);

    T3=O8 *(S3-S1);
    T5=O24*(S5-S1);
    T7=O48*(S7-S1);
    T9=O80*(S9-S1);
    U5=O16*(T5-T3);
    U7=O40*(T7-T3);
    U9=O72*(T9-T3);
    V7=O24*(U7-U5);
    V9=O56*(U9-U5);

    p=O32*(V9-V7);
    g=V7-84.0*p;
    e=U5-35.0*g-996.0*p;
    c=T3-10.0*e-91.0*g-820.0*p;
    a=S1-c-e-g-p;

    S1=0.5000000*(S[5+0]-S[5-1]);
    S3=0.1666667*(S[5+1]-S[5-2]);
    S5=0.1000000*(S[5+2]-S[5-3]);
    S7=0.0714286*(S[5+3]-S[5-4]);
    S9=0.0555556*(S[5+4]-S[5-5]);

    T3=O8 *(S3-S1);
    T5=O24*(S5-S1);
    T7=O48*(S7-S1);
    T9=O80*(S9-S1);
    U5=O16*(T5-T3);
    U7=O40*(T7-T3);
    U9=O72*(T9-T3);
    V7=O24*(U7-U5);
    V9=O56*(U9-U5);

    q=O32*(V9-V7);
    h=V7-84.0*q;
    f=U5-35.0*h-996.0*q;
    d=T3-10.0*f-91.0*h-820.0*q;
    b=S1-d-f-h-q;

    return a+X*(b+X*(c+X*(d+X*(e+X*(f+X*(g+X*(h+X*(p+X*q))))))));
}

/**
 * Interpolate value of f(z) from regularly sampled f defined on axis az using poly9fit.
 */
SE_inline double interp_axa(const axa_t *az, const float *f, double z) {
    int64_t iz;

    if (!axa_x_in_bounds(az, z)) {
        return 0.0;
    }

    iz=mini64( maxi64( axa_closest_indx(az,z), 4), az->n-6);
    return intp_poly9fit(z-axa_x_from_idx(az, iz), f+iz-4, (float)az->d);
}

typedef void* trace_interpolator;

/**
 * input(nz,ntraces) - sampled at zz
 * zz(nz)
 * output(nout,ntraces) - with regular sampling of dzout
 */
trace_interpolator create_trace_interpolator(const float* zz, int nz,
                                                     float dzout, int nout,
                                                     int flag );

void interpolate_trace( trace_interpolator intp,
                            const float* input,
                            float* output,
                            int ntraces );

void destroy_trace_interpolator( trace_interpolator intp );

/**
 * Linear interpolation from array in[] sampled on axis ain to out[] sampled on axis aout
 * between axis values [xmin,xmax]. Does not extrapolate if part of the output axis falls
 * outside the boundaries of the input axis.
 */
void interp1_regular(const axa_t ain, const real *in,
                         const axa_t aout, real *out,
                         real xmin, real xmax);
void interp1_regular_float(const axa_t ain, const reala *in,
                               const axa_t aout, reala *out,
                               real xmin, real xmax);

/**
 * Given a function v(x) defined by an array of length lx via v_i=v(x[i]), i in [0,lx),
 * compute via linear interpolation the values vq_i=v(xq[i]), i in [0,lxq).
 * Assumes that the arrays x,xq are sorted. Does not extrapolate if part of the range
 * defined by xq is outside that defined by x.
 *
 * 2015/02 Ioan Sturzu: Flag out the start index for vq in order to enforce outside "no extrapolation"
 */
int interp1_irregular(const real *x, const real *v, const real *xq, real *vq, int64_t lx, int64_t lxq);
int interp1_irregular_float(const reala *x, const reala *v, const reala *xq, reala *vq, int64_t lx, int64_t lxq);

int check_coeff(const real x);


#endif