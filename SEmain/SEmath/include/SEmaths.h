#ifndef SEmaths_h
#define SEmaths_h
#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEstruct.h"
#include "../../SEbasic/include/SEprint.h"

typedef double real;
typedef float reala;

SE_private double se_mind(double a, double b){
    return a < b ? a : b;
}
SE_private double se_maxd(double a, double b){
    return a > b ? a : b;
}

SE_private float se_minf(float a, float b){
    return a < b ? a : b;
}
SE_private float se_maxf(float a, float b){
    return a > b ? a : b;
}

SE_private real se_max3r(real a, real b, real c){
    real max = a;
    if(b > max) max = b;
    if(c > max) max = c;
    return max;
}


/** A small number, but not zero.
    Useful for floating point comparations */
#define EPSILON ((real)1.0E-20)

#define FNOTEQUAL(a, b) ( fabs((a)-(b)) > EPSILON )
#define FEQUAL(a, b) ( fabs((a)-(b)) < EPSILON )
#define FEQUALEPS(a, b, EPS) ( fabs((a)-(b)) < (EPS) )
#define FISZERO(a) ( fabs((a)) < EPSILON )
#define FISZEROEPS(a, EPS) ( fabs((a)) < (EPS) )

#define VERIFYDEQUAL(a,b) do {                                          \
        if(!dequal((double)(a),(double)(b))) {                      \
            SE_ERROR("Expected equal numbers, but got %f and %f", (a), (b)); \
        }                                                               \
    } while(0)

#define VERIFYFEQUAL(a,b) do {                                          \
        if(!fequal((float)(a),(float)(b))) {                        \
            SE_ERROR("Expected equal numbers, but got %f and %f", (a), (b)); \
        }                                                               \
    } while(0)

SE_inline int fuzzy_compare(double x, double y)
{
    if(FEQUAL(x, y))
        return 0;
    if(x < y) return -1;
    else      return 1;
}

/** uniformly distributed random values inside [0, 1) */
SE_inline real se_rand1()
{
    return ( (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, 1) */
SE_inline real se_rand1_r(unsigned int *seedp)
{
    return ( (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, x) */
SE_inline real se_rand(real x)
{
    return ( x * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, x) */
SE_inline real se_rand_r(real x, unsigned int *seedp)
{
    return ( x * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [x1, x2) */
SE_inline real se_rand2(real x1, real x2)
{
    return x1 + se_rand(x2-x1);
}

/** uniformly distributed random values inside [x1, x2) */
SE_inline real rand2_r(real x1, real x2, unsigned int *seedp)
{
    return x1 + se_rand_r(x2-x1, seedp);
}

/** uniformly distributed integer random values inside [0, n) */
SE_inline int se_rand_int(int n)
{
    int r = (int)(n * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed integer random values inside [0, n) */
SE_inline int se_rand_int_r(int n, unsigned int *seedp)
{
    int r = (int)(n * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed integer random values inside [i1, i2) */
SE_inline int se_rand_int2(int i1, int i2)
{
    return i1 + se_rand_int(i2 - i1);
}

/** uniformly distributed integer random values inside [i1, i2) */
SE_inline int se_rand_int2_r(int i1, int i2, unsigned int *seedp)
{
    return i1 + se_rand_int_r(i2 - i1, seedp);
}

/** uniformly distributed 64-bit integer random values inside [0, n) */
SE_inline int64_t se_rand_int64(int64_t n)
{
    int64_t r = (int64_t)(n * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed 64-bit integer random values inside [0, n) */
SE_inline int64_t se_rand_int64_r(int64_t n, unsigned int *seedp)
{
    int64_t r = (int64_t)(n * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

SE_inline uint64_t se_rand_uint64(uint64_t n)
{
    uint64_t r = (uint64_t)(n * (long double) rand() / ( (long double)RAND_MAX + (long double)1.0 ) );
    if(r >= n) r = n-1;
    return r;
}

SE_inline uint64_t se_rand_uint64_r(uint64_t n, unsigned int *seedp)
{
    uint64_t r = (uint64_t)(n * (long double) rand_r(seedp) / ( (long double)RAND_MAX + (long double)1.0 ) );
    if(r >= n) r = n-1;
    return r;
}

/** Gaussian ("normally") distributed random values with mean 0.0 and
    standard deviation 1.0.

    This uses the polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia,
    as described by Donald E. Knuth in The Art of Computer Programming, 
    Volume 2: Seminumerical Algoztrhms, section 3.4.1, subsection C,
    algoztrhm P. 
*/
real rand_normal();


SE_inline int is_integer(real d)
{
    if (FISZEROEPS(d,1e-10)) {
        return 1;
    }
    return FISZEROEPS(fabs( d - round(d) )/fabs(d),1e-10);
}

/* Next two functions, work only if the numbers are normalized */
/* SE_inline int fequal_mag(float f1, float f2, int nbits) */
/* { */
/*     if( f1 == f2 ) { */
/*         return 1; */
/*     } else { */
/*         union { */
/*             float f; */
/*             int i; */
/*         } u1, u2; */
/*         u1.f = f1; */
/*         u2.f = f2; */
/*         if( u1.i < 0 ) u1.i = 0x80000000 - u1.i; */
/*         if( u2.i < 0 ) u2.i = 0x80000000 - u2.i; */
/*         int di = u1.i - u2.i; */
/*         if( di < 0 ) di = -di; */
/*         return ( di <= nbits ); */
/*     } */
/* } */

/* SE_inline int dequal_mag(double f1, double f2, int nbits) */
/* { */
/*     if( f1 == f2 ) { */
/*         return 1; */
/*     } else { */
/*         union { */
/*             double f; */
/*             int64_t i; */
/*         } u1, u2; */
/*         u1.f = f1; */
/*         u2.f = f2; */
/*         if( u1.i < 0 ) u1.i = (int64_t)0x8000000000000000 - u1.i; */
/*         if( u2.i < 0 ) u2.i = (int64_t)0x8000000000000000 - u2.i; */
/*         int64_t di = u1.i - u2.i; */
/*         if( di < 0 ) di = -di; */
/*         return ( di <= nbits ); */
/*     } */
/* } */

SE_inline int fequal(float f1, float f2)
{
    if(isfinite(f1)) {
        return fabs(f1-f2) <= 0.00001f*se_minf(fabs(f1), fabs(f2));
    } else {
        int inf1 = isinf(f1);
        if(inf1) return ( inf1 == isinf(f2) );
        return isnan(f2);
    }
}

SE_inline int dequal(double f1, double f2)
{
    /*return dequal_mag(f1, f2, 10);*/
    if(isfinite(f1)) {
        return fabs(f1-f2) <= 0.000000000001*se_mind(fabs(f1), fabs(f2));
    } else {
        int inf1 = isinf(f1);
        if(inf1) return ( inf1 == isinf(f2) );
        return isnan(f2);
    }
}

SE_inline int fgreaterequal(float f1, float f2)
{
	return ((f1 > f2) || (fequal(f1,f2)));
}


SE_inline int dgreaterequal(double f1, double f2)

{
	return ((f1 > f2) || (dequal(f1,f2)));
}

SE_inline int sign(real v)
{
    return (signbit(v))?(-1):(1);
}

SE_inline real sinc(real v)
{
    if(fabs(v) < 1.0E-10) {
        return (real)1.0;
    } else {
        return sin(v) / v;
    }
}


reala percentile(real p, reala* samples, size_t ns);


SE_inline reala medianf(reala* samples, size_t ns)
{
    return percentile(0.5, samples, ns);
}


SE_inline int roundi(real v)
{
    return (int)round(v);
}

SE_inline int64_t roundi64(real v)
{
    return (int64_t)round(v);
}

SE_inline int floori(real v)
{
    return (int)floor(v);
}

SE_inline int64_t floori64(real v)
{
    return (int64_t)floor(v);
}
SE_inline int ceili(real v)
{
    return (int)ceil(v);
}

SE_inline int64_t ceili64(real v)
{
    return (int64_t)ceil(v);
}
/** Reproductible uniform distribution between 0 and n-1 */
SE_inline int rrand(unsigned int* seed, int n)
{
    *seed = (*seed) * 1103515245 + 12345;
    int rnd = (int)(((*seed)/65536) % 32768);
    int r = (n * rnd) / 32768;
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/* SE_inline real hypot(const real x, const real y) */
/* { */
/*     real xabs = fabs(x) ; */
/*     real yabs = fabs(y) ; */
/*     real min, max; */

/*     if( xabs < yabs ) { */
/*         min = xabs; */
/*         max = yabs; */
/*     } else { */
/*         min = yabs; */
/*         max = xabs; */
/*     } */

/*     if( min < EPSILON ) { */
/*         return max; */
/*     } else { */
/*         real u = min / max ; */
/*         return max * sqrt(1 + u * u); */
/*     } */
/* } */

SE_inline real hypot3(const real x, const real y, const real z)
{
    real xabs = fabs(x);
    real yabs = fabs(y);
    real zabs = fabs(z);
    real w = se_max3r(xabs, yabs, zabs);

    if( w < EPSILON ) {
        return (real)0.0;
    } else {
        real r = w * sqrt((xabs / w) * (xabs / w) +
                              (yabs / w) * (yabs / w) +
                              (zabs / w) * (zabs / w));
        return r;
    }
}

SE_inline double safe_sqrt(double x){
    if(x > 0.0)
        return sqrt(x);
    else
        return 0.0;
}

SE_inline float safe_sqrtf(float x){
    if(x > 0.0f)
        return sqrtf(x);
    else
        return 0.0f;
}

// defined for x in [0,1)
SE_inline double smooth_cutoff(double x) {
	return 1.0-exp(1.0-1.0/(1.0-x));
}

// Tuckey windows for FFT

SE_inline double tukey_time(double x, double f) {
	double w = 1.0f;
	if      (      x  < f*0.5) w = 0.5 * (1.0 + cos(M_PI/f*(x -       f)));
	else if ((1.0 -x) < f*0.5) w = 0.5 * (1.0 + cos(M_PI/f*(x - 1.0 + f)));
	return w;
}

SE_inline double tukey_freq(double x, double f) {
	double w = 1.0;
	double f1 = f / 5.0;
	if (f1 > 1.0 - f) f1 = 1.0 - f;
	if (x > f) w = 0.5 * (1.0 + cos(M_PI/f1*(x -f)));
	return w;
}

SE_inline void dips2polarangles(const reala dx, const reala dy, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
    *theta = (reala)acos(1.0f / sqrtf(1.0f + dx * dx + dy *dy));
	if (!FISZERO(rho2))  *phi =  (reala)atan2(-dy,-dx);
	else *phi  = 0.0f;
}


SE_inline void polarcoords(const reala dz, const reala dx, const reala dy, reala *rho, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
	*theta = 0.0f;
	*phi  = 0.0f;
	*rho  = 0.0f;
	if (!FISZERO(rho2)) {
		reala r = sqrtf(dz * dz + dx * dx + dy *dy);
	    *theta = (reala)acos(dz / r);
	    *phi =  (reala)atan2(dy,dx);
	    *rho = r;
	}
}

SE_inline void polarangles(const reala dz, const reala dx, const reala dy, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
	*theta = 0.0f;
	*phi  = 0.0f;
	if (!FISZERO(rho2)) {
	    *theta = (reala)acos(dz / sqrtf(dz * dz + dx * dx + dy *dy));
	    *phi =  (reala)atan2(dy,dx);
	}
}
/*
  A is a 2 x 2 matrix in column-major order.
  Returns ret, the dimension of the nullspace of A using relative tolerance rtol.
  The first ret columns of V (also 2 x 2 in column-major order) form a basis for N(A).
*/
SE_inline int null2(const double *A, double *V, double rtol) {
    double B[4];
    double sum=0.0, tol;
    int i;

    memcpy(B,A,4*sizeof(double));
    memset(V,0,4*sizeof(double));

    // compute Frobenius norm of A
    for(i=0;i<4;i++) {
        sum+=B[i]*B[i];
    }
    tol=sqrt(sum)*rtol;

    // row pivot, if needed
    if (fabs(B[0]) < fabs(B[1])) {
        double t;
        t=B[0]; B[0]=B[1]; B[1]=t;
        t=B[2]; B[2]=B[3]; B[3]=t;
    }

    // eliminate entries and construct basis
    if (fabs(B[0]) > tol) {
        double r=B[1]/B[0];

        B[1]-=r*B[0];
        B[3]-=r*B[2];

        if (fabs(B[3]) < tol) {
            double nrmv;
            V[1]=1.0;
            V[0]=-B[2]*V[1]/B[0];
            nrmv=sqrt(V[0]*V[0]+V[1]*V[1]);
            V[0]/=nrmv;
            V[1]/=nrmv;
            return 1;
        }
        else {
            return 0;
        }
    }
    else if (fabs(B[2]) > tol || fabs(B[3]) > tol) {
        V[0] = 1.0;
        return 1;
    }
    else {
        V[0]=1.0;
        V[3]=1.0;
        return 2;
    }
}

#define UPDATE_PINV(lambda)                        \
{                                                  \
    int dim, j;                                    \
    double sigma=sqrt(lambda);                 \
    B[0]=a-lambda;                                 \
    B[1]=b;                                        \
    B[2]=b;                                        \
    B[3]=c-lambda;                                 \
                                                   \
    dim=null2(B,V,eps);                        \
                                                   \
    for(j=0;j<dim;j++) {                           \
        double *v=V+2*j;                           \
        double u[2];                               \
                                                   \
        u[0]=(A[0]*v[0]+A[2]*v[1])/sigma;          \
        u[1]=(A[1]*v[0]+A[3]*v[1])/sigma;          \
                                                   \
        X[0]+=(u[0]*v[0])/sigma;                   \
        X[1]+=(u[0]*v[1])/sigma;                   \
        X[2]+=(u[1]*v[0])/sigma;                   \
        X[3]+=(u[1]*v[1])/sigma;                   \
    }                                              \
}
/*
  Compute the pseudo-inverse X of the 2 x 2 matrix A (A and X both in column-major order),
  using the drop tolerance rtol.
 */
SE_inline void pinv2(const double *A, double *X, double rtol) {
    double B[4], V[4];
    double a,b,c,d;
    double l1, l2; // eigenvalues of B=A^T*A
    const double eps=2.22044604925031e-16;

    memset(X,0,4*sizeof(double));

    a=A[0]*A[0]+A[1]*A[1];
    b=A[0]*A[2]+A[1]*A[3];
    c=A[2]*A[2]+A[3]*A[3];
    d=sqrt((a-c)*(a-c)+4*b*b);

    l1=0.5*(a+c+d);
    l2=0.5*(a+c-d);

    if (FISZERO(l1)) {
        return;
    }
    else {
        UPDATE_PINV(l1);
    }

    if (!FEQUALEPS(l1,l2,eps) && l2 > rtol*l1) {
        UPDATE_PINV(l2);
    }
}
#undef UPDATE_PINV

#endif