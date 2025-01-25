#ifndef SEV4D_M4D_h
#define SEV4D_M4D_h

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEstruct.h"
#include "SEmaths.h"
#include "SEmath.h"



SE_inline void v4d_init(v4d *v, const double v1, const double v2, const double v3, const double v4)
{
    v->v[0]=v1;
    v->v[1]=v2;
    v->v[2]=v3;
    v->v[3]=v4;
}

SE_inline void m4d_sym_fill_lower(m4d *M) {
    M->m[ 4]=M->m[ 1];
    M->m[ 8]=M->m[ 2];  M->m[ 9]=M->m[ 6];
    M->m[12]=M->m[ 3];  M->m[13]=M->m[ 7];  M->m[14]=M->m[11];
}
SE_inline void m4d_sym_init(m4d *M, 
                                 double m11, double m12, double m13,  double m14,
                                 double m22, double m23, double m24,
                                 double m33, double m34,
                                 double m44) 
{
    M->m[ 0]=m11;  M->m[ 1]=m12;  M->m[ 2]=m13;  M->m[3]=m14;
                   M->m[ 5]=m22;  M->m[ 6]=m23;  M->m[7]=m24;
                                  M->m[10]=m33;  M->m[11]=m34;
                                                 M->m[15]=m44;
    m4d_sym_fill_lower(M);
}
SE_inline void m4d_diag_init(m4d *M, 
                                  double m11, double m22, double m33, double m44) 
{
    memset(M->m, 0, 16*sizeof(double));
    M->m[ 0]=m11;
    M->m[ 5]=m22;
    M->m[10]=m33;
    M->m[15]=m44;
}


SE_inline void v4d_accum_scaled(v4d *u, const double d, const v4d *v)
{
    u->v[0] += d*v->v[0];
    u->v[1] += d*v->v[1];
    u->v[2] += d*v->v[2];
    u->v[3] += d*v->v[3];
}


/**
 * Get a matrix column as a vector. First column number is 1.
 */
SE_inline void m4d_v4d_get_col(v4d* v, const m4d *M, const int col)
{
    v->v[0] = M->m[     col-1];
    v->v[1] = M->m[4  + col-1];
    v->v[2] = M->m[8  + col-1];
    v->v[3] = M->m[12 + col-1];
}



SE_inline void m4d_v4d_mult(v4d *u, const m4d *M, const v4d *v) 
{
    u->v[0] = M->m[ 0] * v->v[0] + M->m[ 1] * v->v[1] + M->m[ 2] * v->v[2] + M->m[ 3] * v->v[3];
    u->v[1] = M->m[ 4] * v->v[0] + M->m[ 5] * v->v[1] + M->m[ 6] * v->v[2] + M->m[ 7] * v->v[3];
    u->v[2] = M->m[ 8] * v->v[0] + M->m[ 9] * v->v[1] + M->m[10] * v->v[2] + M->m[11] * v->v[3];
    u->v[3] = M->m[12] * v->v[0] + M->m[13] * v->v[1] + M->m[14] * v->v[2] + M->m[15] * v->v[3];
}


/**
 * Eigenvalues and eigenvectors of a symmetric 4x4 matrix. The matrix
 * E will contain the eigenvectors as columns in the same order as the
 * eigenvalues.  If E is NULL, the eigenvectors will not be
 * computed. ev is asumed to be an array of 4 doubles.
 */
void m4d_sym_eigensystem(const m4d *M, double *ev, m4d *E);

/**
 * Eigenvalues of a symmetric 4x4 matrix, ev is assumed to be array of 4.
 */
SE_inline void m4d_sym_eigenvals(const m4d *M, double *ev)
{
    m4d_sym_eigensystem(M,ev,NULL);
}

SE_inline void m4d_info(const m4d *M){
    SE_MESSAGE("Matrix 4x4:\n");
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",M->m[0], M->m[1], M->m[2], M->m[3]);
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",M->m[4], M->m[5], M->m[6], M->m[7]);
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",M->m[8], M->m[9], M->m[10], M->m[11]);
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",M->m[12], M->m[13], M->m[14], M->m[15]);
}

SE_inline void v4d_info(const v4d *v){
    SE_MESSAGE("Vector 4:\n");
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",v->v[0], v->v[1], v->v[2], v->v[3]);
}


#endif