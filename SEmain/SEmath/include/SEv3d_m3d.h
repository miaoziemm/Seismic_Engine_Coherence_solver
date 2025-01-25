#ifndef SEV3D_M3D_h
#define SEV3D_M3D_h

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEstruct.h"
#include "SEmaths.h"
#include "SEmath.h"
// #include "mpi_info.h"

// Vector init
SE_inline void v3d_init_zero(v3d *v) {
    memset(v->v,0,3*sizeof(double));
}
SE_inline void v3d_init(v3d *v, const double v1, const double v2, const double v3)
{
    v->v[0]=v1;
    v->v[1]=v2;
    v->v[2]=v3;
}

// Vector access 

/**
 * Return the value of a specific vector entry. First entry number is 1.
1开始
 */
SE_inline double v3d_get_elem(const v3d *v,
                                   const int elem)
{
    return v->v[elem-1];
}
SE_inline void v3d_set_elem(v3d *v, double val, int elem)
{
    v->v[elem-1] = val;
}


// vectors -> vector ops
void v3d_cross(v3d *u, const v3d *v, const v3d *w);
void v3d_get_orthogonal(v3d *u, const v3d *v);
void v3d_get_orthogonal_pair(v3d *u, v3d *w, const v3d *v);
SE_inline void v3d_get_2d_orthogonal(v3d *u, const v3d  *v)  //正交
{
    u->v[0] = -v->v[2];
    u->v[1] = 0;
    u->v[2] = v->v[0];
}

SE_inline void v3d_assign_scaled(v3d *u, const double d, const v3d *v)
{
    u->v[0] = d*v->v[0];
    u->v[1] = d*v->v[1];
    u->v[2] = d*v->v[2];
}
SE_inline void v3d_add(v3d *u, const v3d *v, const v3d *w)
{
    u->v[0] = v->v[0] + w->v[0];
    u->v[1] = v->v[1] + w->v[1];
    u->v[2] = v->v[2] + w->v[2];
}
SE_inline void v3d_add_scaled(v3d *u, 
                                   double c, const v3d *v, 
                                   double d, const v3d *w)
{
    u->v[0] = c*v->v[0] + d*w->v[0];
    u->v[1] = c*v->v[1] + d*w->v[1];
    u->v[2] = c*v->v[2] + d*w->v[2];
}
SE_inline void v3d_subtract(v3d *u, const v3d *v, const v3d *w)
{
    u->v[0] = v->v[0] - w->v[0];
    u->v[1] = v->v[1] - w->v[1];
    u->v[2] = v->v[2] - w->v[2];
}
SE_inline void v3d_accum(v3d *u, const v3d *v)
{
    u->v[0] += v->v[0];
    u->v[1] += v->v[1];
    u->v[2] += v->v[2];
}
SE_inline void v3d_accum_scaled(v3d *u, const double d, const v3d *v)
{
    u->v[0] += d*v->v[0];
    u->v[1] += d*v->v[1];
    u->v[2] += d*v->v[2];
}
SE_inline void v3d_accum_doubles(v3d *u, double v1, double v2, double v3)
{
    u->v[0] += v1;
    u->v[1] += v2;
    u->v[2] += v3;
}

// vectors -> scalar ops
SE_inline double v3d_dot(const v3d *v, const v3d *w)  //点积
{
    return v->v[0]*w->v[0] + v->v[1]*w->v[1] + v->v[2]*w->v[2];
}

SE_inline double v3d_dist(const v3d *v, const v3d *w) //距离
{
    return sqrt( (v->v[0] - w->v[0]) * (v->v[0] - w->v[0]) + 
                     (v->v[1] - w->v[1]) * (v->v[1] - w->v[1]) +
                     (v->v[2] - w->v[2]) * (v->v[2] - w->v[2]) ) ;
}


// Vector ops
SE_inline double v3d_norm(const v3d *v)  //向量的范数
{
    return sqrt( v3d_dot(v,v) );
}
SE_inline double v3d_safe_over_norm(const v3d *v)
{
    double n = v3d_norm(v);
    return (FISZERO(n)? 1 : 1/n);
}
SE_inline void v3d_scale(v3d *v, double c)
{
    v->v[0] *= c;
    v->v[1] *= c;
    v->v[2] *= c;
}
/**
 * Normalize a vector: if |v|!=0 then v/=|v|.  //归一化向量
 */
SE_inline void v3d_normalize(v3d *v)
{
    double norm = v3d_norm(v);
    v3d_scale(v,FISZERO(norm)? 0 : 1.0/norm);
}
/**
 * Set the norm of a vector: if |v|!=0 then v*=n/|v|
 */
SE_inline void v3d_set_norm(v3d *v, const double n)
{
    v3d_normalize(v);
    v3d_scale(v,n);
}

SE_inline double v3d_cos(const v3d *v, const v3d *w)
{
    return ( v3d_dot(v,w) 
             *v3d_safe_over_norm(v)
             *v3d_safe_over_norm(w) );
}

/**
 * Return the angle in radiants between two vectors.
 */
SE_inline double v3d_angle(const v3d *v, const v3d *w) {
    return acos(v3d_dot(v,w)/(v3d_norm(v)*v3d_norm(w)));
}



// Matrix init
SE_inline void m3d_init_zero(m3d *M) {
    memset(M->m, 0, 9*sizeof(double));
}
SE_inline void m3d_array_init(m3d *M, const double* vals) {
    memcpy(M->m, vals, 9*sizeof(double));
}
SE_inline void m3d_init(m3d *M, 
                             double m11, double m12, double m13, 
                             double m21, double m22, double m23,
                             double m31, double m32, double m33) 
{
    M->m[0]=m11; M->m[1]=m12; M->m[2]=m13;
    M->m[3]=m21; M->m[4]=m22; M->m[5]=m23;
    M->m[6]=m31; M->m[7]=m32; M->m[8]=m33;
}
SE_inline void m3d_sym_fill_lower(m3d *M) {
    M->m[3]=M->m[1];
    M->m[6]=M->m[2]; M->m[7]=M->m[5];
}
SE_inline void m3d_sym_init(m3d *M,   //对称矩阵初始化
                                 double m11, double m12, double m13, 
                                 double m22, double m23, 
                                 double m33) 
{
    M->m[0]=m11; M->m[1]=m12; M->m[2]=m13;
                 M->m[4]=m22; M->m[5]=m23;
                              M->m[8]=m33;
    m3d_sym_fill_lower(M);
}
SE_inline void m3d_diag_init(m3d *M,   //对角阵
                                  double m11, double m22, double m33) 
{
    M->m[0]=m11; M->m[1]=0;   M->m[2]=0;
    M->m[3]=0;   M->m[4]=m22; M->m[5]=0;
    M->m[6]=0;   M->m[7]=0;   M->m[8]=m33;
}

// Matrix access 

/**
 * Return the value of a specific matrix entry. First row and column
 * numbers are 1.
 */
SE_inline double m3d_get_elem(const m3d *M,     //返回元素
                                   const int row, const int col)
{
    return M->m[3*(row-1) + (col-1)];
}


// Matrix -> scalar   //矩阵的迹
SE_inline double m3drace(const m3d *A) 
{
    return ( A->m[0] + A->m[4] + A->m[8] );
}
SE_inline double m3d_determinant(const m3d *A) 
{

    return ( A->m[0]*A->m[4]*A->m[8] + 
             A->m[1]*A->m[5]*A->m[6] + 
             A->m[2]*A->m[3]*A->m[7] - 
             A->m[2]*A->m[4]*A->m[6] - 
             A->m[0]*A->m[5]*A->m[7] - 
             A->m[1]*A->m[3]*A->m[8] );
}

// matrices -> matrix ops

/**
 * Find the inverse of a matrix, return 0 if the matrix is singular.
 */
int m3d_inverse(m3d *A_inv, const m3d *A);  //反矩阵
/**
 * Eigenvalues of a symmetric 3x3 matrix. ev is assumed to be an
 * allocated array of 3 double.
 */
void m3d_sym_eigenvals(const m3d *M, double *ev);
SE_inline void m3dransp(m3d *M, const m3d *A)  //转置
{
    M->m[0] = A->m[0]; M->m[1] = A->m[3]; M->m[2] = A->m[6];
    M->m[3] = A->m[1]; M->m[4] = A->m[4]; M->m[5] = A->m[7];
    M->m[6] = A->m[2]; M->m[7] = A->m[5]; M->m[8] = A->m[8];
}
SE_inline void m3d_add(m3d *M, const m3d *A, const m3d *B) {
    M->m[0] = A->m[0] + B->m[0]; M->m[1] = A->m[1] + B->m[1]; M->m[2] = A->m[2] + B->m[2];
    M->m[3] = A->m[3] + B->m[3]; M->m[4] = A->m[4] + B->m[4]; M->m[5] = A->m[5] + B->m[5];
    M->m[6] = A->m[6] + B->m[6]; M->m[7] = A->m[7] + B->m[7]; M->m[8] = A->m[8] + B->m[8];
}
SE_inline void m3d_sym_add(m3d *M, const m3d *A, const m3d *B) { //对称加

    M->m[0] = A->m[0] + B->m[0]; M->m[1] = A->m[1] + B->m[1]; M->m[2] = A->m[2] + B->m[2];
                                 M->m[4] = A->m[4] + B->m[4]; M->m[5] = A->m[5] + B->m[5];
                                                              M->m[8] = A->m[8] + B->m[8];
    m3d_sym_fill_lower(M);
}
SE_inline void m3d_add_scaled(m3d *M, 
                                   double a, const m3d *A, 
                                   double b, const m3d *B) {
    M->m[0] = a*A->m[0] + b*B->m[0];  M->m[1] = a*A->m[1] + b*B->m[1];  M->m[2] = a*A->m[2] + b*B->m[2];
    M->m[3] = a*A->m[3] + b*B->m[3];  M->m[4] = a*A->m[4] + b*B->m[4];  M->m[5] = a*A->m[5] + b*B->m[5];
    M->m[6] = a*A->m[6] + b*B->m[6];  M->m[7] = a*A->m[7] + b*B->m[7];  M->m[8] = a*A->m[8] + b*B->m[8];
}
SE_inline void m3d_sym_add_scaled(m3d *M, 
                                       double a, const m3d *A, 
                                       double b, const m3d *B) {

    M->m[0] = a*A->m[0] + b*B->m[0];  M->m[1] = a*A->m[1] + b*B->m[1];  M->m[2] = a*A->m[2] + b*B->m[2];
                                      M->m[4] = a*A->m[4] + b*B->m[4];  M->m[5] = a*A->m[5] + b*B->m[5];  
                                                                        M->m[8] = a*A->m[8] + b*B->m[8];
    m3d_sym_fill_lower(M);
}
SE_inline void m3d_subtract(m3d *M, const m3d *A, const m3d *B) {
    M->m[0] = A->m[0] - B->m[0];  M->m[1] = A->m[1] - B->m[1];  M->m[2] = A->m[2] - B->m[2];
    M->m[3] = A->m[3] - B->m[3];  M->m[4] = A->m[4] - B->m[4];  M->m[5] = A->m[5] - B->m[5];
    M->m[6] = A->m[6] - B->m[6];  M->m[7] = A->m[7] - B->m[7];  M->m[8] = A->m[8] - B->m[8];
}
SE_inline void m3d_sym_subtract(m3d *M, const m3d *A, const m3d *B) {

    M->m[0] = A->m[0] - B->m[0];  M->m[1] = A->m[1] - B->m[1];  M->m[2] = A->m[2] - B->m[2];
                                  M->m[4] = A->m[4] - B->m[4];  M->m[5] = A->m[5] - B->m[5];
                                                                M->m[8] = A->m[8] - B->m[8];
    m3d_sym_fill_lower(M);
}
SE_inline void m3d_assign_scaled(m3d *M, const double c, const m3d *A) {
    M->m[0] = c*A->m[0];  M->m[1] = c*A->m[1];  M->m[2] = c*A->m[2];
    M->m[3] = c*A->m[3];  M->m[4] = c*A->m[4];  M->m[5] = c*A->m[5];
    M->m[6] = c*A->m[6];  M->m[7] = c*A->m[7];  M->m[8] = c*A->m[8];

}
SE_inline void m3d_sym_assign_scaled(m3d *M, const double c, const m3d *A) {
    M->m[0] = c*A->m[0];  M->m[1] = c*A->m[1];  M->m[2] = c*A->m[2];
                          M->m[4] = c*A->m[4];  M->m[5] = c*A->m[5];
                                                M->m[8] = c*A->m[8];
    m3d_sym_fill_lower(M);
}
SE_inline void m3d_accum(m3d *A, const m3d *B) {
    A->m[0] += B->m[0];  A->m[1] += B->m[1];  A->m[2] += B->m[2];
    A->m[3] += B->m[3];  A->m[4] += B->m[4];  A->m[5] += B->m[5];
    A->m[6] += B->m[6];  A->m[7] += B->m[7];  A->m[8] += B->m[8];
}
SE_inline void m3d_sym_accum(m3d *A, const m3d *B) {
    A->m[0] += B->m[0];  A->m[1] += B->m[1];  A->m[2] += B->m[2];
                         A->m[4] += B->m[4];  A->m[5] += B->m[5];
                                              A->m[8] += B->m[8];
    m3d_sym_fill_lower(A);
}
SE_inline void m3d_accum_scaled(m3d *A, double b, const m3d *B) {

    A->m[0] += b*B->m[0];  A->m[1] += b*B->m[1];  A->m[2] += b*B->m[2];
    A->m[3] += b*B->m[3];  A->m[4] += b*B->m[4];  A->m[5] += b*B->m[5];
    A->m[6] += b*B->m[6];  A->m[7] += b*B->m[7];  A->m[8] += b*B->m[8];
}
SE_inline void m3d_sym_accum_scaled(m3d *A, double b, const m3d *B) {
    A->m[0] += b*B->m[0];  A->m[1] += b*B->m[1];  A->m[2] += b*B->m[2];
                           A->m[4] += b*B->m[4];  A->m[5] += b*B->m[5];
                                                  A->m[8] += b*B->m[8];
    m3d_sym_fill_lower(A);
}

SE_inline void m3d_mult(m3d *M, const m3d *A, const m3d *B) { //乘法

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[3] = A->m[3]*B->m[0] + A->m[4]*B->m[3] + A->m[5]*B->m[6];
    M->m[4] = A->m[3]*B->m[1] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[6] = A->m[6]*B->m[0] + A->m[7]*B->m[3] + A->m[8]*B->m[6];
    M->m[7] = A->m[6]*B->m[1] + A->m[7]*B->m[4] + A->m[8]*B->m[7];
    M->m[8] = A->m[6]*B->m[2] + A->m[7]*B->m[5] + A->m[8]*B->m[8];

}
SE_inline void m3d_sym_mult(m3d *M, const m3d *A, const m3d *B) {

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[3] = A->m[3]*B->m[0] + A->m[4]*B->m[3] + A->m[5]*B->m[6];
    M->m[4] = A->m[1]*B->m[3] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[6] = A->m[6]*B->m[0] + A->m[7]*B->m[3] + A->m[8]*B->m[6];
    M->m[7] = A->m[6]*B->m[1] + A->m[7]*B->m[4] + A->m[8]*B->m[7];
    M->m[8] = A->m[2]*B->m[6] + A->m[5]*B->m[7] + A->m[8]*B->m[8];

}
/** 
 * Product of two matrices that is known a propri to be symmetric
 */
SE_inline void m3d_mult_sym_result(m3d *M, const m3d *A, const m3d *B) {

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[4] = A->m[3]*B->m[1] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[8] = A->m[6]*B->m[2] + A->m[7]*B->m[5] + A->m[8]*B->m[8];

    m3d_sym_fill_lower(M);
}
/**
 * Produces the M=H*C*H product of two symmetric matrices H and C. 
 */
SE_inline void m3d_mult_HCH_sym(m3d *M, const m3d *H, const m3d *C) {


    M->m[0] = ( + ( H->m[0]*C->m[0] + 2*H->m[1]*C->m[1] + 2*H->m[2]*C->m[2] )*H->m[0]
                + (                 +   H->m[1]*C->m[4] + 2*H->m[2]*C->m[5] )*H->m[1]
                + (                                     +   H->m[2]*C->m[8] )*H->m[2] );

    M->m[1] = ( + ( H->m[0]*C->m[0] + H->m[1]*C->m[1] + H->m[2]*C->m[2] )*H->m[1]
                + ( H->m[0]*C->m[1] + H->m[1]*C->m[4] + H->m[2]*C->m[5] )*H->m[4]
                + ( H->m[0]*C->m[2] + H->m[1]*C->m[5] + H->m[2]*C->m[8] )*H->m[5] );

    M->m[2] = ( + ( H->m[0]*C->m[0] + H->m[1]*C->m[1] + H->m[2]*C->m[2] )*H->m[2]
                + ( H->m[0]*C->m[1] + H->m[1]*C->m[4] + H->m[2]*C->m[5] )*H->m[5]
                + ( H->m[0]*C->m[2] + H->m[1]*C->m[5] + H->m[2]*C->m[8] )*H->m[8] );

    M->m[4] = ( + ( H->m[1]*C->m[0] + 2*H->m[4]*C->m[1] + 2*H->m[5]*C->m[2] )*H->m[1]
                + (                 +   H->m[4]*C->m[4] + 2*H->m[5]*C->m[5] )*H->m[4]
                + (                                     +   H->m[5]*C->m[8] )*H->m[5] );

    M->m[5] = ( + ( H->m[1]*C->m[0] + H->m[4]*C->m[1] + H->m[5]*C->m[2] )*H->m[2]
                + ( H->m[1]*C->m[1] + H->m[4]*C->m[4] + H->m[5]*C->m[5] )*H->m[5]
                + ( H->m[1]*C->m[2] + H->m[4]*C->m[5] + H->m[5]*C->m[8] )*H->m[8] );


    M->m[8] = ( + ( H->m[2]*C->m[0] + 2*H->m[5]*C->m[1] + 2*H->m[8]*C->m[2] )*H->m[2]
                + (                 +   H->m[5]*C->m[4] + 2*H->m[8]*C->m[5] )*H->m[5]
                + (                                     +   H->m[8]*C->m[8] )*H->m[8] );


    m3d_sym_fill_lower(M);
}

// Matrix ops

/**
 * Test if a matrix is symmetric. Return 0 if not.
 */
int m3d_is_sym(const m3d *A);
SE_inline void m3d_scale(m3d *M, const double c) {
    M->m[0] *= c;  M->m[1] *= c;  M->m[2] *= c;
    M->m[3] *= c;  M->m[4] *= c;  M->m[5] *= c;
    M->m[6] *= c;  M->m[7] *= c;  M->m[8] *= c;
}
SE_inline void m3d_sym_scale(m3d *M, const double c) {
    M->m[0] *= c;  M->m[1] *= c;  M->m[2] *= c;
                   M->m[4] *= c;  M->m[5] *= c;
                                  M->m[8] *= c;
    m3d_sym_fill_lower(M);
}



// Vectors -> Matrix

/**
 * Outer vector product: M = v^t * v .
 */
SE_inline void v3d_outer(m3d *M, const v3d *v, const v3d *w)
{
    M->m[0] = v->v[0] * w->v[0];  M->m[1] = v->v[0] * w->v[1];  M->m[2] = v->v[0] * w->v[2];
    M->m[3] = v->v[1] * w->v[0];  M->m[4] = v->v[1] * w->v[1];  M->m[5] = v->v[1] * w->v[2];
    M->m[6] = v->v[2] * w->v[0];  M->m[7] = v->v[2] * w->v[1];  M->m[8] = v->v[2] * w->v[2];
}


// // info
// SE_inline void v3d_info(int verb, const v3d *v, const char * str) {
//    SE_MESSAGE((verb, "%s -> (%g,%g,%g)",str, v->v[0], v->v[1], v->v[2]));
// } 
// SE_inline void m3d_info(int verb, const m3d *M, char *desc)
// {
//     INFOV((verb,
//            "%s\n"
//            "  [%+1.4e, %+1.4e, %+1.4e]\n"
//            "  [%+1.4e, %+1.4e, %+1.4e]\n"
//            "  [%+1.4e, %+1.4e, %+1.4e]\n",
//            desc,
//            M->m[0], M->m[1], M->m[2],
//            M->m[3], M->m[4], M->m[5],
//            M->m[6], M->m[7], M->m[8]));
// }

// Vector,Matrix -> Vector
SE_inline void m3d_v3d_mult(v3d *u, const m3d *M, const v3d *v) 
{
    u->v[0] = M->m[0] * v->v[0] + M->m[1] * v->v[1] + M->m[2] * v->v[2];
    u->v[1] = M->m[3] * v->v[0] + M->m[4] * v->v[1] + M->m[5] * v->v[2];
    u->v[2] = M->m[6] * v->v[0] + M->m[7] * v->v[1] + M->m[8] * v->v[2];
}

// Vectors, Matrix -> Scalar

/**
 * Return the the scalar product return <v,Mw>.
 */
SE_inline double v3d_m3d_v3d_inner(const v3d *v, const m3d *M, const v3d *w) 
{
    return (v->v[0] * ( M->m[0] * w->v[0] + M->m[1] * w->v[1] + M->m[2] * w->v[2] ) +
            v->v[1] * ( M->m[3] * w->v[0] + M->m[4] * w->v[1] + M->m[5] * w->v[2] ) +
            v->v[2] * ( M->m[6] * w->v[0] + M->m[7] * w->v[1] + M->m[8] * w->v[2] ) );
}
/*
 * Return the the scalar product return <v,Mv> for a symmetric matrix M.
返回标量积
 */
SE_inline double v3d_m3d_sym_inner(const m3d *M, const v3d *v) 
{

    return (v->v[0] * ( M->m[0] * v->v[0] + 2 * M->m[1] * v->v[1] + 2 * M->m[2] * v->v[2] ) +
            v->v[1] * (                         M->m[4] * v->v[1] + 2 * M->m[5] * v->v[2] ) +
            v->v[2] * (                                                 M->m[8] * v->v[2] ) );


}

SE_inline void m3d_info(const m3d *M){
    SE_MESSAGE("Matrix 3x3:\n");
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e]\n",M->m[0], M->m[1], M->m[2]);
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e]\n",M->m[3], M->m[4], M->m[5]);
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e]\n",M->m[6], M->m[7], M->m[8]);

}

SE_inline void v3d_info(const v3d *v){
    SE_MESSAGE("Vector 3:\n");
    SE_MESSAGE("  [%+1.4e, %+1.4e, %+1.4e]\n",v->v[0], v->v[1], v->v[2]);
}



#endif