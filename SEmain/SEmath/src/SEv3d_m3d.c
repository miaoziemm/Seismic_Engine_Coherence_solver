#include "../include/SEv3d_m3d.h"


/*Matrix Vector*/

// vectors -> vector ops
void v3d_cross(v3d *u, const v3d *v, const v3d *w)
{
        
    u->v[0] =   v->v[1] * w->v[2] - v->v[2] * w->v[1];
    u->v[1] = - v->v[0] * w->v[2] + v->v[2] * w->v[0];
    u->v[2] =   v->v[0] * w->v[1] - v->v[1] * w->v[0];

    if (FISZERO(v3d_norm(u))) {

        if ( v3d_norm(v) > v3d_norm(w) ) {
            v3d_get_orthogonal(u,v);
        } else {
            v3d_get_orthogonal(u,w);            
        }
    }
}

void v3d_get_orthogonal(v3d *u, const v3d *v)
{
    int min_ind=1;
    double min=fabs(v->v[0]);
    if (min > fabs(v->v[1])) {
        min = fabs(v->v[1]);
        min_ind=2;
    }
    if (min > fabs(v->v[2])) {
        min = fabs(v->v[2]);
        min_ind=3;
    }
    
    switch(min_ind) {
    case 1:
        u->v[0]=0;
        u->v[1]= v->v[2];
        u->v[2]=-v->v[1];
        break;
    case 2:
        u->v[0]= v->v[2];
        u->v[1]=0;
        u->v[2]=-v->v[0];
        break;
    case 3:
        u->v[0]= v->v[1];
        u->v[1]=-v->v[0];
        u->v[2]=0;
        break;
    default:
        u->v[0]=0;
        u->v[1]=0;
        u->v[2]=0;
        break;
    } 
}

void v3d_get_orthogonal_pair(v3d *u, v3d *w, const v3d *v)
{
    v3d_get_orthogonal(u, v);
    v3d_cross(w, v, u);
}

// Matrix ops
int m3d_is_sym(const m3d *A)
{

    double diff;

    diff = ( FISZERO(fabs(A->m[1])) && FISZERO(fabs(A->m[3])) ) ? 0 : 
        ( A->m[1] - A->m[3] ) / ( fabs(A->m[1]) + fabs(A->m[3]) ) ; 
    if (fabs(diff)  > 1e-6) {
        SE_MESSAGE("Matrix is not symmetric entry (1,2) != (2,1)");
        return 0;
    }

    diff = ( FISZERO(fabs(A->m[2])) && FISZERO(fabs(A->m[6])) ) ? 0 : 
        ( A->m[2] - A->m[6] ) / ( fabs(A->m[2]) + fabs(A->m[6]) ) ; 
    if (fabs(diff)  > 1e-6) {
        SE_MESSAGE("Matrix is not symmetric entry (1,3) != (3,1)");
        return 0;
    }

    diff = ( FISZERO(fabs(A->m[5])) && FISZERO(fabs(A->m[7])) ) ? 0 : 
        ( A->m[5] - A->m[7] ) / ( fabs(A->m[5]) + fabs(A->m[7]) ) ; 
    if (fabs(diff)  > 1e-6) {
        SE_MESSAGE("Matrix is not symmetric entry (2,3) != (3,2)");
        return 0;
    }

    return 1;

}



// Matrix -> Matrix ops

static double get_2d_det(double a11, double a12, 
                              double a21, double a22) {
    return a11*a22 - a12*a21;
}

int m3d_inverse(m3d *A_inv, const m3d *A)
{

    double d=m3d_determinant(A);

    if (FISZERO(d)) {

        m3d_init_zero(A_inv);
        return 1;

    } else {

        d = 1/d;

        A_inv->m[0] = get_2d_det(A->m[4],A->m[5],A->m[7],A->m[8]) * d;
        A_inv->m[1] = get_2d_det(A->m[2],A->m[1],A->m[8],A->m[7]) * d;
        A_inv->m[2] = get_2d_det(A->m[1],A->m[2],A->m[4],A->m[5]) * d;

        A_inv->m[3] = get_2d_det(A->m[5],A->m[3],A->m[8],A->m[6]) * d;
        A_inv->m[4] = get_2d_det(A->m[0],A->m[2],A->m[6],A->m[8]) * d;
        A_inv->m[5] = get_2d_det(A->m[2],A->m[0],A->m[5],A->m[3]) * d;

        A_inv->m[6] = get_2d_det(A->m[3],A->m[4],A->m[6],A->m[7]) * d;
        A_inv->m[7] = get_2d_det(A->m[1],A->m[0],A->m[7],A->m[6]) * d;
        A_inv->m[8] = get_2d_det(A->m[0],A->m[1],A->m[3],A->m[4]) * d;

        return 0;
    }
}


/**
 * for the 3x3 symmetric case use the characteristic polynomial 
 * det(m*I - M) = m^3 - m^2 tr(M) - .5m (tr(M^2)-tr(M)^2) - det(M)
 * use the tranform M=f B + g I, M and B have the same eigenvectors  
 * and eigenvalues of M are eigenvalues of B times p plus q
 * use g = tr(M)/3 and f = tr((A-g I)^2/6)^.5
 * then 
 * det (b*I+B) = b^3 - 3b - det(B)
 * solve using b=2cos(3t)
 * get b = 2cos(1/3 (arccos(det(B)/2) + 2 pi n)) for n=0,1,2
 *
 */
void m3d_sym_eigenvals(const m3d *M, double *ev) 
{
    double f, g, h, r, t;
    m3d B;
    
    h = M->m[1]*M->m[1] + M->m[2]*M->m[2] + M->m[5]*M->m[5];
    g = m3drace(M)/3;
 
    // check if M is diagonal
    if ( (h/fabs(g))<1e-6 ) {
        ev[0] = M->m[0];
        ev[1] = M->m[4];
        ev[2] = M->m[8];
        return;
    }

    f = sqrt( ( (M->m[0]-g)*(M->m[0]-g) + (M->m[4]-g)*(M->m[4]-g) + (M->m[8]-g)*(M->m[8]-g) + 2*h ) / 6 );
    
    B = *M;
    B.m[0] -= g;
    B.m[4] -= g;
    B.m[8] -= g;
    
    r = m3d_determinant(&B) / (2*f*f*f);

    if (r<=-1) {
        t = M_PI/3.0;
    } else if (r>=1) {
        t = 0;
    } else {
        t = acos(r)/3.0;
    }
    
    ev[0] = g + 2*f*cos(t);
    ev[1] = g + 2*f*cos(t + 2.0*M_PI/3.0);
    ev[2] = 3*g - ev[0] - ev[1];

}