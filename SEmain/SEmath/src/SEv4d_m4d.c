#include "../include/SEv4d_m4d.h"


SE_private double m4d_sym_find_largest_off_entry(const m4d *M, int *m_i, int *m_j)
{
    int i,j;
    double ent_norm;

    (*m_i)=0;
    (*m_j)=1;
    ent_norm = fabs(M->m[(*m_i)*4+(*m_j)]);

    for (i=0;i<4;i++) {
        for (j=i+1; j<4; j++) {
            if (ent_norm < fabs(M->m[i*4+j]) ) {
                ent_norm = fabs(M->m[i*4+j]);
                (*m_i) = i;
                (*m_j) = j;
            }
        }
    }
    
    return ent_norm;
}

SE_private double m4d_sym_find_largest_diag_entry(const m4d *M) 
{
    int i;
    double ent_norm;

    ent_norm=0;
    for (i=0;i<4;i++) {
        if (ent_norm < fabs(M->m[4*i+i])) {
            ent_norm = fabs(M->m[4*i+i]);
        }
    }

    return ent_norm;
}

/**
 * Use Jacobi eigenvalue algorithm: Apply Givens rotations to zero
 * out the non-diagonal entries -- not possible to completely zero
 * out the off-diagonal entries, but this algorithm boosts the size
 * of the diagonal entries compared to the off-diagonal entries.
 * The diagonal entries will be close to the eigenvalues, after several
 * iterations. The same rotations will convert the identity matrix
 * to a matrix that contains the eigenvectors as columns.
 * 
 */
void m4d_sym_eigensystem(const m4d *M, double *ev, m4d *E)
{
    int i,j,iter;
    m4d S, Sp;

    if (E) m4d_diag_init(E, 1, 1, 1, 1);

    memcpy(S.m, M->m, 16*sizeof(double));
    iter=0;
    while ( ( m4d_sym_find_largest_off_entry(&S, &i, &j) 
              > 1e-8 * m4d_sym_find_largest_diag_entry(&S)) 
            && ( iter<100 ) ) {

        int k;
        double th, s, c;

        memcpy(Sp.m, S.m, 16*sizeof(double));

        if (FISZERO(S.m[4*i+i]-S.m[4*j+j])) {
            th = M_PI/4;
        } else {
            th = .5*atan(2*S.m[4*i+j]/(S.m[4*j+j]-S.m[4*i+i]));
        }

        sincos(th,&s,&c);

        Sp.m[4*i+i]=c*c*S.m[4*i+i] - 2*s*c*S.m[4*i+j] + s*s*S.m[4*j+j];
        Sp.m[4*j+j]=s*s*S.m[4*i+i] + 2*s*c*S.m[4*i+j] + c*c*S.m[4*j+j];

        Sp.m[4*i+j]=(c*c-s*s)*S.m[4*i+j] + s*c*(S.m[4*i+i] - S.m[4*j+j]);
        Sp.m[4*j+i] = Sp.m[4*i+j];
        
        for (k=0; k<4; k++) {
           
            if (E) {
                double Eik, Ejk;

                Eik = c*E->m[4*k+i] - s*E->m[4*k+j];
                Ejk = s*E->m[4*k+i] + c*E->m[4*k+j];

                E->m[4*k+i] = Eik;
                E->m[4*k+j] = Ejk;
            }

            if ((k==i)||(k==j)) continue;

            Sp.m[4*i+k] = c*S.m[4*i+k] - s*S.m[4*j+k];
            Sp.m[4*k+i] = Sp.m[4*i+k];

            Sp.m[4*j+k] = s*S.m[4*i+k] + c*S.m[4*j+k];
            Sp.m[4*k+j] = Sp.m[4*j+k];
                
        }

        memcpy(S.m, Sp.m, 16*sizeof(double));
        iter++;

    }

    ev[0] = S.m[ 0];
    ev[1] = S.m[ 5];
    ev[2] = S.m[10];
    ev[3] = S.m[15];


}
