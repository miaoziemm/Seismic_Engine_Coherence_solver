#ifndef SEsufft_h
#define SEsufft_h

#include "../../SEbasic/include/SEbasic.h"
#ifndef COMPLEX_H
#define COMPLEX_H
#include "../../SEbasic/include/SEcomplex.h"
#endif

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, se_complex z[]);
void pfarc (int isign, int n, float rz[], se_complex cz[]);
void pfacr (int isign, int n, se_complex cz[], float rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, se_complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, float rz[], se_complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, se_complex cz[], float rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, se_complex z[]);

#endif