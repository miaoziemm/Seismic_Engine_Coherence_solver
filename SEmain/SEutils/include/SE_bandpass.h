#ifndef _SE_BANDPASS_H_
#define _SE_BANDPASS_H_
#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEsys.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEgrid.h"

typedef struct bandpass_s {
    int n;
    double dt;
    double flo, fhi;
    int nplo, nphi;
    int phase;

    double orig_flo, orig_fhi;
    int orig_nplo, orig_nphi;

    double* b;
    double* d;
    int nb, nd;
    int noddhi, noddlo;

    int buff_len;
    double* buff;
}bandpass_t;

bandpass_t* bandpass_init( const double dt, 
                                   const double flo, const double fhi, 
                                   const int nplo, const int nphi, const int phase );

bandpass_t* bandpass_clone(const bandpass_t* src);

void bandpass_destroy( bandpass_t* bp );

void apply_bandpass( bandpass_t* bp, reala* data, int nsamples, int ntraces );

void bandpass_get_pars( const bandpass_t* bp, double* dt, 
                            double* flo, double* fhi, 
                            int* nplo, int* nphi, int* phase);


#endif