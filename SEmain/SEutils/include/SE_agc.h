#ifndef _SE_AGC_H_
#define _SE_AGC_H_
#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEsys.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEgrid.h"

typedef struct agc_s {
    int window;
    int dwind;
    int detect;
    real threshold;

    reala* buff;
    int crt_buff_size;
}agc_t;


agc_t* agc_init(int window, int dwind, int detect, real threshold);

agc_t* agc_clone(const agc_t* src);

void agc_destroy( agc_t* agc );

void apply_agc( agc_t* agc, reala* data, int n );

void agc_get_pars( const agc_t* agc, int *window, int *dwind, int *detect, real *threshold);

#endif