#ifndef SEmath_h
#define SEmath_h

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEstruct.h"
#include "SEmaths.h"
// #include "mpi_info.h"

// #define PI 3.1415926
#define MAX_N_FACTORS 128
typedef unsigned char uint8_t;
typedef uint8_t byte;

// double sinc(double x);
double bessi0(double x);

double kaiser_windowed_sinc(double x, double dx, int r);
void laplace_filter(grid_cfg *sim, float ***in, float ***out);
void hilbert_by_conv(const float* d, const int n, float*  h);
float franuni (void);
float franuni2(float x1, float x2);
void sranuni (int seed);

float distance_2_dim(int sx, int sy, int gx, int gy, int isx, int isy, int igx, int igy);

SE_inline void swap_real(real* a, real* b)
{
    real tmp = *a;
    *a = *b;
    *b = tmp;
}

SE_inline void swapf(float *a, float *b)
{
    float tmp = *a;
    *a = *b;
    *b = tmp;
}

SE_inline void swap_reala(reala* a, reala* b)
{
    reala tmp = *a;
    *a = *b;
    *b = tmp;
}

SE_inline void sincos(real x, real* s, real* c)
{
    *s = sin(x);
    *c = cos(x);
}

int factor(uint64_t n0, uint64_t *factors);
void split_factors(uint64_t number, uint64_t* fct, size_t nfct);
uint64_t next_int_with_max_factor(uint64_t number, uint64_t max_factor, int64_t growth_limit);


#endif