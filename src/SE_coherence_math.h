#ifndef _SE_COHERENCE_MATH_H_
#define _SE_COHERENCE_MATH_H_
#include "../SEmain/include/SEcore.h"
int rand_int(int max);
float rand_float();
void gaussian_filter_1d(const float *input, float *output, int length, float sigma, int order);
void gaussian_filter_3d(const float *input, float *output, int nx, int ny, int nt, float sigma, int order);


#endif  