#ifndef SE_COHERENCE_CUT_H
#define SE_COHERENCE_CUT_H
#include "../SEmain/include/SEcore.h"

typedef struct SE_coherence_cut_cfg_s
{
    char *coherence_cube_name;
    char *coherence_slice_name;
    char *coherence_cut_cube_name;
    int nx;
    int ny;
    int nz;
    int direction; // 0: x, 1: y, 2: z
    int cut_point;

} SE_coherence_cut_cfg_t;

typedef struct SE_coherence_cut_slice_s
{
    int n1;
    int n2;
    int n_cut_direction;
    int direction;
    int cut_point;
    float *data;
} SE_coherence_cut_slice_t;

typedef struct SE_coherence_cut_cube_s
{
    int n1;
    int n2;
    int n3;
    int no1;
    int no2;
    int no3;//偏移量
    float *data;
} SE_coherence_cut_cube_t;




#endif // SE_COHERENCE_CUT_H