#ifndef SE_COHERENCE_H
#define SE_COHERENCE_H

#include "../SEmain/include/SEcore.h"
#include <omp.h>

typedef struct SE_coherence_subcube_s
{
    int ix;
    int iy;
    int it;
    float dx;
    float dy;
    float dt;
    int window_x;
    int window_y;
    int window_t;
    int window_dip; //用于倾角估计的窗口大小
    int type; 

    int *idlist; //存储窗口数据索引
    float *idx;
    float *idy; //周边道索引
    float **data;

    float semblance_dip;
    float px;
    float py;

}SE_coherence_subcube_t;

typedef struct SE_coherence_cube_s
{

    char *data_cube_name;
    char *coherence_cube_name;
    int type; //0:semblance; 1:variation
    int nx;
    int ny;
    int nt;
    float dx;
    float dy;
    float dt;
   
    float *coherence_distribution; //结果

    int total_subcubes; //总共有多少个子快
    SE_coherence_subcube_t *subcubes;
    long max_memory_size;

    int dip_flag; //是否进行倾角矫正
    float dip_max; //倾角最大值
    char *dip_result_fname; //倾角结果文件
    
    int popSize; //种群大小
    int Max_gen; //最大迭代次数

    int flag_along_layer; //是否沿着层面进行计算, 若进行沿层计算，则不输出coherence_cube_name文件
    char *layer_fname; //沿层的坐标文件
    char *layer_result_fname; //沿层的结果文件

}SE_coherence_cube_t;

typedef struct
{
    float position[2];
    float fitness;
} Individual;

void se_coherence_cube_create(SE_coherence_cube_t *coherence_cube,int nx, int ny, int nt, int dip_window, int window_x, int window_y, int window_t);
void se_coherence_cube_compute(int dip_type,SE_coherence_cube_t *coherence_cube);
float se_mar_semblance(SE_coherence_subcube_t *subcube);
float se_var_semblance0(SE_coherence_subcube_t *subcube); //variation_0
float se_var_semblance1(SE_coherence_subcube_t *subcube); //variation_1
float se_var_semblance2(SE_coherence_subcube_t *subcube); //variation_2
float se_eigenstructure(SE_coherence_subcube_t *subcube); //eigenstructure
float se_gradiant_structure_tensor(SE_coherence_subcube_t *subcube); //gradiant_structure_tensor
void se_coherence_cube_destroy(SE_coherence_cube_t *coherence_cube);
void se_coherence_cube_report(SE_coherence_cube_t *coherence_cube);
float se_coherence_no_dip(int type, SE_coherence_subcube_t *subcube);

// dip_calculation
float bundle_dip_residual(float *params, SE_coherence_subcube_t *subcube);
void DE(float *ui, float *popold, Individual *pop, float *bm, int st, float F, float CR, int n, int NP, int j, float *xl, float *xu);
void ncde(float *bounds, int popSize, int Max_gen, SE_coherence_subcube_t *subcube);
void shift_trace_by_interp(float shift, float *data, int n1, int i1, int w1, float *shift_data);

void direct_dip(float *bounds, int popSize, int Max_gen, SE_coherence_subcube_t *subcube);


#endif