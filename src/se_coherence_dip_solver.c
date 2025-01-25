#include "SE_coherence.h"
#include "SE_coherence_math.h"

float bundle_dip_residual_1(float *params, SE_coherence_subcube_t *subcube)
{
    int dip_window = subcube->window_dip;
    int window_x = subcube->window_x;
    int window_y = subcube->window_y;
    int window_t = subcube->window_t;
    int center_it = (window_t - 1) / 2; // i1

    float dx = subcube->dx;
    float dy = subcube->dy;
    float dt = subcube->dt;
    int i, j, k;
    float shift;
    SE_coherence_subcube_t subcube_temp;
    subcube_temp.window_t = dip_window;
    subcube_temp.window_x = window_x;
    subcube_temp.window_y = window_y;
    subcube_temp.data = (float **)malloc(window_x * window_y * sizeof(float *));
    float res;

    float *data_moved;
    data_moved = (float *)malloc(window_x * window_y * dip_window * sizeof(float));
    for (i = 0; i < window_x * window_y; i++)
    {
        shift = params[0] * subcube->idx[i] * dx + params[1] * subcube->idy[i] * dy;
        shift = shift / dt;
        shift_trace_by_interp(shift, subcube->data[i], window_t, center_it, dip_window, data_moved + i * dip_window);
        subcube_temp.data[i] = data_moved + i * dip_window;
    }

    res = 1 - se_mar_semblance(&subcube_temp);
    free(data_moved);
    free(subcube_temp.data);
    return res;
}

// 差分进化算法主函数
void ncde_1(float *bounds, int popSize, int Max_gen, SE_coherence_subcube_t *subcube)
{
    float best_individual[2];
    float best_fitness;
    int i, j, k, l;
    int NP = popSize;
    int D = 2; // 修改维度为 2
    float F = 0.9;
    float CR = 0.1;
    int st = 2;
    int Max_FES = NP * Max_gen;
    int nfeval = 0, iter = 0, min_idx;
    float min_dist, dist, temp_val, *xl, *xu, *pop, *val, *tempval, *ui, bm[D], popold[D], newpop1[NP * D];

    xl = (float *)malloc(D * sizeof(float));
    xu = (float *)malloc(D * sizeof(float));
    for (i = 0; i < D; i++)
    {
        xl[i] = bounds[i * 2];
        xu[i] = bounds[i * 2 + 1];
    }

    pop = (float *)malloc(NP * D * sizeof(float));
    val = (float *)malloc(NP * sizeof(float));
    tempval = (float *)malloc(NP * sizeof(float));
    ui = (float *)malloc(NP * D * sizeof(float));

    // 初始化种群
    for (i = 0; i < NP; i++)
    {
        for (j = 0; j < D; j++)
        {
            pop[i * D + j] = xl[j] + (xu[j] - xl[j]) * rand_float();
        }
    }

    // 评估初始种群
    for (i = 0; i < NP; i++)
    {
        val[i] = bundle_dip_residual_1(&pop[i * D], subcube);
        nfeval++;
    }

    while (nfeval < Max_FES)
    {
        // 根据适应度进行排序
        for (i = 0; i < NP - 1; i++)
        {
            for (j = i + 1; j < NP; j++)
            {
                if (val[i] < val[j])
                {
                    temp_val = val[i];
                    val[i] = val[j];
                    val[j] = temp_val;

                    for (k = 0; k < D; k++)
                    {
                        float temp_pop = pop[i * D + k];
                        pop[i * D + k] = pop[j * D + k];
                        pop[j * D + k] = temp_pop;
                    }
                }
            }
        }

        for (j = 0; j < NP; j++)
        {
            for (k = 0; k < D; k++)
            {
                bm[k] = pop[j * D + k];
                popold[k] = pop[j * D + k];
            }

            for (k = 0; k < NP; k++)
            {
                for (l = 0; l < D; l++)
                {
                    newpop1[k * D + l] = pop[k * D + l];
                }
            }

            DE(&ui[j * D], popold, newpop1, bm, st, F, CR, D, NP, j, xl, xu);
            nfeval++;
        }

        // 评估新种群
        for (i = 0; i < NP; i++)
        {
            tempval[i] = bundle_dip_residual_1(&ui[i * D], subcube);
            nfeval++;
        }

        // 更新种群
        for (i = 0; i < NP; i++)
        {
            min_idx = 0;
            min_dist = INFINITY;
            for (j = 0; j < NP; j++)
            {
                dist = 0.0;
                for (k = 0; k < D; k++)
                {
                    dist += (ui[i * D + k] - pop[j * D + k]) * (ui[i * D + k] - pop[j * D + k]);
                }
                dist = sqrt(dist);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_idx = j;
                }
            }

            if (val[min_idx] > tempval[i])
            {
                for (k = 0; k < D; k++)
                {
                    pop[min_idx * D + k] = ui[i * D + k];
                }
                val[min_idx] = tempval[i];
            }
        }
        iter++;
    }



    // 找到最好的个体
    int best_idx = 0;
    for (i = 1; i < NP; i++)
    {
        if (val[i] > val[best_idx])
        {
            best_idx = i;
        }
    }

    for (j = 0; j < D; j++)
    {
        best_individual[j] = pop[best_idx * D + j];
    }
    best_fitness = val[best_idx];

    subcube->semblance_dip = 1 - best_fitness; // 将最佳适应度保存到subcube中
    
    subcube->px = best_individual[0];
    subcube->py = best_individual[1]; //将个体最佳参数保存到subcube中

    free(xl);
    free(xu);
    free(pop);
    free(val);
    free(tempval);
    free(ui);
}

void se_coherence_cube_compute_1(int dip_type, SE_coherence_cube_t *coherence_cube)
{
     int i, j, k;
    int ii, jj, kk;
    int ix, iy, it;
    int window_x, window_y, window_t;
    int id;
    float params[2];

    float dip_max = coherence_cube->dip_max;
    float bounds[4] = {-dip_max, dip_max, -dip_max, dip_max};
    float best_individual[2];
    float best_fitness;
    float *zero_temp;
    zero_temp = alloc1float(coherence_cube->nt);
    FILE *data_cube_fp;
    data_cube_fp = fopen(coherence_cube->data_cube_name, "rb");
    if (data_cube_fp == NULL)
    {
        SE_ERROR("Cannot open data cube file %s\n", coherence_cube->data_cube_name);
        exit(1);
    }

    float *buff;
    buff = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
    fread(buff, sizeof(float), coherence_cube->nx * coherence_cube->ny * coherence_cube->nt, data_cube_fp);

    float *px, *py;
    px = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        py = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);

    int complete_jobs_count = 0;

#pragma omp parallel for private(id, window_x, window_y, window_t, ii) schedule(dynamic)
    for (id = 0; id < coherence_cube->total_subcubes; id++)
    {
        coherence_cube->subcubes[id].data = (float **)malloc(coherence_cube->subcubes[id].window_x * coherence_cube->subcubes[id].window_y * sizeof(float *));
        window_t = coherence_cube->subcubes[id].window_t;
        window_x = coherence_cube->subcubes[id].window_x;
        window_y = coherence_cube->subcubes[id].window_y;
        for (ii = 0; ii < window_x * window_y; ii++)
        {
            if (coherence_cube->subcubes[id].idlist[ii] >= 0)
            {
                coherence_cube->subcubes[id].data[ii] = &buff[coherence_cube->subcubes[id].idlist[ii]];
            }
            else
            {
                coherence_cube->subcubes[id].data[ii] = &zero_temp[0];
            }
        }

        ncde_1(bounds, coherence_cube->popSize, coherence_cube->Max_gen, &coherence_cube->subcubes[id]);
            coherence_cube->coherence_distribution[id] = coherence_cube->subcubes[id].semblance_dip;
            if (coherence_cube->type == 0)
            {
                px[id] = coherence_cube->subcubes[id].px;
                py[id] = coherence_cube->subcubes[id].py;
            } // 只有semblance需要保存倾角

        free(coherence_cube->subcubes[id].data);

        #pragma omp atomic
        complete_jobs_count++;

        if (complete_jobs_count % 20000 == 0)
        {
            SE_MESSAGE("Complete %d/%d jobs\n", complete_jobs_count, coherence_cube->total_subcubes);
        }
    }
#pragma omp barrier
SE_MESSAGE("Complete %d/%d jobs\n", coherence_cube->total_subcubes, coherence_cube->total_subcubes);

 FILE *dip_fp;
        dip_fp = fopen("dip.temp", "w");
        for (i = 0; i < coherence_cube->total_subcubes; i++)
        {
            fprintf(dip_fp, "%f,%f\n", px[i], py[i]);
        }
        fclose(dip_fp);
    free(buff);
    free(px);
        free(py);
    free1float(zero_temp);
    fclose(data_cube_fp);

}

int main(int argc, char **argv)
{

    int nx,ny,nt;
    int dip_window;
    int window_x,window_y,window_t;
    float dx,dy,dt;
    int dip_type;
    SE_coherence_cube_t coherence_cube;

    if(!getparint("nx",&nx)) nx=100;
    if(!getparint("ny",&ny)) ny=100;
    if(!getparint("nt",&nt)) nt=100;
    
    if(!getparfloat("dx",&dx)) dx=10.0;
    if(!getparfloat("dy",&dy)) dy=10.0;
    if(!getparfloat("dt",&dt)) dt=0.004;
    if(!getparfloat("dip_max",&coherence_cube.dip_max)) coherence_cube.dip_max=0.2;
    coherence_cube.dx=dx;
    coherence_cube.dy=dy;
    coherence_cube.dt=dt;

    if(!getparint("pop_size",&coherence_cube.popSize)) coherence_cube.popSize=100;
    if(!getparint("max_gen",&coherence_cube.Max_gen)) coherence_cube.Max_gen=100;

    if(!getparint("half_window_x",&window_x)) window_x=1;
    if(!getparint("half_window_y",&window_y)) window_y=1;
    if(!getparint("half_window_t",&window_t)) window_t=10;
    if(!getparint("half_dip_window",&dip_window)) dip_window=window_t/2;
    if(dip_window>window_t){
        SE_WARNING("half_dip_window is larger than window_t, set dip_window = window_t/2\n");
        dip_window = window_t/2;
    }

    window_t = 2 * window_t + 1;
    window_x = 2 * window_x + 1;
    window_y = 2 * window_y + 1;
    dip_window = 2 * dip_window + 1;


    if(!getparstring("data_cube_file",&coherence_cube.data_cube_name)) coherence_cube.data_cube_name="inputfile";
    if(!getparstring("coherence_cube_file",&coherence_cube.coherence_cube_name)) coherence_cube.coherence_cube_name="coherencefile";

    if(!getparint("dip_type",&dip_type)) dip_type=0; //0:no_dip 1:ncde 

    if(!getparint("along_layer",&coherence_cube.flag_along_layer)) coherence_cube.flag_along_layer=0;
    if(coherence_cube.flag_along_layer==1)
    {
        if(!getparstring("layer_file",&coherence_cube.layer_fname)) coherence_cube.layer_fname="layer.txt";
        if(!getparstring("layer_result_file",&coherence_cube.layer_result_fname)) coherence_cube.layer_result_fname="layer_result.txt";
    }

    SE_MESSAGE("Create coherence cube\n");
    se_coherence_cube_create(&coherence_cube,nx,ny,nt,dip_window,window_x,window_y,window_t);
    SE_MESSAGE("Create coherence cube done\n");
    se_coherence_cube_report(&coherence_cube);

    SE_MESSAGE("Compute coherence cube\n");
    se_coherence_cube_compute_1(dip_type,&coherence_cube);
    SE_MESSAGE("Compute coherence cube done\n");

    SE_MESSAGE("Destroy coherence cube\n");
    se_coherence_cube_destroy(&coherence_cube);
    SE_MESSAGE("Destroy coherence cube done\n");
 
}
   