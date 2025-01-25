#include "SE_coherence.h"
#include "SE_coherence_math.h"

#define MAX_POP_SIZE 1000
#define MAX_DIM 100

typedef struct sort_double_s
{
    double p;
    int k;
} sort_double_t;


static int compare_sort_double(const void *a, const void *b) {
    if (((Individual *)a)->fitness > ((Individual *)b)->fitness) {
        return -1;
    } else if (((Individual *)a)->fitness < ((Individual *)b)->fitness) {
        return 1;
    } else {
        return 0;
    }
}

// 差分进化算法函数
void DE(float *ui, float *popold, Individual *pop, float *bm, int st, float F, float CR, int n, int NP, int j, float *xl, float *xu)
{
    int r1, r2, r3;
    float pm1[n], pm2[n], pm3[n];
    int mui[n], mpo[n];
    int sum_mui = 0, nn, i;
    unsigned int seed = 0;

    do
    {
        r1 = rand_int(NP);
    } while (r1 == j);
    do
    {
        r2 = rand_int(NP);
    } while (r2 == j || r2 == r1);
    do
    {
        r3 = rand_int(NP);
    } while (r3 == j || r3 == r1 || r3 == r2);

    for (i = 0; i < n; i++)
    {
        pm1[i] = pop[r1].position[i];
        pm2[i] = pop[r2].position[i];
        pm3[i] = pop[r3].position[i];
    }

    for (i = 0; i < n; i++)
    {
        mui[i] = rand_float() < CR ? 1 : 0;
    }

    for (i = 0; i < n; i++)
    {
        sum_mui += mui[i];
    }

    if (sum_mui == 0)
    {
        nn = rand_int(n);
        mui[nn] = 1;
    }

    for (i = 0; i < n; i++)
    {
        mpo[i] = mui[i] < 0.5 ? 1 : 0;
    }

    for (i = 0; i < n; i++)
    {
        ui[i] = pm3[i] + F * (pm1[i] - pm2[i]);
        ui[i] = popold[i] * mpo[i] + ui[i] * mui[i];
    }

    for (i = 0; i < n; i++)
    {
        if (ui[i] < xl[i])
            ui[i] = xl[i];
        if (ui[i] > xu[i])
            ui[i] = xu[i];
    }
}

float bundle_dip_residual(float *params, SE_coherence_subcube_t *subcube)
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

    if (subcube->type == 0)
    {
        res = 1 - se_mar_semblance(&subcube_temp);
    }
    else if (subcube->type == 1)
    {
        res = 1 - se_var_semblance0(&subcube_temp);
    }
    else if (subcube->type == 2)
    {
        res = 1 - se_var_semblance1(&subcube_temp);
    }
    else if (subcube->type == 3)
    {
        res = 1 - se_var_semblance2(&subcube_temp);
    }
    else if (subcube->type == 4)
    {
        res = 1 - se_eigenstructure(&subcube_temp);
    }
    else if (subcube->type == 5)
    {
        res = 1 - se_gradiant_structure_tensor(&subcube_temp);
    }
    else
    {
        res = 1 - se_mar_semblance(&subcube_temp);
    }

    free(data_moved);
    free(subcube_temp.data);
    return res;
}



// 差分进化算法主函数
void ncde(float *bounds, int popSize, int Max_gen, SE_coherence_subcube_t *subcube)
{
    Individual *pop;
    float best_individual[2];
    float best_fitness;
    int i, j, k, l;
    int NP = popSize;
    int D = 2; // 修改维度为 2
    float F = 0.9;
    float CR = 0.1;
    int Max_FES = NP * Max_gen;
    int nfeval = 0, iter = 0, min_idx;
    float min_dist, dist, temp_val, *xl, *xu, *tempval, *ui, bm[D], popold[D], newpop1[NP * D];

    xl = (float *)malloc(D * sizeof(float));
    xu = (float *)malloc(D * sizeof(float));
    for (i = 0; i < D; i++)
    {
        xl[i] = bounds[i * 2];
        xu[i] = bounds[i * 2 + 1];
    }

    pop = (Individual *)malloc(NP * sizeof(Individual));
    tempval = (float *)malloc(NP * sizeof(float));
    ui = (float *)malloc(NP * D * sizeof(float));
    memset(ui, 0, NP * D * sizeof(float));

    // 初始化种群
    for (i = 0; i < NP; i++)
    {
        for (j = 0; j < D; j++)
        {
            pop[i].position[j] = xl[j] + (xu[j] - xl[j]) * rand_float();
        }
    }

    // 评估初始种群
    for (i = 0; i < NP; i++)
    {
        pop[i].fitness = bundle_dip_residual(pop[i].position, subcube);
        nfeval++;
    }

    while (nfeval < Max_FES)
    {
        // 根据适应度进行排序
        qsort(pop, NP, sizeof(Individual), compare_sort_double);

        for (j = 0; j < NP; j++)
        {
            for (k = 0; k < D; k++)
            {
                bm[k] = pop[j].position[k];
                popold[k] = pop[j].position[k];
            }

            for (k = 0; k < NP; k++)
            {
                for (l = 0; l < D; l++)
                {
                    newpop1[k * D + l] = pop[k].position[l];
                }
            }

            DE(&ui[j * D], popold, pop, bm, 2, F, CR, D, NP, j, xl, xu);
            nfeval++;
        }

        // 评估新种群
        for (i = 0; i < NP; i++)
        {
            tempval[i] = bundle_dip_residual(&ui[i * D], subcube);
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
                    dist += (ui[i * D + k] - pop[j].position[k]) * (ui[i * D + k] - pop[j].position[k]);
                }
                dist = sqrt(dist);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_idx = j;
                }
            }

            if (pop[min_idx].fitness > tempval[i])
            {
                for (k = 0; k < D; k++)
                {
                    pop[min_idx].position[k] = ui[i * D + k];
                }
                pop[min_idx].fitness = tempval[i];
            }
        }
        iter++;
    }

    // 找到最好的个体
    int best_idx = 0;
    for (i = 1; i < NP; i++)
    {
        if (pop[i].fitness > pop[best_idx].fitness)
        {
            best_idx = i;
        }
    }

    for (j = 0; j < D; j++)
    {
        best_individual[j] = pop[best_idx].position[j];
    }
    best_fitness = pop[best_idx].fitness;

    subcube->semblance_dip = 1 - best_fitness; // 将最佳适应度保存到subcube中
    
    subcube->px = best_individual[0];
    subcube->py = best_individual[1]; //将个体最佳参数保存到subcube中

    free(xl);
    free(xu);
    free(pop);
    free(tempval);
    free(ui);
}

void shift_trace_by_interp(float shift, float *data, int n1, int i1, int w1, float *shift_data)
// w1 窗口
// i1 窗口中心
// shift 位移
// data 输入数据
{
    int hw1 = (w1 - 1) / 2;
    int shift_int = round(shift);
    int shift_dir = (shift - shift_int) >= 0 ? 1 : -1;
    shift_dir = (shift - shift_int) == 0 ? 0 : shift_dir;
    float shift_dec = fabs(shift - shift_int);

    int st_i = 0, j1;
    for (j1 = i1 - hw1 + shift_int; j1 <= i1 + hw1 + shift_int; j1++)
    {
        if (j1 < 0 || j1 >= n1)
        {
            shift_data[st_i] = 0;
        }
        else
        {
            if (j1 + shift_dir < 0 || j1 + shift_dir >= n1)
            {
                shift_data[st_i] = data[j1];
            }
            else
            {
                shift_data[st_i] = (1 - shift_dec) * data[j1] + shift_dec * data[j1 + shift_dir];
            }
        }
        st_i++;
    }
}



void direct_dip(float *bounds, int popSize, int Max_gen, SE_coherence_subcube_t *subcube)
{
 //将bounds的两个边界平均分成popSize份,分别计算对应的值，保留最好的一组
    int i, j, k;
    float step_x = (bounds[1] - bounds[0]) / popSize;
    float step_y = (bounds[3] - bounds[2]) / popSize;
    float best_fitness = INFINITY;
    float best_individual[2];
    float fitness;

    float pop[2];

    for (i = 0; i < popSize; i++)
    {
        for (j = 0; j < popSize; j++)
        {
            pop[0] = bounds[0] + i * step_x;
            pop[1] = bounds[2] + j * step_y;
            fitness = bundle_dip_residual(pop, subcube);
            if (fitness < best_fitness)
            {
                best_fitness = fitness;
                best_individual[0] = pop[0];
                best_individual[1] = pop[1];
            }
        }
    }

    subcube->semblance_dip = 1 - best_fitness;
    subcube->px = best_individual[0];
    subcube->py = best_individual[1];
    
}




