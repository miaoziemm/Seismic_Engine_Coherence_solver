#include "SE_coherence.h"

#define POWER_ITERATION_MAX 100
#define POWER_ITERATION_TOL 1e-6

float se_coherence_no_dip(int type, SE_coherence_subcube_t *subcube)
{
    float res;
    if (type == 0)
    {
        res = se_mar_semblance(subcube);
    }
    else if (type  == 1)
    {
        res  = se_var_semblance0(subcube);
    }
    else if (type  == 2)
    {
        res = se_var_semblance1(subcube);
    }
    else if (type  == 3)
    {
        res  = se_var_semblance2(subcube);
    }
    else if (type  == 4)
    {
        res  = se_eigenstructure(subcube);
    }
    else if (type == 5)
    {
        res  = se_gradiant_structure_tensor(subcube);
    }
    else
    {
        res  = se_mar_semblance(subcube);
    }
    return res; 
}

float se_mar_semblance(SE_coherence_subcube_t *subcube)
{
    int nt = subcube->window_t;
    int nx = subcube->window_x;
    int ny = subcube->window_y;
    float square_of_sums = 0.0;
    float sum_of_squares = 0.0;
    int i, j, k;

    float temp1 = 0;
    float temp2 = 0;
    for (k = 0; k < nt; k++)
    {
        square_of_sums = 0.0;
        sum_of_squares = 0.0;
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                square_of_sums += subcube->data[i * ny + j][k];
                sum_of_squares += subcube->data[i * ny + j][k] * subcube->data[i * ny + j][k];
            }
        }

        square_of_sums = square_of_sums * square_of_sums;

        temp1 += square_of_sums;
        temp2 += sum_of_squares;
    }
    if (temp2 == 0)
    {
        return 0;
    }
    float semblance = temp1 / temp2;
    return semblance / (nx * ny);
}

float se_var_semblance0(SE_coherence_subcube_t *subcube)
{
    float square_of_sums = 0.0;
    float sum_of_squares = 0.0;
    int i, j, k;
    int nt = subcube->window_t;
    int nx = subcube->window_x;
    int ny = subcube->window_y;

    float *mean_trace;
    mean_trace = (float *)malloc(nt * sizeof(float));

    float *a;
    a = alloc1float(nx * ny);

    for (i = 0; i < nt; i++)
    {
        mean_trace[i] = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            mean_trace[i] += subcube->data[j][i] / (nx * ny);
        } // 计算标准道
    }
    float temp1, temp2, temp3, temp4;

    for (i = 0; i < nx * ny; i++)
    {
        a[i] = 0.0;
        temp1 = 0;
        temp2 = 0;
        for (j = 0; j < nt; j++)
        {
            temp1 += mean_trace[j] * mean_trace[j];
            temp2 += subcube->data[i][j] * mean_trace[j];
        }
        a[i] = temp2 / temp1;
    }

    temp3 = 0;
    temp2 = 0;
    temp1 = 0;
    for (i = 0; i < nt; i++)
    {
        square_of_sums = 0.0;
        sum_of_squares = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            square_of_sums += (subcube->data[j][i] - a[j] * mean_trace[i]) * (subcube->data[j][i] - a[j] * mean_trace[i]);
            sum_of_squares += subcube->data[j][i] * subcube->data[j][i];
        }

        temp1 += square_of_sums;
        temp2 += sum_of_squares;
    }

    temp3 = temp3 / (nx * ny);
    float semblance = 1 - temp1 / temp2;
    free1float(a);
    free1float(mean_trace);

    if (temp2 == 0)
    {
        return 0;
    }
    return semblance;
}

float se_var_semblance1(SE_coherence_subcube_t *subcube)
{
    float square_of_sums = 0.0;
    float sum_of_squares = 0.0;
    int i, j, k;
    int nt = subcube->window_t;
    int nx = subcube->window_x;
    int ny = subcube->window_y;

    float *mean_trace;
    mean_trace = (float *)malloc(nt * sizeof(float));

    float *a;
    a = alloc1float(nx * ny);

    for (i = 0; i < nt; i++)
    {
        mean_trace[i] = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            mean_trace[i] += subcube->data[j][i] / (nx * ny);
        } // 计算标准道
    }
    float temp1, temp2, temp3, temp4;

    for (i = 0; i < nx * ny; i++)
    {
        a[i] = 0.0;
        temp1 = 0;
        temp2 = 0;
        for (j = 0; j < nt; j++)
        {
            temp1 += mean_trace[j] * mean_trace[j];
            temp2 += subcube->data[i][j] * mean_trace[j];
        }
        a[i] = temp2 / temp1;
    }

    temp3 = 0;
    temp2 = 0;
    temp1 = 0;
    for (i = 0; i < nt; i++)
    {
        square_of_sums = 0.0;
        sum_of_squares = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            square_of_sums += (subcube->data[j][i] - a[j] * mean_trace[i]) * (subcube->data[j][i] - a[j] * mean_trace[i]);
            sum_of_squares += subcube->data[j][i] * subcube->data[j][i] + a[j] * a[j] * mean_trace[i] * mean_trace[i]; // 分母
        }

        temp1 += square_of_sums;
        temp2 += sum_of_squares;
    }

    temp3 = temp3 / (nx * ny);
    float semblance = 1 - temp1 / temp2;
    free1float(a);
    free1float(mean_trace);

    if (isnan(semblance))
    {
        return 0;
    }
    return semblance;
}

float se_var_semblance2(SE_coherence_subcube_t *subcube)
{
    float square_of_sums = 0.0;
    float sum_of_squares = 0.0;
    int i, j, k;
    int nt = subcube->window_t;
    int nx = subcube->window_x;
    int ny = subcube->window_y;

    float *mean_trace;
    mean_trace = (float *)malloc(nt * sizeof(float));

    float *a;
    a = alloc1float(nx * ny);

    for (i = 0; i < nt; i++)
    {
        mean_trace[i] = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            mean_trace[i] += subcube->data[j][i] / (nx * ny);
        } // 计算标准道
    }
    float temp1, temp2, temp3, temp4;

    for (i = 0; i < nx * ny; i++)
    {
        a[i] = 0.0;
        temp1 = 0;
        temp2 = 0;
        for (j = 0; j < nt; j++)
        {
            temp1 += mean_trace[j] * mean_trace[j];
            temp2 += subcube->data[i][j] * mean_trace[j];
        }
        a[i] = fabs(temp2 / temp1);
    }

    temp3 = 0;
    temp2 = 0;
    temp1 = 0;
    for (i = 0; i < nt; i++)
    {
        square_of_sums = 0.0;
        sum_of_squares = 0.0;
        for (j = 0; j < nx * ny; j++)
        {
            square_of_sums += (subcube->data[j][i] - a[j] * mean_trace[i]) * (subcube->data[j][i] - a[j] * mean_trace[i]);
            sum_of_squares += subcube->data[j][i] * subcube->data[j][i] + a[j] * a[j] * mean_trace[i] * mean_trace[i]; // 分母
        }

        temp1 += square_of_sums;
        temp2 += sum_of_squares;
    }

    temp3 = temp3 / (nx * ny);
    float semblance = 1 - temp1 / temp2;
    free1float(a);
    free1float(mean_trace);

    if (isnan(semblance))
    {
        return 0;
    }
    return semblance;
}

float se_eigenstructure(SE_coherence_subcube_t *subcube)
{
    int n = subcube->window_x * subcube->window_y;
    int m = subcube->window_t; // Assuming the last dimension is the number of traces

    // 计算协方差矩阵
    float *cov = (float *)calloc(m * m, sizeof(float));
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            for (k = 0; k < m; k++)
            {
                cov[j * m + k] += subcube->data[i][j] * subcube->data[i][k];
            }
        }
    }

    // 使用幂法计算最大特征值
    float *b = (float *)malloc(m * sizeof(float));
    float *b_next = (float *)malloc(m * sizeof(float));
    for (i = 0; i < m; i++) b[i] = 1.0; // 初始向量

    float max_eigval = 0.0;
    int iter;
    for (iter = 0; iter < POWER_ITERATION_MAX; iter++)
    {
        // 计算 b_next = cov * b
        for (i = 0; i < m; i++)
        {
            b_next[i] = 0.0;
            for (j = 0; j < m; j++)
            {
                b_next[i] += cov[i * m + j] * b[j];
            }
        }

        // 归一化 b_next
        float norm = 0.0;
        for (i = 0; i < m; i++) norm += b_next[i] * b_next[i];
        norm = sqrt(norm);
        for (i = 0; i < m; i++) b_next[i] /= norm;

        // 检查收敛性
        float diff = 0.0;
        for (i = 0; i < m; i++) diff += fabs(b_next[i] - b[i]);
        if (diff < POWER_ITERATION_TOL) break;

        // 更新 b
        for (i = 0; i < m; i++) b[i] = b_next[i];
    }

    // 计算最大特征值
    for (i = 0; i < m; i++)
    {
        float temp = 0.0;
        for (j = 0; j < m; j++)
        {
            temp += cov[i * m + j] * b[j];
        }
        max_eigval += b[i] * temp;
    }

    // 计算特征值之和
    float sum_eigvals = 0.0;
    for (i = 0; i < m; i++)
    {
        sum_eigvals += cov[i * m + i];
    }

    float result = max_eigval / sum_eigvals;

    free(cov);
    free(b);
    free(b_next);

    if (isnan(result))
    {
        return 0;
    }

    return result;
}

float se_gradiant_structure_tensor(SE_coherence_subcube_t *subcube)
{

    return 0;
}