#include "SE_coherence.h"
#include "SE_coherence_math.h"

#define CHECK_TIME 1000

void se_coherence_cube_create(SE_coherence_cube_t *coherence_cube, int nx, int ny, int nt, int dip_window, int window_x, int window_y, int window_t)
{
    int i, j, k;
    int is;
    int ii, jj, kk;
    int trace_center_x, trace_center_y, trace_center_t;
    int trace_save_x, trace_save_y, trace_save_t;
    int center_id, id_temp;

    coherence_cube->nx = nx;
    coherence_cube->ny = ny;
    coherence_cube->nt = nt;

    coherence_cube->max_memory_size = get_total_memory() - 104857600; // 预留100Mb

    coherence_cube->coherence_distribution = alloc1float(nt * nx * ny);
    memset(coherence_cube->coherence_distribution, 0, nt * nx * ny * sizeof(float));

    if (coherence_cube->flag_along_layer == 0)
    {
        coherence_cube->total_subcubes = nx * ny * nt; // 理论上每个点都要处理
        coherence_cube->subcubes = (SE_coherence_subcube_t *)malloc(coherence_cube->total_subcubes * sizeof(SE_coherence_subcube_t));
        // printf("total_subcubes = %d\n", coherence_cube->total_subcubes);
// #pragma omp parallel for private(i, j, k, ii, jj, trace_center_x, trace_center_y, trace_center_t, trace_save_x, trace_save_y, id_temp, center_id) schedule(dynamic)
        for (i = 0; i < coherence_cube->ny; i++)
        {
            for (j = 0; j < coherence_cube->nx; j++)
            {
                for (k = 0; k < coherence_cube->nt; k++)
                {
                    // printf("i = %d, j = %d, k = %d\n", i, j, k);
                    center_id = i * coherence_cube->nx * coherence_cube->nt + j * coherence_cube->nt + k;
                    coherence_cube->subcubes[center_id].ix = j;
                    coherence_cube->subcubes[center_id].iy = i;
                    coherence_cube->subcubes[center_id].it = k;
                    coherence_cube->subcubes[center_id].window_x = window_x;
                    coherence_cube->subcubes[center_id].window_y = window_y;
                    coherence_cube->subcubes[center_id].window_t = window_t;
                    coherence_cube->subcubes[center_id].window_dip = dip_window;
                    coherence_cube->subcubes[center_id].type = coherence_cube->type;
                    coherence_cube->subcubes[center_id].dx = coherence_cube->dx;
                    coherence_cube->subcubes[center_id].dy = coherence_cube->dy;
                    coherence_cube->subcubes[center_id].dt = coherence_cube->dt;
                    // coherence_cube->subcubes[center_id].data = (float **)malloc(window_x * window_y * sizeof(float *));
                    coherence_cube->subcubes[center_id].idlist = (int *)malloc(window_x * window_y * sizeof(int));
                    coherence_cube->subcubes[center_id].idx = (float *)malloc(window_x * window_y * sizeof(float));
                    coherence_cube->subcubes[center_id].idy = (float *)malloc(window_x * window_y * sizeof(float));

                    trace_center_t = k - (window_t - 1) / 2;
                    trace_center_x = j;
                    trace_center_y = i;

                    for (ii = 0; ii < window_x; ii++)
                    {
                        for (jj = 0; jj < window_y; jj++)
                        {
                            trace_save_x = trace_center_x + ii - (window_x - 1) / 2;
                            trace_save_y = trace_center_y + jj - (window_y - 1) / 2;

                            if (trace_save_x < 0 || trace_save_x >= coherence_cube->nx || trace_save_y < 0 || trace_save_y >= coherence_cube->ny || trace_center_t < 0 || trace_center_t >= coherence_cube->nt)
                            {
                                id_temp = -1;
                            }
                            else
                            {
                                id_temp = trace_save_y * coherence_cube->nx * coherence_cube->nt + trace_save_x * coherence_cube->nt + trace_center_t;
                            }

                            coherence_cube->subcubes[center_id].idlist[ii + jj * window_x] = id_temp;
                            coherence_cube->subcubes[center_id].idx[ii + jj * window_x] = ii - (window_x - 1) / 2;
                            coherence_cube->subcubes[center_id].idy[ii + jj * window_x] = jj - (window_y - 1) / 2;
                        }
                    }
                }
            }
        }
    }
    else if (coherence_cube->flag_along_layer == 1)
    {
        int layer_count = 0;
        FILE *layer_fp;

        layer_fp = fopen(coherence_cube->layer_fname, "r");
        if (layer_fp == NULL)
        {
            SE_ERROR("Cannot open layer file %s\n", coherence_cube->layer_fname);
            exit(1);
        }

        fscanf(layer_fp, "%d\n", &layer_count);
        coherence_cube->total_subcubes = layer_count;

        coherence_cube->subcubes = (SE_coherence_subcube_t *)malloc(coherence_cube->total_subcubes * sizeof(SE_coherence_subcube_t));

        for (is = 0; is < coherence_cube->total_subcubes; is++)
        {
            fscanf(layer_fp, "%d,%d,%d\n", &j, &i, &k);
            // center_id = i * coherence_cube->nx * coherence_cube->nt + j * coherence_cube->nt + k;
            center_id = is;
            coherence_cube->subcubes[center_id].ix = j;
            coherence_cube->subcubes[center_id].iy = i;
            coherence_cube->subcubes[center_id].it = k;
            coherence_cube->subcubes[center_id].window_x = window_x;
            coherence_cube->subcubes[center_id].window_y = window_y;
            coherence_cube->subcubes[center_id].window_t = window_t;
            coherence_cube->subcubes[center_id].window_dip = dip_window;
            coherence_cube->subcubes[center_id].type = coherence_cube->type;
            coherence_cube->subcubes[center_id].dx = coherence_cube->dx;
            coherence_cube->subcubes[center_id].dy = coherence_cube->dy;
            coherence_cube->subcubes[center_id].dt = coherence_cube->dt;
            // coherence_cube->subcubes[center_id].data = (float **)malloc(window_x * window_y * sizeof(float *));
            coherence_cube->subcubes[center_id].idlist = (int *)malloc(window_x * window_y * sizeof(int));
            coherence_cube->subcubes[center_id].idx = (float *)malloc(window_x * window_y * sizeof(float));
            coherence_cube->subcubes[center_id].idy = (float *)malloc(window_x * window_y * sizeof(float));

            trace_center_t = k - (window_t - 1) / 2;
            trace_center_x = j;
            trace_center_y = i;

            for (ii = 0; ii < window_x; ii++)
            {
                for (jj = 0; jj < window_y; jj++)
                {
                    trace_save_x = trace_center_x + ii - (window_x - 1) / 2;
                    trace_save_y = trace_center_y + jj - (window_y - 1) / 2;

                    if (trace_save_x < 0 || trace_save_x >= coherence_cube->nx || trace_save_y < 0 || trace_save_y >= coherence_cube->ny || trace_center_t < 0 || trace_center_t >= coherence_cube->nt)
                    {
                        id_temp = -1;
                    }
                    else
                    {
                        id_temp = trace_save_y * coherence_cube->nx * coherence_cube->nt + trace_save_x * coherence_cube->nt + trace_center_t;
                    }

                    coherence_cube->subcubes[center_id].idlist[ii + jj * window_x] = id_temp;
                    coherence_cube->subcubes[center_id].idx[ii + jj * window_x] = ii - (window_x - 1) / 2;
                    coherence_cube->subcubes[center_id].idy[ii + jj * window_x] = jj - (window_y - 1) / 2;
                }
            }
        }
        fclose(layer_fp);
    }
    else
    {
        SE_WARNING("flag_along_layer is not defined\n");
        exit(1);
    }

    // se_coherence_cube_report(coherence_cube);
}

void se_coherence_cube_destroy(SE_coherence_cube_t *coherence_cube)
{
    int i;
    free1float(coherence_cube->coherence_distribution);

    for (i = 0; i < coherence_cube->total_subcubes; i++)
    {
        free(coherence_cube->subcubes[i].idlist);
        free(coherence_cube->subcubes[i].idx);
        free(coherence_cube->subcubes[i].idy);
    }
    free(coherence_cube->subcubes);
    // free(coherence_cube);
}

// TODO : add version report function
// void se_coherence_cube_version_report(SE_coherence_cube_t *coherence_cube)
// {

// }
// dip_type 是否进行倾角矫正
void se_coherence_cube_compute(int dip_type, SE_coherence_cube_t *coherence_cube)
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
    zero_temp = alloc1float(coherence_cube->nt); //
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

    float *px = NULL, *py = NULL;

        px = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        py = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        if (dip_type == 2)
        {
            FILE *dip_temp;
            dip_temp = fopen(coherence_cube->dip_result_fname, "r");
            for (i = 0; i < coherence_cube->nx * coherence_cube->ny * coherence_cube->nt; i++)
            {
                fscanf(dip_temp, "%f,%f\n", &px[i], &py[i]);
            }
            fclose(dip_temp);
        }
    
    // 判断是否为GST,如果是GST则需要计算梯度
    if (coherence_cube->type == 5)
    {
        SE_MESSAGE("GST need gradient\n");
        SE_MESSAGE("Calculating gradient\n");
        float *buff_gradient;
        buff_gradient = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        gaussian_filter_3d(buff, buff_gradient, coherence_cube->nx, coherence_cube->ny, coherence_cube->nt, 1.0, 1);
        memcpy(buff, buff_gradient, coherence_cube->nx * coherence_cube->ny * coherence_cube->nt * sizeof(float));
        free1float(buff_gradient);
        SE_MESSAGE("Gradient has been calculated\n");
    }

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
                // coherence_cube->subcubes[id].data[ii] = (float *)malloc(coherence_cube->nt * sizeof(float));
                // memcpy(coherence_cube->subcubes[id].data[ii], &buff[coherence_cube->subcubes[id].idlist[ii]], coherence_cube->nt * sizeof(float));
                coherence_cube->subcubes[id].data[ii] = &buff[coherence_cube->subcubes[id].idlist[ii]];
            }
            else
            {
                coherence_cube->subcubes[id].data[ii] = zero_temp;
                // coherence_cube->subcubes[id].data[ii] = (float *)malloc(coherence_cube->nt * sizeof(float));
                // memcpy(coherence_cube->subcubes[id].data[ii], zero_temp, coherence_cube->nt * sizeof(float));
            }
        }

        if (dip_type == 0) // 无需倾角矫正
        {
            coherence_cube->coherence_distribution[id] = se_coherence_no_dip(coherence_cube->type, &coherence_cube->subcubes[id]);
        }
        else if (dip_type == 1) // 需要倾角矫正
        {
            // ncde(bounds, coherence_cube->popSize, coherence_cube->Max_gen, &coherence_cube->subcubes[id]);
            direct_dip(bounds, coherence_cube->popSize, coherence_cube->Max_gen, &coherence_cube->subcubes[id]);
            coherence_cube->coherence_distribution[id] = coherence_cube->subcubes[id].semblance_dip;
           
                px[id] = coherence_cube->subcubes[id].px;
                py[id] = coherence_cube->subcubes[id].py;
          
        }else if(dip_type==2 && coherence_cube->type !=0)
        {
            params[0] = px[id];
            params[1] = py[id];
            coherence_cube->coherence_distribution[id] = 1-bundle_dip_residual(params, &coherence_cube->subcubes[id]);
        }
        else
        {
            SE_WARNING("Dip type is not defined\n");
            exit(1);
        }

        // for(ii = 0; ii < window_x * window_y; ii++)
        //     free(coherence_cube->subcubes[id].data[ii]);
        free(coherence_cube->subcubes[id].data);

        #pragma omp atomic
        complete_jobs_count++;

        if (complete_jobs_count % CHECK_TIME == 0)
        {
            SE_MESSAGE("Complete %d/%d jobs\n", complete_jobs_count, coherence_cube->total_subcubes);
        }
    }
#pragma omp barrier
SE_MESSAGE("Complete %d/%d jobs\n", coherence_cube->total_subcubes, coherence_cube->total_subcubes);


    if (coherence_cube->flag_along_layer == 0)
    {
        FILE *coherence_cube_fp;
        coherence_cube_fp = fopen(coherence_cube->coherence_cube_name, "wb");
        // write to file
        fwrite(coherence_cube->coherence_distribution, sizeof(float), coherence_cube->nx * coherence_cube->ny * coherence_cube->nt, coherence_cube_fp);
        fclose(coherence_cube_fp);

        SE_MESSAGE("Coherence cube has been written to %s\n", coherence_cube->coherence_cube_name);
    }
    else if (coherence_cube->flag_along_layer == 1)
    {
        FILE *layer_result_fp;
        layer_result_fp = fopen(coherence_cube->layer_result_fname, "w");
        for (i = 0; i < coherence_cube->total_subcubes; i++)
        {
            fprintf(layer_result_fp, "%d,%d,%d,%f\n", coherence_cube->subcubes[i].ix, coherence_cube->subcubes[i].iy, coherence_cube->subcubes[i].it, coherence_cube->coherence_distribution[i]);
        }
        fclose(layer_result_fp);
        SE_MESSAGE("Coherence layer has been written to %s\n", coherence_cube->layer_result_fname);
    }
    else
    {
        SE_WARNING("flag_along_layer is not defined\n");
        exit(1);
    }

    if(dip_type==1)
    {
        FILE *dip_fp;
        dip_fp = fopen(coherence_cube->dip_result_fname, "w");
        for (i = 0; i < coherence_cube->total_subcubes; i++)
        {
            fprintf(dip_fp, "%f,%f\n", px[i], py[i]);
        }
        fclose(dip_fp);
    }

    free(buff);

        free(px);
        free(py);
    
    free1float(zero_temp);
    fclose(data_cube_fp);
}

void se_coherence_cube_pre(int dip_type, SE_coherence_cube_t *coherence_cube)
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
    zero_temp = alloc1float(coherence_cube->nt); //
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
    if (coherence_cube->type == 0)
    {
        px = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        py = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
    }
    else if (coherence_cube->type != 0 && dip_type == 2)
    {
        px = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        py = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        FILE *dip_temp;
        dip_temp = fopen("dip.temp", "r");
        for (i = 0; i < coherence_cube->nx * coherence_cube->ny * coherence_cube->nt; i++)
        {
            fscanf(dip_temp, "%f,%f\n", &px[i], &py[i]);
        }
        fclose(dip_temp);
    }
    // 判断是否为GST,如果是GST则需要计算梯度
    if (coherence_cube->type == 5)
    {
        SE_MESSAGE("GST need gradient\n");
        SE_MESSAGE("Calculating gradient\n");
        float *buff_gradient;
        buff_gradient = alloc1float(coherence_cube->nx * coherence_cube->ny * coherence_cube->nt);
        gaussian_filter_3d(buff, buff_gradient, coherence_cube->nx, coherence_cube->ny, coherence_cube->nt, 1.0, 1);
        memcpy(buff, buff_gradient, coherence_cube->nx * coherence_cube->ny * coherence_cube->nt * sizeof(float));
        free1float(buff_gradient);
        SE_MESSAGE("Gradient has been calculated\n");
    }

    int complete_jobs_count = 0;

#pragma omp parallel for private(id) schedule(dynamic)
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
                coherence_cube->subcubes[id].data[ii] = (float *)malloc(coherence_cube->nt * sizeof(float));
                memcpy(coherence_cube->subcubes[id].data[ii], &buff[coherence_cube->subcubes[id].idlist[ii]], coherence_cube->nt * sizeof(float));
                // coherence_cube->subcubes[id].data[ii] = &buff[coherence_cube->subcubes[id].idlist[ii]];
            }
            else
            {
                // coherence_cube->subcubes[id].data[ii] = &zero_temp[0];
                coherence_cube->subcubes[id].data[ii] = (float *)malloc(coherence_cube->nt * sizeof(float));
                memcpy(coherence_cube->subcubes[id].data[ii], zero_temp, coherence_cube->nt * sizeof(float));
            }
        }

        if (dip_type == 0) // 无需倾角矫正
        {
            coherence_cube->coherence_distribution[id] = se_coherence_no_dip(coherence_cube->type, &coherence_cube->subcubes[id]);
        }
        else if (dip_type == 1) // 需要倾角矫正
        {
            ncde(bounds, coherence_cube->popSize, coherence_cube->Max_gen, &coherence_cube->subcubes[id]);
            coherence_cube->coherence_distribution[id] = coherence_cube->subcubes[id].semblance_dip;
            if (coherence_cube->type == 0)
            {
                px[id] = coherence_cube->subcubes[id].px;
                py[id] = coherence_cube->subcubes[id].py;
            } // 只有semblance需要保存倾角
        }else if(dip_type==2 && coherence_cube->type !=0)
        {
            params[0] = px[id];
            params[1] = py[id];
            coherence_cube->coherence_distribution[id] = 1-bundle_dip_residual(params, &coherence_cube->subcubes[id]);
        }
        else
        {
            SE_WARNING("Dip type is not defined\n");
            exit(1);
        }

        for(ii = 0; ii < window_x * window_y; ii++)
            free(coherence_cube->subcubes[id].data[ii]);
        free(coherence_cube->subcubes[id].data);

        // #pragma omp atomic
        // complete_jobs_count++;

        // if (complete_jobs_count % CHECK_TIME == 0)
        // {
        //     SE_MESSAGE("Complete %d/%d jobs\n", complete_jobs_count, coherence_cube->total_subcubes);
        // }
    }
#pragma omp barrier
SE_MESSAGE("Complete %d/%d jobs\n", coherence_cube->total_subcubes, coherence_cube->total_subcubes);


    if (coherence_cube->flag_along_layer == 0)
    {
        FILE *coherence_cube_fp;
        coherence_cube_fp = fopen(coherence_cube->coherence_cube_name, "wb");
        // write to file
        fwrite(coherence_cube->coherence_distribution, sizeof(float), coherence_cube->nx * coherence_cube->ny * coherence_cube->nt, coherence_cube_fp);
        fclose(coherence_cube_fp);

        SE_MESSAGE("Coherence cube has been written to %s\n", coherence_cube->coherence_cube_name);
    }
    else if (coherence_cube->flag_along_layer == 1)
    {
        FILE *layer_result_fp;
        layer_result_fp = fopen(coherence_cube->layer_result_fname, "w");
        for (i = 0; i < coherence_cube->total_subcubes; i++)
        {
            fprintf(layer_result_fp, "%d,%d,%d,%f\n", coherence_cube->subcubes[i].ix, coherence_cube->subcubes[i].iy, coherence_cube->subcubes[i].it, coherence_cube->coherence_distribution[i]);
        }
        fclose(layer_result_fp);
        SE_MESSAGE("Coherence layer has been written to %s\n", coherence_cube->layer_result_fname);
    }
    else
    {
        SE_WARNING("flag_along_layer is not defined\n");
        exit(1);
    }

    if(coherence_cube->type==0&&dip_type==1)
    {
        FILE *dip_fp;
        dip_fp = fopen("dip.temp", "w");
        for (i = 0; i < coherence_cube->total_subcubes; i++)
        {
            fprintf(dip_fp, "%f,%f\n", px[i], py[i]);
        }
        fclose(dip_fp);
    }

    free(buff);
    if (coherence_cube->type == 0)
    {
        free(px);
        free(py);
    }
    free1float(zero_temp);
    fclose(data_cube_fp);
}



void se_coherence_cube_report(SE_coherence_cube_t *coherence_cube)
{
    SE_MESSAGE("Coherence cube created\n");
    SE_MESSAGE("nx = %d, ny = %d, nt = %d\n", coherence_cube->nx, coherence_cube->ny, coherence_cube->nt);
    SE_MESSAGE("window_x = %d, window_y = %d, window_t = %d\n", coherence_cube->subcubes[0].window_x, coherence_cube->subcubes[0].window_y, coherence_cube->subcubes[0].window_t);
    SE_MESSAGE("Available memory size = %ld Mb\n", coherence_cube->max_memory_size / 1048576);
    if (coherence_cube->type == 0)
    {
        SE_MESSAGE("Coherence cube type = Semblance\n");
    }
    else if (coherence_cube->type == 1)
    {
        SE_MESSAGE("Coherence cube type = Variation0\n");
    }
    else if (coherence_cube->type == 2)
    {
        SE_MESSAGE("Coherence cube type = Variation1\n");
    }
    else if (coherence_cube->type == 3)
    {
        SE_MESSAGE("Coherence cube type = Variation2\n");
    }
    else if (coherence_cube->type == 4)
    {
        SE_MESSAGE("Coherence cube type = Eigenstructure\n");
    }
    else if (coherence_cube->type == 5)
    {
        SE_MESSAGE("Coherence cube type = Gradiant Structure Tensor\n");
    }
    else
    {
        SE_WARNING("Coherence cube type is not defined, use default Semblance\n");
        coherence_cube->type = 0;
    }
}
