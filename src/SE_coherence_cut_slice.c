/*
Auther: Chang Zhimiao
Date: 2024/12/1
Version: 1.0
*/

#include "SE_coherence_cut.h"
#include "./help/se_coherence_help.h"

void se_coherence_cut_cfg_init(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_slice_t *slice)
{
    int nx = cfg->nx;
    int ny = cfg->ny;
    int nz = cfg->nz;
    int direction = cfg->direction;
    int cut_point = cfg->cut_point;

    slice->direction = direction;
    slice->cut_point = cut_point;

    cfg->coherence_slice_name = (char *)malloc(512 * sizeof(char));

    switch (direction)
    {
    case 0:
        SE_MESSAGE("Cutting coherence cube along x direction\n");
        slice->n1 = nz;
        slice->n2 = ny;
        slice->n_cut_direction = nx;
        sprintf(cfg->coherence_slice_name, "%s_%d_%d_%d_x_%d", cfg->coherence_cube_name, nx, ny, nz, cut_point);
        SE_MESSAGE("Slice will saved in %s\n", cfg->coherence_slice_name);
        break;
    case 1:
        SE_MESSAGE("Cutting coherence cube along y direction\n");
        slice->n1 = nz;
        slice->n2 = nx;
        slice->n_cut_direction = ny;
        sprintf(cfg->coherence_slice_name, "%s_%d_%d_%d_y_%d", cfg->coherence_cube_name, nx, ny, nz, cut_point);
        SE_MESSAGE("Slice will saved in %s\n", cfg->coherence_slice_name);
        break;
    case 2:
        SE_MESSAGE("Cutting coherence cube along z direction\n");
        slice->n1 = ny;
        slice->n2 = nx;
        slice->n_cut_direction = nz;
        sprintf(cfg->coherence_slice_name, "%s_%d_%d_%d_z_%d", cfg->coherence_cube_name, nx, ny, nz, cut_point);
        SE_MESSAGE("Slice will saved in %s\n", cfg->coherence_slice_name);
        break;
    default:
        SE_ERROR("Invalid direction\n");
        SE_WARNING("Automatic set direction to z\n");
        slice->n1 = ny;
        slice->n2 = nx;
        slice->n_cut_direction = nz;
        direction = 2;
        slice->direction = direction;
        sprintf(cfg->coherence_slice_name, "%s_%d_%d_%d_z_%d", cfg->coherence_cube_name, nx, ny, nz, cut_point);
        SE_MESSAGE("Slice will saved in %s\n", cfg->coherence_slice_name);
        break;
    }

    if (cut_point > slice->n_cut_direction || cut_point < 0)
    {
        SE_ERROR("Invalid cut point\n");
        SE_ERROR("Cutpoint should be in range [0, %d]\n", slice->n_cut_direction);
        exit(1);
    }

    slice->data = (float *)malloc(sizeof(float) * slice->n1 * slice->n2);
}

void se_coherence_cut(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_slice_t *slice)
{
    int i, j, k;
    int direction = slice->direction;
    float coherence_value;
    float ***coherecee_cube;
    coherecee_cube = alloc3float(cfg->nz, cfg->nx, cfg->ny);
    FILE *coherence_cube_fp;
    FILE *coherence_slice_fp;

    coherence_cube_fp = fopen(cfg->coherence_cube_name, "rb");
    if (coherence_cube_fp == NULL)
    {
        SE_ERROR("Cannot open coherence cube file %s\n", cfg->coherence_cube_name);
        exit(1);
    }
    coherence_slice_fp = fopen(cfg->coherence_slice_name, "wb");
    for (i = 0; i < cfg->ny; i++)
    {
        for (j = 0; j < cfg->nx; j++)
        {
            for (k = 0; k < cfg->nz; k++)
            {
                fread(&coherecee_cube[i][j][k], sizeof(float), 1, coherence_cube_fp);
            }
        }
    }

    switch (direction)
    {
    case 0:
        for (i = 0; i < slice->n2; i++)
        {
            for (j = 0; j < slice->n1; j++)
            {
                coherence_value = coherecee_cube[i][slice->cut_point][j];
                slice->data[i * slice->n1 + j] = coherence_value;
            }
        }
        break;
    case 1:
        for (i = 0; i < slice->n2; i++)
        {
            for (j = 0; j < slice->n1; j++)
            {
                coherence_value = coherecee_cube[slice->cut_point][i][j];
                slice->data[i * slice->n1 + j] = coherence_value;
            }
        }
        break;
    case 2:
        for (i = 0; i < slice->n2; i++)
        {
            for (j = 0; j < slice->n1; j++)
            {
                coherence_value = coherecee_cube[j][i][slice->cut_point];
                slice->data[i * slice->n1 + j] = coherence_value;
            }
        }
        break;
    default:
        for (i = 0; i < slice->n2; i++)
        {
            for (j = 0; j < slice->n1; j++)
            {
                coherence_value = coherecee_cube[j][i][slice->cut_point];
                slice->data[i * slice->n1 + j] = coherence_value;
            }
        }
        break;
    }

    // 判断是否有nan
    for(i=0;i<slice->n1*slice->n2;i++){
        if(isnan(slice->data[i])){
            slice->data[i] = 0;
        }
    }

    fwrite(slice->data, sizeof(float), slice->n1 * slice->n2, coherence_slice_fp);

    SE_MESSAGE("Slice done!\n");
    SE_MESSAGE("Slice saved in %s\n", cfg->coherence_slice_name);
    fclose(coherence_cube_fp);
    fclose(coherence_slice_fp);

    free3float(coherecee_cube);
}

void se_coherence_cut_cfg_free(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_slice_t *slice)
{
    free(slice->data);
}

int main(int argc, char *argv[])
{
    if(argc==1)
    {
        printf("%s",se_coherence_cut_slice_help);
        return 0;
    }
    initargs(argc, argv);
    SE_coherence_cut_cfg_t cut_cfg;
    SE_coherence_cut_slice_t cut_slice;

    if (!getparstring("coherence_cube_file", &cut_cfg.coherence_cube_name))
        cut_cfg.coherence_cube_name = "coherence_cube.bin";
    if (!getparint("nx", &cut_cfg.nx))
        cut_cfg.nx = 100;
    if (!getparint("ny", &cut_cfg.ny))
        cut_cfg.ny = 100;
    if (!getparint("nz", &cut_cfg.nz))
        cut_cfg.nz = 100;
    if (!getparint("direction", &cut_cfg.direction))
        cut_cfg.direction = 2;
    if (!getparint("cutpoint", &cut_cfg.cut_point))
        cut_cfg.cut_point = 50;

    se_coherence_cut_cfg_init(&cut_cfg, &cut_slice);
    se_coherence_cut(&cut_cfg, &cut_slice);
    se_coherence_cut_cfg_free(&cut_cfg, &cut_slice);

    return 0;
}