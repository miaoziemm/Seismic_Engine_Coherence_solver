/*
Auther: Chang Zhimiao
Date: 2024/12/1
Version: 1.0
*/

#include "SE_coherence_cut.h"
#include "./help/se_coherence_help.h"


void se_coherence_cut_cube_cfg_init(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_cube_t *cut_cube)
{
    int nx = cfg->nx;
    int ny = cfg->ny;
    int nz = cfg->nz;

    cfg->coherence_cut_cube_name = (char *)malloc(512 * sizeof(char));

    sprintf(cfg->coherence_cut_cube_name, "%s_%d_%d_%d_%d_%d_%d_o_%d_%d_%d", cfg->coherence_cube_name, nx, ny, nz, cut_cube->n1, cut_cube->n2, cut_cube->n3, cut_cube->no1, cut_cube->no2, cut_cube->no3);

}

void se_coherence_cut_cube_cfg_free(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_cube_t *cut_cube)
{
    free(cut_cube->data);
    free(cfg->coherence_cut_cube_name);
}

void se_coherence_cut_cube(SE_coherence_cut_cfg_t *cfg, SE_coherence_cut_cube_t *cut_cube)
{
    int i, j, k;
    int nx = cfg->nx;
    int ny = cfg->ny;
    int nz = cfg->nz;
    int n1 = cut_cube->n1; //x
    int n2 = cut_cube->n2; //y
    int n3 = cut_cube->n3; //z
    int no1 = cut_cube->no1;
    int no2 = cut_cube->no2;
    int no3 = cut_cube->no3;

    cut_cube->data = (float *)malloc(sizeof(float) * n1 * n2 * n3);

    float ***coherecee_cube;
    coherecee_cube = alloc3float(nz, nx, ny);
    FILE *coherence_cube_fp;
    FILE *coherence_cut_cube_fp;

    coherence_cube_fp = fopen(cfg->coherence_cube_name, "rb");
    if (coherence_cube_fp == NULL)
    {
        SE_ERROR("Cannot open coherence cube file %s\n", cfg->coherence_cube_name);
        exit(1);
    }
    coherence_cut_cube_fp = fopen(cfg->coherence_cut_cube_name, "wb");

    for (i = 0; i < ny; i++)
    {
        for (j = 0; j < nx; j++)
        {
            for (k = 0; k < nz; k++)
            {
                fread(&coherecee_cube[i][j][k], sizeof(float), 1, coherence_cube_fp);
            }
        }
    }

    for (i = 0; i < n2; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n3; k++)
            {
                cut_cube->data[i * n1 * n3 + j * n3 + k] = coherecee_cube[i + no2][j + no1][k + no3];
            }
        }
    }

    fwrite(cut_cube->data, sizeof(float), n1 * n2 * n3, coherence_cut_cube_fp);

    SE_MESSAGE("Cut cube done!\n");
    SE_MESSAGE("Cut cube saved in %s\n", cfg->coherence_cut_cube_name);
    fclose(coherence_cube_fp);
    fclose(coherence_cut_cube_fp);

    free3float(coherecee_cube);
}



int main(int argc, char *argv[])
{
    if(argc==1)
    {
        printf("%s",se_coherence_cut_cube_help);
        return 0;
    }

    initargs(argc, argv);
    SE_coherence_cut_cfg_t cut_cfg;
    SE_coherence_cut_cube_t cut_cube;

    if (!getparstring("coherence_cube_file", &cut_cfg.coherence_cube_name))
        cut_cfg.coherence_cube_name = "coherence_cube.bin";
    if (!getparint("nx", &cut_cfg.nx))
        cut_cfg.nx = 100;
    if (!getparint("ny", &cut_cfg.ny))
        cut_cfg.ny = 100;
    if (!getparint("nz", &cut_cfg.nz))
        cut_cfg.nz = 100;
    
    if(!getparint("n1",&cut_cube.n1))
        cut_cube.n1 = 10;
    if(!getparint("n2",&cut_cube.n2))
        cut_cube.n2 = 10;
    if(!getparint("n3",&cut_cube.n3))
        cut_cube.n3 = 10;
    
    if(!getparint("no1",&cut_cube.no1))
        cut_cube.no1 = 0;
    if(!getparint("no2",&cut_cube.no2))
        cut_cube.no2 = 0;
    if(!getparint("no3",&cut_cube.no3))
        cut_cube.no3 = 0;
        
    se_coherence_cut_cube_cfg_init(&cut_cfg, &cut_cube);
    
    se_coherence_cut_cube(&cut_cfg, &cut_cube);

    se_coherence_cut_cube_cfg_free(&cut_cfg, &cut_cube);

    return 0;
}