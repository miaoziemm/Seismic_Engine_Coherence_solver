#include "../SEmain/include/SEcore.h"
#include "./help/se_coherence_help.h"

int main(int argc, char *argv[])
{
    initargs(argc, argv);
    if(argc==1)
    {
        printf("%s",se_layer_generator_help);
        return 0;
    }
   
    int i,j,k;
    int nx,ny,nt;
    int direction;
    int cut_point;
    int total_points;
    char *layer_fname;

    if(!getparint("nx",&nx)) nx = 100;
    if(!getparint("ny",&ny)) ny = 100;
    if(!getparint("nt",&nt)) nt = 100;
    if(!getparint("direction",&direction)) direction = 2;
    if(!getparint("cut_point",&cut_point)) cut_point = 50;
    if(!getparstring("layer_fname",&layer_fname)) layer_fname = "layer.txt";


    FILE *fp = fopen(layer_fname,"w");

    if(direction==0) //沿x方向切割
    {
        if(cut_point>=nx||cut_point<0)
        {
            SE_ERROR("Error: cut_point must be less than nx %d and greater than 0\n",nx);
            return -1;
        }

        total_points = ny * nt;
        fprintf(fp,"%d\n",total_points);
        for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                for(k=0;k<nt;k++)
                {
                    if(i==cut_point)
                    {
                        fprintf(fp,"%d,%d,%d\n",i,j,k);
                    }
                }
            }
        }
    }
    else if(direction==1) //沿y方向切割
    {
        if(cut_point>=ny||cut_point<0)
        {
            SE_ERROR("Error: cut_point must be less than ny %d and greater than 0\n",ny);
            return -1;
        }
        total_points = nx * nt;
        fprintf(fp,"%d\n",total_points);
        for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                for(k=0;k<nt;k++)
                {
                    if(j==cut_point)
                    {
                        fprintf(fp,"%d,%d,%d\n",i,j,k);
                    }
                }
            }
        }
    }
    else if(direction==2) //沿z方向切割
    {
        if(cut_point>=nt||cut_point<0)
        {
            SE_ERROR("Error: cut_point must be less than nt %d and greater than 0\n",nt);
            return -1;
        }
        total_points = nx * ny;
        fprintf(fp,"%d\n",total_points);
        for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                for(k=0;k<nt;k++)
                {
                    if(k==cut_point)
                    {
                        fprintf(fp,"%d,%d,%d\n",i,j,k);
                    }
                }
            }
        }
    }
    else
    {
        SE_ERROR("Error: direction must be 0, 1 or 2\n");
        return -1;
    }

    SE_MESSAGE("Create layer file done\n");
    SE_MESSAGE("layer file name: %s\n",layer_fname);
    
    fclose(fp);
    return 0;
}