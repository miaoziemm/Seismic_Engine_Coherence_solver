#include "../SEmain/include/SEcore.h"
#include "./help/se_coherence_help.h"

typedef struct SE_interplatation_2d_s
{
    int ix;
    int iy;
    int it;
    float coherence_value;
} SE_interplatation_2d_t;

int main(int argc, char *argv[])
{
    // 将数据插值为一个二维剖面

    initargs(argc, argv);
    if(argc==1)
    {
        printf("%s",se_interplatation_2d_help);
        return 0;
    }
    int i, j, k, is;
    int n1, n2, nt;
    int axis1, axis2, axis;
    int cut_point;
    int total_points;
    char *layer_fname;
    char *layer_interplatation_fname;
    float **data;
    SE_interplatation_2d_t interplatation_data;

    if (!getparint("n1", &n1))
        n1 = 100;
    if (!getparint("n2", &n2))
        n2 = 100;
    if (!getparint("axis1", &axis1))
        axis1 = 1;
    if (!getparint("axis2", &axis2))
        axis2 = 2;
    if (!getparstring("layer_fname", &layer_fname))
        layer_fname = "layer.txt";
    if (!getparstring("layer_interplatation_fname", &layer_interplatation_fname))
        layer_interplatation_fname = "layer_interplatation.bin";

    data = alloc2float(n2, n1);
    memset(*data, 0, n1 * n2 * sizeof(float));
    FILE *fp = fopen(layer_fname, "r");
    if (fp == NULL)
    {
        SE_ERROR("Can't open file %s\n", layer_fname);
        return -1;
    }
    total_points = n1*n2;
    axis = axis1 * 10 + axis2;

    for (is = 0; is < total_points; is++)
    {
        fscanf(fp, "%d,%d,%d,%f\n", &i, &j, &k, &interplatation_data.coherence_value);
        switch (axis)
        {
        case 12:
            data[i][j] = interplatation_data.coherence_value;
            break;
        case 21:
            data[j][i] = interplatation_data.coherence_value;
            break;
        case 13:
            data[i][k] = interplatation_data.coherence_value;
            break;
        case 31:
            data[k][i] = interplatation_data.coherence_value;
            break;
        case 23:
            data[j][k] = interplatation_data.coherence_value;
            break;
        case 32:
            data[k][j] = interplatation_data.coherence_value;
            break;
        default:
            SE_ERROR("axis1 and axis2 must be 1,2,3\n");
        }
    }

    fclose(fp);
    FILE *fp2 = fopen(layer_interplatation_fname, "wb");
    fwrite(data[0], sizeof(float), n1 * n2, fp2);
    fclose(fp2);

    SE_MESSAGE("Interplatation done\n");
    SE_MESSAGE("Interplatation data saved in %s\n", layer_interplatation_fname);
    
    free2float(data);
    return 0;
}
