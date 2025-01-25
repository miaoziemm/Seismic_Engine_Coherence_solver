/*
Auther: Chang Zhimiao
Date: 2024/12/1
Version: 1.0
*/


#include "SE_coherence.h"
#include "./help/se_coherence_help.h"

int main(int argc, char **argv)
{
    if(argc==1)
    {
        printf("%s",se_coherence_solver_help);
        return 0;
    }
    initargs(argc, argv);

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
    
    if(!getparint("type",&coherence_cube.type)) coherence_cube.type=0;

    if(!getparstring("data_cube_file",&coherence_cube.data_cube_name)) coherence_cube.data_cube_name="inputfile";
    if(!getparstring("coherence_cube_file",&coherence_cube.coherence_cube_name)) coherence_cube.coherence_cube_name="coherencefile";

    if(!getparint("dip_type",&dip_type)) dip_type=0; //0:no_dip 1:ncde 

    if(!getparstring("dip_result_file",&coherence_cube.dip_result_fname)) coherence_cube.dip_result_fname="dip.temp";

    if(!getparint("along_layer",&coherence_cube.flag_along_layer)) coherence_cube.flag_along_layer=0;
    if(coherence_cube.flag_along_layer==1)
    {
        if(!getparstring("layer_file",&coherence_cube.layer_fname)) coherence_cube.layer_fname="layer.txt";
        if(!getparstring("layer_result_file",&coherence_cube.layer_result_fname)) coherence_cube.layer_result_fname="layer_result.txt";
    }

    double start_time = omp_get_wtime();

    SE_MESSAGE("Create coherence cube\n");
    se_coherence_cube_create(&coherence_cube,nx,ny,nt,dip_window,window_x,window_y,window_t);
    SE_MESSAGE("Create coherence cube done\n");
    se_coherence_cube_report(&coherence_cube);

    SE_MESSAGE("Compute coherence cube\n");
    se_coherence_cube_compute(dip_type,&coherence_cube);
    SE_MESSAGE("Compute coherence cube done\n");

    SE_MESSAGE("Destroy coherence cube\n");
    se_coherence_cube_destroy(&coherence_cube);
    SE_MESSAGE("Destroy coherence cube done\n");

    double end_time = omp_get_wtime();
    printf("Total execution time: %f seconds\n", end_time - start_time);

    
    return 0;
}