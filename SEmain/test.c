#include "./include/SEcore.h"

int main(int argc, char* argv[])
{
    initargs(argc, argv);
    int e;
    if(!getparint("e", &e)) e = 0;
    SE_MESSAGE("e = %d\n", e);
    m3d m;
    m.m[0]=1;
    m.m[1]=2;
    m.m[2]=3;
    m.m[3]=4;
    m.m[4]=5;
    m.m[5]=6;
    m.m[6]=4;
    m.m[7]=5;
    m.m[8]=6;

    m4d m4;
    m4d_sym_init(&m4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);

    m4d_info(&m4);


    SE_VERSION();
    SE_MESSAGE("Hello, World!\n");
    SE_ERROR("Hello, World!\n");
    SE_WARNING("Hello, World!\n");
    printf("%d\n", rand());
    int a=1, b=2;
    double d=rand_normal();
    printf("rand_normal() = %f\n", d);
    int c=mini(a, b);
    printf("mini(%d, %d) = %d\n", a, b, c);

    return 0;
}