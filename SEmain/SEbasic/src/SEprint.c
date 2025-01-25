#include "../include/SEprint.h"
#include "../include/SEversion.h"

void SE_VERSION()
{
    SE_MESSAGE("************************************************************\n");
    SE_MESSAGE("%s\n", VERSION);
    SE_MESSAGE("************************************************************\n");
    printf("\n");
}


void SE_MESSAGE(const char *fmt, ...)
{
    //获取当前时间，年月日时分秒
    time_t now;
    struct tm *timenow;
    time(&now);
    timenow = localtime(&now);
    printf("SE_MESSAGE :: %d-%02d-%02d %02d:%02d:%02d :: ", timenow->tm_year + 1900, timenow->tm_mon + 1, timenow->tm_mday, timenow->tm_hour, timenow->tm_min, timenow->tm_sec);
    //打印信息
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

void SE_ERROR(const char *fmt, ...)
{
    //获取当前时间，年月日时分秒
    time_t now;
    struct tm *timenow;
    time(&now);
    timenow = localtime(&now);
    printf("SE_ERROR :: %d-%02d-%02d %02d:%02d:%02d :: ", timenow->tm_year + 1900, timenow->tm_mon + 1, timenow->tm_mday, timenow->tm_hour, timenow->tm_min, timenow->tm_sec);
    //打印信息
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

void SE_WARNING(const char *fmt, ...)
{
    //获取当前时间，年月日时分秒
    time_t now;
    struct tm *timenow;
    time(&now);
    timenow = localtime(&now);
    printf("SE_WARNING :: %d-%02d-%02d %02d:%02d:%02d :: ", timenow->tm_year + 1900, timenow->tm_mon + 1, timenow->tm_mday, timenow->tm_hour, timenow->tm_min, timenow->tm_sec);
    //打印信息
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

void SE_REPORT_P(int rank,int thread,const char *fmt, ...)
{
    if(rank==0 && thread==0){
        //获取当前时间，年月日时分秒
        time_t now;
        struct tm *timenow;
        time(&now);
        timenow = localtime(&now);
        printf("SE_REPORT :: %d-%02d-%02d %02d:%02d:%02d :: ", timenow->tm_year + 1900, timenow->tm_mon + 1, timenow->tm_mday, timenow->tm_hour, timenow->tm_min, timenow->tm_sec);
        //打印信息
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

