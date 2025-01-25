#include "SEbasic.h"


void SE_VERSION();
void SE_MESSAGE(const char *fmt, ...);
void SE_ERROR(const char *fmt, ...);
void SE_WARNING(const char *fmt, ...);

//用于并行报告
void SE_REPORT_P(int rank,int thread,const char *fmt, ...);
