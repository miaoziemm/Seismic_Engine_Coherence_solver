#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#elif defined(__APPLE__) || defined(__MACH__)
#include <sys/types.h>
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

size_t get_total_memory();
char *link_chars(const char *a, const char *b);
size_t remove_nans(float* a, size_t n, int copy_last, float replacement);