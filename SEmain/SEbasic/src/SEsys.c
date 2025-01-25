#include "../include/SEsys.h"

size_t remove_nans(float* a, size_t n, int copy_last, float replacement)
{
    size_t i, cnt;
    for(cnt = i = 0; i < n; ++i) {
        float v = a[i];
        if( ! isfinite(v) || v > 1.0E+20 || v < -1.0E20 ) {
            a[i] = replacement;
            ++cnt;
        } else if( copy_last ) {
            replacement = v;
        }
    }

    return cnt;
}

char *link_chars(const char *a, const char *b) {
    size_t len_a = strlen(a);
    size_t len_b = strlen(b);
    char *result = malloc(len_a + len_b + 1);
    if (result == NULL) {
        perror("malloc");
        return NULL;
    }

    memcpy(result, a, len_a);
    memcpy(result + len_a, b, len_b);
    result[len_a + len_b] = '\0';

    return result;
}


size_t get_total_memory() {
#if defined(_WIN32) || defined(_WIN64)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return status.ullTotalPhys;
    } else {
        fprintf(stderr, "GlobalMemoryStatusEx failed\n");
        return 0;
    }
#elif defined(__APPLE__) || defined(__MACH__)
    int mib[2];
    size_t mem;
    size_t len = sizeof(mem);

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;

    if (sysctl(mib, 2, &mem, &len, NULL, 0) == -1) {
        perror("sysctl");
        return 0;
    }

    return mem;
#elif defined(__linux__)
    FILE *file = fopen("/proc/meminfo", "r");
    if (file == NULL) {
        perror("fopen");
        return 0;
    }

    char line[256];
    size_t mem_available = 0;

    while (fgets(line, sizeof(line), file)) {
        // 解析 /proc/meminfo 文件中的内存信息
        if (sscanf(line, "MemTotal: %zu kB", &mem_available) == 1) {
            mem_available *= 1024; // 将 kB 转换为字节
            break;
        }
    }

    fclose(file);
    return mem_available;
#else
    fprintf(stderr, "Unsupported platform\n");
    return 0;
#endif
}