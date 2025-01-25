#include "SE_coherence_math.h"

// 随机数生成函数
int rand_int(int max) {
    return rand() % max;
}

float rand_float() {
    return (float)rand() / RAND_MAX;
}

void gaussian_filter_1d(const float *input, float *output, int length, float sigma, int order) {
    int i,j;
    int kernel_radius = (int)(3.0 * sigma);
    int kernel_size = 2 * kernel_radius + 1;
    float *kernel = (float *)malloc(kernel_size * sizeof(float));
    float sum = 0.0;

    // Create Gaussian kernel
    for (i = -kernel_radius; i <= kernel_radius; ++i) {
        kernel[i + kernel_radius] = expf(-0.5f * (i * i) / (sigma * sigma));
        sum += kernel[i + kernel_radius];
    }

    // Normalize the kernel
    for (i = 0; i < kernel_size; ++i) {
        kernel[i] /= sum;
    }

    // Apply Gaussian filter
    for (i = 0; i < length; ++i) {
        output[i] = 0.0;
        for (j = -kernel_radius; j <= kernel_radius; ++j) {
            int index = i + j;
            if (index < 0) index = -index; // Reflect boundary
            if (index >= length) index = 2 * length - index - 2; // Reflect boundary
            output[i] += input[index] * kernel[j + kernel_radius];
        }
    }

    // Compute gradient if order is 1
    if (order == 1) {
        float *gradient = (float *)malloc(length * sizeof(float));
        for (i = 1; i < length - 1; ++i) {
            gradient[i] = (output[i + 1] - output[i - 1]) / 2.0f;
        }
        gradient[0] = output[1] - output[0];
        gradient[length - 1] = output[length - 1] - output[length - 2];

        for (i = 0; i < length; ++i) {
            output[i] = gradient[i];
        }
        free(gradient);
    }

    free(kernel);
}

void gaussian_filter_3d(const float *input, float *output, int nx, int ny, int nt, float sigma, int order) {
    int i, j, k;
    float *temp_input = (float *)malloc(nx * sizeof(float));
    float *temp_output = (float *)malloc(nx * sizeof(float));

    // Filter along nx direction
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            for (k = 0; k < nt; k++) {
                temp_input[j] = input[k + j * nt + i * nx * nt];
            }
            gaussian_filter_1d(temp_input, temp_output, nx, sigma, order);
            for (j = 0; j < nx; j++) {
                output[k + j * nt + i * nx * nt] = temp_output[j];
            }
        }
    }

    // Filter along ny direction
    temp_input = (float *)realloc(temp_input, ny * sizeof(float));
    temp_output = (float *)realloc(temp_output, ny * sizeof(float));
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nt; k++) {
                temp_input[j] = output[k + j * nt + i * nx * nt];
            }
            gaussian_filter_1d(temp_input, temp_output, ny, sigma, order);
            for (j = 0; j < ny; j++) {
                output[k + j * nt + i * nx * nt] = temp_output[j];
            }
        }
    }

    // Filter along nt direction
    temp_input = (float *)realloc(temp_input, nt * sizeof(float));
    temp_output = (float *)realloc(temp_output, nt * sizeof(float));
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nt; k++) {
                temp_input[k] = output[k + j * nt + i * nx * nt];
            }
            gaussian_filter_1d(temp_input, temp_output, nt, sigma, order);
            for (k = 0; k < nt; k++) {
                output[k + j * nt + i * nx * nt] = temp_output[k];
            }
        }
    }

    free(temp_input);
    free(temp_output);
}



