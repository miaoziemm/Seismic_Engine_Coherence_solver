
#ifndef SESORT_H
#define SESORT_H

#include "../include/SEmath.h"
#include "../include/SEmaths.h"
#include "../include/SEminmax.h"


/** Generic signature for a function that can be used to sort float arrays
 *  \param x   - the float array to be sorted
 *  \param off - the offset in the array to start sorting
 *  \param len - the number of elements to sort
 */
typedef void (*sortf_fc)(float* x, size_t off, size_t len);
void insert_sortf(float* x, size_t off, size_t len);
void quick_sortf(float* x, size_t off, size_t len);
void bubble_sortf(float* x, size_t off, size_t len);
ssize_t binary_searchf(float key, const float* a, size_t from, size_t to);

typedef void (*sortd_fc)(double* x, size_t off, size_t len);
void insert_sortd(double* x, size_t off, size_t len);
void quick_sortd(double* x, size_t off, size_t len);
void bubble_sortd(double* x, size_t off, size_t len);
ssize_t binary_searchd(double key, const double* a, size_t from, size_t to);

typedef void (*sortr_fc)(real* x, int off, int len);
void insert_sortr(real* x, size_t off, size_t len);
void quick_sortr(real* x, size_t off, size_t len);
void bubble_sortr(real* x, size_t off, size_t len);
ssize_t binary_searchr(real key, const real* a, size_t from, size_t to);

typedef void (*sortra_fc)(reala* x, size_t off, size_t len);
void insert_sortra(reala* x, size_t off, size_t len);
void quick_sortra(reala* x, size_t off, size_t len);
void bubble_sortra(reala* x, size_t off, size_t len);
ssize_t binary_searchra(reala key, const reala* a, size_t from, size_t to);

typedef void (*sorti_fc)(int* x, size_t off, size_t len);
void insert_sorti(int* x, size_t off, size_t len);
void quick_sorti(int* x, size_t off, size_t len);
void bubble_sorti(int* x, size_t off, size_t len);
ssize_t binary_searchi(int key, const int* a, size_t from, size_t to);

typedef void (*sorti64_fc)(int64_t* x, size_t off, size_t len);
void insert_sorti64(int64_t* x, size_t off, size_t len);
void quick_sorti64(int64_t* x, size_t off, size_t len);
void bubble_sorti64(int64_t* x, size_t off, size_t len);
ssize_t binary_searchi64(int64_t key, const int64_t* a, size_t from, size_t to);

typedef void (*sortsz_fc)(size_t* x, size_t off, size_t len);
void insert_sortsz(size_t* x, size_t off, size_t len);
void quick_sortsz(size_t* x, size_t off, size_t len);
void bubble_sortsz(size_t* x, size_t off, size_t len);
ssize_t binary_searchsz(size_t key, const size_t* a, size_t from, size_t to);

typedef void (*sortb_fc)(byte* x, size_t off, size_t len);
void insert_sortb(byte* x, size_t off, size_t len);
void quick_sortb(byte* x, size_t off, size_t len);
void bubble_sortb(byte* x, size_t off, size_t len);
ssize_t binary_searchb(byte key, const byte* a, size_t from, size_t to);

/*some specialization routines */
typedef struct int_pairs_s {
    int key;
    int val;
} int_keyval_t;

typedef void (*sort_int_keyval_fc)(int_keyval_t* x, size_t off, size_t len);
void insert_sort_int_keyval(int_keyval_t* x, size_t off, size_t len);
void quick_sort_int_keyval(int_keyval_t* x, size_t off, size_t len);
void bubble_sort_int_keyval(int_keyval_t* x, size_t off, size_t len);
ssize_t binary_search_int_keyval(int_keyval_t key, const int_keyval_t* a, size_t from, size_t to);

typedef struct int64_pairs_s {
    int64_t key;
    int64_t val;
} int64_keyval_t;

typedef void (*sort_int64_keyval_fc)(int64_keyval_t* x, size_t off, size_t len);
void insert_sort_int64_keyval(int64_keyval_t* x, size_t off, size_t len);
void quick_sort_int64_keyval(int64_keyval_t* x, size_t off, size_t len);
void bubble_sort_int64_keyval(int64_keyval_t* x, size_t off, size_t len);
ssize_t binary_search_int64_keyval(int64_keyval_t key, const int64_keyval_t* a, 
                                       size_t from, size_t to);

typedef struct int64_double_s {
    double key;
    int64_t val;
} int64_double_keyval_t;

typedef void (*sort_int64_double_keyval_fc)(int64_double_keyval_t* x, size_t off, size_t len);
void insert_sort_int64_double_keyval(int64_double_keyval_t* x, size_t off, size_t len);
void quick_sort_int64_double_keyval(int64_double_keyval_t* x, size_t off, size_t len);
void bubble_sort_int64_double_keyval(int64_double_keyval_t* x, size_t off, size_t len);
ssize_t binary_search_int64_double_keyval(int64_double_keyval_t key, const int64_double_keyval_t* a,
                                              size_t from, size_t to);

typedef struct double_int64_s {
    int64_t key;
    double val;
} double_int64_keyval_t;

typedef void (*sort_double_int64_keyval_fc)(double_int64_keyval_t* x, size_t off, size_t len);
void insert_sort_double_int64_keyval(double_int64_keyval_t* x, size_t off, size_t len);
void quick_sort_double_int64_keyval(double_int64_keyval_t* x, size_t off, size_t len);
void bubble_sort_double_int64_keyval(double_int64_keyval_t* x, size_t off, size_t len);
ssize_t binary_search_double_int64_keyval(double_int64_keyval_t key, const double_int64_keyval_t* a,
                                              size_t from, size_t to);

#endif