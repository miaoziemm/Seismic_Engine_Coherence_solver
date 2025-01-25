#include "../include/SEsort.h"

#define __SESORT_INTER_H__

#define sort_lt(a, b) ((a) < (b))
#define sort_gt(a, b) ((a) > (b))
#define sort_eq(a, b) FEQUAL((a), (b))

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type float
#define sort_name(name) name ## f
#include "SEsort_inter.c"

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type double
#define sort_name(name) name ## d
#include "SEsort_inter.c"

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type real
#define sort_name(name) name ## r
#include "SEsort_inter.c"

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type reala
#define sort_name(name) name ## ra
#include "SEsort_inter.c"


#ifdef sort_eq
#undef sort_eq
#endif
#define sort_eq(a, b) ((a) == (b))


#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type int
#define sort_name(name) name ## i
#include "SEsort_inter.c"


#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type int64_t
#define sort_name(name) name ## i64
#include "SEsort_inter.c"


#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type size_t
#define sort_name(name) name ## sz
#include "SEsort_inter.c"


#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type byte
#define sort_name(name) name ## b
#include "SEsort_inter.c"

// int64_keyval_t

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#ifdef sort_lt
#undef sort_lt
#endif
#ifdef sort_gt
#undef sort_gt
#endif
#ifdef sort_eq
#undef sort_eq
#endif

#define sort_lt(a, b) ((a).key < (b).key)
#define sort_gt(a, b) ((a).key > (b).key)
#define sort_eq(a, b) ((a).key == (b).key)

#define sort_type int_keyval_t
#define sort_name(name) name ## _int_keyval
#include "SEsort_inter.c"

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#define sort_type int64_keyval_t
#define sort_name(name) name ## _int64_keyval
#include "SEsort_inter.c"

// int64_double_keyval_t

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#ifdef sort_lt
#undef sort_lt
#endif
#ifdef sort_gt
#undef sort_gt
#endif
#ifdef sort_eq
#undef sort_eq
#endif

#define sort_lt(a, b) ((a).key < (b).key)
#define sort_gt(a, b) ((a).key > (b).key)
#define sort_eq(a, b) (FEQUAL((a).key,(b).key))

#define sort_type int64_double_keyval_t
#define sort_name(name) name ## _int64_double_keyval
#include "SEsort_inter.c"

// double_int64_keyval_t

#ifdef sort_type
#undef sort_type
#endif
#ifdef sort_name
#undef sort_name
#endif
#ifdef sort_lt
#undef sort_lt
#endif
#ifdef sort_gt
#undef sort_gt
#endif
#ifdef sort_eq
#undef sort_eq
#endif

#define sort_lt(a, b) ((a).key < (b).key)
#define sort_gt(a, b) ((a).key > (b).key)
#define sort_eq(a, b) ((a).key == (b).key)

#define sort_type double_int64_keyval_t
#define sort_name(name) name ## _double_int64_keyval
#include "SEsort_inter.c"
