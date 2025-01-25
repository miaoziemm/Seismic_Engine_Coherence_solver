#ifndef __SEMINMAX_H__
#define __SEMINMAX_H__

#include "../include/SEmath.h"
#include "../include/SEmaths.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type double
#define minmax_name(name) name ## d
#include "SEminmax_inter.h"


#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type real
#define minmax_name(name) name ## r
#include "SEminmax_inter.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type float
#define minmax_name(name) name ## f
#include "SEminmax_inter.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type int
#define minmax_name(name) name ## i
#include "SEminmax_inter.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type int64_t
#define minmax_name(name) name ## i64
#include "SEminmax_inter.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type size_t
#define minmax_name(name) name ## sz
#include "SEminmax_inter.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type ssize_t
#define minmax_name(name) name ## ssz
#include "SEminmax_inter.h"


#endif