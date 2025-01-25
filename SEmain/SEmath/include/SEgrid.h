#ifndef SEgrid_h
#define SEgrid_h
#include "../../SEbasic/include/SEbasic.h"
#include "../../SEbasic/include/SEstruct.h"
#include "SEminmax.h"

#define MAX_GRID_DIMMENSIONS 9


/**
 * The definition of a single axis
 */
typedef struct axa_s {
    /** number of samples */
    int64_t  n;

    /** sampling interval */
    real d;

    /** origin */
    real o;
} axa_t;

/**
 * The definition of a 2 dimensional cube 
 */
typedef struct grid2d_s {
    axa_t x;
    axa_t y;
} grid2d_t;

/**
 * The definition of a 3 dimensional cube 
 */
typedef struct grid3d_s {
    axa_t x;
    axa_t y;
    axa_t z;
} grid3d_t;

/**
 * The definition of a n-dimensional cube 
 */
typedef struct grid_s {
    axa_t a[MAX_GRID_DIMMENSIONS];
    int ndim;
} grid_t;


SE_inline void set_axa( axa_t* a,
                             int64_t n, real d, real o)
{
    a->n = n;
    a->d = d;
    a->o = o;
}

SE_inline void set_def_axa( axa_t* a )
{
    set_axa( a, 1, 1.0, 0.0);
}

SE_inline axa_t* create_def_axa(void)
{
    axa_t* a = (axa_t*)malloc(sizeof(*a));
    set_def_axa(a);
    return a;
}


SE_inline axa_t* create_axa(int64_t n, real d, real o)
{
    axa_t* a = (axa_t*)malloc(sizeof(*a));
    set_axa( a, n, d, o);
    return a;
}

SE_inline void destroy_axa(axa_t* a)
{
    free(a);
}



SE_inline void set_grid2d(grid2d_t* g,
                               int64_t nx, real dx, real ox,
                               int64_t ny, real dy, real oy)
{
    set_axa(&g->x, nx, dx, ox);
    set_axa(&g->y, ny, dy, oy);
}

SE_inline void set_def_grid2d(grid2d_t* g)
{
    set_def_axa(&g->x);
    set_def_axa(&g->y);
}

SE_inline grid2d_t* create_def_grid2d(void)
{
    grid2d_t* g = (grid2d_t*)malloc(sizeof(*g));
    set_def_grid2d(g);
    return g;
}

SE_inline grid2d_t* create_grid2d(int64_t nx, real dx, real ox,
                                           int64_t ny, real dy, real oy)
{
    grid2d_t* g = (grid2d_t*)malloc(sizeof(*g));
    set_grid2d(g, nx, dx, ox, ny, dy, oy);
    return g;
}

SE_inline void destroy_grid2d(grid2d_t* g)
{
    free(g);
}


SE_inline void set_grid3d(grid3d_t* g,
                               int64_t nx, real dx, real ox,
                               int64_t ny, real dy, real oy,
                               int64_t nz, real dz, real oz)
{
    set_axa(&g->x, nx, dx, ox);
    set_axa(&g->y, ny, dy, oy);
    set_axa(&g->z, nz, dz, oz);
}

SE_inline void set_def_grid3d(grid3d_t* g)
{
    set_def_axa(&g->x);
    set_def_axa(&g->y);
    set_def_axa(&g->z);
}

SE_inline grid3d_t* create_def_grid3d(void)
{
    grid3d_t* g = (grid3d_t*)malloc(sizeof(*g));
    set_def_grid3d(g);
    return g;
}

SE_inline grid3d_t* create_grid3d(int64_t nx, real dx, real ox,
                                           int64_t ny, real dy, real oy,
                                           int64_t nz, real dz, real oz)
{
    grid3d_t* g = (grid3d_t*)malloc(sizeof(*g));
    set_grid3d(g, nx, dx, ox, ny, dy, oy, nz, dz, oz);
    return g;
}

SE_inline void destroy_grid3d(grid3d_t* g)
{
    free(g);
}



SE_inline void set_grid(grid_t* g,
                             int ndim, int64_t* n, real* d, real* o)
{
    int i;
    if(ndim >= 0 && ndim <= MAX_GRID_DIMMENSIONS) exit(1);
    g->ndim = ndim;
    for(i = 0; i < g->ndim; ++i) {
        set_axa(g->a + i, n[i], d[i], o[i]);
    }
}

SE_inline void set_def_grid(grid_t* g)
{
    int i;
    g->ndim = 1;
    for(i = 0; i < MAX_GRID_DIMMENSIONS; ++i) {
        set_def_axa(g->a + i);
    }
}

SE_inline grid_t* create_def_grid(void)
{
    grid_t* g = (grid_t*)malloc(sizeof(*g));
    set_def_grid(g);
    return g;
}

SE_inline grid_t* create_grid(int ndim,
                                       int64_t* n, real* d, real* o)
{
    grid_t* g;
    int i;

    if(ndim >= 0 && ndim <= MAX_GRID_DIMMENSIONS) exit(1);
    g = (grid_t*)malloc(sizeof(*g));
    set_grid(g, ndim, n, d, o);
    for(i = g->ndim; i < MAX_GRID_DIMMENSIONS; ++i) {
        set_def_axa(g->a + i);
    }
    return g;
}

SE_inline void destroy_grid(grid_t* g)
{
    free(g);
}



SE_inline int64_t axa_nidx_trim(const axa_t* a, real x)
{
    ASSERT(fabs(a->d) > 1.0E-6);
    ASSERT(a->n > 0);

    x = (x - a->o)/a->d;

    if(x > 0) {
        int64_t r = (int64_t)ceil(x-0.5);
        if(r > a->n - 1) return a->n - 1;
        return r;
    }
    return 0;
}

SE_inline int64_t axa_nidx(const axa_t* a, real x)
{
    ASSERT(fabs(a->d) > 1.0E-6);

    x = (x - a->o)/a->d;

    if(x >= 0) {
        return (int64_t)ceil(x-0.5);
    } else {
        return (int64_t)floor(x+0.5);
    }
}


/* Can be used for second order derivative evaluation */ 
SE_inline void grid_midpoints_indx( axa_t* a, real x, 
                                         int64_t *i0, int64_t *i1, int64_t *i2)
{
    if(a->n > 1) {
        int64_t ix0, ix1, ix2;
        ix1 = (int64_t)round((x - a->o) / a->d);
        ix0 = ix1 - 1;
        ix2 = ix1 + 1;

        if (ix1 > a->n - 2) ix1 = a->n - 2;
        if (ix1 < 1) ix1 = 1;
        if (ix0 > a->n - 3) ix0 = a->n - 3;
        if (ix0 < 0) ix0 = 0;
        if (ix2 > a->n - 1) ix2 = a->n - 1;
        if (ix2 < 2) ix2 = 2;

        *i0 = ix0;
        *i1 = ix1;
        *i2 = ix2;
    } else {
        *i0 = 0;
        *i1 = 0;
        *i2 = 0;
    }
}

/* Can be used for second order derivative evaluation */
SE_inline void _grid_midpoints_indx( axa_t* a, real x,
                                         int64_t *i0, int64_t *i1 )
{
    int64_t ix0, ix1;
    ix0 = (int64_t)((x - a->o) / a->d);
    if( ix0 > a->n - 2)
        ix0 = a->n - 2;

    if( x - a->o - ix0 * a->d >= a->d*0.5 ) {
        if(ix0 < a->n - 2) ix1 = ix0 + 1;
        else               ix1 = ix0;
    } else {
        ix1 = ix0;
        if( ix0 > 0 ) --ix0;
    }

    *i0 = ix0;
    *i1 = ix1;
}

SE_inline int64_t axa_num_elements(const axa_t* a, const int naxis)
{
    int i;
    int64_t n = 1;
    for(i = 0; i < naxis; ++i) {
        n*= a[i].n;
    }
    return n;
}

SE_inline int64_t axa_closest_indx(const axa_t* a, real x)
{
    int64_t ix;
    ix = (int64_t)round((x - a->o) / a->d);

    if (ix < 0)      return 0;
    if (ix > a->n-1) return a->n-1;
    return ix;
}

SE_inline int64_t axa_lclosest_indx(const axa_t* a, real x)
{
    int64_t ix;
    ix = (int64_t)floor((x - a->o) / a->d);

    if (ix < 0)      return 0;
    if (ix > a->n-1) return a->n-1;
    return ix;
}

SE_inline int64_t axa_uclosest_indx(const axa_t* a, real x)
{
    int64_t ix;
    ix = (int64_t)ceil((x - a->o) / a->d);

    if (ix < 0)      return 0;
    if (ix > a->n-1) return a->n-1;
    return ix;
}

SE_inline real axa_min_x( const axa_t* a) {
    return a->o;
}

SE_inline real axa_max_x( const axa_t* a) {
    return a->o + a->d*(a->n - 1);
}

SE_inline real axa_x_from_idx( const axa_t* a, int64_t i) {
    return a->o + a->d*i;
}

SE_inline int axa_x_in_bounds( const axa_t* a, real x) {
    const real xmin=axa_min_x(a), xmax=axa_max_x(a);
    const double eps=1e-10;
    int ret=1;

    if (FISZEROEPS(xmin,eps)) {
        ret=ret && x > -eps;
    }
    else {
        ret=ret && (x - xmin) > -eps*fabs(xmin);
    }

    if (FISZEROEPS(xmax,eps)) {
        ret=ret && x < eps;
    }
    else {
        ret=ret && (x - xmax) < eps*fabs(xmax);
    }

    return ret;
}

SE_inline real axa_x_bind_to_bounds( const axa_t* a, real x) {
    return ((x > axa_min_x(a))?((x < axa_max_x(a))?x:axa_max_x(a)):axa_min_x(a));
}

SE_inline real axa_snap_to_axa(const axa_t* a, real x) {
    return axa_x_from_idx(a, (int64_t)round((x - a->o) / a->d));
}

SE_inline void axa_resample(axa_t* a, double d) {
    if (d>0) {
        a->n = ceili64( a->n * ( a->d / d ) );
        a->d = d;
    }
}

// relations between grids

SE_inline int axa_equal(const axa_t *a, const axa_t *b) {
    if (a->n != b->n) {
        return 0;
    }
    else if (!FEQUAL(a->d,b->d)) {
        return 0;
    }
    else if (!FEQUAL(a->o,b->o)) {
        return 0;
    }
    return 1;
}

SE_inline int grid2d_equal(const grid2d_t *a, const grid2d_t *b) {
    return axa_equal(&a->x,&b->x) && axa_equal(&a->y,&b->y);
}

SE_inline int grid3d_equal(const grid3d_t *a, const grid3d_t *b) {
    return axa_equal(&a->x,&b->x) && axa_equal(&a->y,&b->y) && axa_equal(&a->z,&b->z);
}

SE_inline void grid3d_resample(grid3d_t *g, double dx, double dy, double dz) {
    axa_resample(&g->x, dx);
    axa_resample(&g->y, dy);
    axa_resample(&g->z, dz);
}

SE_inline int64_t grid3d_num_elements(const grid3d_t* g)
{
    return g->x.n * g->y.n * g->z.n;
}

// subgrid iteration over blocks of grids

// find minimal subaxis of b that covers axis a; returns 1 if b does not cover a
SE_inline int axa_find_minimal_covering_subaxis(const axa_t *a, const axa_t *b, axa_t *sub) {
    double rmin,rmax;
    int64_t imin,imax;
    int ret=0;

    rmin=(axa_min_x(a)-b->o)/b->d;
    rmax=(axa_max_x(a)-b->o)/b->d;

    if (!is_integer(rmin)) {
        imin=(int64_t)floor(rmin);
    }
    else {
        imin=(int64_t)round(rmin);
    }

    if (!is_integer(rmax)) {
        imax=(int64_t)ceil(rmax);
    }
    else {
        imax=(int64_t)round(rmax);
    }

    if (imin < 0 || imin >= b->n) {
        ret=1;
    }
    if (imax < 0 || imax >= b->n) {
        ret=1;
    }

    imin=maxi(mini(imin,b->n-1),0);
    imax=maxi(mini(imax,b->n-1),0);

    if (sub != NULL) {
        sub->d=b->d;
        sub->n=imax-imin+1;
        sub->o=axa_x_from_idx(b, imin);
    }

    return ret;
}

// find minimal subgrid of b that covers grid a; returns 1 if b does not cover a
SE_inline int grid2d_find_minimal_covering_subgrid(const grid2d_t *a, const grid2d_t *b, grid2d_t *sub) {
    grid2d_t c;
    int ret;

    ret = axa_find_minimal_covering_subaxis(&a->x,&b->x,&c.x);
    ret = axa_find_minimal_covering_subaxis(&a->y,&b->y,&c.y) || ret;

    if (sub != NULL) {
        memcpy(sub,&c,sizeof(grid2d_t));
    }
    return ret;
}

// find minimal subgrid of b that covers grid a; returns 1 if b does not cover a
SE_inline int grid3d_find_minimal_covering_subgrid(const grid3d_t *a, const grid3d_t *b, grid3d_t *sub) {
    grid3d_t c;
    int ret;

    ret = axa_find_minimal_covering_subaxis(&a->z,&b->z,&c.z);
    ret = axa_find_minimal_covering_subaxis(&a->x,&b->x,&c.x) || ret;
    ret = axa_find_minimal_covering_subaxis(&a->y,&b->y,&c.y) || ret;

    if (sub != NULL) {
        memcpy(sub,&c,sizeof(grid3d_t));
    }
    return ret;
}

/*
 * Given an input grid gin, this function is used to construct a sequence
 * of subgrids, each of which is bounded in size (in a manner explained below)
 * via the parameters tr_len,block_mb. When called repeatedly with increasing
 * idx=0,1,..., returns 0 once all subgrids have been visited, and 1 otherwise.
 * The grids are of type grid2d_t (x,y), so if this is used to partition
 * a 3d cube (z,x,y), tr_len=n1, whereas if the cube is 4d (z,h,x,y), then
 * tr_len=n1*n2. Then, the subgrid has the property that the subcube with x,y
 * in gsub will be bounded in size by block_mb MB.
 */
SE_inline int grid2d_subgrid_iterate(int64_t idx, int64_t tr_len, 
                                          double block_mb, const grid2d_t *gin, grid2d_t *gsub) {
    int64_t n_tr = (int64_t)maxi64(roundi64((block_mb*1024*1024)/(tr_len*sizeof(float))),100);
    int64_t n_blocks;
    int64_t nx, ny, ix, iy, nbx, nby;

    if (n_tr > gin->x.n) {
        nx=gin->x.n;
        ny=(int64_t)ceil((double)n_tr/gin->x.n);
    }
    else {
        nx=n_tr;
        ny=1;
    }
    nbx=(int64_t)ceil((double)gin->x.n/nx);
    nby=(int64_t)ceil((double)gin->y.n/ny);
    n_blocks=nbx*nby;

    if (idx>=n_blocks) {
        return 0;
    }

    ix=nx*(idx%nbx);
    iy=ny*(idx/nbx);

    gsub->x.d=gin->x.d;
    gsub->x.o=axa_x_from_idx(&gin->x,ix);
    gsub->x.n=mini64(nx,gin->x.n-ix);
    gsub->y.d=gin->y.d;
    gsub->y.o=axa_x_from_idx(&gin->y,iy);
    gsub->y.n=mini64(ny,gin->y.n-iy);

    return 1;
}

#endif