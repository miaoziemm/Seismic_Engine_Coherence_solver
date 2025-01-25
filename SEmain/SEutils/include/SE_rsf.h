/**
 * \file rsf.h
 *
 * Definitions for working with regularly sampled functions.
 *关于用常规采样函数工作的定义
 * Copyright (c) 2015 Z-Terra Inc. All Rights Reserved.
 * THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF Z-Terra Inc.
 *
 * \author Iulian Musat
 */

#ifndef __SE_RSF_H__
#define __SE_RSF_H__

#include <string.h>

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEgrid.h"


/*#define RSF_CAN_USE_NONDEFAULT_LAYOUT*/  //非默认布局 非非默认格式布局

typedef reala rsf_value;  //reala = float

typedef struct rsf2d_s rsf2d_t;
struct rsf2d_s {
    grid2d_t g;       //grid2d_t 2D值 变长数组
    rsf_value* (*v)(const rsf2d_t* ths, int64_t ix, int64_t iy);
    void* prvt;
};

typedef struct rsf3d_s rsf3d_t;
struct rsf3d_s {
    grid3d_t g;
    rsf_value* (*v)(const rsf3d_t* ths, int64_t ix, int64_t iy, int64_t iz);
    void* prvt;
};

typedef struct rsf_s rsf_t;
struct rsf_s {
    grid_t g;    //n维数组
    rsf_value* (*v)(const rsf_t* ths, ...);
    void* prvt;
};



typedef struct rsf2d_blocking_private_s {
    int64_t bsx;
    int64_t bsy;
    rsf_value* data;
}rsf2d_blocking_private_t;

typedef struct rsf3d_blocking_private_s {
    int64_t bsx;
    int64_t bsy;
    int64_t bsz;
    rsf_value* data;
}rsf3d_blocking_private_t;



SE_inline rsf_value* rsf_idx2_default(const rsf2d_t* rsf,
                                               int64_t ix, int64_t iy)
{
#ifdef DEBUG
    if(ix < 0 || ix >= rsf->g.x.n)
        ERROR(("Internal error: accessing RSF2D ix=%Ld out of range [0,%Ld)", ix, rsf->g.x.n));
    if(iy < 0 || iy >= rsf->g.y.n)
        ERROR(("Internal error: accessing RSF2D iy=%Ld out of range [0,%Ld)", iy, rsf->g.y.n));
#endif
    return (rsf_value*)(rsf->prvt) + 
        iy*rsf->g.x.n + ix;
}

SE_inline rsf_value* rsf_idx3_default(const rsf3d_t* rsf,
                                               int64_t ix, int64_t iy, int64_t iz)
{
#ifdef DEBUG
    if(ix < 0 || ix >= rsf->g.x.n)
        ERROR(("Internal error: accessing RSF3D ix=%Ld out of range [0,%Ld)", ix, rsf->g.x.n));
    if(iy < 0 || iy >= rsf->g.y.n)
        ERROR(("Internal error: accessing RSF3D iy=%Ld out of range [0,%Ld)", iy, rsf->g.y.n));
    if(iz < 0 || iz >= rsf->g.z.n)
        ERROR(("Internal error: accessing RSF3D iz=%Ld out of range [0,%Ld)", iz, rsf->g.z.n));
#endif
    return (rsf_value*)rsf->prvt + 
        iy*rsf->g.x.n*rsf->g.z.n + ix*rsf->g.z.n + iz;
}


SE_inline rsf_value* rsf_idx2_blocks(rsf2d_t* rsf,
                                              int64_t ix, int64_t iy)
{
    rsf2d_blocking_private_t* rp = (rsf2d_blocking_private_t*)rsf->prvt;

    /* indexes of the current block */
    int64_t ibx = ix / rp->bsx;
    int64_t iby = iy / rp->bsy;

    /* indexes inside the current block */
    int64_t iix = ix % rp->bsx;
    int64_t iiy = iy % rp->bsy;
    
    /* size of the current block */ 
    int64_t bsx = ((ibx+1)*rp->bsx <= rsf->g.x.n) ? 
        rp->bsx : rsf->g.x.n%rp->bsx;
    int64_t bsy = ((iby+1)*rp->bsy <= rsf->g.y.n) ? 
        rp->bsy : rsf->g.y.n%rp->bsy;
    
    return rp->data + 
        /* size of full columns of blocks until current column */
        iby*rp->bsy*rsf->g.x.n + 
        /* size of all blocks from current column until current block */
        ibx*bsy*rp->bsx +
        /* offset inside current block */
        iiy*bsx + iix;
}

SE_inline rsf_value* rsf_idx3_blocks(rsf3d_t* rsf,
                                              int64_t ix, int64_t iy, int64_t iz)
{
    rsf3d_blocking_private_t* rp = (rsf3d_blocking_private_t*)rsf->prvt;

    /* indexes of the current block */
    int64_t ibx = ix / rp->bsx;
    int64_t iby = iy / rp->bsy;
    int64_t ibz = iz / rp->bsz;
    /* indexes inside the current block */
    int64_t iix = ix % rp->bsx;
    int64_t iiy = iy % rp->bsy;
    int64_t iiz = iz % rp->bsz;
    
    /* size of the current block */ 
    int64_t bsx = ((ibx+1)*rp->bsx <= rsf->g.x.n) ? 
        rp->bsx : rsf->g.x.n%rp->bsx;
    int64_t bsy = ((iby+1)*rp->bsy <= rsf->g.y.n) ? 
        rp->bsy : rsf->g.y.n%rp->bsy;
    int64_t bsz = ((ibz+1)*rp->bsz <= rsf->g.z.n) ? 
        rp->bsz : rsf->g.z.n%rp->bsz;
    
    return rp->data + 
        iby*rsf->g.z.n*rsf->g.x.n*rp->bsy + 
        ibx*rsf->g.z.n*rp->bsx*bsy +
        ibz*rp->bsz*bsx*bsy +
        iiy*bsx*bsz + iix*bsz + iiz;
}

SE_inline void rsf_set2d_defmem(int64_t nx, double dx, double ox,
                                     int64_t ny, double dy, double oy,
                                     rsf_value* data, rsf2d_t* rsf)
{
    set_grid2d(&rsf->g, nx, dx, ox, ny, dy, oy);
    rsf->v = rsf_idx2_default;
    if(data == NULL) {
        data = alloc1float(nx*ny);
    }
    rsf->prvt = data;
}

SE_inline rsf2d_t* rsf_create2d_defmem(int64_t nx, double dx, double ox,
                                                int64_t ny, double dy, double oy,
                                                rsf_value* data)
{
    rsf2d_t* rsf = (rsf2d_t*)malloc(sizeof(*rsf));
    rsf_set2d_defmem(nx, dx, ox, ny, dy, oy, data, rsf);
    return rsf;
}

SE_inline void rsf_set2d_defmem_from_grid(const grid2d_t* g, rsf_value* data, rsf2d_t* rsf) {
    rsf_set2d_defmem(g->x.n, g->x.d, g->x.o,
                         g->y.n, g->y.d, g->y.o,
                         data, rsf);
}

SE_inline rsf2d_t*
rsf_create2d_defmem_from_grid(const grid2d_t* g, rsf_value* data)
{
    return rsf_create2d_defmem(g->x.n, g->x.d, g->x.o,
                                   g->y.n, g->y.d, g->y.o,
                                   data);
}

SE_inline void rsf_destroy2d_defmem(rsf2d_t* rsf)
{
    free(rsf);
}

SE_inline void rsf_destroy2d_and_data_defmem(rsf2d_t* rsf)
{
    free(rsf->prvt);
    free(rsf);
}

SE_inline void rsf_set3d_defmem(int64_t nx, double dx, double ox,
                                     int64_t ny, double dy, double oy,
                                     int64_t nz, double dz, double oz,
                                     rsf_value* data, rsf3d_t* rsf) {
    set_grid3d(&rsf->g, nx, dx, ox, ny, dy, oy, nz, dz, oz);
    rsf->v = rsf_idx3_default;
    if(data == NULL) {
        data = alloc1float(nx*ny*nz);
    }
    rsf->prvt = data;
}

SE_inline rsf3d_t* rsf_create3d_defmem(int64_t nx, double dx, double ox,
                                                int64_t ny, double dy, double oy,
                                                int64_t nz, double dz, double oz,
                                                rsf_value* data)
{
    rsf3d_t* rsf = (rsf3d_t*)malloc(sizeof(*rsf));
    rsf_set3d_defmem(nx, dx, ox, ny, dy, oy, nz, dz, oz, data, rsf);
    return rsf;
}

SE_inline void rsf_set3d_defmem_from_grid(const grid3d_t* g, rsf_value* data, rsf3d_t *rsf) {
    rsf_set3d_defmem(g->x.n, g->x.d, g->x.o,
                         g->y.n, g->y.d, g->y.o,
                         g->z.n, g->z.d, g->z.o,
                         data, rsf);
}

SE_inline rsf3d_t*
rsf_create3d_defmem_from_grid(const grid3d_t* g,
                                  rsf_value* data)
{
    return rsf_create3d_defmem(g->x.n, g->x.d, g->x.o,
                                   g->y.n, g->y.d, g->y.o,
                                   g->z.n, g->z.d, g->z.o,
                                   data);
}

SE_inline rsf3d_t*
rsf_copy3d_defmem(const rsf3d_t* org)
{
    rsf3d_t* rsf = rsf_create3d_defmem_from_grid(&org->g, NULL);
    memcpy(rsf->prvt, org->prvt, org->g.x.n*org->g.y.n*org->g.z.n*sizeof(rsf_value));
    return rsf;
}

SE_inline rsf2d_t*
rsf_copy2d_defmem(const rsf2d_t* org)
{
    rsf2d_t* rsf = rsf_create2d_defmem_from_grid(&org->g, NULL);
    memcpy(rsf->prvt, org->prvt, org->g.x.n*org->g.y.n*sizeof(rsf_value));
    return rsf;
}

SE_inline reala* rsf_getinternal_buffer3d_defmem(const rsf3d_t* rsf)
{
    return (reala*)rsf->prvt;
}
SE_inline reala* rsf_getinternal_buffer2d_defmem(const rsf2d_t* rsf)
{
    return (reala*)rsf->prvt;
}

SE_inline void rsf_zero3d_defmem(rsf3d_t* rsf)
{
    memset(rsf->prvt, 0, rsf->g.x.n*rsf->g.y.n*rsf->g.z.n*sizeof(rsf_value));
}

SE_inline void rsf_zero2d_defmem(rsf2d_t* rsf)
{
    memset(rsf->prvt, 0, rsf->g.x.n*rsf->g.y.n*sizeof(rsf_value));
}

SE_inline void rsf_destroy3d_defmem(rsf3d_t* rsf)
{
    free(rsf);
}

SE_inline void rsf_destroy3d_and_data_defmem(rsf3d_t* rsf)
{
    free(rsf->prvt);
    free(rsf);
}

SE_inline rsf2d_t*
rsf_create2d_defmem_from_3d(rsf3d_t* rsf3d)
{
    grid2d_t g;
    g.x = rsf3d->g.z;
    g.y = rsf3d->g.x;
    return rsf_create2d_defmem_from_grid(&g,
                                             (rsf_value*)rsf3d->prvt);
}

SE_inline rsf2d_t*
rsf_create2d_defmem_from_3d_transp(rsf3d_t* rsf3d)
{
    grid2d_t g;
    g.x = rsf3d->g.z;
    g.y = rsf3d->g.y;
    return rsf_create2d_defmem_from_grid(&g,
                                             (rsf_value*)rsf3d->prvt);
}
#ifdef RSF_CAN_USE_NONDEFAULT_LAYOUT
SE_inline rsf2d_t*
rsf_create2d_blockmem(int64_t nx, double dx, double ox,
                          int64_t ny, double dy, double oy,
                          int64_t bsx, int64_t bsy,
                          rsf_value* data)
{
    rsf2d_t* rsf = (rsf2d_t*)
        alloc(sizeof(*rsf) + sizeof(rsf2d_blocking_private_t));
    rsf2d_blocking_private_t* rp;

    set_grid2d(&rsf->g, nx, dx, ox, ny, dy, oy);
    rsf->v = rsf_idx2_blocks;
    rsf->prvt = &rsf->prvt + 1;
    rp = rsf->prvt;
    rp->bsx = bsx;
    rp->bsy = bsy;
    if(data == NULL) {
        data = alloc_float(nx*ny);
    }
    rp->data = data;
    return rsf;
}

SE_inline void rsf_destroy2d_blockmem(rsf2d_t* rsf)
{
    free(rsf);
}

SE_inline rsf3d_t*
rsf_create3d_blockmem(int64_t nx, double dx, double ox,
                          int64_t ny, double dy, double oy,
                          int64_t nz, double dz, double oz,
                          int64_t bsx, int64_t bsy, int64_t bsz,
                          rsf_value* data)
{
    rsf3d_t* rsf = (rsf3d_t*)
        alloc(sizeof(*rsf) + sizeof(rsf3d_blocking_private_t));
    rsf3d_blocking_private_t* rp;

    set_grid3d(&rsf->g, nx, dx, ox, ny, dy, oy, nz, dz, oz);
    rsf->v = rsf_idx3_blocks;
    rsf->prvt = &rsf->prvt + 1;
    rp = rsf->prvt;
    rp->bsx = bsx;
    rp->bsy = bsy;
    rp->bsz = bsz;
    if(data == NULL) {
        data = alloc_float(nx*ny*nz);
    }
    rp->data = data;
    return rsf;
}

SE_inline void rsf_destroy3d_blockmem(rsf3d_t* rsf)
{
    free(rsf);
}

#define I2(rsf,ix,iy) (rsf)->v((rsf), (ix), (iy))
#define I3(rsf,ix,iy,iz) (rsf)->v((rsf), (ix), (iy), (iz))

#else /*!RSF_CAN_USE_NONDEFAULT_LAYOUT*/

SE_inline rsf2d_t*
rsf_create2d_blockmem(int64_t nx, double dx, double ox,
                          int64_t ny, double dy, double oy,
                          int64_t bsx, int64_t bsy,
                          rsf_value* data)
{
    return rsf_create2d_defmem(nx, dx, ox, ny, dy, oy, data);
}

SE_inline void rsf_destroy2d_blockmem(rsf2d_t* rsf)
{
    rsf_destroy2d_defmem(rsf);
}

SE_inline rsf3d_t*
rsf_create3d_blockmem(int64_t nx, double dx, double ox,
                          int64_t ny, double dy, double oy,
                          int64_t nz, double dz, double oz,
                          int64_t bsx, int64_t bsy, int64_t bsz,
                          rsf_value* data)
{
    return rsf_create3d_defmem(nx, dx, ox, ny, dy, oy, nz, dz, oz, data);
}

SE_inline void rsf_destroy3d_blockmem(rsf3d_t* rsf)
{
    rsf_destroy3d_defmem(rsf);
}

#define I2(rsf,ix,iy) rsf_idx2_default((rsf), (ix), (iy))
#define I3(rsf,ix,iy,iz) rsf_idx3_default((rsf), (ix), (iy), (iz))

#endif /*RSF_CAN_USE_NONDEFAULT_LAYOUT*/


#endif /* __RSF_H__ */
