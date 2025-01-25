#ifndef SE_VEL_CFG_H
#define SE_VEL_CFG_H
#include "../../SEbasic/include/SEbasic.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEgrid.h"
/**
 * Beam depth migration configuration.
 */
typedef struct vel_cfg_s {

    char *vel_file; /**< Velocity filename. */
    char *eps_file; /**< Epsilon filename. */
    char *del_file; /**< Delta filename. */
    char *tetx_file; /**< Theta X filename (for TTI). */ //Theta
    char *tety_file; /**< Theta Y filename (for TTI). */  //Theta

    int aniso; /**< Flag for anisotropic migration. */  
    int tti; /**< Flag for TTI anisotropic migration. */  
    int tet_dip; /**< Flag for TTI Theta files contain dips or angles (degrees). */  

    int force_3d; /**< Flag for forcing 3D. */

    int vel_stencil; /**< Velocity interpolation stencil half size. */ 
    int detect_vel_bdry; /**< Flag to detect large velocity changes. */
    double vel_bdry_frac; /**< Fractional threshold for detecting boundaries (max_vel-min_vel)/max_vel */ //分数阈值 探测边界

} vel_cfg_t;

vel_cfg_t *vel_cfg_create(int argc, char **argv);

SE_inline void vel_cfg_destroy(vel_cfg_t* p) 
{
    if (p) {
        free(p->vel_file);
        free(p->eps_file);
        free(p->del_file);
        free(p->tetx_file);
        free(p->tety_file);    
        free(p);
    }
}



#endif