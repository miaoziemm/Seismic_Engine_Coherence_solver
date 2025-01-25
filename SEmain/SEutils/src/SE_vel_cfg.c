#include "../include/SE_vel_cfg.h"

vel_cfg_t *vel_cfg_create(int argc, char **argv)
{
    initargs(argc, argv);

    vel_cfg_t *cfg;
    cfg = (vel_cfg_t *)malloc(sizeof(vel_cfg_t));

    if(!getparstring("vel", &cfg->vel_file)) {
        SE_ERROR("Velocity file not specified.");
    }
    if(!getparint("anisotropic", &cfg->aniso)) {
        cfg->aniso = 0;
    }
    if(!getparstring("epsilon", &cfg->eps_file)) {
        cfg->eps_file = NULL;
    }
    if(!getparstring("delta", &cfg->del_file)) {
        cfg->del_file = NULL;
    }
    if(!getparint("tti", &cfg->tti)) {
        cfg->tti = 0;
    }
    if(!getparint("theta_dip", &cfg->tet_dip)) {
        cfg->tet_dip = 1;
    }
    if(!getparstring("theta_x", &cfg->tetx_file)) {
        cfg->tetx_file = NULL;
    }
    if(!getparstring("theta_y", &cfg->tety_file)) {
        cfg->tety_file = NULL;
    }
    if(!getparint("force_3d", &cfg->force_3d)) {
        cfg->force_3d = 0;
    }
    if(!getparint("interpolation_stencil", &cfg->vel_stencil)) {
        cfg->vel_stencil = 3;
    }
    if(!getparint("detect_velocity_boundaries", &cfg->detect_vel_bdry)) {
        cfg->detect_vel_bdry = 1;
    }
    if(!getpardouble("velocity_boundary_fraction", &cfg->vel_bdry_frac)) {
        cfg->vel_bdry_frac = 0.25;
    }

    return cfg;

}
