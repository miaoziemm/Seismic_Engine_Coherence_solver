#ifndef SE_COORDINATE_SYSTEM_H
#define SE_COORDINATE_SYSTEM_H

#include "../../SEbasic/include/SEbasic.h"
#include "../../SEmath/include/SEmath.h"
#include "../../SEmath/include/SEgrid.h"
typedef enum
{
    SURVEY_UNITS_METERS = 0,
    SURVEY_UNITS_FEET = 1,

    SURVEY_UNITS_UNKNOWN = -1
} survey_units_t; // 测量单位

/**
 * This structure contains the necessary information to define an
 * arbitrarily positioned and rotated rectangular grid with an
 * associated labeling convention (in-line and cross-line numbers).
 *
 * Usually, this structure is created from parameters available to the
 * program at run-time, either directly specified or read from some
 * common place (ex. project db). See get_survey_description and
 * get_named_survey_description functions.
 *
 * Each member of the structure has an associated flag (has_XXX) which
 * is non-zero if the corresponding member was found among the
 * availabe parameters, and zero if it just contains a default value.
 */
typedef struct survey_description_s
{
    /**
     * Optional name of this survey. Not read from parameters, but
     * given by the application.
     */
    char name[256];
    int has_name;

    /**
     * Origin of the coordinate system along the X direction (usually
     * West->East direction
     */
    double ox_survey;
    int has_ox_survey;

    /**
     * Origin of the coordinate system along the Y direction (usually
     * South->North direction
     */
    double oy_survey;
    int has_oy_survey;

    /**
     * The angle, in degrees of the local cross-line axis from the
     * global Y axis. It increases clockwise.
     *<code>
     *          ^     Y axis
     *          | usually South->North direction
     *          |
     *          |        / Line axis
     *          |       /    lines number increase/decrease
     *          |      /     along this axis - n_xline axis
     *          |     /
     *          |    /
     *          |   /
     *          |  /
     *          | /
     *          |/
     *           -----------------------------------> X axis,
     *           |                       usually West->East direction
     *            |
     *             |
     *              |   Cross-line axis
     *               |    cross-line numbers increase/decrease
     *                |   along this axis - n_iline axis
     *</code>
     */
    double survey_azimuth;
    int has_survey_azimuth;
    /** cached values for azimuth sin and cos */
    double sina, cosa;

    /**
     * Flag to specify a left- or right-handed survey.
     *<code>
     * A normal survey is left handed :
     *
     *     ^ In-Line increase
     *     |
     *     |
     *     |-------> Cross-Line increase
     *
     * The thumb points in the direction of the cross-line numbers (or
     * CDP numbers) increase while the index finger points in the
     * direction of the in-line numbers increase.
     *
     * A right handed survey :
     *
     *                      ^ In-Line
     *                      | increase
     *                      |
     *                      |
     * Cross-Line  <--------|
     *  increase
     *
     * The thumb points in the direction of the cross-line numbers
     * increase while the index finger points in the direction of the
     * in-line numbers increase.
     *
     *</code>
     */
    int left_handed;
    int has_left_handed;

    /**
     * The first in-line from this survey. Note that the
     * in-/cross-lines can have fractional values.
     */
    double first_inline;
    int has_first_inline;

    /**
     * The first cross-line from this survey. Note that the
     * in-/cross-lines can have fractional values.
     */
    double first_crossline;
    int has_first_crossline;

    /**
     * The last in-line from this survey. Note that the
     * in-/cross-lines can have fractional values.
     */
    double last_inline;
    int has_last_inline;

    /**
     * The last cross-line from this survey. Note that the
     * in-/cross-lines can have fractional values.
     */
    double last_crossline;
    int has_last_crossline;

    /**
     * The in-line number corresponding to the survey origin.
     *
     * Usually this will be the same as the first_inline. However, for
     * defining sub-surveys, this will be set some other value.
     */
    double inline_at_origin;
    int has_inline_at_origin;

    /**
     * The cross-line number corresponding to the survey origin.
     *
     * Usually this will be the same as the first_crossline. However,
     * for defining sub-surveys, this will be set some other value.
     */
    double crossline_at_origin;
    int has_crossline_at_origin;

    /**
     * The increment for in-line numbering.
     *
     * The second in-line number will be
     * first_inline + inline_increment and so on.
     *
     * The effective distance between in-lines is
     * inline_increment * inline_spacing.
     */
    double inline_increment;
    int has_inline_increment;

    /**
     * The increment for cross-line numbering.
     *
     * The second cross-line number will be
     * first_crossline + crossline_increment and so on.
     *
     * The effective distance between cross-lines is
     * crossline_increment * crossline_spacing.
     */
    double crossline_increment;
    int has_crossline_increment;

    /**
     * In-line spacing - the distance between two in-line numbers.
     *
     * If the in-line increment is one, then this is the effective
     * distance between two in-lines in this survey.
     */
    double inline_spacing;
    int has_inline_spacing;

    /**
     * Cross-line spacing - the distance between two cross-line
     * numbers.
     *
     * If the cross-line increment is one, then this is the effective
     * distance between two cross-lines in this survey.
     */
    double crossline_spacing;
    int has_crossline_spacing;

    /**
     * Units - the survey is in meters or feet.
     */
    survey_units_t survey_units;
    int has_survey_units;

    /**
     * Measure the angle against the inline direction.
     */
    int inline_azimuth;
    int has_inline_azimuth;
} survey_description_t;

SE_inline void get_survey_description(survey_description_t *sd, int argc, char **argv, char *name)
{
    initargs(argc, argv);


    if (!getpardouble(link_chars(name,".ox_survey"), &sd->ox_survey)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".ox_survey"));
    }
    if (!getpardouble(link_chars(name,".oy_survey"), &sd->oy_survey)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".oy_survey"));
    }
    if (!getpardouble(link_chars(name,".survey_azimuth"), &sd->survey_azimuth)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".survey_azimuth"));
    }
    if (!getparint(link_chars(name,".left_handed"), &sd->left_handed)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".left_handed"));
    }
    if (!getpardouble(link_chars(name,".first_inline"), &sd->first_inline)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".first_inline"));
    }
    if (!getpardouble(link_chars(name,".first_crossline"), &sd->first_crossline)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".first_crossline"));
    }
    if (!getpardouble(link_chars(name,".last_inline"), &sd->last_inline)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".last_inline"));
    }
    if (!getpardouble(link_chars(name,".last_crossline"), &sd->last_crossline)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".last_crossline"));
    }
    if (!getpardouble(link_chars(name,".inline_at_origin"), &sd->inline_at_origin)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".inline_at_origin"));
    }
    if (!getpardouble(link_chars(name,".crossline_at_origin"), &sd->crossline_at_origin)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".crossline_at_origin"));
    }
    if (!getpardouble(link_chars(name,".inline_increment"), &sd->inline_increment)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".inline_increment"));
    }
    if (!getpardouble(link_chars(name,".crossline_increment"), &sd->crossline_increment)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".crossline_increment"));
    }
    if (!getpardouble(link_chars(name,".inline_spacing"), &sd->inline_spacing)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".inline_spacing"));
    }
    if (!getpardouble(link_chars(name,".crossline_spacing"), &sd->crossline_spacing)){
        SE_ERROR("Error reading parameter::%s\n",link_chars(name,".crossline_spacing"));
    }
    
    if (!getparint(link_chars(name,".survey_units"), (int *)&sd->survey_units)){
        sd->survey_units = SURVEY_UNITS_METERS;
    }
    if (!getparint(link_chars(name,".inline_azimuth"), &sd->inline_azimuth)){
        sd->inline_azimuth = 0;
    }

}

SE_inline double
_get_yazimuth(const survey_description_t *sd)
{
    return sd->inline_azimuth ? (sd->survey_azimuth - 90.0) : sd->survey_azimuth;
}

SE_inline void
survey_description_update_cached_values(survey_description_t *sd)
{
    if (sd)
    {
        double a = _get_yazimuth(sd);
        a *= (M_PI / 180.0);
        sd->sina = sin(a);
        sd->cosa = cos(a);
    }
}

SE_inline void init_survey_description(survey_description_t *sd)
{
    if (sd == NULL)
        return;
    /* name will be empty string, all has_* will be false,
       origin will be (0,0) */
    memset(sd, 0, sizeof(*sd));
    sd->left_handed = 1;
    sd->first_inline = 1.0;
    sd->first_crossline = 1.0;
    sd->last_inline = 1.0;
    sd->last_crossline = 1.0;
    sd->inline_at_origin = 1.0;
    sd->crossline_at_origin = 1.0;
    sd->inline_increment = 1.0;
    sd->crossline_increment = 1.0;
    sd->inline_spacing = 1.0;
    sd->crossline_spacing = 1.0;
    sd->survey_units = SURVEY_UNITS_METERS;
    survey_description_update_cached_values(sd);
}

SE_inline const char *survey_get_units_name(survey_units_t survey_units)
{
    switch (survey_units)
    {
    case SURVEY_UNITS_METERS:
        return "meters";
    case SURVEY_UNITS_FEET:
        return "feet";
    default:
        return "unknown";
    }
}

SE_inline const char *survey_get_units_small_name(survey_units_t survey_units)
{
    switch (survey_units)
    {
    case SURVEY_UNITS_METERS:
        return "m";
    case SURVEY_UNITS_FEET:
        return "f";
    default:
        return "?";
    }
}

SE_inline survey_units_t survey_get_units_from_int(int v)
{
    switch (v)
    {
    case SURVEY_UNITS_METERS:
        return SURVEY_UNITS_METERS;
    case SURVEY_UNITS_FEET:
        return SURVEY_UNITS_FEET;
    default:
        return SURVEY_UNITS_UNKNOWN;
    }
}

SE_inline void survey_description_set_name(survey_description_t *sd,
                                           const char *name)
{
    if (name)
    {
        strncpy(sd->name, name, sizeof(sd->name) - 1);
        sd->name[sizeof(sd->name) - 1] = 0;
    }
    else
    {
        sd->name[0] = 0;
    }
}

/**
 * Computes the local coordinates of the point (gx,gy) relative to the
 * survey description defined by the sd argument, and stores the
 * result to lx and ly arguments.
 *
 * A NULL value for sd, means the identity affine transform - that is
 * the global and local coordinates are the same.
 */
SE_inline void to_local_cs(const survey_description_t *sd, // 坐标变换
                           double gx, double gy,
                           double *lx, double *ly)
{
    if (sd == NULL)
    {
        *lx = gx;
        *ly = gy;
    }
    else
    {
        gy -= sd->oy_survey;
        gx -= sd->ox_survey;

        if (sd->left_handed)
        {
            *ly = gy * sd->cosa + gx * sd->sina;
            *lx = gx * sd->cosa - gy * sd->sina;
        }
        else
        {
            *ly = gy * sd->cosa + gx * sd->sina;
            *lx = -gx * sd->cosa + gy * sd->sina;
        }
    }
}

SE_inline void to_local_cs_float(const survey_description_t *sd,
                                 double gx, double gy,
                                 float *lx, float *ly)
{
    double dlx, dly;
    to_local_cs(sd, gx, gy, &dlx, &dly);
    *lx = (float)dlx;
    *ly = (float)dly;
}

SE_inline void to_local_cs_real(const survey_description_t *sd,
                                double gx, double gy,
                                real *lx, real *ly)
{
    double dlx, dly;
    to_local_cs(sd, gx, gy, &dlx, &dly);
    *lx = (real)dlx;
    *ly = (real)dly;
}

SE_inline void to_local_cs_reala(const survey_description_t *sd,
                                 double gx, double gy,
                                 reala *lx, reala *ly)
{
    double dlx, dly;
    to_local_cs(sd, gx, gy, &dlx, &dly);
    *lx = (reala)dlx;
    *ly = (reala)dly;
}

/**
 * Computes the global coordinates of the point (lx,ly) specified
 * relative to the survey description defined by the sd argument, and
 * stores the result to gx and gy arguments.
 *
 * A NULL value for sd, means the identity affine transform - that is
 * the global and local coordinates are the same.
 */
SE_inline void to_global_cs(const survey_description_t *sd,
                            double lx, double ly,
                            double *gx, double *gy)
{
    if (sd == NULL)
    {
        *gx = lx;
        *gy = ly;
    }
    else
    {
        if (sd->left_handed)
        {
            *gy = ly * sd->cosa - lx * sd->sina + sd->oy_survey;
            *gx = lx * sd->cosa + ly * sd->sina + sd->ox_survey;
        }
        else
        {
            *gy = ly * sd->cosa + lx * sd->sina + sd->oy_survey;
            *gx = -lx * sd->cosa + ly * sd->sina + sd->ox_survey;
        }
    }
}

SE_inline void to_global_cs_float(const survey_description_t *sd,
                                  double lx, double ly,
                                  float *gx, float *gy)
{
    double dgx, dgy;
    to_global_cs(sd, lx, ly, &dgx, &dgy);
    *gx = (float)dgx;
    *gy = (float)dgy;
}

SE_inline void to_global_cs_real(const survey_description_t *sd,
                                 double lx, double ly,
                                 real *gx, real *gy)
{
    double dgx, dgy;
    to_global_cs(sd, lx, ly, &dgx, &dgy);
    *gx = (real)dgx;
    *gy = (real)dgy;
}

SE_inline void to_global_cs_reala(const survey_description_t *sd,
                                  double lx, double ly,
                                  reala *gx, reala *gy)
{
    double dgx, dgy;
    to_global_cs(sd, lx, ly, &dgx, &dgy);
    *gx = (reala)dgx;
    *gy = (reala)dgy;
}

/**
 * Returns true (non-zero) if the two surveys described by sd1 and sd2
 * arguments are rotated relative to each other.
 *
 * A NULL value for any of the surveys, means the identity affine
 * transform and a grid of one point, starting from 1 in both
 * direction with spacing of 1.
 */
SE_inline int survey_is_rotated(const survey_description_t *sd1,
                                const survey_description_t *sd2)
{
    double a1 = sd1 != NULL ? sd1->survey_azimuth : 0.0;
    double a2 = sd2 != NULL ? sd2->survey_azimuth : 0.0;
    a1 -= 360 * floor(a1 / 360.0);
    a2 -= 360 * floor(a2 / 360.0);
    return (fabs(a1 - a2) > 1.0E-20) ? 1 : 0;
}

/**
 * Returns true (non-zero) if one of the two surveys and grids
 * described by sd1 and sd2 arguments have all the grid nodes included
 * in the other and the in- and cross-lines are parallel (I guess that
 * mathematically this means that the 2-D rectangluar latice described
 * by one survey is a sub-latice of the other).
 *
 * This means that they are not rotated with respect to each other,
 * their effective distance between in-lines and cross-lines differ by
 * at most an integer multiplier, and the translation vector between
 * the two points to a grid node from one of the two surveys.
 *
 * A NULL value for any of the surveys, means the identity affine
 * transform and a grid of one point, starting from 1 in both
 * direction with spacing of 1.
 */
// int survey_is_subgrid(const survey_description_t *sd1,
//                       const survey_description_t *sd2)
// {
//     (void)sd1;
//     (void)sd2;
//     SE_ERROR("survey_is_subgrid not yet implemented");
//     return 0;
// }

/**
 * Returns true (non-zero) if the two surveys and grids described by
 * sd1 and sd2 arguments have the same nodes.
 *
 * This means that they are not rotated with respect to each other,
 * their effective distance between in-lines and cross-lines are the
 * same, and the translation vector between the two points to a grid
 * node.
 *
 * A NULL value for any of the surveys, means the identity affine
 * transform and a grid of one point, starting from 1 in both
 * direction with spacing of 1.
 */
// int survey_is_equivalent(const survey_description_t *sd1,
//                          const survey_description_t *sd2)
// {
//     (void)sd1;
//     (void)sd2;
//     SE_ERROR("survey_is_equivalent not yet implemented");
//     return 0;
// }

/**
 * Computes the (o,d,n) sets (origin, distance, number of samples) for
 * the two spatial directions based on the survey defined by the
 * the sd argument.
 *
 * The ox,dx,nx refers to the X direction in the local coordinate
 * system, corresponding to the direction along an in-line, or the
 * direction the cross-line numbers increase.
 *
 * The oy,dy,ny refers to the Y direction in the local coordinate
 * system, corresponding to the direction along a cross-line, or the
 * direction the in-line numbers increase.
 */
SE_inline void survey_computation_pars(const survey_description_t *sd,
                                       real *ox, real *dx, int *nx,
                                       real *oy, real *dy, int *ny)
{
    if (sd == NULL)
    {
        *ox = *oy = (real)0.0;
        *dx = *dy = (real)1.0;
        *nx = *ny = 1;
    }
    else
    {
        *dy = (real)(sd->inline_spacing * sd->inline_increment);
        *dx = (real)(sd->crossline_spacing * sd->crossline_increment);

        *oy = (real)((sd->first_inline - sd->inline_at_origin) *
                     sd->inline_spacing);
        *ox = (real)((sd->first_crossline - sd->crossline_at_origin) *
                     sd->crossline_spacing);

        *ny = (int)ceil((sd->last_inline - sd->first_inline) /
                        sd->inline_increment) +
              1;
        *nx = (int)ceil((sd->last_crossline - sd->first_crossline) /
                        sd->crossline_increment) +
              1;
    }
}

SE_inline void survey_to_axa(const survey_description_t *sd,
                             axa_t *x, axa_t *y)
{
    int nx, ny;
    survey_computation_pars(sd,
                            &(x->o), &(x->d), &nx,
                            &(y->o), &(y->d), &ny);
    x->n = nx;
    y->n = ny;
}

/**
 * Adjusts the first, and last and increment for the inlines and
 * crosslines to match the computational parameters. All other
 * parameters are not modified.
 */
SE_inline void survey_adjust_from_computation_pars(real ox, real dx, int nx,
                                                   real oy, real dy, int ny,
                                                   survey_description_t *sd)
{
    if (sd)
    {

        sd->inline_increment = dy / sd->inline_spacing;
        sd->first_inline = oy / sd->inline_spacing + sd->inline_at_origin;
        sd->last_inline = sd->first_inline + (ny - 1) * sd->inline_increment;

        sd->crossline_increment = dx / sd->crossline_spacing;
        sd->first_crossline = ox / sd->crossline_spacing + sd->crossline_at_origin;
        sd->last_crossline = sd->first_crossline + (nx - 1) * sd->crossline_increment;
    }
}

SE_inline void survey_get_range(const survey_description_t *sd,
                                real *fx, real *lx,
                                real *fy, real *ly)
{
    if (sd == NULL)
    {
        *fx = *fy = *lx = *ly = (real)0.0;
    }
    else
    {
        *fy = (real)((sd->first_inline - sd->inline_at_origin) * sd->inline_spacing);
        *ly = (real)((sd->last_inline - sd->inline_at_origin) * sd->inline_spacing);
        *fx = (real)((sd->first_crossline - sd->crossline_at_origin) * sd->crossline_spacing);
        *lx = (real)((sd->last_crossline - sd->crossline_at_origin) * sd->crossline_spacing);
    }
}

SE_inline void survey_get_minmax(const survey_description_t *sd,
                                 real *minx, real *maxx,
                                 real *miny, real *maxy)
{
    survey_get_range(sd, minx, maxx, miny, maxy);
    if (*miny > *maxy)
    {
        real tmp = *miny;
        *miny = *maxy;
        *maxy = tmp;
    }

    if (*minx > *maxx)
    {
        real tmp = *minx;
        *minx = *maxx;
        *maxx = tmp;
    }
}

SE_inline void
survey_computation_pars_float(const survey_description_t *sd,
                              float *_ox, float *_dx, int *nx,
                              float *_oy, float *_dy, int *ny)
{
    real ox, dx, oy, dy;
    survey_computation_pars(sd, &ox, &dx, nx, &oy, &dy, ny);
    *_ox = (float)ox;
    *_dx = (float)dx;
    *_oy = (float)oy;
    *_dy = (float)dy;
}

SE_inline void survey_line_numbers_from_global(const survey_description_t *sd,
                                               double gx, double gy,
                                               double *xl, double *il)
{
    if (sd == NULL)
    {
        *xl = gx;
        *il = gy;
    }
    else
    {
        double lx, ly;
        to_local_cs(sd, gx, gy, &lx, &ly);

        *il = sd->inline_at_origin + ly / sd->inline_spacing;
        *xl = sd->crossline_at_origin + lx / sd->crossline_spacing;
    }
}

SE_inline void survey_line_numbers_to_global(const survey_description_t *sd,
                                             double xl, double il,
                                             double *gx, double *gy)
{
    if (sd == NULL)
    {
        *gx = xl;
        *gy = il;
    }
    else
    {
        double ly = sd->inline_spacing * (il - sd->inline_at_origin);
        double lx = sd->crossline_spacing * (xl - sd->crossline_at_origin);

        to_global_cs(sd, lx, ly, gx, gy);
    }
}

SE_inline void survey_line_numbers_from_local(const survey_description_t *sd,
                                              double lx, double ly,
                                              double *xl, double *il)
{
    if (sd == NULL)
    {
        *xl = lx;
        *il = ly;
    }
    else
    {
        *il = sd->inline_at_origin + ly / sd->inline_spacing;
        *xl = sd->crossline_at_origin + lx / sd->crossline_spacing;
    }
}

SE_inline void survey_line_numbers_to_local(const survey_description_t *sd,
                                            double xl, double il,
                                            double *lx, double *ly)
{
    if (sd == NULL)
    {
        *lx = xl;
        *ly = il;
    }
    else
    {
        *ly = sd->inline_spacing * (il - sd->inline_at_origin);
        *lx = sd->crossline_spacing * (xl - sd->crossline_at_origin);
    }
}

/**
 * If the effective verbosity level is bigger or equal than verb, it
 * prints out (to the current log device) the current values for the
 * sd survey.
 */
void report_survey_description(const survey_description_t *sd,
                               int verb);

void print_survey_description(FILE *f, const survey_description_t *sd);

/**
 * Returns true (non-zero) if the two surveys belong to the same project,
 * of false (zero) otherwise.
 */
SE_inline int same_project_surveys(const survey_description_t *sd1,
                                   const survey_description_t *sd2)
{
    if (sd1 == NULL)
    {
        return sd2 == NULL;
    }
    else if (sd1 == sd2)
    {
        return 1;
    }
    else
    {
        const int l1 = (sd1->left_handed != 0) ? 1 : 0;
        const int l2 = (sd2->left_handed != 0) ? 1 : 0;

        if (l1 != l2)
            return 0;
        if (sd1->survey_units != sd2->survey_units)
            return 0;

        if (!dequal(sd1->ox_survey, sd2->ox_survey))
            return 0;
        if (!dequal(sd1->oy_survey, sd2->oy_survey))
            return 0;
        if (!dequal(sd1->crossline_at_origin, sd2->crossline_at_origin))
            return 0;
        if (!dequal(sd1->inline_at_origin, sd2->inline_at_origin))
            return 0;
        if (!dequal(sd1->crossline_spacing, sd2->crossline_spacing))
            return 0;
        if (!dequal(sd1->inline_spacing, sd2->inline_spacing))
            return 0;

        return is_integer((_get_yazimuth(sd1) - _get_yazimuth(sd2)) / 360.0);
    }
}

/**
 * Returns true (non-zero) if the destination survey description was changed,
 * of false (zero) otherwise.
 */
SE_inline int merge_surveys(survey_description_t *dst,
                            const survey_description_t *src)
{
    if (dst == NULL)
        return 0;
    if (!same_project_surveys(dst, src))
    {
        SE_WARNING("Trying to merge information from incompatible surveys\n");
        return 0;
    }
    else
    {
        int ret = 0;
        real sfx, sfy, slx, sly;
        real dfx, dfy, dlx, dly;
        real sminx, smaxx, sminy, smaxy;

        if (fabs(dst->crossline_increment) > fabs(src->crossline_increment))
        {
            double factor = ceil(fabs(dst->crossline_increment) / fabs(src->crossline_increment));
            dst->crossline_increment /= factor;
            ret = 1;
        }
        if (fabs(dst->inline_increment) > fabs(src->inline_increment))
        {
            double factor = ceil(fabs(dst->inline_increment) / fabs(src->inline_increment));
            dst->inline_increment /= factor;
            ret = 1;
        }

        survey_get_minmax(src, &sminx, &smaxx, &sminy, &smaxy);
        survey_get_range(src, &sfx, &slx, &sfy, &sly);
        survey_get_range(dst, &dfx, &dlx, &dfy, &dly);

        if (dfx < dlx)
        {
            if (dfx > sminx)
            {
                dst->first_crossline = dst->crossline_at_origin + sminx / dst->crossline_spacing;
                ret = 1;
            }
            if (dlx < smaxx)
            {
                dst->last_crossline = dst->crossline_at_origin + smaxx / dst->crossline_spacing;
                ret = 1;
            }
        }
        else
        {
            if (dlx > sminx)
            {
                dst->last_crossline = dst->crossline_at_origin + sminx / dst->crossline_spacing;
                ret = 1;
            }
            if (dfx < smaxx)
            {
                dst->first_crossline = dst->crossline_at_origin + smaxx / dst->crossline_spacing;
                ret = 1;
            }
        }

        if (dfy < dly)
        {
            if (dfy > sminy)
            {
                dst->first_inline = dst->inline_at_origin + sminy / dst->inline_spacing;
                ret = 1;
            }
            if (dly < smaxy)
            {
                dst->last_inline = dst->inline_at_origin + smaxy / dst->inline_spacing;
                ret = 1;
            }
        }
        else
        {
            if (dly > sminy)
            {
                dst->last_inline = dst->inline_at_origin + sminy / dst->inline_spacing;
                ret = 1;
            }
            if (dfy < smaxy)
            {
                dst->first_inline = dst->inline_at_origin + smaxy / dst->inline_spacing;
                ret = 1;
            }
        }

        return ret;
    }
}

#endif
