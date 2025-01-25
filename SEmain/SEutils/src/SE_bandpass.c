#include "../include/SE_bandpass.h"


SE_private void _bandpass_init_bp_lo(bandpass_t* bp);
SE_private void _bandpass_init_bp_hi(bandpass_t* bp);
SE_private void _bandpass_do_cut( bandpass_t* bp, reala* data, int n, double* coef, int ncoef );

bandpass_t* bandpass_init( const double dt,
                                   const double flo, const double fhi,
                                   const int nplo, const int nphi, const int phase )
{
    bandpass_t* bp;
    bp = (bandpass_t*)malloc(sizeof(bandpass_t));
    memset(bp, 0, sizeof(bandpass_t));


    bp->phase  = phase;
    bp->dt = dt;

    bp->orig_nphi = nphi;
    bp->orig_nplo = nplo;
    bp->orig_flo = flo;
    bp->orig_fhi = fhi;

    bp->buff_len = 0;
    bp->buff = NULL;

    bp->nphi   = nphi;
    bp->noddhi = nphi % 2;
    if( ! bp->phase ) {
        if( bp->noddhi != 0) bp->nphi = ( nphi + 1 ) / 2;
        else                 bp->nphi = bp->nphi / 2;
        bp->noddhi = bp->nphi%2;
    }
    bp->nb = (bp->nphi + 1) / 2;
    bp->b = alloc1double(5 * bp->nb);

    bp->nplo   = nplo;
    bp->noddlo = nplo % 2;
    if( ! bp->phase ) {
        if( bp->noddlo != 0) bp->nplo = ( nplo + 1 ) / 2;
        else                 bp->nplo = nplo / 2;
        bp->noddlo = bp->nplo%2;
    }
    bp->nd = (bp->nplo + 1) / 2;
    bp->d = alloc1double(5 * bp->nd);
    memset(bp->d, 0, 5 * bp->nd * sizeof(double));


    if( flo <= 0 ) bp->flo = 0.0;
    else           bp->flo = flo;

    if( fhi <= 0 ) bp->fhi = 0.5 / dt;
    else           bp->fhi = fhi;

    bp->flo *= dt;
    bp->fhi *= dt;
    if( bp->flo > 0.000001 ) _bandpass_init_bp_lo(bp);
    if( bp->fhi < 4.999999 ) _bandpass_init_bp_hi(bp);

    return bp;
}

bandpass_t* bandpass_clone(const bandpass_t* src) 
{
    return bandpass_init(src->dt, src->orig_flo, src->orig_fhi, src->orig_nplo, src->orig_nphi, src->phase);
}

void bandpass_destroy( bandpass_t* bp )
{
    if(bp) {
        free(bp->b);
        free(bp->d);
        if(bp->buff_len > 0)
            free(bp->buff);
        free(bp);
    }
}

void apply_bandpass( bandpass_t* bp, reala* data, int nsamples, int ntraces )
{
    /* if we have to apply any filter, cleanup the data */
    if( bp->flo > 0.000001 || bp->fhi < 4.999999 ) {
        remove_nans(data, nsamples*ntraces, 1, 0.0f);
    }

    if( bp->flo > 0.000001 ) {
        int i;
        for(i = 0; i < ntraces; ++i) {
            _bandpass_do_cut( bp, data + i*nsamples, nsamples, bp->d, bp->nd );
        }
    }
    if( bp->fhi < 4.999999 ) {
        int i;
        for(i = 0; i < ntraces; ++i) {
            _bandpass_do_cut( bp, data + i*nsamples, nsamples, bp->b, bp->nb );
        }
    }
}

void bandpass_get_pars( const bandpass_t* bp, double* dt,
                            double* flo, double* fhi,
                            int* nplo, int* nphi, int* phase)
{
    *dt = bp->dt;
    *phase = bp->phase;

    *flo = bp->orig_flo;
    *fhi = bp->orig_fhi;
    *nplo = bp->orig_nplo;
    *nphi = bp->orig_nphi;
}

SE_private void _bandpass_apply_internal( double* src, double* dst, const int n,
                                               const double* coef, const int ncoef )
{
    int k;
    for( k = 0; k < ncoef; k++ ) {
        const double c0 = coef[k + 0*ncoef];
        const double c1 = coef[k + 1*ncoef];
        const double c2 = coef[k + 2*ncoef];
        const double c3 = coef[k + 3*ncoef];
        const double c4 = coef[k + 4*ncoef];
        int i;
        for( i = 2; i < n; i++ ) {
            dst[i] =
                c0 * src[i]   +
                c1 * src[i-1] +
                c2 * src[i-2] -
                c3 * dst[i-1] -
                c4 * dst[i-2];
        }

        for( i = 0; i < n; i++ ) src[i] = dst[i];
    }
}

SE_private void _bandpass_do_cut( bandpass_t* bp, reala* data, int n, double* coef, int ncoef )
{
    double* newdata;
    double* tmpdata;
    int i, npad = n + 2;

    if( npad > bp->buff_len ) {
        bp->buff = alloc1double(2*npad);
        bp->buff_len = npad;
    }

    newdata = bp->buff;
    tmpdata = bp->buff + npad;

    tmpdata[0] = tmpdata[1] = 0.0;
    for(i = 0; i < n; i++ )
        tmpdata[i+2] = (double)data[i];

    newdata[0] = newdata[1] = 0.0;
    _bandpass_apply_internal( tmpdata, newdata, npad, coef, ncoef );

    if( ! bp->phase ) { /* again in reverse */
        newdata[0] = newdata[1] = 0.0;
        for( i = 2; i < npad; i++ )
            newdata[i] = tmpdata[npad+1-i];

        _bandpass_apply_internal( newdata, tmpdata, npad, coef, ncoef );

        for( i = 0; i < n; i++ )
            data[i] = (float)newdata[npad-1-i];
    }
    else {
        for( i = 0; i < n; i++ )
            data[i] = (float)tmpdata[i+2];
    }
}


SE_private void _bandpass_init_bp_lo(bandpass_t* bp)
{
    double a,aa,aap4,dtheta,theta0,e,ee,fno2;
    int j;

    fno2 = 0.25;   /*   Nyquist frequency over two?? */
    a = 2.0 * sin(M_PI * fno2) / cos(M_PI * fno2);
    aa = a*a;
    aap4 = aa + 4;

    e = -cos( M_PI * (bp->flo + fno2) ) / cos( M_PI * (bp->flo - fno2) );
    ee = e*e;

    dtheta = M_PI / bp->nplo;	          /*  angular separation of poles  */
    if( bp->noddlo != 0 ) theta0 = 0;     /*  pole closest to real s axis */
    else                  theta0 = dtheta / 2.0;

    if( bp->noddlo != 0 ) {
        double b1 = a / (a + 2);
        double b2 = (a - 2) / (a + 2);
        double den = 1.0 - b2*e;
        bp->d[0+0*bp->nd] = b1 * (1.0 - e) / den;
        bp->d[0+1*bp->nd] = -bp->d[0];
        bp->d[0+2*bp->nd] = 0.0;
        bp->d[0+3*bp->nd] = (e - b2) / den;
        bp->d[0+4*bp->nd] = 0.0;
    }

    for(j = bp->noddlo; j < bp->nd; ++j) {
        double c = 4.0 * a * cos(theta0 + j*dtheta);
        double b1 = aa / (aap4 + c);
        double b2 = (2.0*aa - 8.0) / (aap4 + c);
        double b3 = (aap4 - c) / (aap4 + c);
        double den = 1.0 - b2*e + b3*ee;
        bp->d[j+0*bp->nd] = b1 * (1.0-e) * (1.0-e) / den;
        bp->d[j+1*bp->nd] = -2.0 * bp->d[j];
        bp->d[j+2*bp->nd] = bp->d[j];
        bp->d[j+3*bp->nd] = ( 2.0 * e * (1.0 + b3) - b2 * (1.0 + ee) ) / den;
        bp->d[j+4*bp->nd] = (ee - b2*e + b3) / den;
    }
}

SE_private void _bandpass_init_bp_hi(bandpass_t* bp)
{
    double a,aa,aap4,dtheta,theta0;
    int j;

    a = 2.0 * tan(M_PI * bp->fhi);  /* radius of poles in s-plane */
    aa = a*a;
    aap4 = aa + 4.0;

    dtheta = M_PI / bp->nphi;            /*  angular separation of poles  */
    if( bp->noddhi != 0 ) theta0 = 0.0;  /*  pole closest to real s axis */
    else                  theta0 = dtheta / 2.0;

    if( bp->noddhi != 0 ) {
        bp->b[0*bp->nb] = a / (a+2.0);        /* 0,0 */
        bp->b[1*bp->nb] = bp->b[0];           /* 0,1 */
        bp->b[2*bp->nb] = 0.0;                /* 0,2 */
        bp->b[3*bp->nb] = (a-2.0) / (a+2.0);  /* 0,3 */
        bp->b[4*bp->nb] = 0.0;                /* 0,4 */
    }

    for(j = bp->noddhi; j < bp->nb; ++j) {
        double c = 4.0 * a * cos(theta0 + j*dtheta);
        bp->b[j+0*bp->nb] = aa / (aap4 + c);                  /* j,0 */
        bp->b[j+1*bp->nb] = 2.0 * bp->b[j];                   /* j,1 */
        bp->b[j+2*bp->nb] = bp->b[j];                         /* j,2 */
        bp->b[j+3*bp->nb] = (2.0 * aa - 8.0) / (aap4 + c);    /* j,3 */
        bp->b[j+4*bp->nb] = (aap4-c) / (aap4+c);              /* j,4 */
    }
}
