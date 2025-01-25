#include "../include/SEmath.h"
#include "../include/SEmaths.h"


real rand_normal()
{
    static unsigned char have_next = 0;
    static real next;
    real v1, v2, s, tmp;

    if ( have_next ) {
        have_next = 0;
        return next;
    }
    
    do { 
        v1 = 2 * se_rand1() - 1;
        v2 = 2 * se_rand1() - 1;
        s = v1*v1 + v2*v2;
    } while ( s >= 1 );

    if( ! FISZERO(s) ) {
        tmp = sqrt( -2 * log(s) / s );
        next = v2 * tmp;
        have_next = 1;
        return v1 * tmp;
    }
    else {
        return 0;
    }
}


reala percentile(real p, reala* samples, size_t ns)
{
    ssize_t k, lo, hi;

    k = (ssize_t) (ns * p / 100.0);
    if(k < 0)   k = 0;
    if(k >= ns) k = ns - 1;

    for (lo = 0, hi = ns - 1; lo < hi; ) {
        reala xk = samples[k];
        ssize_t i = lo;
        ssize_t j = hi;
        do {
            while (samples[i] < xk) {
                i++;
            }
            while (samples[j] > xk) {
                j--;
            }
            if (i <= j) {
                reala tmp = samples[i];
                samples[i] = samples[j];
                samples[j] = tmp;
                i++;
                j--;
            } else { // spare a comparison
                break;
            }
        } while (i <= j);

        if (j < k) {
            lo = i;
        }
        if (k < i) {
            hi = j;
        }
    }

    return samples[k];
}
