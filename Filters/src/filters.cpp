#include <math.h>
#include <malloc.h>
#include <stdint.h>
#include <string.h>

#include "filters.h"

/**********************************************************************
  binomial_mult - multiplies a series of binomials together and returns
  the coefficients of the resulting polynomial.

  The multiplication has the following form:

  (x+p[0])*(x+p[1])*...*(x+p[n-1])

  The p[i] coefficients are assumed to be complex and are passed to the
  function as a pointer to an array of doubles of length 2n.

  The resulting polynomial has the following form:

  x^n + a[0]*x^n-1 + a[1]*x^n-2 + ... +a[n-2]*x + a[n-1]

  The a[i] coefficients can in general be complex but should in most
  cases turn out to be real. The a[i] coefficients are returned by the
  function as a pointer to an array of doubles of length 2n. Storage
  for the array is allocated by the function and should be freed by the
  calling program when no longer needed.

  Function arguments:

  n  -  The number of binomials to multiply
  p  -  Pointer to an array of doubles where p[2i] (i=0...n-1) is
        assumed to be the real part of the coefficient of the ith binomial
        and p[2i+1] is assumed to be the imaginary part. The overall size
        of the array is then 2n.
*/

static double *binomial_mult( int n, double *p )
{
    int i, j;
    double *a;

    a = (double *) calloc(2 * n, sizeof(double));
    if (a == NULL) return( NULL );

    for( i = 0; i < n; ++i )
    {
        for( j = i; j > 0; --j )
        {
            a[2*j] += p[2*i] * a[2*(j-1)] - p[2*i+1] * a[2*(j-1)+1];
            a[2*j+1] += p[2*i] * a[2*(j-1)+1] + p[2*i+1] * a[2*(j-1)];
        }
        a[0] += p[2*i];
        a[1] += p[2*i+1];
    }
    return( a );
}

/**********************************************************************
  dcof_bwlp - calculates the d coefficients for a butterworth lowpass
  filter. The coefficients are returned as an array of doubles.
*/

static double *dcof_bwlp( int n, double fcf )
{
    int k;            // loop variables
    double theta;     // M_PI * fcf / 2.0
    double st;        // sine of theta
    double ct;        // cosine of theta
    double parg;      // pole angle
    double sparg;     // sine of the pole angle
    double cparg;     // cosine of the pole angle
    double a;         // workspace variable
    double *rcof;     // binomial coefficients
    double *dcof;     // dk coefficients

    rcof = (double *) calloc(2 * n, sizeof(double));
    if (rcof == NULL) return (NULL);

    theta = M_PI * fcf;
    st = sin(theta);
    ct = cos(theta);

    for( k = 0; k < n; ++k )
    {
        parg = M_PI * (double)(2*k+1)/(double)(2*n);
        sparg = sin(parg);
        cparg = cos(parg);
        a = 1.0 + st * sparg;
        rcof[2*k] = -ct / a;
        rcof[2*k+1] = -st * cparg / a;
    }

    dcof = binomial_mult(n, rcof);
    free(rcof);

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for(k = 3; k <= n; ++k) {
        dcof[k] = dcof[2 * k - 2];
    }
    return dcof;
}

/**********************************************************************
  dcof_bwhp - calculates the d coefficients for a butterworth highpass
  filter. The coefficients are returned as an array of doubles.

*/

static double *dcof_bwhp( int n, double fcf )
{
    return dcof_bwlp(n, fcf);
}
/**********************************************************************
  ccof_bwlp - calculates the c coefficients for a butterworth lowpass
  filter. The coefficients are returned as an array of integers.

*/

static double *ccof_bwlp(int n)
{
    double *ccof;
    int m;
    int i;

    ccof = (double *) calloc(n+1, sizeof(double));
    if (ccof == NULL) return (NULL);

    ccof[0] = 1;
    ccof[1] = n;
    m = n / 2;
    for(i = 2; i <= m; ++i)
    {
        ccof[i] = (n - i + 1) * ccof[i-1] / i;
        ccof[n-i] = ccof[i];
    }
    ccof[n-1] = n;
    ccof[n] = 1;

    return ccof;
}

/**********************************************************************
  ccof_bwhp - calculates the c coefficients for a butterworth highpass
  filter. The coefficients are returned as an array of integers.

*/

static double *ccof_bwhp(int n)
{
    double *ccof;
    int i;

    ccof = ccof_bwlp( n );
    if (ccof == NULL) return(NULL);

    for(i = 0; i <= n; ++i)
        if (i % 2) ccof[i] = -ccof[i];

    return(ccof);
}

static uint32_t upper_power_of_two(uint32_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

static double filter_sample(struct FilterInfo *filterInfo, const double x) {
    double accC = x, accD = 0.0;
    uint32_t idx = filterInfo->queue_head;
    uint32_t idx0 = idx;
    filterInfo->x[idx0] = x;
    uint32_t mask = filterInfo->queue_mask;
    filterInfo->queue_head = (filterInfo->queue_head + 1) & mask;
    double *c = filterInfo->c, *d = filterInfo->d, *px = filterInfo->x, *py = filterInfo->y;
    for (unsigned int i = 0; i < filterInfo->nc; ++i) {
        accC = accC + c[i] * px[idx];
        accD = accD + d[i] * py[idx];
        idx = (idx - 1) & mask;
    }
    double y = accD - accC;
    filterInfo->y[idx0] = y;
    return y;
}

static void filter_batch(FilterInfo *filterInfo, const double *x, const int n, double *y) {
    for (int i = 0; i < n; ++i) {
        y[i] = filter_sample(filterInfo, x[i]);
    }
}

static double filter_amplitude(FilterInfo *filterInfo, const double f, const unsigned int runTime) {
    /* use of filter is invasive. make individual copies for each real and imaginary component */
    FilterInfo filterInfoR, filterInfoI;
    filterInfoR = *filterInfo;
    filterInfoI = *filterInfo;

    double dt = 1.0 / filterInfo->fs;
    double phase0 = 0.25 * M_PI;
    double aR = cos(phase0), aI = sin(phase0);
    double accR = aR * aR, accI = aI * aI;
    unsigned int i, startUpTime = runTime / 10;
    for (i = 1; i < startUpTime; ++i) {
        const double theta = i * dt;
        const double c = cos(theta + phase0);
        const double s = sin(theta + phase0);
        aR = filterInfoR.sample(&filterInfoR, c);
        aI = filterInfoI.sample(&filterInfoI, s);
    }
    for (; i < runTime; ++i) {
        const double theta = i * dt;
        const double c = cos(theta + phase0);
        const double s = sin(theta + phase0);
        aR = filterInfoR.sample(&filterInfoR, c);
        aI = filterInfoI.sample(&filterInfoI, s);
        accR += (aR * aR);
        accI += (aI * aI);
    }
    double phase = atan2(aI, aR);
    return (phase - phase0);
}

static double filter_phase(FilterInfo *filterInfo, const double f, unsigned int runTime) {
    /* use of filter is invasive. make individual copies for each real and imaginary component */
    FilterInfo filterInfoR, filterInfoI;
    filterInfoR = *filterInfo;
    filterInfoI = *filterInfo;

    double dt = 1.0 / filterInfo->fs;
    double phase0 = 0.25 * M_PI;
    double omega = M_PI * f;
    double aR = cos(phase0), aI = sin(phase0);
    for (unsigned int i = 1; i < runTime; ++i) {
        const double theta = omega * i * dt;
        const double c = cos(theta + phase0);
        const double s = sin(theta + phase0);
        aR = filterInfoR.sample(&filterInfoR, c);
        aI = filterInfoI.sample(&filterInfoI, s);
    }
    double phase = atan2(aI, aR);
    return (phase - phase0);
}

static int filter_close(FilterInfo *filterInfo) {
    free (filterInfo->x);
    free (filterInfo->y);
    free (filterInfo->c);
    free (filterInfo->d);
}

int filter_init(FilterInfo *filterInfo) {
    double r_cutoff = M_PI * filterInfo->f1 / filterInfo->fs;
    switch (filterInfo->options & FilterTypeMask) {
    case FilterTypeLowPass:
        filterInfo->c = ccof_bwlp(filterInfo->order);
        filterInfo->d = dcof_bwlp(filterInfo->order, r_cutoff);
        filterInfo->ni = 0;
        filterInfo->queue_size = upper_power_of_two(filterInfo->order + 1) * 2;
        filterInfo->x = (double *) calloc(filterInfo->queue_size, sizeof(double));
        filterInfo->y = (double *) calloc(filterInfo->queue_size, sizeof(double));
        filterInfo->queue_mask = filterInfo->queue_size - 1;
        filterInfo->queue_head = 0;
        filterInfo->batch = filter_batch;
        filterInfo->sample = filter_sample;
        filterInfo->amplitude = filter_amplitude;
        filterInfo->phase = filter_phase;
        filterInfo->close = filter_close;
        break;
    default:
        break;
    }
    return 0;
}
