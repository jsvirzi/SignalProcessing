#include <math.h>
#include <malloc.h>

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

double *binomial_mult( int n, double *p )
{
    int i, j;
    double *a;

    a = (double *)calloc( 2 * n, sizeof(double) );
    if( a == NULL ) return( NULL );

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

double *dcof_bwlp( int n, double fcf )
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

    rcof = (double *)calloc( 2 * n, sizeof(double) );
    if( rcof == NULL ) return( NULL );

    theta = M_PI * fcf;
    st = sin(theta);
    ct = cos(theta);

    for( k = 0; k < n; ++k )
    {
        parg = M_PI * (double)(2*k+1)/(double)(2*n);
        sparg = sin(parg);
        cparg = cos(parg);
        a = 1.0 + st*sparg;
        rcof[2*k] = -ct/a;
        rcof[2*k+1] = -st*cparg/a;
    }

    dcof = binomial_mult( n, rcof );
    free( rcof );

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for( k = 3; k <= n; ++k )
        dcof[k] = dcof[2*k-2];
    return( dcof );
}

/**********************************************************************
  dcof_bwhp - calculates the d coefficients for a butterworth highpass
  filter. The coefficients are returned as an array of doubles.

*/

double *dcof_bwhp( int n, double fcf )
{
    return( dcof_bwlp( n, fcf ) );
}
/**********************************************************************
  ccof_bwlp - calculates the c coefficients for a butterworth lowpass
  filter. The coefficients are returned as an array of integers.

*/

int *ccof_bwlp( int n )
{
    int *ccof;
    int m;
    int i;

    ccof = (int *)calloc( n+1, sizeof(int) );
    if( ccof == NULL ) return( NULL );

    ccof[0] = 1;
    ccof[1] = n;
    m = n/2;
    for( i=2; i <= m; ++i)
    {
        ccof[i] = (n-i+1)*ccof[i-1]/i;
        ccof[n-i]= ccof[i];
    }
    ccof[n-1] = n;
    ccof[n] = 1;

    return( ccof );
}

/**********************************************************************
  ccof_bwhp - calculates the c coefficients for a butterworth highpass
  filter. The coefficients are returned as an array of integers.

*/

int *ccof_bwhp( int n )
{
    int *ccof;
    int i;

    ccof = ccof_bwlp( n );
    if( ccof == NULL ) return( NULL );

    for( i = 0; i <= n; ++i)
        if( i % 2 ) ccof[i] = -ccof[i];

    return( ccof );
}

int init(FilterInfo *filterInfo, double fs, double f1, double f2, int type, int order) {

    double a = tan( M_PI * f1 / fs);
    double a2 = a * a;
    double r;
    double *A = (double *)malloc(order * sizeof(double));
    double *d1 = (double *)malloc(order * sizeof(double));
    double *d2 = (double *)malloc(order * sizeof(double));

    for(int i = 0; i < order; ++i){
        r = sin(M_PI*(2.0*i+1.0)/(4.0*n));
        s = a2 + 2.0*a*r + 1.0;
        A[i] = a2/s;
        d1[i] = 2.0*(1-a2)/s;
        d2[i] = -(a2 - 2.0*a*r + 1.0)/s;}

}

double filter(struct FilterInfo *filterInfo, double *y) {
}