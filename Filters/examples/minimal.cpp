#include <math.h>
#include <stdio.h>

#include "filters.h"

/* pretty much the absolute minimum skeleton to run a filter */

int main(int argc, char **argv) {
    FilterInfo filterInfo;
    filterInfo.options = FilterTypeLowPass | FilterVariantButterworth | FilterImplementationIIR;
    filterInfo.order = 2;
    filterInfo.fs = 1000.0;
    filterInfo.f1 = 100.0;
    filter_init(&filterInfo);

    double x;
    while (scanf("%lf", &x)) {
        double y = filterInfo.sample(&filterInfo, x); /* calculate next value */
        printf("%lf %lf\n", x, y);
    }
}
