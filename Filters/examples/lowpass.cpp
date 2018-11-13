#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "filters.h"

/*
 * this serves as an example of usage of filters library
 *
 */

int main(int argc, char **argv) {

    double cornerFrequency = 100.0;
    double samplingFrequency = 1000.0;
    int filterOrder = 2;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-fc") == 0) {
            cornerFrequency = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "-fs") == 0) {
            samplingFrequency = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "-order") == 0) {
            filterOrder = atoi(argv[++i]);
        }
    }

    /* TODO how large can arguments (to sin() and cos()) get before bad things happen? */
    const double nyquistFrequency = samplingFrequency * 0.5;
    const double stopCycles = 100.0; /* run for 100 cycles */
    const double steadyStateCycles = 10.0; /* how many cycles to achieve steady state */
    const double stopAngle = (stopCycles + steadyStateCycles) * 2.0 * M_PI;
    const double steadyStateAngle = steadyStateCycles * 2.0 * M_PI; /* let filter run for 10 cycles before analysis */

    /* initialize filter(s). running two instances simultaneously */
    FilterInfo filterInfo;
    memset(&filterInfo, 0, sizeof(FilterInfo));
    filterInfo.options = FilterTypeLowPass | FilterVariantButterworth | FilterImplementationIIR;
    filterInfo.order = filterOrder;
    filterInfo.fs = samplingFrequency;
    filterInfo.f1 = cornerFrequency;
    filter_init(&filterInfo);
    filterInfo.print(&filterInfo);

    FILE *fp = fopen("filter_response.txt", "w");

    /* generate input samples, at different frequencies, and run filter on them */
    for (double frequency = 1.0; frequency < nyquistFrequency; frequency += 1.0) {
        double accSR = 0.0, accSI = 0.0, accFR = 0.0, accFI = 0.0;
        double acc2SR = 0.0, acc2SI = 0.0, acc2FR = 0.0, acc2FI = 0.0;
        double acc2R = 0.0, acc2I = 0.0;
        double dp = 2.0 * M_PI * frequency / samplingFrequency;
        for (double arg = 0.0; arg < stopAngle; arg += dp) {
            double sr = cos(arg), si = sin(arg);
            double fr = filterInfo.sample(&filterInfo, sr);
            if (arg >= steadyStateAngle) {
                accSR = accSR + sr;
                acc2SR = acc2SR + sr * sr;
                accFR = accFR + fr;
                acc2FR = acc2FR + fr * fr;
                acc2R = acc2R + fr * sr;
                acc2I = acc2I + fr * si;
            }
        }

        double amp = sqrt(acc2FR / acc2SR);
        double phase = atan2(acc2I, acc2R) / M_PI;
        printf("frequency=%8.f response: amplitude=%10.7f/phase=%10.7f\n", frequency, amp, phase);
        fprintf(fp, "%f,%f,%f\n", frequency, amp, phase);
        // getchar();
    }
}
