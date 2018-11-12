#include <string.h>
#include <stdio.h>
#include <math.h>

#include "filters.h"

/*
 * this serves as an example of usage of filters library
 *
 *
 */

int main(int argc, char **argv) {

    /* initialize filter(s). running two instances simultaneously */
    FilterInfo filterInfoR, filterInfoI;
    memset(&filterInfoR, 0, sizeof(FilterInfo));
    memset(&filterInfoI, 0, sizeof(FilterInfo));
    filterInfoR.options = FilterTypeLowPass | FilterVariantButterworth | FilterImplementationIIR;
    filterInfoI.options = FilterTypeLowPass | FilterVariantButterworth | FilterImplementationIIR;
    filterInfoR.order = filterInfoI.order = 2;
    double cornerFrequency = 200.0;
    double samplingFrequency = 1000.0;
    filterInfoR.fs = filterInfoI.fs = samplingFrequency;
    filterInfoR.f1 = filterInfoI.f1 = cornerFrequency;
    filter_init(&filterInfoR);
    filter_init(&filterInfoI);

    filterInfoR.print(&filterInfoR);
    filterInfoI.print(&filterInfoI);

    /* generate input samples, at different frequencies, and run filter on them */

    unsigned int runLength = 1000; /* let filter run for 5 seconds at 200Hz */
    unsigned int steadyStatePoint = 100; /* let filter run for a little bit before start to analyze */
    for (double frequency = 0.0; frequency < samplingFrequency; frequency += 1.0) {

        double accSR = 0.0, accSI = 0.0, accFR = 0.0, accFI = 0.0;
        double acc2SR = 0.0, acc2SI = 0.0, acc2FR = 0.0, acc2FI = 0.0;
        double dp = 2.0 * M_PI * frequency / samplingFrequency;
        for (unsigned int i = 0; i < runLength; ++i) {
            double arg = i * dp;

            double so, si = cos(arg); /* TODO how large do arguments get before bad things happen? */
            so = filterInfoR.sample(&filterInfoR, si);
            if (i >= steadyStatePoint) {
                accSR = accSR + si;
                acc2SR = acc2SR + si * si;
                accFR = accFR + so;
                acc2FR = acc2FR + so * so;
            }

//            si = sin(arg);
//            so = filterInfoI.sample(&filterInfoI, si);
//            if (i > steadyStatePoint) {
//                accSI = accSI + si;
//                acc2SI = acc2SI + si * si;
//                accFI = accFI + so;
//                acc2FI = acc2FI + so * so;
//            }
        }

//        double phaseS = atan2(accSI, accSR);
//        double phaseF = atan2(accFI, accFR);
//        double phase = phaseS - phaseF;
        double ampR = sqrt(acc2FR / acc2SR);
//        double ampI = sqrt(acc2FI / acc2SI);
        // double phase = filterInfoR.phase(&filterInfoR, frequency, runLength);
//        printf("phase response at frequency = %f is %f radians. amp = %f, %f\n", frequency, phase, ampR, ampI);
        printf("amplitude response = %f at frequency=%f\n", ampR, frequency);
        getchar();
    }
}

