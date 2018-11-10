#include <string.h>
#include <stdio.h>

#include "filters.h"

int main(int argc, char **argv) {
    FilterInfo filterInfo;
    memset(&filterInfo, 0, sizeof(FilterInfo));
    filterInfo.options = FilterTypeLowPass | FilterVariantButterworth | FilterImplementationIIR;
    filterInfo.order = 2;
    float samplingFrequency = 1000.0;
    filterInfo.fs = samplingFrequency;
    filterInfo.f1 = 200.0;
    // filter_init(&filterInfo, samplingFrequency, cornerFrequency, cornerFrequency, filterType, filterOrder);
    filter_init(&filterInfo);
    unsigned int runLength = 1000; /* let filter run for 5 seconds at 200Hz */
    for (double frequency = 0.0; frequency < samplingFrequency; frequency += 1.0) {
        double phase = filterInfo.phase(&filterInfo, frequency, runLength);
        printf("phase response at frequency = %f is %f radians\n", frequency, phase);
    }
}

