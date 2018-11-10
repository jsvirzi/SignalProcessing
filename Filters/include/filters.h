#ifndef FILTER_H
#define FILTER_H

#include <stdint.h>

struct FilterInfo;

typedef struct FilterInfo {

    /* user specified inputs */
    double fs; /* sampling frequency */
    double f1; /* 1st corner */
    double f2; /* 2nd corner (if bandpass/notch) */
    int order; /* filter order */
    unsigned int options; /* specify filter options - Butterworth, IIR vs FIR, etc */

    /* internal */
	int nx, ny;
	double *x, *y;
	int nc, nd;
	double *c, *d;
	int ni;
	uint32_t queue_size, queue_mask, queue_head;

	/* user methods */
	int (*init)(struct FilterInfo *filterInfo, double fs, double f1, double f2, int type, int order);
	double (*sample)(struct FilterInfo *filterInfo, const double x);
	void (*batch)(struct FilterInfo *filterInfo, const double *x, const int n, double *y);
	double (*amplitude)(struct FilterInfo *filterInfo, const double f, const unsigned int run_length);
	double (*phase)(struct FilterInfo *filterInfo, const double f, unsigned int run_length);
	int (*close)(struct FilterInfo *filterInfo);

} FilterInfo;

enum {
	FilterTypeLowPass = 0,
	FilterTypeHighPass,
	FilterTypeBandPass,
	FilterTypeBandStop,
	FilterTypes
};

const uint32_t FilterTypeMask = 0xf;

enum {
	FilterImplementationFIR = 0 * 0x100,
	FilterImplementationIIR = 1 * 0x100
};

const uint32_t FilterImplementationMask = 0x300;

enum {
	FilterVariantButterworth = 0 * 0x10000,
	FilterVariantChebyshev = 1 * 0x10000,
	FilterVariants = 2 * 0x10000
};

const uint32_t FilterVariantsMask = 0x70000;

int filter_init(FilterInfo *filterInfo);
// int filter_init(FilterInfo *filterInfo, double fs, double f1, double f2, int type, int order);
// double filter_sample(FilterInfo *filterInfo, const double x);
// void filter_batch(FilterInfo *filterInfo, const double *x, const int n, double *y);
// int filter_close(FilterInfo *filterInfo);

#endif
