#ifndef FILTER_H
#define FILTER_H

typedef struct FilterInfo {
	int nx, ny;
	double *x, *y;
	int nc, nd;
	double *c, *d;
} FilterInfo;

enum {
	FilterTypeLowPass,
	FilterTypeHighPass,
	FilterTypeBandPass,
	FilterTypeBandStop,
	FilterTypes
};

int init(FilterInfo *filterInfo, double fs, double f1, double f2, int type, int order);
double filter_sample(FilterInfo *filterInfo, const double *y);
int filter_batch(FilterInfo *filterInfo, const double *y, const int n, double *y_out);
int close(FilterInfo *filterInfo);

#endif
