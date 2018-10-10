#ifndef FILTER_H
#define FILTER_H

typedef struct FilterInfo {
	int nx, ny;
	double *x, *y;
} FilterInfo;

enum {
	FilterTypeLowPass,
	FilterTypeHighPass,
	FilterTypeBandPass,
	FilterTypeBandStop,
	FilterTypes
};

int init(FilterInfo *filterInfo, double fs, double f1, double f2, int type, int order);
double filter(FilterInfo *filterInfo, double *y);

#endif
