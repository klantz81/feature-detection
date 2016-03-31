#ifndef KERNEL_H
#define KERNEL_H

#include <math.h>
#include <iostream>

class cKernel {
private:
protected:
public:
	int kernel_radius;
	double mu, sigma; // mean and standard deviation
	double *kernel, *kernel_derivative;

	cKernel(int kernel_radius, double mu, double sigma);
	~cKernel();
};

#endif