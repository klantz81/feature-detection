#include "kernel.h"

cKernel::cKernel(int kernel_radius, double mu, double sigma) : kernel_radius(kernel_radius), mu(mu), sigma(sigma) {
	this->kernel            = new double[this->kernel_radius * 2 + 1];
	this->kernel_derivative = new double[this->kernel_radius * 2 + 1];

	for (int _i = -this->kernel_radius; _i <= this->kernel_radius; _i++) {
		this->kernel[_i + this->kernel_radius] = 1.0 / (this->sigma * sqrt(2.0 * M_PI)) * exp(-0.5 * pow((_i - this->mu) / this->sigma, 2.0));
		this->kernel_derivative[_i + this->kernel_radius] = - (_i - mu) / (pow(this->sigma, 2.0) * sqrt(2.0 * M_PI)) * exp(-0.5 * pow((_i - this->mu) / this->sigma, 2.0));
	}

	double sum0 = 0.0, sum1 = 0.0;
	for (int _i = -this->kernel_radius; _i <= this->kernel_radius; _i++) {
		sum0 += this->kernel[_i + this->kernel_radius];
		sum1 -= _i * this->kernel_derivative[_i + this->kernel_radius];
	}
	for (int _i = -this->kernel_radius; _i <= this->kernel_radius; _i++) {
		this->kernel[_i + this->kernel_radius] /= sum0;
		this->kernel_derivative[_i + this->kernel_radius] /= sum1;
	}

	sum0 = 0.0; sum1 = 0.0;
	for (int _i = -this->kernel_radius; _i <= this->kernel_radius; _i++) {
		sum0 += this->kernel[_i + this->kernel_radius];
		sum1 -= _i * this->kernel_derivative[_i + this->kernel_radius];
	}
	std::cout << "kernel summation: " << sum0 << std::endl;
	std::cout << "kernel derivative summation: " << sum1 << std::endl;
}

cKernel::~cKernel() {
	delete [] this->kernel;
	delete [] this->kernel_derivative;
}
