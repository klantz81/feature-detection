#ifndef BUFFERS_H
#define BUFFERS_H

#include <SDL/SDL.h>
#include <SDL/SDL_image.h>
#include <fstream>
#include <sstream>
#include "matrix.h"
#include "vector.h"
#include "kernel.h"
#include "defines.h"

class cFloatBuffer {
private:
protected:
public:
	int width, height;
	float* buffer;
	
	cFloatBuffer();
	cFloatBuffer(int width, int height);

	~cFloatBuffer();

	cFloatBuffer& operator=(const cFloatBuffer& buffer);
};

class cDoubleBuffer {
private:
protected:
public:
	static int save_index;

	int width, height;
	double* buffer;
	
	cDoubleBuffer();
	cDoubleBuffer(int width, int height);
	cDoubleBuffer(const cDoubleBuffer& buffer);

	~cDoubleBuffer();

	cDoubleBuffer& operator=(const cDoubleBuffer& buffer);
	cFloatBuffer toFloat(float scale_factor = 32.0f);
	void save();
	void save(const char* file);
	double operator[](int index) const;
	double at(int index) const;
	double at(int x, int y) const;
	double interpolate(double x, double y) const;
	double sum(int x, int y, const cKernel& kernel);

	cDoubleBuffer operator*(const cDoubleBuffer& buffer);
	cDoubleBuffer magnitude(const cDoubleBuffer& buffer);

	cDoubleBuffer sub(int x, int y, const cKernel& kernel) const;
	cDoubleBuffer subInterpolate(double x, double y, const cKernel& kernel) const;
	cDoubleBuffer subInterpolate(matrix22 A, vector2 b, const cKernel& kernel) const;

	cDoubleBuffer downsample();

	cDoubleBuffer scharrOperatorX();
	cDoubleBuffer scharrOperatorY();
	cDoubleBuffer scharrOperator();

	cDoubleBuffer gaussianConvolutionX(const cKernel& kernel);
	cDoubleBuffer gaussianConvolutionY(const cKernel& kernel);
	cDoubleBuffer gaussianConvolution(const cKernel& kernel);

	cDoubleBuffer gaussianConvolutionDX(const cKernel& kernel);
	cDoubleBuffer gaussianConvolutionDY(const cKernel& kernel);

	cDoubleBuffer gaussianDerivativeX(const cKernel& kernel);
	cDoubleBuffer gaussianDerivativeY(const cKernel& kernel);
};

struct _rgb{
	unsigned char r, g, b;
};

class cColorBuffer {
private:
protected:
public:
	static int save_index;

	int width, height, bpp;
	unsigned char* buffer;
	
	cColorBuffer();
	cColorBuffer(int width, int height, int bpp);
	cColorBuffer(const char* image);
	void save();
	void save(const char* file);

	_rgb at(int x, int y) const;
	_rgb interpolate(double x, double y) const;

	cColorBuffer sub(int x, int y, const cKernel& kernel) const;
	cColorBuffer subInterpolate(double x, double y, const cKernel& kernel) const;
	cColorBuffer subInterpolate(matrix22 A, vector2 b, const cKernel& kernel) const;

	~cColorBuffer();

	void load(const char* image);
	cDoubleBuffer grayscale();
};

#endif