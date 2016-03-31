#ifndef PATCH_H
#define PATCH_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "defines.h"
#include "kernel.h"
#include "buffers.h"

class cPatch {
private:
protected:
public:
	static int save_index;
	int kernel_radius;
	double *patch;

	cPatch();
	cPatch(const cPatch& B);
	cPatch& operator=(const cPatch& B);
	~cPatch();
	void setup(double* buffer, int width, int height, int x, int y, const cKernel& kernel);
	void setup(const cDoubleBuffer& buffer, int x, int y, const cKernel& kernel);
	void save();
	void save(const char* file);
};

#endif