#include "patch.h"

int cPatch::save_index = 0;

cPatch::cPatch() : kernel_radius(0) {
	this->patch = 0;
}

cPatch::cPatch(const cPatch& B) : kernel_radius(B.kernel_radius), patch(0) {
	if (this->kernel_radius > 0) {
		int size = (this->kernel_radius * 2 + 1) * (this->kernel_radius * 2 + 1);
		this->patch = new double[size];
		for (int _i = 0; _i < size; _i++) this->patch[_i] = B.patch[_i];
	}
}

cPatch& cPatch::operator=(const cPatch& B) {
	int size = (B.kernel_radius * 2 + 1) * (B.kernel_radius * 2 + 1);
	if (this->kernel_radius != B.kernel_radius) {
		if (this->patch) delete [] this->patch;
		this->patch = 0;
		this->kernel_radius = B.kernel_radius;
		if (this->kernel_radius > 0) this->patch = new double[size];
	}
	if (this->kernel_radius > 0) for (int _i = 0; _i < size; _i++) this->patch[_i] = B.patch[_i];
	return *this;
}

cPatch::~cPatch() {
	if (this->patch) delete [] this->patch;
}

void cPatch::setup(double* buffer, int width, int height, int x, int y, const cKernel& kernel) {
	if (this->kernel_radius != kernel.kernel_radius) {
		int size = (kernel.kernel_radius * 2 + 1) * (kernel.kernel_radius * 2 + 1);
		if (this->patch) delete [] this->patch;
		this->patch = 0;
		this->kernel_radius = kernel.kernel_radius;
		if (this->kernel_radius > 0) this->patch = new double[size];
	}
	int index = 0;
	if (this->kernel_radius > 0)
		for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++)
			for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++)
				this->patch[index++] = buffer[(y + _j) * width + (x + _i)];
}

void cPatch::setup(const cDoubleBuffer& buffer, int x, int y, const cKernel& kernel) {
	if (this->kernel_radius != kernel.kernel_radius) {
		int size = (kernel.kernel_radius * 2 + 1) * (kernel.kernel_radius * 2 + 1);
		if (this->patch) delete [] this->patch;
		this->patch = 0;
		this->kernel_radius = kernel.kernel_radius;
		if (this->kernel_radius > 0) this->patch = new double[size];
	}
	int index = 0;
	if (this->kernel_radius > 0)
		for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++)
			for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++)
				this->patch[index++] = buffer.at(x + _i, y + _j);
}

void cPatch::save() {
	cPatch::save_index++;

	int width = this->kernel_radius * 2 + 1;

	if (width > 1) {
		std::stringstream s;
		s << "media_save/patch" << (save_index < 10 ? "000" : (save_index < 100 ? "00" : "0")) << save_index << ".pgm";
		std::fstream of(s.str().c_str(), std::ios_base::out | std::ios_base::trunc);

		
		of << "P2" << std::endl;
		of << width << " " << width << std::endl;
		of << "255" << std::endl;
		
		for (int j = 0; j < width; j++) {
			for (int i = 0; i < width; i++) {
				int p = (int)(this->patch[j * width + i] * 255.0);
				p = CLAMP(p, 0, 255);
				of << p << " ";
			}
			of << std::endl;
		}
	}
}

void cPatch::save(const char* file) {
	int width = this->kernel_radius * 2 + 1;

	if (width > 1) {
		std::fstream of(file, std::ios_base::out | std::ios_base::trunc);

		of << "P2" << std::endl;
		of << width << " " << width << std::endl;
		of << "255" << std::endl;
		
		for (int j = 0; j < width; j++) {
			for (int i = 0; i < width; i++) {
				int p = (int)(this->patch[j * width + i] * 255.0);
				p = CLAMP(p, 0, 255);
				of << p << " ";
			}
			of << std::endl;
		}
	}
}
