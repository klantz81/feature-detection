#include "buffers.h"

// cFloatBuffer /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cFloatBuffer::cFloatBuffer() : width(0), height(0), buffer(0) { }
cFloatBuffer::cFloatBuffer(int width, int height) : width(width), height(height) { this->buffer = new float[this->width*this->height]; for (int i = 0; i < this->width * this->height; i++) this->buffer[i] = 0.0f; }

cFloatBuffer::~cFloatBuffer() { if (this->buffer) delete [] buffer; }

cFloatBuffer& cFloatBuffer::operator=(const cFloatBuffer& buffer) {
	if (this->width != buffer.width || this->height != buffer.height) {
		this->width = buffer.width;
		this->height = buffer.height;
		if (this->buffer) delete [] this->buffer;
		this->buffer = new float[this->width*this->height];
	}
	for (int i = 0; i < this->width * this->height; i++)
		this->buffer[i] = buffer.buffer[i];
	return *this;
}



// cDoubleBuffer /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int cDoubleBuffer::save_index = 0;
cDoubleBuffer::cDoubleBuffer() : width(0), height(0), buffer(0) { }
cDoubleBuffer::cDoubleBuffer(int width, int height) : width(width), height(height) { this->buffer = new double[this->width*this->height]; for (int i = 0; i < this->width * this->height; i++) this->buffer[i] = 0.0; }
cDoubleBuffer::cDoubleBuffer(const cDoubleBuffer& buffer) : width(buffer.width), height(buffer.height), buffer(0) {
	this->buffer = new double[this->width*this->height];
	for (int i = 0; i < this->width*this->height; i++) this->buffer[i] = buffer.buffer[i];
}

cDoubleBuffer::~cDoubleBuffer() { if (this->buffer) delete [] buffer; }

cDoubleBuffer& cDoubleBuffer::operator=(const cDoubleBuffer& buffer) {
	if (this->width != buffer.width || this->height != buffer.height) {
		this->width = buffer.width;
		this->height = buffer.height;
		if (this->buffer) delete [] this->buffer;
		this->buffer = new double[this->width*this->height];
	}
	for (int i = 0; i < this->width * this->height; i++)
		this->buffer[i] = buffer.buffer[i];
	return *this;
}

cFloatBuffer cDoubleBuffer::toFloat(float scale_factor) {
	cFloatBuffer buf(this->width, this->height);
	for (int i = 0; i < this->width * this->height; i++)
		buf.buffer[i] = (float)this->buffer[i] * scale_factor;
	return buf;
}

void cDoubleBuffer::save() {
	cDoubleBuffer::save_index++;

	std::stringstream s;
	s << "media_save/buffer" << (save_index < 10 ? "000" : (save_index < 100 ? "00" : (save_index < 1000 ? "0" : ""))) << save_index << ".pgm";
	std::fstream of(s.str().c_str(), std::ios_base::out | std::ios_base::trunc);

	of << "P2" << std::endl;
	of << this->width << " " << this->height << std::endl;
	of << "255" << std::endl;

	for (int j = 0; j < this->height; j++) {
		for (int i = 0; i < this->width; i++) {
			int p = (int)(this->buffer[j * this->width + i] * 255.0);
			p = CLAMP(p, 0, 255);
			of << p << " ";
		}
		of << std::endl;
	}
}

void cDoubleBuffer::save(const char* file) {
	std::fstream of(file, std::ios_base::out | std::ios_base::trunc);

	of << "P2" << std::endl;
	of << this->width << " " << this->height << std::endl;
	of << "255" << std::endl;

	for (int j = 0; j < this->height; j++) {
		for (int i = 0; i < this->width; i++) {
			int p = (int)(this->buffer[j * this->width + i] * 255.0);
			p = CLAMP(p, 0, 255);
			of << p << " ";
		}
		of << std::endl;
	}
}

double cDoubleBuffer::operator[](int index) const {
	return (index >= 0 && index < this->width * this->height) ? this->buffer[index] : 0.0;
}
double cDoubleBuffer::at(int index) const {
	return (index >= 0 && index < this->width * this->height) ? this->buffer[index] : 0.0;
}
double cDoubleBuffer::at(int x, int y) const {
	int index = y * this->width + x;
	return (index >= 0 && index < this->width * this->height) ? this->buffer[index] : 0.0;
}

double cDoubleBuffer::interpolate(double x, double y) const {
	x = CLAMP(x, 0.000001, this->width  - 1.000001);
	y = CLAMP(y, 0.000001, this->height - 1.000001);
	int   x0 = (int)x,
	      y0 = (int)y;
	double s0 = x - x0,
	       s1 = y - y0;

	double t = this->buffer[(y0 + 0) * this->width + x0 + 0]*(1 - s0)*(1 - s1) +
		   this->buffer[(y0 + 0) * this->width + x0 + 1]*(    s0)*(1 - s1) +
		   this->buffer[(y0 + 1) * this->width + x0 + 0]*(1 - s0)*(    s1) +
		   this->buffer[(y0 + 1) * this->width + x0 + 1]*(    s0)*(    s1);
	return t;
}

double cDoubleBuffer::sum(int x, int y, const cKernel& kernel) {
	double sum = 0.0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++)
			sum += this->at(x + i, y + j);
	return sum;
}

cDoubleBuffer cDoubleBuffer::operator*(const cDoubleBuffer& buffer) {
	cDoubleBuffer buf(this->width, this->height);
	for (int i = 0; i < this->width * this->height; i++)
		buf.buffer[i] = this->buffer[i] * buffer.buffer[i];
	return buf;
}

cDoubleBuffer cDoubleBuffer::magnitude(const cDoubleBuffer& buffer) {
	cDoubleBuffer buf(this->width, this->height);
	for (int i = 0; i < this->width * this->height; i++)
		buf.buffer[i] = sqrt(pow(this->buffer[i], 2.0) + pow(buffer.buffer[i], 2.0));
	return buf;
}

cDoubleBuffer cDoubleBuffer::sub(int x, int y, const cKernel& kernel) const {
	cDoubleBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++)
			buf.buffer[index++] = this->at(x + i, y + j);

	return buf;
}

cDoubleBuffer cDoubleBuffer::subInterpolate(double x, double y, const cKernel& kernel) const {
	cDoubleBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++)
			buf.buffer[index++] = this->interpolate(x + i, y + j);

	return buf;
}

cDoubleBuffer cDoubleBuffer::subInterpolate(matrix22 A, vector2 b, const cKernel& kernel) const {
	cDoubleBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++) {
			vector2 x = A * vector2(i, j) + b;
			buf.buffer[index++] = this->interpolate(x.x, x.y);
		}
	return buf;
}

cDoubleBuffer cDoubleBuffer::downsample() {
	cDoubleBuffer buf(this->width/2, this->height/2);
	for (int j = 0; j < this->height / 2; j++)
		for (int i = 0; i < this->width / 2; i++)
			buf.buffer[j * this->width / 2 + i] = 0.25   *  this->buffer[j * 2 * this->width + i * 2] +
							      0.125  * (
									(j > 0				? this->buffer[(j - 1) * 2 * this->width + i * 2] : 0.0) +
									(j < this->height / 2 - 1	? this->buffer[(j + 1) * 2 * this->width + i * 2] : 0.0) +
									(i > 0				? this->buffer[j * 2 * this->width + (i - 1) * 2] : 0.0) +
									(i < this->width / 2 - 1	? this->buffer[j * 2 * this->width + (i + 1) * 2] : 0.0)
									) +
							      0.0625 * (
									(j > 0 && i > 0						? this->buffer[(j - 1) * 2 * this->width + (i - 1) * 2] : 0.0) +
									(j > 0 && i < this->width / 2 - 1			? this->buffer[(j - 1) * 2 * this->width + (i + 1) * 2] : 0.0) +
									(j < this->height / 2 - 1 && i > 0			? this->buffer[(j + 1) * 2 * this->width + (i - 1) * 2] : 0.0) +
									(j < this->height / 2 - 1 && i < this->width / 2 - 1	? this->buffer[(j + 1) * 2 * this->width + (i + 1) * 2] : 0.0)
									);
	return buf;
}


cDoubleBuffer cDoubleBuffer::scharrOperatorX() {
	cDoubleBuffer Ix(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			Ix.buffer[j * this->width + i] =  (i > 0 ?
							      (
								      (j > 0 ?		      this->buffer[(j - 1) * this->width + (i - 1)] *  -3.0 : 0.0) +
											      this->buffer[ j      * this->width + (i - 1)] * -10.0 +
								      (j < this->height - 1 ? this->buffer[(j + 1) * this->width + (i - 1)] *  -3.0 : 0.0)
							      ) : 0.0) +
							  (i < this->width - 1 ?
							      (
								      (j > 0 ?		      this->buffer[(j - 1) * this->width + (i + 1)] *  3.0 : 0.0) +
											      this->buffer[ j      * this->width + (i + 1)] * 10.0 +
								      (j < this->height - 1 ? this->buffer[(j + 1) * this->width + (i + 1)] *  3.0 : 0.0)
							      ) : 0.0);
	return Ix;
}

cDoubleBuffer cDoubleBuffer::scharrOperatorY() {
	cDoubleBuffer Iy(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			Iy.buffer[j * this->width + i] =  (j > 0 ?
							      (
								      (i > 0 ?		     this->buffer[(j - 1) * this->width + (i - 1)] *  -3.0 : 0.0) +
											     this->buffer[(j - 1) * this->width +  i     ] * -10.0 +
								      (i < this->width - 1 ? this->buffer[(j - 1) * this->width + (i + 1)] *  -3.0 : 0.0)
							      ) : 0.0) +
							  (j < this->height - 1 ?
							      (
								      (i > 0 ?		     this->buffer[(j + 1) * this->width + (i - 1)] *  3.0 : 0.0) +
											     this->buffer[(j + 1) * this->width +  i     ] * 10.0 +
								      (i < this->width - 1 ? this->buffer[(j + 1) * this->width + (i + 1)] *  3.0 : 0.0)
							      ) : 0.0);
	return Iy;
}

cDoubleBuffer cDoubleBuffer::scharrOperator() {
	cDoubleBuffer Ix = this->scharrOperatorX();
	cDoubleBuffer Iy = this->scharrOperatorY();
	cDoubleBuffer I(this->width, this->height);
	for (int i = 0; i < this->width * this->height; i++)
		I.buffer[i] = sqrt(Ix.buffer[i]*Ix.buffer[i] + Iy.buffer[i]*Iy.buffer[i]);
	return I;
}

cDoubleBuffer cDoubleBuffer::gaussianConvolutionX(const cKernel& kernel) {
	cDoubleBuffer buf(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			for (int k = -kernel.kernel_radius; k <= kernel.kernel_radius; k++)
				if (i + k >= 0 && i + k < this->width)
					buf.buffer[j * this->width + i] += this->buffer[j * this->width + i + k] * kernel.kernel[k + kernel.kernel_radius];
	return buf;
}
cDoubleBuffer cDoubleBuffer::gaussianConvolutionY(const cKernel& kernel) {
	cDoubleBuffer buf(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			for (int k = -kernel.kernel_radius; k <= kernel.kernel_radius; k++)
				if (j + k >= 0 && j + k < this->height)
					buf.buffer[j * this->width + i] += this->buffer[(j + k) * this->width + i] * kernel.kernel[k + kernel.kernel_radius];
	return buf;
}
cDoubleBuffer cDoubleBuffer::gaussianConvolution(const cKernel& kernel) {
	return (this->gaussianConvolutionX(kernel)).gaussianConvolutionY(kernel);
}

cDoubleBuffer cDoubleBuffer::gaussianConvolutionDX(const cKernel& kernel) {
	cDoubleBuffer buf(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			for (int k = -kernel.kernel_radius; k <= kernel.kernel_radius; k++)
				if (i + k >= 0 && i + k < this->width)
					buf.buffer[j * this->width + i] += this->buffer[j * this->width + i + k] * kernel.kernel_derivative[k + kernel.kernel_radius];
	return buf;
}
cDoubleBuffer cDoubleBuffer::gaussianConvolutionDY(const cKernel& kernel) {
	cDoubleBuffer buf(this->width, this->height);
	for (int j = 0; j < this->height; j++)
		for (int i = 0; i < this->width; i++)
			for (int k = -kernel.kernel_radius; k <= kernel.kernel_radius; k++)
				if (j + k >= 0 && j + k < this->height)
					buf.buffer[j * this->width + i] += this->buffer[(j + k) * this->width + i] * kernel.kernel_derivative[k + kernel.kernel_radius];
	return buf;
}

cDoubleBuffer cDoubleBuffer::gaussianDerivativeX(const cKernel& kernel) {
	return (this->gaussianConvolutionDX(kernel)).gaussianConvolutionY(kernel);
}
cDoubleBuffer cDoubleBuffer::gaussianDerivativeY(const cKernel& kernel) {
	return (this->gaussianConvolutionX(kernel)).gaussianConvolutionDY(kernel);
}



// cColorBuffer /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int cColorBuffer::save_index = 0;
cColorBuffer::cColorBuffer() : width(0), height(0), bpp(0), buffer(0) { }
cColorBuffer::cColorBuffer(int width, int height, int bpp) : width(width), height(height), bpp(bpp) { this->buffer = new unsigned char[this->width*this->height*this->bpp]; for (int i = 0; i < this->width * this->height * this->bpp; i++) this->buffer[i] = 0; }
cColorBuffer::cColorBuffer(const char* image) {
	SDL_Surface *s = IMG_Load(image);
	this->width = s->w;
	this->height = s->h;
	this->bpp = s->format->BitsPerPixel;
	this->buffer = new unsigned char[width * height * bpp / 8];
	memcpy(buffer, s->pixels, width * height * bpp / 8);
	SDL_FreeSurface(s);
}

void cColorBuffer::save() {
	cColorBuffer::save_index++;

	std::stringstream s;
	s << "media_save/cbuffer" << (save_index < 10 ? "000" : (save_index < 100 ? "00" : (save_index < 1000 ? "0" : ""))) << save_index << ".ppm";
	std::fstream of(s.str().c_str(), std::ios_base::out | std::ios_base::trunc);

	of << "P3" << std::endl;
	of << this->width << " " << this->height << std::endl;
	of << "255" << std::endl;

	for (int j = 0; j < this->height; j++) {
		for (int i = 0; i < this->width; i++) {
			int r = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 0]; r = CLAMP(r, 0, 255);
			int g = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 1]; g = CLAMP(g, 0, 255);
			int b = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 2]; b = CLAMP(b, 0, 255);
			of << r << " " << g << " " << b << " ";
		}
		of << std::endl;
	}
}

void cColorBuffer::save(const char* file) {
	std::fstream of(file, std::ios_base::out | std::ios_base::trunc);

	of << "P3" << std::endl;
	of << this->width << " " << this->height << std::endl;
	of << "255" << std::endl;

	for (int j = 0; j < this->height; j++) {
		for (int i = 0; i < this->width; i++) {
			int r = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 0]; r = CLAMP(r, 0, 255);
			int g = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 1]; g = CLAMP(g, 0, 255);
			int b = this->buffer[j * this->width * this->bpp/8 + i * this->bpp/8 + 2]; b = CLAMP(b, 0, 255);
			of << r << " " << g << " " << b << " ";
		}
		of << std::endl;
	}
}

_rgb cColorBuffer::at(int x, int y) const {
	int index = y * this->width * this->bpp/8 + x * this->bpp/8;
	_rgb rgb = {0,0,0};

	if (index >= 0 && index < this->width * this->bpp / 8 * this->height - 2) {
		rgb.r = this->buffer[index + 0];
		rgb.g = this->buffer[index + 1];
		rgb.b = this->buffer[index + 2];
	}
	
	return rgb;
}

_rgb cColorBuffer::interpolate(double x, double y) const {
	x = CLAMP(x, 0.000001, this->width  - 1.000001);
	y = CLAMP(y, 0.000001, this->height - 1.000001);
	int   x0 = (int)x,
	      y0 = (int)y;
	double s0 = x - x0,
	       s1 = y - y0;

	_rgb rgb = {0,0,0};

	rgb.r = this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 0]*(1 - s0)*(1 - s1) +
		this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 0]*(    s0)*(1 - s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 0]*(1 - s0)*(    s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 0]*(    s0)*(    s1);
	rgb.g = this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 1]*(1 - s0)*(1 - s1) +
		this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 1]*(    s0)*(1 - s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 1]*(1 - s0)*(    s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 1]*(    s0)*(    s1);
	rgb.b = this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 2]*(1 - s0)*(1 - s1) +
		this->buffer[(y0 + 0) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 2]*(    s0)*(1 - s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 0) * this->bpp/8 + 2]*(1 - s0)*(    s1) +
		this->buffer[(y0 + 1) * this->width * this->bpp/8 + (x0 + 1) * this->bpp/8 + 2]*(    s0)*(    s1);
	return rgb;
}

cColorBuffer cColorBuffer::sub(int x, int y, const cKernel& kernel) const {
	cColorBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1, this->bpp);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++) {
			_rgb rgb = this->at(x + i, y + j);
			buf.buffer[index + 0] = rgb.r;
			buf.buffer[index + 1] = rgb.g;
			buf.buffer[index + 2] = rgb.b;
			index += buf.bpp / 8;
		}

	return buf;
}

cColorBuffer cColorBuffer::subInterpolate(double x, double y, const cKernel& kernel) const {
	cColorBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1, this->bpp);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++) {
			_rgb rgb = this->interpolate(x + i, y + j);
			buf.buffer[index + 0] = rgb.r;
			buf.buffer[index + 1] = rgb.g;
			buf.buffer[index + 2] = rgb.b;
			index += buf.bpp / 8;
		}

	return buf;
}

cColorBuffer cColorBuffer::subInterpolate(matrix22 A, vector2 b, const cKernel& kernel) const {
	cColorBuffer buf(kernel.kernel_radius * 2 + 1, kernel.kernel_radius * 2 + 1, this->bpp);

	int index = 0;
	for (int j = -kernel.kernel_radius; j <= kernel.kernel_radius; j++)
		for (int i = -kernel.kernel_radius; i <= kernel.kernel_radius; i++) {
			vector2 x = A * vector2(i, j) + b;
			_rgb rgb = this->interpolate(x.x, x.y);
			buf.buffer[index + 0] = rgb.r;
			buf.buffer[index + 1] = rgb.g;
			buf.buffer[index + 2] = rgb.b;
			index += buf.bpp / 8;
		}
	return buf;
}

cColorBuffer::~cColorBuffer() { if (this->buffer) delete [] buffer; }

void cColorBuffer::load(const char* image) {
	SDL_Surface *s = IMG_Load(image);
	if (this->width != s->w || this->height != s->h || this->bpp != s->format->BitsPerPixel) {
		this->width = s->w;
		this->height = s->h;
		this->bpp = s->format->BitsPerPixel;
		if (this->buffer) delete [] this->buffer;
		this->buffer = new unsigned char[this->width*this->height*this->bpp];
	}
	memcpy(this->buffer, s->pixels, this->width * this->height * this->bpp / 8);
	SDL_FreeSurface(s);
}

cDoubleBuffer cColorBuffer::grayscale() {
	cDoubleBuffer buf(this->width, this->height);
	for (int i = 0; i < this->width * this->height; i++)
		buf.buffer[i] = (this->buffer[i * this->bpp / 8 + 0] / 255.0 * 0.30 +
			         this->buffer[i * this->bpp / 8 + 1] / 255.0 * 0.59 +
			         this->buffer[i * this->bpp / 8 + 2] / 255.0 * 0.11);
	return buf;
}
