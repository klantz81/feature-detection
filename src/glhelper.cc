#include "glhelper.h"

#define CLAMP(x,a,b) ( ((x)<(a))?(a):(((x)>(b))?(b):(x)) )

unsigned char* loadImage(const char* image, int& width, int& height, int& bpp) {
	SDL_Surface *s = IMG_Load(image);
	width = s->w; height = s->h; bpp = s->format->BitsPerPixel;
	unsigned char *buffer = new unsigned char[width * height * bpp / 8];
	memcpy(buffer, s->pixels, width * height * bpp / 8);
	SDL_FreeSurface(s);
	return buffer;
}

double* grayscaleConversion(unsigned char* buffer, int width, int height, int bpp) {
	double *buffer_ = new double[width * height * sizeof(double)];
	double scale = 255.0;
	for (int i = 0; i < width * height; i++)
		buffer_[i] = (buffer[i * bpp / 8 + 0] / scale * 0.30 +
			      buffer[i * bpp / 8 + 1] / scale * 0.59 +
			      buffer[i * bpp / 8 + 2] / scale * 0.11);
	return buffer_;
}

float* doubleToFloat(double* buffer, int width, int height) {
	float *buffer_ = new float[width * height * sizeof(double)];
	for (int i = 0; i < width * height; i++)
		buffer_[i] = buffer[i];
	return buffer_;
}

void scharrOperator(double*& Ix, double*& Iy, double*& I, double* buffer, int width, int height) {
	Ix = new double[width * height];
	Iy = new double[width * height];
	I  = new double[width * height];
	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++) {
			Ix[j * width + i] = (i > 0 ?
						(
							(j > 0 ?          buffer[(j - 1) * width + (i - 1)] *  -3.0 : 0.0) +
							                  buffer[ j      * width + (i - 1)] * -10.0 +
							(j < height - 1 ? buffer[(j + 1) * width + (i - 1)] *  -3.0 : 0.0)
						) : 0.0) +
					    (i < width - 1 ?
						(
							(j > 0 ?          buffer[(j - 1) * width + (i + 1)] *  3.0 : 0.0) +
							                  buffer[ j      * width + (i + 1)] * 10.0 +
							(j < height - 1 ? buffer[(j + 1) * width + (i + 1)] *  3.0 : 0.0)
						) : 0.0);
			Iy[j * width + i] = (j > 0 ?
						(
							(i > 0 ?         buffer[(j - 1) * width + (i - 1)] *  -3.0 : 0.0) +
							                 buffer[(j - 1) * width +  i     ] * -10.0 +
							(i < width - 1 ? buffer[(j - 1) * width + (i + 1)] *  -3.0 : 0.0)
						) : 0.0) +
					    (j < height - 1 ?
						(
							(i > 0 ?         buffer[(j + 1) * width + (i - 1)] *  3.0 : 0.0) +
							                 buffer[(j + 1) * width +  i     ] * 10.0 +
							(i < width - 1 ? buffer[(j + 1) * width + (i + 1)] *  3.0 : 0.0)
						) : 0.0);
			I[j * width + i] = sqrt(pow(Ix[j * width + i], 2.0) + pow(Iy[j * width + i], 2.0));
		}
}

double* gaussianConvolution(double* buffer, int width, int height, int kernel_radius, double mean, double std_deviation) {
	double *kernel = new double[kernel_radius * 2 + 1];
	for (int i = -kernel_radius; i <= kernel_radius; i++)
		kernel[i + kernel_radius] = 1.0 / (std_deviation * sqrt(2.0 * M_PI)) * exp(-0.5 * pow((i - mean) / std_deviation, 2.0));

	double *_buffer = new double[width * height * sizeof(double)];
	double *buffer_ = new double[width * height * sizeof(double)];
	memset(_buffer, 0, width * height * sizeof(double));
	memset(buffer_, 0, width * height * sizeof(double));

	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			for (int k = -kernel_radius; k <= kernel_radius; k++)
				if (i + k >= 0 && i + k < width)
					_buffer[j * width + i] += buffer[j * width + i + k] * kernel[k + kernel_radius];

	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			for (int k = -kernel_radius; k <= kernel_radius; k++)
				if (j + k >= 0 && j + k < height)
					buffer_[j * width + i] += _buffer[(j + k) * width + i] * kernel[k + kernel_radius];

	delete [] kernel;
	delete [] _buffer;
	return buffer_;
}

void structureTensor(double*& Ixx, double*& Ixy, double*& Iyy, double* Ix, double* Iy, int width, int height, int kernel_radius, double mean, double std_deviation) {
	double *Ixx_ = new double[width * height];
	double *Ixy_ = new double[width * height];
	double *Iyy_ = new double[width * height];
	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++) {
			Ixx_[j * width + i] = Ix[j * width + i] * Ix[j * width + i];
			Ixy_[j * width + i] = Ix[j * width + i] * Iy[j * width + i];
			Iyy_[j * width + i] = Iy[j * width + i] * Iy[j * width + i];
		}
	Ixx = gaussianConvolution(Ixx_, width, height, kernel_radius, mean, std_deviation);
	Ixy = gaussianConvolution(Ixy_, width, height, kernel_radius, mean, std_deviation);
	Iyy = gaussianConvolution(Iyy_, width, height, kernel_radius, mean, std_deviation);
	delete [] Ixx_;
	delete [] Ixy_;
	delete [] Iyy_;
}

unsigned char* downsample(unsigned char* buffer, int width, int height, int bpp) {
	unsigned char *buffer_ = new unsigned char[width * height / 4 * sizeof(unsigned char) * bpp / 8];
	for (int j = 0; j < height / 2; j++) {
		for (int i = 0; i < width / 2; i++) {
			for (int k = 0; k < 3; k++) {
				buffer_[j * width / 2 * bpp / 8 + i * bpp / 8 + k] =
							    0.25   *  buffer[j * 2 * width * bpp / 8 + i * 2 * bpp / 8 + k] +
							    0.125  * (
								      (j > 0              ? buffer[(j - 1) * 2 * width * bpp / 8 + i       * 2 * bpp / 8 + k] : 0.0) +
								      (j < height / 2 - 1 ? buffer[(j + 1) * 2 * width * bpp / 8 + i       * 2 * bpp / 8 + k] : 0.0) +
								      (i > 0              ? buffer[j       * 2 * width * bpp / 8 + (i - 1) * 2 * bpp / 8 + k] : 0.0) +
								      (i < width / 2 - 1  ? buffer[j       * 2 * width * bpp / 8 + (i + 1) * 2 * bpp / 8 + k] : 0.0)
								      ) +
							    0.0625 * (
								      (j > 0 && i > 0                          ? buffer[(j - 1) * 2 * width * bpp / 8 + (i - 1) * 2 * bpp / 8 + k] : 0.0) +
								      (j > 0 && i < width / 2 - 1              ? buffer[(j - 1) * 2 * width * bpp / 8 + (i + 1) * 2 * bpp / 8 + k] : 0.0) +
								      (j < height / 2 - 1 && i > 0             ? buffer[(j + 1) * 2 * width * bpp / 8 + (i - 1) * 2 * bpp / 8 + k] : 0.0) +
								      (j < height / 2 - 1 && i < width / 2 - 1 ? buffer[(j + 1) * 2 * width * bpp / 8 + (i + 1) * 2 * bpp / 8 + k] : 0.0)
								     );
			}
		}
	}
	return buffer_;
}

double* downsample(double* buffer, int width, int height) {
	double *buffer_ = new double[width * height / 4 * sizeof(double)];
	for (int j = 0; j < height / 2; j++)
		for (int i = 0; i < width / 2; i++) {
			buffer_[j * width / 2 + i] = 0.25   *  buffer[j * 2 * width + i * 2] +
						     0.125  * (
							       (j > 0              ? buffer[(j - 1) * 2 * width + i * 2] : 0.0) +
							       (j < height / 2 - 1 ? buffer[(j + 1) * 2 * width + i * 2] : 0.0) +
							       (i > 0              ? buffer[j * 2 * width + (i - 1) * 2] : 0.0) +
							       (i < width / 2 - 1  ? buffer[j * 2 * width + (i + 1) * 2] : 0.0)
							      ) +
						     0.0625 * (
							       (j > 0 && i > 0                          ? buffer[(j - 1) * 2 * width + (i - 1) * 2] : 0.0) +
							       (j > 0 && i < width / 2 - 1              ? buffer[(j - 1) * 2 * width + (i + 1) * 2] : 0.0) +
							       (j < height / 2 - 1 && i > 0             ? buffer[(j + 1) * 2 * width + (i - 1) * 2] : 0.0) +
							       (j < height / 2 - 1 && i < width / 2 - 1 ? buffer[(j + 1) * 2 * width + (i + 1) * 2] : 0.0)
							      );
		}
	return buffer_;
}

float* downsample(float* buffer, int width, int height) {
	float *buffer_ = new float[width * height / 4 * sizeof(float)];
	for (int j = 0; j < height / 2; j++)
		for (int i = 0; i < width / 2; i++) {
			buffer_[j * width / 2 + i] = 0.25   *  buffer[j * 2 * width + i * 2] +
						     0.125  * (
							       (j > 0              ? buffer[(j - 1) * 2 * width + i * 2] : 0.0) +
							       (j < height / 2 - 1 ? buffer[(j + 1) * 2 * width + i * 2] : 0.0) +
							       (i > 0              ? buffer[j * 2 * width + (i - 1) * 2] : 0.0) +
							       (i < width / 2 - 1  ? buffer[j * 2 * width + (i + 1) * 2] : 0.0)
							      ) +
						     0.0625 * (
							       (j > 0 && i > 0                          ? buffer[(j - 1) * 2 * width + (i - 1) * 2] : 0.0) +
							       (j > 0 && i < width / 2 - 1              ? buffer[(j - 1) * 2 * width + (i + 1) * 2] : 0.0) +
							       (j < height / 2 - 1 && i > 0             ? buffer[(j + 1) * 2 * width + (i - 1) * 2] : 0.0) +
							       (j < height / 2 - 1 && i < width / 2 - 1 ? buffer[(j + 1) * 2 * width + (i + 1) * 2] : 0.0)
							      );
		}
	return buffer_;
}

double interpolate(double* buffer, int width, int height, double x, double y) {
	x = CLAMP(x, 0.000001, width  - 1.000001);
	y = CLAMP(y, 0.000001, height - 1.000001);
	int   x0 = (int)x,
	      y0 = (int)y;
	double s0 = x - x0,
	       s1 = y - y0;

	double t = buffer[(y0 + 0) * width + x0 + 0]*(1 - s0)*(1 - s1) +
		   buffer[(y0 + 0) * width + x0 + 1]*(    s0)*(1 - s1) +
		   buffer[(y0 + 1) * width + x0 + 0]*(1 - s0)*(    s1) +
		   buffer[(y0 + 1) * width + x0 + 1]*(    s0)*(    s1);
	return t;
}

void drawBox(unsigned char* buffer, int width, int height, int bpp, int x, int y, int half, unsigned char r, unsigned char g, unsigned char b) {
	int delta = 32;
	for (int j = y - half; j <= y + half; j++) {
		for (int i = x - half; i <= x + half; i++) {
			if (j >= 0 && j < height && i >= 0 && i < width) {
				buffer[j * width * bpp / 8 + i * bpp / 8 + 0] = (buffer[j * width * bpp / 8 + i * bpp / 8 + 0] < (255 - delta)) ? (buffer[j * width * bpp / 8 + i * bpp / 8 + 0] + delta) : 255;
				buffer[j * width * bpp / 8 + i * bpp / 8 + 1] = (buffer[j * width * bpp / 8 + i * bpp / 8 + 1] < (255 - delta)) ? (buffer[j * width * bpp / 8 + i * bpp / 8 + 1] + delta) : 255;
				buffer[j * width * bpp / 8 + i * bpp / 8 + 2] = (buffer[j * width * bpp / 8 + i * bpp / 8 + 2] < (255 - delta)) ? (buffer[j * width * bpp / 8 + i * bpp / 8 + 2] + delta) : 255;
			}
		}
	}
	for (int j = y - half; j <= y + half; j++) {
		int i = x;// - half;
		if (j >= 0 && j < height && i >= 0 && i < width) {
			buffer[j * width * bpp / 8 + i * bpp / 8 + 0] = r;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 1] = g;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 2] = b;
		}
/*		i = x + half;
		if (j >= 0 && j < height && i >= 0 && i < width) {
			buffer[j * width * bpp / 8 + i * bpp / 8 + 0] = r;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 1] = g;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 2] = b;
		}
*/	}
	for (int i = x - half; i <= x + half; i++) {
		int j = y;// - half;
		if (j >= 0 && j < height && i >= 0 && i < width) {
			buffer[j * width * bpp / 8 + i * bpp / 8 + 0] = r;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 1] = g;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 2] = b;
		}
/*		j = y + half;
		if (j >= 0 && j < height && i >= 0 && i < width) {
			buffer[j * width * bpp / 8 + i * bpp / 8 + 0] = r;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 1] = g;
			buffer[j * width * bpp / 8 + i * bpp / 8 + 2] = b;
		}
*/	}
}

void line(unsigned char* buffer, int width, int height, int bpp,
	  int x0, int y0, int x1, int y1,
	  int r, int g, int b, int a) {
	int dx = abs(x1 - x0);
	int dy = abs(y1 - y0);
	int sx = x0 < x1 ? 1 : -1;
	int sy = y0 < y1 ? 1 : -1;
	int err = dx - dy;
	int e2;

	do {
		if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) {
			buffer[y0 * width * bpp / 8 + x0 * bpp / 8 + 0] = r;
			buffer[y0 * width * bpp / 8 + x0 * bpp / 8 + 1] = g;
			buffer[y0 * width * bpp / 8 + x0 * bpp / 8 + 2] = b;
			if (bpp == 32) buffer[y0 * width * bpp / 8 + x0 * bpp / 8 + 3] = a;
		}
		if (x0 == x1 && y0 == y1) break;
		e2 = 2 * err;
		if (e2 > -dy) {
			err -= dy;
			x0 += sx;
		}
		if (e2 < dx) {
			err += dx;
			y0 += sy;
		}
	} while (true);
}

void setupOrtho(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, width, height, 0.0f, -1.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

SDL_Surface* mySDLInit(const int WIDTH, const int HEIGHT, const int BPP, const bool fullscreen) {
	SDL_Init(SDL_INIT_EVERYTHING);

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE,     8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,   8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,    8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,   8);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,  16);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
//	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
//	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);

	SDL_Surface *s = SDL_SetVideoMode(WIDTH, HEIGHT, BPP, (fullscreen ? SDL_FULLSCREEN : 0) | SDL_HWSURFACE | SDL_OPENGL);

	glewInit();

	glViewport(0, 0, WIDTH, HEIGHT);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0f);

	//glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	//glEnable(GL_MULTISAMPLE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//std::cout << glGetString(GL_VERSION)<< std::endl;
	//std::cout << glGetString(GL_SHADING_LANGUAGE_VERSION)<< std::endl;
	//std::cout << glewGetString(GLEW_VERSION)<< std::endl;
	//std::cout << glGetString(GL_EXTENSIONS)<< std::endl;

	return s;
}

void setupTexture(GLuint& texture) {
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void setupTexture(GLuint& texture, SDL_Surface *s) {
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, s->w, s->h, 0, s->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, s->pixels);
}

void setupTextureFloat(GLuint& texture, int w, int h, const float texturef[]) {
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_FLOAT, texturef);
}

void setupTextureFloat32(GLuint& texture, int w, int h, const float texturef[]) {
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_LUMINANCE, GL_FLOAT, texturef);
}

void setupTextureRGBA(GLuint& texture, int w, int h, const unsigned char texture_[]) {
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, texture_);
}

void setupTextureRGB(GLuint& texture, int w, int h, const unsigned char texture_[]) {
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_BGR, GL_UNSIGNED_BYTE, texture_);
}

void setupTextureImage(GLuint& texture, int w, int h, int bpp, const unsigned char texture_[]) {
	if      (bpp == 32) setupTextureRGBA(texture, w, h, texture_);
	else if (bpp == 24) setupTextureRGB(texture, w, h, texture_);
}

void setupTextureTGA(GLuint& texture, const char* tga_file, unsigned char*& buffer, int& width, int& height, int& bpp) {
	std::fstream inf(tga_file, std::ios_base::in | std::ios_base::binary);
	unsigned char header[18];
	inf.read((char *)header, 18);
	width  = header[12] + (header[13] << 8);
	height = header[14] + (header[15] << 8);
	bpp    = header[16];
	buffer = new unsigned char[width*height*bpp/8];
	std::cout<<(long long)(buffer)<<std::endl;
	inf.read((char *)buffer, width*height*bpp/8);
	std::cout << width << " " << height << " " << bpp << std::endl;
	setupTexture(texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, bpp == 32 ? GL_BGRA : GL_BGR, GL_UNSIGNED_BYTE, buffer);
}

void deleteTextureTGA(GLuint& texture, unsigned char*& buffer) {
	delete [] buffer;
	glDeleteTextures(1, &texture);
}

void deleteTexture(GLuint& texture) {
	glDeleteTextures(1, &texture);
}

void setupCubeMap(GLuint& texture) {
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_CUBE_MAP);
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_CUBE_MAP, texture);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}

void setupCubeMap(GLuint& texture, SDL_Surface *xpos, SDL_Surface *xneg, SDL_Surface *ypos, SDL_Surface *yneg, SDL_Surface *zpos, SDL_Surface *zneg) {
	setupCubeMap(texture);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGBA, xpos->w, xpos->h, 0, xpos->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, xpos->pixels);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGBA, xneg->w, xneg->h, 0, xneg->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, xneg->pixels);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGBA, ypos->w, ypos->h, 0, ypos->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, ypos->pixels);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGBA, yneg->w, yneg->h, 0, yneg->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, yneg->pixels);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGBA, zpos->w, zpos->h, 0, zpos->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, zpos->pixels);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGBA, zneg->w, zneg->h, 0, zneg->format->BytesPerPixel == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, zneg->pixels);
}

void deleteCubeMap(GLuint& texture) {
	glDeleteTextures(1, &texture);
}

void createProgram(GLuint& glProgram, GLuint& glShaderV, GLuint& glShaderF, const char* vertex_shader, const char* fragment_shader) {
	glShaderV = glCreateShader(GL_VERTEX_SHADER);
	glShaderF = glCreateShader(GL_FRAGMENT_SHADER);
	const GLchar* vShaderSource = loadFile(vertex_shader);
	const GLchar* fShaderSource = loadFile(fragment_shader);
	glShaderSource(glShaderV, 1, &vShaderSource, NULL);
	glShaderSource(glShaderF, 1, &fShaderSource, NULL);
	delete [] vShaderSource;
	delete [] fShaderSource;
	glCompileShader(glShaderV);
	glCompileShader(glShaderF);
	glProgram = glCreateProgram();
	glAttachShader(glProgram, glShaderV);
	glAttachShader(glProgram, glShaderF);
	glLinkProgram(glProgram);
	glUseProgram(glProgram);

	int  vlength,    flength,    plength;
	char vlog[2048], flog[2048], plog[2048];
	glGetShaderInfoLog(glShaderV, 2048, &vlength, vlog);
	glGetShaderInfoLog(glShaderF, 2048, &flength, flog);
	glGetProgramInfoLog(glProgram, 2048, &flength, plog);
	std::cout << vlog << std::endl << std::endl << flog << std::endl << std::endl << plog << std::endl << std::endl;
}

void releaseProgram(GLuint& glProgram, GLuint glShaderV, GLuint glShaderF) {
	glDetachShader(glProgram, glShaderF);
	glDetachShader(glProgram, glShaderV);
	glDeleteShader(glShaderF);
	glDeleteShader(glShaderV);
	glDeleteProgram(glProgram);
}

void saveTGA(unsigned char* buffer, int width, int height, bool video) {
	static int i = 0;
	std::stringstream out;
	if (video) {
		if      (i < 10)
			out << "video000" << (i++) << ".tga";
		else if (i < 100)
			out << "video00" << (i++) << ".tga";
		else if (i < 1000)
			out << "video0" << (i++) << ".tga";
		else if (i < 10000)
			out << "video" << (i++) << ".tga";
	} else {
		if      (i < 10)
			out << "media/debug/capture000" << (i++) << ".tga";
		else if (i < 100)
			out << "media/debug/capture00" << (i++) << ".tga";
		else if (i < 1000)
			out << "media/debug/capture0" << (i++) << ".tga";
		else if (i < 10000)
			out << "media/debug/capture" << (i++) << ".tga";
	}
	std::string s = out.str();
	
	glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, buffer);
	std::fstream of(s.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	char header[18] = { 0 };
	header[2] = 2;
	header[12] = width & 0xff;
	header[13] = width >> 8;
	header[14] = height & 0xff;
	header[15] = height >> 8;
	header[16] = 32;
	of.write(header, 18);
	of.write((char *)buffer, width * height * 4);
}

void saveTGARGBA(unsigned char* buffer, int width, int height, bool video) {
	static int i = 0;
	std::stringstream out;
	if (video) {
		if      (i < 10)
			out << "video000" << (i++) << ".tga";
		else if (i < 100)
			out << "video00" << (i++) << ".tga";
		else if (i < 1000)
			out << "video0" << (i++) << ".tga";
		else if (i < 10000)
			out << "video" << (i++) << ".tga";
	} else {
		if      (i < 10)
			out << "media/debug/capture000" << (i++) << ".tga";
		else if (i < 100)
			out << "media/debug/capture00" << (i++) << ".tga";
		else if (i < 1000)
			out << "media/debug/capture0" << (i++) << ".tga";
		else if (i < 10000)
			out << "media/debug/capture" << (i++) << ".tga";
	}
	std::string s = out.str();
	
	std::fstream of(s.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	char header[18] = { 0 };
	header[2] = 2;
	header[12] = width & 0xff;
	header[13] = width >> 8;
	header[14] = height & 0xff;
	header[15] = height >> 8;
	header[16] = 32;
	of.write(header, 18);
	of.write((char *)buffer, width * height * 4);
}

void saveTGADouble(double* buffer, int width, int height) {
	static int i = 0;
	std::stringstream out;
	
		if      (i < 10)
			out << "media/debug/image000" << (i++) << ".tga";
		else if (i < 100)
			out << "media/debug/image00" << (i++) << ".tga";
		else if (i < 1000)
			out << "media/debug/image0" << (i++) << ".tga";
		else if (i < 10000)
			out << "media/debug/image" << (i++) << ".tga";

	std::string s = out.str();

	std::fstream of(s.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

	unsigned char *buffer_ = new unsigned char[width * height * 4];
	for (int j = 0; j < width; j++) {
		for (int i = 0; i < width; i++) {
			buffer_[j * (width * 4) + i * 4 + 0] = (unsigned char)(buffer[j * width + i]);
			buffer_[j * (width * 4) + i * 4 + 1] = (unsigned char)(buffer[j * width + i]);
			buffer_[j * (width * 4) + i * 4 + 2] = (unsigned char)(buffer[j * width + i]);
			buffer_[j * (width * 4) + i * 4 + 3] = 255;
		}
	}
	
	char header[18] = { 0 };
	header[2] = 2;
	header[12] = width & 0xff;
	header[13] = width >> 8;
	header[14] = height & 0xff;
	header[15] = height >> 8;
	header[16] = 32;
	of.write(header, 18);
	of.write((char *)buffer_, width * height * 4);

	delete [] buffer_;
}

// saves the texture int the portable pixmap format
void savePPM(unsigned char* buffer, int width, int height, int bpp, const char* image_out) {
	std::fstream of(image_out, std::ios_base::out | std::ios_base::trunc);
	of << "P3" << std::endl;
	of << width << " " << height << std::endl;
	of << "255" << std::endl;
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			of << (int)buffer[j * width * bpp / 8 + i * bpp / 8 + 0] << " "
			   << (int)buffer[j * width * bpp / 8 + i * bpp / 8 + 1] << " "
			   << (int)buffer[j * width * bpp / 8 + i * bpp / 8 + 2] << " ";
		}
		of << std::endl;
	}
}

void glerror(const char* prepend) {
	std::cout << prepend << " -- " << gluErrorString(glGetError()) << std::endl;
}