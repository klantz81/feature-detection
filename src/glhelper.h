#ifndef GLHELPER_H
#define GLHELPER_H

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>
#include <SDL/SDL_image.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "misc.h"

SDL_Surface* mySDLInit(const int WIDTH, const int HEIGHT, const int BPP, const bool fullscreen);

unsigned char* loadImage(const char* image, int& width, int& height, int& bpp);
double* grayscaleConversion(unsigned char* buffer, int width, int height, int bpp);
float* doubleToFloat(double* buffer, int width, int height);
void scharrOperator(double*& Ix, double*& Iy, double*& I, double* buffer, int width, int height);
double* gaussianConvolution(double* buffer, int width, int height, int kernel_radius, double mean, double std_deviation);
void structureTensor(double*& Ixx, double*& Ixy, double*& Iyy, double* Ix, double* Iy, int width, int height, int kernel_radius, double mean, double std_deviation);
unsigned char* downsample(unsigned char* buffer, int width, int height, int bpp);
double* downsample(double* buffer, int width, int height);
float* downsample(float* buffer, int width, int height);
double interpolate(double* buffer, int width, int height, double x, double y);

void drawBox(unsigned char* buffer, int width, int height, int bpp, int x, int y, int half, unsigned char r, unsigned char g, unsigned char b);
void line(unsigned char* buffer, int width, int height, int bpp, int x0, int y0, int x1, int y1, int r, int g, int b, int a);

void setupOrtho(int width, int height);

void setupTexture(GLuint& texture);
void setupTexture(GLuint& texture, SDL_Surface *s);
void setupTextureFloat(GLuint& texture, int w, int h, const float texturef[]);
void setupTextureFloat32(GLuint& texture, int w, int h, const float texturef[]);
void setupTextureRGBA(GLuint& texture, int w, int h, const unsigned char texture_[]);
void setupTextureRGB(GLuint& texture, int w, int h, const unsigned char texture_[]);
void setupTextureImage(GLuint& texture, int w, int h, int bpp, const unsigned char texture_[]);
void setupTextureTGA(GLuint& texture, const char* tga_file, unsigned char*& buffer, int& width, int& height, int& bpp);
void deleteTexture(GLuint& texture);
void deleteTextureTGA(GLuint& texture, unsigned char*& buffer);

void setupCubeMap(GLuint& texture);
void setupCubeMap(GLuint& texture, SDL_Surface *xpos, SDL_Surface *xneg, SDL_Surface *ypos, SDL_Surface *yneg, SDL_Surface *zpos, SDL_Surface *zneg);
void deleteCubeMap(GLuint& texture);

void createProgram(GLuint& glProgram, GLuint& glShaderV, GLuint& glShaderF, const char* vertex_shader, const char* fragment_shader);
void releaseProgram(GLuint& glProgram, GLuint glShaderV, GLuint glShaderF);

void saveTGA(unsigned char* buffer, int width, int height, bool video = false);
void saveTGARGBA(unsigned char* buffer, int width, int height, bool video);
void saveTGADouble(double* buffer, int width, int height);
void savePPM(unsigned char* buffer, int width, int height, int bpp, const char* image_out);

void glerror(const char* prepend);

#endif
