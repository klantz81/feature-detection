#include <sstream>
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
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "glm-0.9.2.6/glm/glm.hpp"
#include "glm-0.9.2.6/glm/gtc/matrix_transform.hpp"
#include "glm-0.9.2.6/glm/gtc/type_ptr.hpp"

#include "src/timer.h"
#include "src/keyboard.h"
#include "src/interrupt.h"
#include "src/buffer.h"
#include "src/glhelper.h"
#include "src/misc.h"
#include "src/vector.h"
#include "src/matrix.h"
#include "src/complex.h"
#include "src/kernel.h"
#include "src/patch.h"
#include "src/defines.h"
#include "src/buffers.h"



struct feature_ {
	bool alive;
	int first_frame;

	double lambda0, lambda1;
  
	std::vector<vector2> locations;

	// translation
	matrix22 H;
	matrix22 Hinv;
	double Hdet;

	// affine
	cMatrix H_;
	cMatrix H_inv;
	double H_det;
	cMatrix o; // last estimate
	bool H_singular;

	cDoubleBuffer T_buffer, Tx_buffer, Ty_buffer;

	bool operator()(feature_ a, feature_ b) { return MIN(a.lambda0, a.lambda1) > MIN(b.lambda0, b.lambda1); }
};

class buffer_ {
private:
protected:
public:
	int width, height, bpp;

	cColorBuffer color_buffer;
	cDoubleBuffer grayscale_buffer;
	cFloatBuffer grayscale_buffer_float;
	cDoubleBuffer Ix_buffer, Iy_buffer, I_buffer;
	cDoubleBuffer Ixx_buffer, Ixy_buffer, Iyy_buffer;

	buffer_(const char* image, const int levels, const int extract_features, const int frame, std::vector<feature_>& features_, const cKernel& kernel);
	~buffer_();
};

buffer_::buffer_(const char* image, const int levels, const int extract_features, const int frame, std::vector<feature_>& features_, const cKernel& kernel) {

	this->color_buffer.load(image);
	this->width  = this->color_buffer.width;
	this->height = this->color_buffer.height;
	this->bpp    = this->color_buffer.bpp;

	this->grayscale_buffer       = this->color_buffer.grayscale();
	this->grayscale_buffer_float = this->grayscale_buffer.toFloat(1.0);

	this->Ix_buffer              = this->grayscale_buffer.scharrOperatorX();
	this->Iy_buffer              = this->grayscale_buffer.scharrOperatorY();

	if (extract_features > 0) {
		double min_feature_distance = kernel.kernel_radius;// / 2.0;
		double min_eigenvalue = 1000.0;// 150.0

		std::vector<feature_> features;

		features.reserve(this->width * this->height);

		for (int j = (kernel.kernel_radius + 1) * 2; j < this->height - (kernel.kernel_radius + 1) * 2; j++) {
			for (int i = (kernel.kernel_radius + 1) * 2; i < this->width - (kernel.kernel_radius + 1) * 2; i++) {
				int index = j * this->width + i;

				vector2 l;
				l.x = i;
				l.y = j;

				feature_ f;
				f.locations.push_back(l);
				f.alive = false;
				f.first_frame = frame;

				double Hxx = 0, Hxy = 0, Hyy = 0;
				for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++) {
					for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++) {
						double factor = 1.0;
						double Tx = this->Ix_buffer.at(i + _i, j + _j);
						double Ty = this->Iy_buffer.at(i + _i, j + _j);
						Hxx += Tx * Tx * factor;
						Hxy += Tx * Ty * factor;
						Hyy += Ty * Ty * factor;
					}
				}
				matrix22 m(Hxx, Hxy, Hxy, Hyy);

				m.eigenvalues(f.lambda0, f.lambda1);

				if (MIN(f.lambda0, f.lambda1) > min_eigenvalue) {
					f.H    = m;
					f.Hdet = m.det();
					f.Hinv = m.inv();

					features.push_back(f);
				}
			}
		}

		std::cout << features.size() << " potential features" << std::endl;
		if (features.size() > 0) std::sort(features.begin(), features.end(), features[0]);

		int feature_count = 0;
		for (int i = 0; i < features.size(); i++) {
			bool flag = true;

			for (int j = 0; j < features_.size(); j++) {
				if (features_[j].alive && sqrt(pow(features[i].locations[ features[i].locations.size() - 1 ].x - features_[j].locations[ features_[j].locations.size() - 1 ].x,2.0) +
								pow(features[i].locations[ features[i].locations.size() - 1 ].y - features_[j].locations[ features_[j].locations.size() - 1 ].y,2.0)) < min_feature_distance) {
					flag = false;
					break;
				}
			}


			if (flag) {
				features[i].alive = true;

				features[i].T_buffer  = this->grayscale_buffer.sub(features[i].locations[0].x, features[i].locations[0].y, kernel);
				features[i].Tx_buffer = this->Ix_buffer.sub       (features[i].locations[0].x, features[i].locations[0].y, kernel);
				features[i].Ty_buffer = this->Iy_buffer.sub       (features[i].locations[0].x, features[i].locations[0].y, kernel);
//				features[i].T_buffer.save();

				features[i].o = cMatrix(6,1);
				features[i].o.zero();

				// evalute affine hessian
				cMatrix H_(6,6);
				H_.zero();

				for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++) {
					for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++) {
						cMatrix J_(1,6);
						J_.zero();
						
						double factor = 1.0;//kernel.kernel[_j + kernel.kernel_radius] * kernel.kernel[_i + kernel.kernel_radius];
						double x = _i;
						double y = _j;
						double Tx = features[i].Tx_buffer.at(_i + kernel.kernel_radius, _j + kernel.kernel_radius);//[(_j + kernel.kernel_radius) * (kernel.kernel_radius * 2 + 1) + (_i + kernel.kernel_radius)];
						double Ty = features[i].Ty_buffer.at(_i + kernel.kernel_radius, _j + kernel.kernel_radius);//[(_j + kernel.kernel_radius) * (kernel.kernel_radius * 2 + 1) + (_i + kernel.kernel_radius)];
//						double T  = features[i].T_buffer.buffer[(_j + kernel.kernel_radius) * (kernel.kernel_radius * 2 + 1) + (_i + kernel.kernel_radius)];
						
						J_.A[0] = x * Tx;
						J_.A[1] = y * Tx;
						J_.A[2] = x * Ty;
						J_.A[3] = y * Ty;
						J_.A[4] =     Tx;
						J_.A[5] =     Ty;
//						J_.A[6] =     T;
//						J_.A[7] =     1;

						H_ = H_ + ((J_.transpose() * J_) * factor);
					}
				}

				features[i].H_    = H_;
				features[i].H_det = H_.det();
				features[i].H_inv = H_.gaussJordanElimination(features[i].H_singular);

				features_.push_back(features[i]);

				feature_count++;
			}
			if (feature_count >= extract_features) break;
		}
	}
}

buffer_::~buffer_() { }

int trackFeatures(const buffer_& b0, const buffer_& b1, const int frame, std::vector<feature_>& features, const cKernel& kernel, cInterrupt& interrupt) {

	// stopping criteria: minimum determinant (Hessian not invertible), out of bounds, delta p too large after maximum iterations, residue too large
	const double min_determinant = 0.0001;
	const double min_delta_p = 0.0001;
	const int max_iterations = 2000;
	const double max_residue = (kernel.kernel_radius * 2 + 1) * (kernel.kernel_radius * 2 + 1) * 0.1;//5;

	int count = 0;
	for (int i = 0; i < features.size() && !interrupt.interrupted; i++) {
		if (!features[i].alive || features[i].first_frame == frame) continue; // exclude dead features or features that were detected this frame

		bool feature_lost = false;

		if (fabs(features[i].Hdet) < min_determinant) {
			std::cout << "feature lost: MIN DETERMINANT" << std::endl;
			feature_lost = true;
		}

		vector2  x = features[i].locations[ features[i].locations.size() - 1]; // last location
		vector2 _x = x;

		if (!feature_lost) {

			for (int j = 0; j < max_iterations && !feature_lost && !interrupt.interrupted; j++) {
				if (_x.x - kernel.kernel_radius >= 0 && _x.x + kernel.kernel_radius < b0.width && _x.y - kernel.kernel_radius >= 0 && _x.y + kernel.kernel_radius < b0.height) {
					vector2 s(0.0, 0.0);
					double residue = 0.0;
					for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++) {
						for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++) {
							double factor = 1.0;//(kernel.kernel[_j + kernel.kernel_radius] * kernel.kernel[_i + kernel.kernel_radius]);
							double Tx = b0.Ix_buffer.interpolate(x.x + _i, x.y + _j);
							double Ty = b0.Iy_buffer.interpolate(x.x + _i, x.y + _j);
							double T =  b0.grayscale_buffer.interpolate( x.x + _i,  x.y + _j);
							double I =  b1.grayscale_buffer.interpolate(_x.x + _i, _x.y + _j);

							vector2 _s = vector2(Tx, Ty) * (I - T) * factor;
							s = s + _s;

							residue += fabs(I - T);
						}
					}
					vector2 delta_p = features[i].Hinv * s;
					_x = _x - delta_p;
					if (delta_p.length() < min_delta_p) {
						if (fabs(residue) > max_residue) {
							std::cout << "feature lost: MAX RESIDUE" << std::endl;
							feature_lost = true;
						}
						break;
					} else if (j == max_iterations - 1) {
						std::cout << "feature lost: MAX ITERATIONS" << std::endl;
						feature_lost = true;
					}
				} else {
					std::cout << "feature lost: OUT OF BOUNDS" << std::endl;
					feature_lost = true;
				}
			}
		}

		if (!feature_lost) {//&& i < 2) {
			// affine consistency check
			vector2 x = features[i].locations[0]; // first location
			cMatrix o = features[i].o;
			if (features[i].locations.size() == 1) {
				o.A[4] = x.x;
				o.A[5] = x.y;
			}

			if (features[i].locations.size() > 0) { // only compare a feature with a past frame
				if (fabs(features[i].H_det) < min_determinant) {
					std::cout << "feature lost: MIN DETERMINANT (affine)" << std::endl;
					feature_lost = true;
				}
				if (!feature_lost) {
					for (int j = 0; j < max_iterations && !feature_lost && !interrupt.interrupted; j++) {
						cMatrix s(1,6);
						s.zero();
						double residue = 0;

						matrix22 A(1.0 + o.A[0],       o.A[1],
							         o.A[2], 1.0 + o.A[3]);
						vector2  b(o.A[4],
							   o.A[5]);
//						double alpha = o.A[6];
//						double beta  = o.A[7];

						for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius && !feature_lost; _j++) {
							for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius && !feature_lost; _i++) {

								vector2 x__(_i, _j);
								vector2 x_ = A * x__ + b;

								if (x_.x < 0 || x_.x > b1.width - 1 || x_.y < 0 || x_.y > b1.height - 1) {
									std::cout << "feature lost: OUT OF BOUNDS (affine) : " << x_.x << "," << x_.y << std::endl;
									feature_lost = true;
									break;
								}

								double factor = 1.0;//kernel.kernel[_j + kernel.kernel_radius] * kernel.kernel[_i + kernel.kernel_radius];

								double Tx = features[i].Tx_buffer.at(_i + kernel.kernel_radius, _j + kernel.kernel_radius);
								double Ty = features[i].Ty_buffer.at(_i + kernel.kernel_radius, _j + kernel.kernel_radius);
								double T =  features[i].T_buffer.at (_i + kernel.kernel_radius, _j + kernel.kernel_radius);
								double I =  b1.grayscale_buffer.interpolate(x_.x, x_.y);
//								double I          = (1.0 + alpha) * interpolate(b1.grayscale, b1.width, b1.height, x_.x, x_.y) + beta;

								s.A[0] +=  x__.x * Tx * (I - T) * factor;
								s.A[1] +=  x__.y * Tx * (I - T) * factor;
								s.A[2] +=  x__.x * Ty * (I - T) * factor;
								s.A[3] +=  x__.y * Ty * (I - T) * factor;
								s.A[4] +=          Tx * (I - T) * factor;
								s.A[5] +=          Ty * (I - T) * factor;
								//s.A[6] +=          T  * (I - T) * factor;
								//s.A[7] +=         1.0 * (I - T) * factor;

								residue += fabs(I - T);
							}
						}
							
						cMatrix delta_o = features[i].H_inv * s.transpose();

						matrix22 dA = matrix22(1.0 + delta_o.A[0],       delta_o.A[1],
								             delta_o.A[2], 1.0 + delta_o.A[3]);

						if (fabs(dA.det()) < min_determinant) {
							std::cout << "feature lost: MIN DETERMINANT (affine)" << std::endl;
							feature_lost = true;
						}
						
						if (!feature_lost) {

							matrix22 dAinv = dA.inv();

							matrix22 AdAinv = A * dAinv;

							vector2 db(delta_o.A[4], delta_o.A[5]);
							vector2 AdAinvdb = AdAinv * db;

							cMatrix tmp = o;

							o.A[0]  = AdAinv.a - 1.0;
							o.A[1]  = AdAinv.b;
							o.A[2]  = AdAinv.c;
							o.A[3]  = AdAinv.d - 1.0;
							o.A[4] -= AdAinvdb.x;
							o.A[5] -= AdAinvdb.y;
//							double f = (1 + o.A[6]) / (1 + delta_o.A[6]); 
//							o.A[6]  = f - 1.0;
//							o.A[7] -= f * delta_o.A[7];
							
							delta_o = tmp - o;
							
							double delta =
							   sqrt(delta_o.A[0]*delta_o.A[0] +
								delta_o.A[1]*delta_o.A[1] +
								delta_o.A[2]*delta_o.A[2] +
								delta_o.A[3]*delta_o.A[3] +
								delta_o.A[4]*delta_o.A[4] +
								delta_o.A[5]*delta_o.A[5]);// +
//								delta_o.A[6]*delta_o.A[6] +
//								delta_o.A[7]*delta_o.A[7]);

							if (fabs(delta) > 1e7 || delta!=delta) {
								std::cout << "feature lost: MAX DELTA (affine)" << std::endl;
								feature_lost = true;
							}

							if (fabs(delta) < min_delta_p) {
								if (fabs(residue) > max_residue) {
									std::cout << "feature lost: MAX RESIDUE (affine)" << std::endl;
									feature_lost = true;
								}
								break;
							} else if (j == max_iterations - 1) {
								std::cout << "feature lost: MAX ITERATIONS (affine)" << std::endl;
								feature_lost = true;
							}
						}
					}
				}
			}

			if (!feature_lost) {

//std::cout << "----------------------------------------- feature kept -- locations: " << features[i].locations.size() << " + 1 " << std::endl;
std::stringstream s;
int k0 = i;
int k1 = features[i].locations.size();

if (false) {
	s.str("");
	s << "media_save/p" << (k0 < 10 ? "000" : (k0 < 100 ? "00" : "0")) << k0 << ".pgm";
	features[i].T_buffer.save(s.str().c_str());
}

if (false) {
	s.str("");
	s << "media_save/p" << (k0 < 10 ? "000" : (k0 < 100 ? "00" : "0")) << k0 << "_translation_" << (k1 < 10 ? "000" : (k1 < 100 ? "00" : "0")) << k1 << ".ppm";
	//(b1.grayscale_buffer.subInterpolate(_x.x, _x.y, kernel)).save(s.str().c_str());
	(b1.color_buffer.subInterpolate(_x.x, _x.y, kernel)).save(s.str().c_str());
}

if (false) {
	s.str("");
	s << "media_save/p" << (k0 < 10 ? "000" : (k0 < 100 ? "00" : "0")) << k0 << "_affine_" << (k1 < 10 ? "000" : (k1 < 100 ? "00" : "0")) << k1 << ".ppm";
	matrix22 A(1.0 + o.A[0], o.A[1], o.A[2], 1.0 + o.A[3]);
	vector2 b(o.A[4], o.A[5]);
	//(b1.grayscale_buffer.subInterpolate(A, b, kernel)).save(s.str().c_str());
	(b1.color_buffer.subInterpolate(A, b, kernel)).save(s.str().c_str());
}

				features[i].locations.push_back(_x);

				// evaluate Hessian using b1 centered at _x
				double Hxx = 0, Hxy = 0, Hyy = 0;
				for (int _j = -kernel.kernel_radius; _j <= kernel.kernel_radius; _j++) {
					for (int _i = -kernel.kernel_radius; _i <= kernel.kernel_radius; _i++) {
						double factor = 1.0;//(kernel.kernel[_j + kernel.kernel_radius] * kernel.kernel[_i + kernel.kernel_radius]);
						double Tx = b1.Ix_buffer.interpolate(_x.x + _i, _x.y + _j);
						double Ty = b1.Iy_buffer.interpolate(_x.x + _i, _x.y + _j);
						Hxx += Tx * Tx * factor;
						Hxy += Tx * Ty * factor;
						Hyy += Ty * Ty * factor;
					}
				}

				matrix22 m(Hxx, Hxy, Hxy, Hyy);
				features[i].H    = m;
				features[i].Hdet = m.det();
				features[i].Hinv = m.inv();
				features[i].o = o;
			} else features[i].alive = false;
		} else features[i].alive = false;
	}

	int features_alive = 0;
	for (int i = 0; i < features.size(); i++) if (features[i].alive) features_alive++;
	return features_alive;
}

struct vertex_ {
	double x,  y,   z;
	double tx, ty, tz;
	double t[6];
};

struct dimensions_ {
	int width, height, bpp;
};

std::string getImage() {
	static int frame = 0;
	frame++;
	std::stringstream s;
	s << "media/img" << (frame < 10 ? "000" : (frame < 100 ? "00" : (frame < 1000 ? "0" : ""))) << frame << ".png";
	std::cout << "getImage(): " << s.str() << std::endl;
	return s.str();
}

int main(int argc, char* argv[]) {
  
	dimensions_ dim = {1920, 1080, 32};
//	dimensions_ dim = {640, 360, 32};
	cKeyboard kb;
	cInterrupt interrupt;

	cKernel kernel(12, 0.0, 4);

	buffer_ *buffer0, *buffer1;
	std::vector<feature_> features;
	int frame = 0, alive = 0, min = 16, max = 32, max_frame = 24;//21;//21;//10;//356;

	buffer0 = new buffer_(getImage().c_str(), 1, max, ++frame, features, kernel);
	alive = features.size();
	
	dim.width = buffer0->width;
	dim.height = buffer0->height;
	cBuffer cbuffer(dim.width, dim.height);

for (int i = 0; i < features.size(); i++) {
  //std::cout << "min eigenvalue: " << MIN(features[i].lambda0, features[i].lambda1) << std::endl;
  //features[i].T.save();
  //features[i].Tx.save();
  //features[i].Ty.save();
}

	cTimer t0, t1;
	double elapsed0, elapsed1;
  
	timespec current;
	clock_gettime(CLOCK_REALTIME, &current);
	srand(current.tv_sec + current.tv_nsec / 1000000000.0);
	
	bool active = true;

	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_Surface *screen = mySDLInit(dim.width, dim.height, dim.bpp, false);
	SDL_Event event;

	// shaders
	GLuint glProgram, glShaderV, glShaderF;
	GLint vertex, texCoord, Projection, View, Model, textureSample;
	createProgram(glProgram, glShaderV, glShaderF, "src/vertex.sh", "src/fragment.sh");
	vertex        = glGetAttribLocation(glProgram, "vertex");
	texCoord      = glGetAttribLocation(glProgram, "texCoord");
	Projection    = glGetUniformLocation(glProgram, "Projection");
	View          = glGetUniformLocation(glProgram, "View");
	Model         = glGetUniformLocation(glProgram, "Model");
	textureSample = glGetUniformLocation(glProgram, "textureSample");

	// vertices
	const int vertex_count = 4;
	vertex_ vertices[vertex_count];
	vertices[0].x  = 0;         vertices[0].y  = 0;          vertices[0].z = 0;     vertices[0].tx = 0;     vertices[0].ty = 0;	
	vertices[1].x  = dim.width; vertices[1].y  = 0;          vertices[1].z = 0;     vertices[1].tx = 1;     vertices[1].ty = 0;	
	vertices[2].x  = dim.width; vertices[2].y  = dim.height; vertices[2].z = 0;     vertices[2].tx = 1;     vertices[2].ty = 1;
	vertices[3].x  = 0;         vertices[3].y  = dim.height; vertices[3].z = 0;     vertices[3].tx = 0;     vertices[3].ty = 1;
	GLuint vbo_vertices;
	glGenBuffers(1, &vbo_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_) * vertex_count, vertices, GL_DYNAMIC_DRAW);

	// indices
	const int indices_count = 6;
	unsigned int indices[indices_count] = {0, 1, 2, 2, 3, 0};
	GLuint vbo_indices;
	glGenBuffers(1, &vbo_indices);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indices);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices_count, indices, GL_STATIC_DRAW);

	// create a texture for rendering the result
	GLuint texture;
	setupTexture(texture);

	glm::mat4 projection = glm::ortho(0.0f, (float)dim.width, (float)dim.height, 0.0f, -5.0f, 5.0f); 
	glm::mat4 view       = glm::mat4(1.0f);
	glm::mat4 model      = glm::mat4(1.0f);

//	int index = 0;

	while (active) {

		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			    case SDL_QUIT: active = false; break;
//			    case SDL_KEYDOWN:
//				switch (event.key.keysym.sym) {
//				    case SDLK_a: index = 0; break;
//				    case SDLK_s: index = 1; break;
//				    case SDLK_d: index = 2; break;
//				    case SDLK_f: index = 3; break;
//				}
//				break;
			}
		}
		if (interrupt.interrupted) active = false;
//		if (kb.getKeyState(KEY_A)) index = 0;
//		if (kb.getKeyState(KEY_S)) index = 1;
//		if (kb.getKeyState(KEY_D)) index = 2;
//		if (kb.getKeyState(KEY_F)) index = 3;

		elapsed0 = t0.elapsed(true);
		elapsed1 = t1.elapsed(false);

		glDisable(GL_DEPTH_TEST);
		glUseProgram(glProgram);

		// Project, View, and Model uniforms
		glUniformMatrix4fv(Projection, 1, GL_FALSE, glm::value_ptr(projection));
		glUniformMatrix4fv(View, 1, GL_FALSE, glm::value_ptr(view));
		glUniformMatrix4fv(Model, 1, GL_FALSE, glm::value_ptr(model));

		// vertex and texture attributes
		glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
		glEnableVertexAttribArray(vertex);
		glVertexAttribPointer(vertex, 3, GL_DOUBLE, GL_FALSE, sizeof(vertex_), 0);
		glEnableVertexAttribArray(texture);
		glVertexAttribPointer(texCoord, 2, GL_DOUBLE, GL_FALSE, sizeof(vertex_), (char *)NULL + 24);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indices);

		// texture
		glUniform1i(textureSample, 0);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);


		if (frame == 1) {
			for (int i = 0; i < features.size(); i++)
				if (features[i].alive)
					drawBox(buffer0->color_buffer.buffer, buffer0->width, buffer0->height, buffer0->bpp, features[i].locations[features[i].locations.size()-1].x, features[i].locations[features[i].locations.size()-1].y, kernel.kernel_radius, 0, 255, 0);

			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, buffer0->color_buffer.width, buffer0->color_buffer.height, 0, GL_RGB, GL_UNSIGNED_BYTE, buffer0->color_buffer.buffer);
			//this->grayscale_buffer_float = this->grayscale_buffer.toFloat();
			//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, buffer0->grayscale_buffer_float.width,
			//					   buffer0->grayscale_buffer_float.height, 0, GL_LUMINANCE, GL_FLOAT,
			//					   buffer0->grayscale_buffer_float.buffer);

			// draw it
			glDrawElements(GL_TRIANGLES, indices_count, GL_UNSIGNED_INT, 0);
			SDL_GL_SwapBuffers();
			
			cbuffer.save(frame);
		}

		if (frame < max_frame) {
			buffer1 = new buffer_(getImage().c_str(), 1, alive < min ? (max - alive) : 0, ++frame, features, kernel);
			alive = trackFeatures(*buffer0, *buffer1, frame, features, kernel, interrupt);
			std::cout << "features alive: " << alive << std::endl;

			delete buffer0;
			buffer0 = buffer1;

			for (int i = 0; i < features.size(); i++) {
				if (features[i].alive) {
					drawBox(buffer0->color_buffer.buffer, buffer0->width, buffer0->height, buffer0->bpp, features[i].locations[features[i].locations.size()-1].x, features[i].locations[features[i].locations.size()-1].y, kernel.kernel_radius, 0, 255, 0);
					for (int j = 0; j < features[i].locations.size() - 1; j++) {
						line(buffer0->color_buffer.buffer, buffer0->width, buffer0->height, buffer0->bpp,
						    features[i].locations[j].x,   features[i].locations[j].y,
						    features[i].locations[j+1].x, features[i].locations[j+1].y,
						    255, 0, 0, 255);
					}
				}
			}
		}


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, buffer0->width, buffer0->height, 0, GL_RGB, GL_UNSIGNED_BYTE, buffer0->color_buffer.buffer);

		// draw it
		glDrawElements(GL_TRIANGLES, indices_count, GL_UNSIGNED_INT, 0);
		SDL_GL_SwapBuffers();

		cbuffer.save(frame);
	}

	// release the texture, vertex buffer objects, and shader program.. also shut down SDL
	deleteTexture(texture);
	glDeleteBuffers(1, &vbo_indices);
	glDeleteBuffers(1, &vbo_vertices);
	releaseProgram(glProgram, glShaderV, glShaderF);
	SDL_Quit();

	delete buffer0;
//	delete buffer1;

	return 0;
}