#include "complex.h"

unsigned int cComplex::additions = 0;
unsigned int cComplex::multiplications = 0;

cComplex::cComplex() : a(0.0f), b(0.0f) { }
cComplex::cComplex(float a, float b) : a(a), b(b) { }
cComplex cComplex::conj() { return cComplex(this->a, -this->b); }

cComplex cComplex::operator*(const cComplex& c) const {
	cComplex::multiplications++;
	return cComplex(this->a*c.a - this->b*c.b, this->a*c.b + this->b*c.a);
}

cComplex cComplex::operator+(const cComplex& c) const {
	cComplex::additions++;
	return cComplex(this->a + c.a, this->b + c.b);
}

cComplex cComplex::operator-(const cComplex& c) const {
	cComplex::additions++;
	return cComplex(this->a - c.a, this->b - c.b);
}

cComplex cComplex::operator-() const {
	return cComplex(-this->a, -this->b);
}

cComplex cComplex::operator*(const float c) const {
	return cComplex(this->a*c, this->b*c);
}

cComplex& cComplex::operator=(const cComplex& c) {
	this->a = c.a; this->b = c.b;
	return *this;
}

void cComplex::reset() {
	cComplex::additions = 0;
	cComplex::multiplications = 0;
}