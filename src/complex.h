#ifndef COMPLEX_H
#define COMPLEX_H

class cComplex {
  private:
  protected:
  public:
    float a, b;
    static unsigned int additions, multiplications;
    cComplex();
    cComplex(float a, float b);
    cComplex conj();
    cComplex operator*(const cComplex& c) const;
    cComplex operator+(const cComplex& c) const;
    cComplex operator-(const cComplex& c) const;
    cComplex operator-() const;
    cComplex operator*(const float c) const;
    cComplex& operator=(const cComplex& c);
    static void reset();
};

#endif