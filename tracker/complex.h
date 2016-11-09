#pragma once
#include "float_type.h"

struct Complex {
  Tfloat real;
  Tfloat imag;  

  __device__ __host__
  Complex(Tfloat real, Tfloat imag = 0.): real(real), imag(imag) {}

  __device__ __host__
  Complex & operator=(const Complex & rhs) {
    real = rhs.real;
    imag = rhs.imag;
    return *this;
  }
  
  __device__ __host__
  Complex & operator+=(const Complex & rhs) {
    real += rhs.real;
    imag += rhs.imag;
    return *this;
  }

  __device__ __host__
  Complex & operator-=(const Complex & rhs) {
    real -= rhs.real;
    imag -= rhs.imag;
    return *this;
  }

  __device__ __host__
  Complex operator*=(const Complex & rhs) {
    const Tfloat R = real;
    const Tfloat I = imag;
    real = R*rhs.real - I*rhs.imag;
    imag = R*rhs.imag + I*rhs.real;
    return *this;
  }

  __device__ __host__
  Complex & operator/=(const Tfloat & rhs) {
    real /= rhs;
    imag /= rhs;
    return *this;
  }
};

/////////////////

__device__ __host__
inline Complex operator+(Complex lhs, const Complex rhs) {
  lhs += rhs;
  return lhs;
}

__device__ __host__
inline Complex operator-(Complex lhs, const Complex rhs) {
  lhs -= rhs;
  return lhs;
}

__device__ __host__
inline Complex operator*(Complex lhs, const Complex rhs) {
  lhs *= rhs;
  return lhs;
}

__device__ __host__
inline Complex operator/(Complex lhs, const Tfloat & rhs) {
  lhs /= rhs;
  return lhs;
}

__device__ __host__
Complex pow(const Complex & rhs, int t) {
//  if ( t == 0 ) return Complex(1.,0.);
  Complex res(rhs);
  while ( t > 1 ) {
    res *= rhs;
    --t;
  }
  return res;
}

