#pragma once
#include "float_type.h"
#include "complex.h"

namespace cst {
  const double pi = 3.14159265;
  const double c  = 299792458.0;
  const double r0 = 2.8179402894e-15;
}

struct Particle {
  Tfloat& x;
  Tfloat& xp;
  Tfloat& y;
  Tfloat& yp;
  Tfloat& z;
  Tfloat& d;
};

struct Drift {
  Tfloat L;

  __host__ __device__ 
  Drift(const Tfloat L): L(L) {}
  
  __host__ __device__ 
  void operator()(Particle & p) {
    const Tfloat cxp = p.xp;
    const Tfloat cyp = p.yp;
    p.x += cxp * L;
    p.y += cyp * L;
    p.z += L * 0.5 * (cxp*cxp + cyp*cyp);
  }
};

// A thin dipole with weak focussing
struct Dipole {
  Tfloat angle;
  Tfloat aa_l; //used for the weak focussing

  __host__ __device__
  Dipole(const Tfloat A, const Tfloat L=0): angle(A), aa_l(A*A/L) {}

  __host__ __device__
  void operator()(Particle & p) {
    p.xp += angle * p.d / (p.d + 1.); //energy effect
    p.xp -= aa_l * p.x; // weak focussing
    p.z -= p.x * angle ;
  }
};

struct Quad {
  Tfloat K1L;

  __host__ __device__
  Quad(const Tfloat K1L): K1L(K1L) {}

  __host__ __device__
  void operator()(Particle & p) {
    const Tfloat S = K1L/(p.d + 1.);
    p.xp -= S * p.x;
    p.yp += S * p.y;
  }
};

struct Multipole {
  int order; 	
  Complex strength;

  __host__ __device__
  Multipole(int order, Tfloat KL, Tfloat KLs = 0.): order(order), strength(KL, KLs) {}

  __host__ __device__
  void operator()(Particle & p) {
    Complex k = strength / (p.d + 1.);
    if ( order == 0 ) { // no weak focussing here
      k -= strength;  
    } else {
      int factor=1;
      for (int i=1;i<order+1;i++)
           factor*=i;	      
      k *= pow(Complex(p.x, p.y), order);
      k=k/factor;
    }
    p.xp -= k.real;
    p.yp += k.imag;
  }
};

struct VKicker {
  double kick;

  __host__ __device__
  VKicker(double kick): kick(kick) {}

  __host__ __device__
  void operator()(Particle & p) {
    p.yp += kick/(p.d + 1.);
  }
};

struct HKicker {
  double kick;

  __host__ __device__
  HKicker(double kick): kick(kick) {}

  __host__ __device__
  void operator()(Particle & p) {
    p.xp += kick/(p.d + 1.);
  }
};

struct RF {
  Tfloat f, VE;

  __host__ __device__
  RF(const Tfloat f, const Tfloat VE): f(f),VE(VE) {}

  __host__ __device__
  void operator()(Particle & p) {
    const Tfloat K = (2 * cst::pi * f * (1e6)) /(cst::c);  
    p.d += (VE /1e3) * sin(K * p.z );
  }
};

struct BeamBeam {
  double n, gamma, sigma;

  __host__ __device__
  BeamBeam(size_t n, double gamma, double sigma): n(n), gamma(gamma), sigma(sigma) {}

  __host__ __device__
  void operator()(Particle & p) {  
    double r2 = p.x*p.x + p.y*p.y;
    p.xp += 2.*n*cst::r0/(gamma)* p.x/r2 * (1.-exp(-0.5*r2/(sigma*sigma)));
    p.yp += 2.*n*cst::r0/(gamma)* p.y/r2 * (1.-exp(-0.5*r2/(sigma*sigma)));
  }
};

/*struct Noise_element {
  double kick;
  double mean, sigma;
  __host__ __device__
  Noise_element(double mean, double sigma): mean(mean), sigma(sigma) {}

  __host__ __device__
  void operator()(Particle & p) {
    p.xp += kick;
  }
};*/

//__host__ __device__
//void Drift(const Tfloat L, Particles & p, const size_t t) {
//  const Tfloat cxp = p.xp[t];
//  const Tfloat cyp = p.yp[t];
//  p.x[t] += cxp * L;
//  p.y[t] += cyp * L;
//  p.z[t] += L * 0.5 * (cxp*cxp + cyp*cyp);
//}
//
//
//__host__ __device__
//void Multipole(const int order, const Tfloat KL, const Tfloat KSL, Particles & p, const size_t t) {
//  Complex k(KL,KSL);
//  k /= p.d[t] + 1.;
//  if ( order == 0 ) {
//    k *= -p.d[t];
//  } else {
//    k *= pow(Complex(p.x[t], p.y[t]), order);
//  }
//  p.xp[t] -= k.real;
//  p.yp[t] += k.imag;
//}
//
//__host__ __device__
//void HKicker(const int kick, Particles & p, const size_t t) {
//  p.xp[t] += kick/(p.d[t] + 1.);
//}
//
//__host__ __device__
//void VKicker(const int kick, Particles & p, const size_t t) {
//  p.yp[t] += kick/(p.d[t] + 1.);
//}

