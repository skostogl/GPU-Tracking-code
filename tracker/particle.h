#ifndef PARTICLE_H
#define PARTICLE_H

#include "float_type.h"

struct Particle {
  Tfloat x;
  Tfloat xp;
  Tfloat y;
  Tfloat yp;
  Tfloat z;
  Tfloat d;
};

/*
struct Particles {
  size_t size;
  Tfloat * x;
  Tfloat * xp;
  Tfloat * y;
  Tfloat * yp;
  Tfloat * z;
  Tfloat * d;

  __host__ __device__
  Particle operator[](size_t i) const {
    return {x[i], xp[i], y[i], yp[i], z[i], d[i]};
  }
};
*/

#endif
