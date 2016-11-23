#ifndef BUNCH_H
#define BUNCH_H

#include <vector>
#include <random>

#include "nvrtc_wrap.h"
#include "float_type.h"

struct HostBunch;

std::default_random_engine bunch_rnd_gen(56789);

class DeviceBunch {
  friend class HostBunch;

  CUdeviceptr x,xp,y,yp,z,d;
  size_t bufferSize = 0;
  size_t np; //used only to pass out the void* in the args method

 public: 
  DeviceBunch() {}
  DeviceBunch(const size_t N) { alloc(N); }
  DeviceBunch(const HostBunch & hb);
  ~DeviceBunch() { dealloc(); }

  size_t n() const {
    return bufferSize/sizeof(Tfloat);
  }

  void alloc(const size_t N) {
    if ( bufferSize != 0 ) dealloc();
    bufferSize = N*sizeof(Tfloat);
    CUDA_SAFE_CALL(cuMemAlloc(&x,  bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&xp, bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&y,  bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&yp, bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&z,  bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&d,  bufferSize));
  }

  void dealloc() {
    CUDA_SAFE_CALL(cuMemFree(x));
    CUDA_SAFE_CALL(cuMemFree(xp));
    CUDA_SAFE_CALL(cuMemFree(y));
    CUDA_SAFE_CALL(cuMemFree(yp));
    CUDA_SAFE_CALL(cuMemFree(z));
    CUDA_SAFE_CALL(cuMemFree(d));
    bufferSize = 0;
  }

  std::vector<void*> get_args() {
    np = n();
    std::vector<void*> args = { &np, &x, &xp, &y, &yp, &z, &d };
    return args;
  }

  void operator=(const HostBunch & hb);
};

struct HostBunch {
  std::vector<Tfloat> x;
  std::vector<Tfloat> xp;
  std::vector<Tfloat> y;
  std::vector<Tfloat> yp;
  std::vector<Tfloat> z;
  std::vector<Tfloat> d;

  HostBunch() {}
  HostBunch(const size_t N):
    x (N,0.),
    xp(N,0.),
    y (N,0.),
    yp(N,0.),
    z (N,0.),
    d (N,0.)
  {std::cout << "HostBunch::HostBunch " << N << " "<< size()<< std::endl; }

  void set_z(const double sigma_z, const double mean_z = 0.) {
    std::normal_distribution<Tfloat> dist(mean_z,sigma_z);
    z.clear();
    for (size_t i = 0; i < x.size(); ++i) {
      z.push_back(dist(bunch_rnd_gen));
    }
  }

  void set_d(const double sigma_d, const double mean_d = 0.) {
    d.clear();
    std::normal_distribution<Tfloat> dist(mean_d,sigma_d);
    for (size_t i = 0; i < x.size(); ++i) {
      d.push_back(dist(bunch_rnd_gen));
    }
  }


  HostBunch(const HostBunch & o) = default;
  HostBunch(HostBunch && o) = default;
  HostBunch & operator=(const HostBunch & o) = default;
  HostBunch & operator=(HostBunch && o) = default;
 
  bool operator==(const HostBunch & o) {
    return ( x == o.x and xp == o.xp and y == o.y and yp == o.yp and z == o.z and d == o.d );
  }

  HostBunch(const DeviceBunch & db) { copyFromDeviceBunch(db); }

  void copyFromDeviceBunch(const DeviceBunch & db) {
    x .resize(db.n());
    xp.resize(db.n());
    y .resize(db.n());
    yp.resize(db.n());
    z .resize(db.n());
    d .resize(db.n());
    CUDA_SAFE_CALL(cuMemcpyDtoH(x .data(), db.x , db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyDtoH(xp.data(), db.xp, db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyDtoH(y .data(), db.y , db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyDtoH(yp.data(), db.yp, db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyDtoH(z .data(), db.z , db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyDtoH(d .data(), db.d , db.bufferSize));
  }

  void copyToDeviceBunch(DeviceBunch & db) const {
    db.alloc(size());
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.x ,  x.data(), db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.xp, xp.data(), db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.y ,  y.data(), db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.yp, yp.data(), db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.z ,  z.data(), db.bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(db.d ,  d.data(), db.bufferSize));
  }
  size_t size() const {
    const auto s = x.size();
    if ( s!=xp.size() or s!=y.size() or s!=yp.size() or s!=z.size() or s!=d.size() ) {
      throw std::runtime_error("HostBunch particle vector sizes do not match");
    }
    return s;
  }

  void operator=(const DeviceBunch & db) {
    copyFromDeviceBunch(db);
  }
};

void DeviceBunch::operator=(const HostBunch & hb) {
  hb.copyToDeviceBunch(*this);
}

DeviceBunch::DeviceBunch(const HostBunch & hb) {
  hb.copyToDeviceBunch(*this);
}

std::vector<Tfloat> get_x_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> xs;
  for (const auto & bunch:b) {
    xs.push_back(bunch.x[particle_id]);
  }
  return xs;
}

std::vector<Tfloat> get_xp_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> xp;
  for (const auto & bunch:b) {
    xp.push_back(bunch.xp[particle_id]);
  }
  return xp;
}

std::vector<Tfloat> get_y_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> ys;
  for (const auto & bunch:b) {
    ys.push_back(bunch.y[particle_id]);
  }
  return ys;
}

std::vector<Tfloat> get_yp_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> yp;
  for (const auto & bunch:b) {
    yp.push_back(bunch.yp[particle_id]);
  }
  return yp;
}

std::vector<Tfloat> get_d_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> d;
  for (const auto & bunch:b) {
    d.push_back(bunch.d[particle_id]);
  }
  return d;
}

std::vector<Tfloat> get_z_from_bunches (const std::vector<HostBunch> &b, const size_t particle_id) {
  std::vector <Tfloat> z;
  for (const auto & bunch:b) {
    z.push_back(bunch.z[particle_id]);
  }
  return z;
}
#endif
