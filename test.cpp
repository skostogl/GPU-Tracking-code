// CUDA_PATH=/opt/cuda/ && g++ -I $CUDA_PATH/include -L $CUDA_PATH/lib64 -lnvrtc -lcuda -lcudart -Wl,-rpath,$CUDA_PATH/lib64 nvrtc.cpp -o nvrtc

#include <iostream>
#include <sstream>
#include <iomanip>
#include <future>

#include "include/bunch.h"
#include "include/lattice.h"

int main() {
  CUDA_info();

  Lattice lattice("LHC/lhc.twi");
//  Lattice lattice("small_ring/lattice.twi");
//  Lattice lattice;

//  lattice.add("Drift(1.)");
//  lattice.add("Multipole(1, 0.1, 0.)");
//  lattice.add("Drift(1.)");
//  lattice.add("Multipole(1, -0.1, 0.)");

  lattice.n_turns = 1000;
  std::cout << "N_ELE: " << lattice.get_n_elements() << std::endl;
  lattice.optimize();
  std::cout << "N_ELE: " << lattice.get_n_elements() << std::endl;

  size_t n = 100000;
  //size_t n = 512*512;
  HostBunch b(n);
  DeviceBunch db(b);

  // Execute.
  for (int k = 0; k < 1 ; k++) {
//    auto future = std::async( std::launch::async, [&](){
//      //offline analysis on a separated thread
//      std::ofstream f0( std::to_string(k)+".dat" );
//      for (size_t i = 0; i < 1000; ++i) {
//        f0 << std::setprecision(12) << b.x[i] << " " << b.xp[i] << " " << b.d[i] << " " << b.z[i] << "\n";
//      }
//    });

    lattice.track(db);
    b = db;
    std::ofstream f0( "1000.dat" );
    for (size_t i = 0; i < 1000; ++i) {
      f0 << std::setprecision(12) << b.x[i] << " " << b.xp[i] << " " << b.d[i] << " " << b.z[i] << "\n";
    }
  }

  return 0;
};

