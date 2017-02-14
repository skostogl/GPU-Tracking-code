#include "3Dvec.h"

/*int main() {
  Signal s;
  return 0;
}*/

int main() {
  std::vector<std::complex<double>> signal;
  double i=0.0;
  while (i<=0.1) {
    signal.push_back(std::complex<double>(100.0*sin(2.*3.141592*i*30.),0.0));
    i+=1e-4;
  }
  Signal s(signal);
  Component c(300.0,signal.size());
  Component proj = projection(s,c);
  //signal=signal-proj;
  return 0;
}

