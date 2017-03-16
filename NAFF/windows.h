#pragma once

#include <vector>
#include "consts.h"

std::vector<std::complex<double>> cheb_window(const size_t N, const double a = 5.0) {
  std::vector<double> out;
  std::vector<std::complex<double>> complex_out;
  out.resize(N);
  std::cout << "Chebyshev window calculated with a="<<a<<std::endl;
  auto cheb_poly = [](int n, double x) {
    double res;
    if (fabs(x) <= 1) res = cos(n*acos(x));
    else  res = cosh(n*acosh(x));
    return res;
  };
  double max=0;
  double tg = pow(10, a);
  double x0 = cosh((1.0/(N-1))*acosh(tg));
  double M = (N-1)/2;
  if(N%2==0) M = M + 0.5; 
  for(size_t nn=0; nn<(N/2+1); nn++){
      double n = nn-M;
      double sum = 0;
      for(size_t i=1; i<=M; i++){
          sum += cheb_poly(N-1,x0*cos(pi*i/N))*cos(2.0*n*pi*i/N);
      }
      out[nn] = tg + 2*sum;
      out[N-nn-1] = out[nn];
      if (out[nn]>max)
        max=out[nn];
  }
  for(size_t nn=0; nn<N; nn++) 
    out[nn] /= max; 
  for (const auto i:out) {
    complex_out.emplace_back(i,0.);
  } 
   /*std::ofstream myfile;
   myfile.open("window_cheb.dat");
   for (size_t i= 0; i< complex_out.size(); i++) {
     myfile<<complex_out[i].real()<<" "<<complex_out[i].imag()<<std::endl;
   }*/
  return complex_out;

}

std::vector<std::complex<double>> hann_harm_window(const size_t N, const double n = 3.0) {
  std::vector<std::complex<double>> out;
  double T1 = 0;
  double T2 = N;
  double TM = (T2-T1)/2.0;
  double PIST = pi/TM;
  int factorial_1 = 1, factorial_2 = 1;
  std::cout << "Hann window used with h="<<n<<std::endl;
  for (size_t j=1;j<n+1;j++) {
    factorial_1 *= j;
  }
  for (size_t j=1;j<2*n+1;j++) {
    factorial_2 *= j;
  }
  double cn = pow(2, n) * pow(factorial_1, 2)/(factorial_2*1.0);
  for (size_t i=1; i<N+1; i++ ) {
    double T = (i-1) - TM;
    out.emplace_back(cn*pow((1.0 + cos(T*PIST)), n));
  }
   /*std::ofstream myfile;
   myfile.open("window_hann.dat");
   for (size_t i= 0; i< out.size(); i++) {
     myfile<<abs(out[i])<<std::endl;
   }*/
  return out;
}

class Signal_window {
  public: 
    std::vector<std::complex<double>> data;
    Signal_window() {}
    Signal_window(const std::vector<std::complex<double>> & v){
      for (const auto i:v) data.emplace_back(i);
    }
    ~Signal_window() {}
    
    std::complex<double> operator[] (size_t t) const {
      return data[t];
    }    
    size_t size() const{
      return data.size();
    }
};

class WindowFunc {
  private:
    Signal_window window;
  public:
    double parameter = 3.0;
    char type = 'h';

    WindowFunc(): window() {}
    ~WindowFunc() {}
    
   void compute(const size_t N) {
     if (type == 'c') { 
       window = cheb_window(N, parameter);
     }
     else if (type == 'h') {
       window = hann_harm_window(N, parameter);
     }
     else {
       throw std::runtime_error("window type is not defined");
     }
    }

   void compute(const size_t N, const double param, const char tp) {
     if (type == 'c') { 
       window = cheb_window(N, param);
     }
     else if (type == 'h') {
       window = hann_harm_window(N, param);
     }
     else {
       throw std::runtime_error("window type is not defined");
     }
      parameter = param;
      type = tp;
    }

   std::complex<double> operator()(size_t t, size_t N) {
     if (window.size()!=N) {
       compute(N);
     }
     return window[t];
   }
    
   std::complex<double> operator()(size_t t, size_t N) const {
     if (window.size()!=N) {
       throw std::runtime_error("window and signal sizes do not match"); 
     }
     return window[t];
   } 
};
