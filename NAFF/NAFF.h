#include <fftw3.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <utility>

const double PI = boost::math::constants::pi<double>();
typedef double Tfloat;

namespace NAFF {

  double find_peak(const std::vector<Tfloat> &amps, const std::vector<Tfloat> &freqs)  {
    size_t maxIndex = 0;
    double maxVal = amps[0];
    const size_t N = amps.size();
    for (size_t j = 1; j < N; j++) {
      const double v = amps[j];
      if( v > maxVal ) {
        maxVal   = v;
        maxIndex = j;
      }
    }
    return freqs[maxIndex];
  }
 
  double stat_exp_value(const double f, const std::vector<std::complex<Tfloat>> &data) {
    std::vector<std::complex<double>> sum_vec;
    const std::complex<double> z (0,-1);
    for (size_t t=0; t<data.size(); t++) {
      sum_vec.push_back(std::conj(data[t])*std::exp(2.0*PI*f*t*z));
    }  
    const double sum_of_elems=-abs(std::accumulate(sum_vec.begin(),sum_vec.end(), std::complex<double>(0.0)));
    return sum_of_elems/data.size(); 
  }

  std::vector<std::complex<Tfloat>> apply_window(const std::vector<Tfloat> &init_data_x, const std::vector<Tfloat> &init_data_xp) {
    const size_t windowSize = init_data_x.size();
    std::vector<std::complex<Tfloat>> data;
    for (int i=0; i<windowSize; i++) {
      const double multiplier = 0.5*(1-cos(2*PI*i/(windowSize-1)));
      data.emplace_back(init_data_x[i]* multiplier, init_data_xp[i]*multiplier);
    }
    return data;
  }
}

std::pair<std::vector<double>,std::vector<double>> FFT(const std::vector<std::complex<Tfloat>> & data) {
  const size_t N = data.size();
  fftw_complex in[N], out[N];
  fftw_plan p;
  for (int i = 0; i < N; i++) {
    in[i][0] = data[i].real();
    in[i][1] = data[i].imag();
  }
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  std::vector<double> amps;
  std::vector<double> freqs;
  for (size_t i = 0; 2*i < N; i++) {	
    amps.push_back(sqrt(out[i][0]*out[i][0]+ out[i][1]*out[i][1]));
    freqs.push_back(i/(N-1.0));
  }
  fftw_destroy_plan(p);
  return std::make_pair(amps,freqs);
}

double NAFF_f1(const std::vector<Tfloat> & init_data_x,const std::vector<Tfloat> & init_data_xp)  {
  std::vector<std::complex<Tfloat>> data( NAFF::apply_window(init_data_x,init_data_xp) );
  auto fft = FFT(data);
  std::vector<double> & amps = fft.first;
  std::vector<double> & freqs = fft.second;
  double peak_frequency=NAFF::find_peak(amps, freqs);
  //std::cout<<"Tune from FFT: "<<peak_frequency<<std::endl;
  auto y = [&data](double f){return NAFF::stat_exp_value(f, data);};
  double small_step = 0.001; //This should be derived from init_data.size
  double a = peak_frequency-small_step;
  double b = peak_frequency+small_step;
  std::pair<double, double> r = boost::math::tools::brent_find_minima(y, a, b, 40); //one should foresee the possibility to adjust this 40, a global variable?
  //std::cout << "Tune from NAFF: "<<r.first <<std::endl;
  double scale = 1e-7;
  r.first = (int)(r.first / scale) * scale; // ?
  return r.first;
}


   
  
  
  
  


