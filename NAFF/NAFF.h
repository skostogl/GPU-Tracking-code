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
    double maxVal=-1, maxIndex=-1;
    int N=amps.size();
    for (int j=0; j<N; j++) {
      double v=amps[j];
      if( v > maxVal )  {
        maxVal   = v;
        maxIndex = j;
      }
    }
    return freqs[maxIndex];
  }
 
  double stat_exp_value(double f, const std::vector<std::complex<Tfloat>> &data)  {
    std::vector<std::complex<double>> sum_vec;
    const std::complex<double> z (0,-1);
    for (int t=0; t<data.size(); t++)  {
      sum_vec.push_back(std::conj(data[t])*std::exp(z*2.0*PI*f*double(t)));  
    }  
    double sum_of_elems=-abs(std::accumulate(sum_vec.begin(),sum_vec.end(), std::complex<double>(0.0)));
    return sum_of_elems/data.size(); 
  }

  std::vector<std::complex<Tfloat>>  apply_window(const std::vector<Tfloat> &init_data_x, const std::vector<Tfloat> &init_data_xp)  {
    int windowSize =init_data_x.size();
    std::vector<std::complex<Tfloat>> data;
    for (int i=0; i<windowSize; i++)  {
      double multiplier = 0.5*(1-cos(2*PI*i/(windowSize-1)));
      data.push_back(std::complex<Tfloat>(init_data_x[i]* multiplier, init_data_xp[i]*multiplier));
    } 
    return data;
  }   
}

std::pair<const std::vector<double>,const std::vector<double>> FFT(const std::vector<std::complex<Tfloat>> & data) {
  int N=data.size();
  std::pair<std::vector<double>, std::vector<double> > results;
  fftw_complex in[N], out[N];
  fftw_plan p;
  for (int i = 0; i < N; i++)  {
    in[i][0] = data[i].real();
    in[i][1] = data[i].imag();
  }
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  std::vector<double> amps;
  std::vector<double> freqs;
  int i=0;
  while (i/double(N)<0.5){	
    amps.push_back(sqrt(out[i][0]*out[i][0]+ out[i][1]*out[i][1]));
    freqs.push_back(i/double(N-1));
    i++;
  }
  fftw_destroy_plan(p);
  results=std::make_pair(amps,freqs);
  return results;
}

double NAFF_f1(const std::vector<Tfloat> & init_data_x,const std::vector<Tfloat> & init_data_xp)  {
  std::vector<std::complex<Tfloat>> data;
  data=NAFF::apply_window(init_data_x,init_data_xp);
  std::vector<double> amps;
  std::vector<double> freqs;
  amps=FFT(data).first;
  freqs=FFT(data).second;
  double peak_frequency=NAFF::find_peak(amps, freqs);
  //std::cout<<"Tune from FFT: "<<peak_frequency<<std::endl;
  auto y=[data](double f) {return ((NAFF::stat_exp_value(f, data)));};
  double small_step = 0.001;
  double a = peak_frequency-small_step;
  double b = peak_frequency+small_step;
  std::pair<double, double> r = boost::math::tools::brent_find_minima(y, a, b, 40);
  //std::cout << "Tune from NAFF: "<<r.first <<std::endl;
  double scale = 1e-7;
  r.first = (int)(r.first / scale) * scale;
  return r.first;
}


   
  
  
  
  


