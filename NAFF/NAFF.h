#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <utility>
#include <algorithm>
#include <tuple>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <boost/math/tools/minima.hpp>

#include "3Dvec.h"

typedef double Tfloat;
typedef std::vector<double> double_vec;
typedef std::vector<std::complex<double>> complex_vec;
const std::complex<double> z (0,1);

class NAFF {
  private:
    
    size_t fft_size, f_counter;
    double fft_frequency, fft_phase, RMS, max_index ; 
    fftw_plan fftw_plan_;
    Signal signal;
    std::vector<double> frequencies;
    std::vector<Signal> norm_vectors;
    std::pair<double_vec, double_vec> amps_freqs;
  //////// Initialize signal
  void input(double_vec &init_data_x, double_vec &init_data_xp) {
    for (size_t i=0; i<init_data_x.size(); i++) { 
      signal.data.emplace_back(std::complex<double>(init_data_x[i],init_data_xp[i])); 
     }
  } 
  //////// Fast Fourier Transform
  void FFTw () {
    std::vector<fftw_complex> fftw_(signal.size());
    fftw_plan_ = fftw_plan_dft_1d(signal.size(), reinterpret_cast<fftw_complex*>(&signal.data[0]), fftw_.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fftw_plan_);
    fft_size = fftw_.size();
    max_fft_frequency(fftw_);
    size_t N = signal.size();
    for (size_t i=0; 2*i<N; i++) {
      amps_freqs.first .push_back(sqrt(fftw_[i][0]*fftw_[i][0]+ fftw_[i][1]*fftw_[i][1]));
      amps_freqs.second.push_back(i/(N-1.0));
    }
  }
  //////// Hann filter
  void hann_window () {
    signal.hann_filter();
  }
  //////// First estimation of the peak frequency from FFT 
  void max_fft_frequency (std::vector<fftw_complex> &fftw_) {
    double max_amplitude = sqrt(fftw_[0][0]*fftw_[0][0]+fftw_[0][1]*fftw_[0][1]);
    max_index = 0;
    for (size_t i=1; 2*i<fft_size; i++ ) {
      const double current_amplitude = sqrt(fftw_[i][0]*fftw_[i][0]+fftw_[i][1]*fftw_[i][1]);
      if (current_amplitude > max_amplitude) {
        max_amplitude = current_amplitude;
	max_index = i;
      }      
    }  
    fft_frequency = (max_index / fft_size);
    fft_phase = atan2(fftw_[max_index][1], fftw_[max_index][0]);   
  }
  //////// Maximization of <f(t),exp(i*2*pi*f*t)> function for refined frequency f
  ///////////////////// First method: Golden Section Search
  template <typename FuncMin>
  double cmp_min (const FuncMin& f, double min_value, double med_value,double max_value, double precision) {
    double phi = (1+sqrt(5))/2;
    double resphi=2-phi;
    if (std::abs(min_value-max_value)<precision) {
      return (min_value+max_value)/2;
    }
    double d=med_value+resphi*(max_value-med_value);
    if (f(d) < f(med_value))
          return cmp_min(f,med_value,d,max_value,precision);
    else
          return cmp_min(f,d,med_value,min_value,precision);
   }
  ///////////////////// Second method: BOOST Brent
  double maximize_fourier_integral () {
    /*auto y = [this](double f) { 
        complex_vec sums; 
        sums.reserve(signal.size());
        for (size_t t=0; t<signal.size(); t++) {
          //sums.push_back(std::conj(signal.data[t])*std::exp(2.0*pi*f*t*z));
          sums.push_back(std::conj(std::exp(2.0*pi*f*t*z))*signal.data[t]);
        }  
        return (-abs(std::accumulate(sums.begin(), sums.end(), std::complex<double>(0.0)))/signal.size()); };*/
    
    auto y = [this](double f) { 
	Component c(f,signal.size());
	return (-abs(signal%c)); };

    double step = 1./signal.size();
    double min = fft_frequency - step;
    double max = fft_frequency + step;
    ///
      std::ofstream myfile2;    
      std::string name = "/home/skostogl/cuTrack/lam"+std::to_string(f_counter+1)+".dat";
      myfile2.open(name);
      double i=min-0.3;
      while (i<=max+0.3){  
        myfile2<<i<<" "<<y(i)<<std::endl;
	i=i+0.0001;
      }
      myfile2.close();
    ///
    double r_gold = cmp_min (y, min,fft_frequency, max, 1e-12);
    std::cout<<" frequency from golden section method "<<r_gold<<std::endl;
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = boost::math::tools::brent_find_minima(y, min, max, bits);
    std::cout<<" frequency from BOOST brent "<<r.first<<std::endl;
    return r.first;     
  }
  
  void subtract_contribution() {
    if (f_counter == 0) {
      Component v1(frequencies,signal.size());
      //Component proj = projection (signal,v1);
      Component proj = projection (signal,v1);
      /*std::ofstream myfile2;    
      std::string name = "/home/skostogl/cuTrack/signal.dat";
      myfile2.open(name);
      ///
      for (auto i:signal.data) {
        myfile2<<i.real()<<" "<<i.imag()<<std::endl;
      }
      myfile2.close();*/
      signal=signal-proj.ampl*proj.exp_();
      std::cout<<proj.ampl<<" Amplitude "<<std::endl;
      ///
      
      Signal u1(v1.exp_());
      norm_vectors.push_back(u1);

    }
    else {
      Component v_i(frequencies.back(),signal.size());
      Signal u_i(v_i.exp_());
      std::ofstream myfile2;    
      std::string name = "/home/skostogl/cuTrack/signal.dat";
      /*myfile2.open(name);
      ///
      for (auto i:signal.data) {
      //for (auto i:proj.exp_()) {
        //myfile2<<(proj.ampl*i).real()<<" "<<(proj.ampl*i).imag()<<std::endl;
        myfile2<<i.real()<<" "<<i.imag()<<std::endl;
      }
      myfile2.close();*/
      for (size_t i=1;i<f_counter;i++) {
	Signal proj = projection (v_i,norm_vectors[i]);
	u_i = u_i - proj.ampl*proj.data_();
      } 

      /*myfile2.open(name);
      ///
      //for (auto i:signal.data) {
      for (auto i:u_i.data) {
        //myfile2<<(proj.ampl*i).real()<<" "<<(proj.ampl*i).imag()<<std::endl;
        myfile2<<i.real()<<" "<<i.imag()<<std::endl;
      }
      myfile2.close();*/
      Signal proj = projection (signal, u_i);

      /*myfile2.open(name);
      ///
      //for (auto i:signal.data) {
      for (auto i:proj.data) {
        myfile2<<(proj.ampl*i).real()<<" "<<(proj.ampl*i).imag()<<std::endl;
        //myfile2<<(proj.ampl*i).real()<<" "<<i.imag()<<std::endl;
      }
      myfile2.close();*/


      signal=signal-proj.ampl*proj.data_();
      myfile2.open(name);
      ///
      for (auto i:signal.data) {
      //for (auto i:proj.data) {
        //myfile2<<(proj.ampl*i).real()<<" "<<(proj.ampl*i).imag()<<std::endl;
        myfile2<<i.real()<<" "<<i.imag()<<std::endl;
      }
      myfile2.close();
      norm_vectors.push_back(u_i);
    }
  }

  public:
 
  ~NAFF() {  
    fftw_destroy_plan(fftw_plan_); 
  } 

  double_vec get_f1 (double_vec &init_data_x,double_vec &init_data_xp) {
    if (frequencies.size() == 0) {
      input(init_data_x, init_data_xp);
    }
    FFTw();
    //hann_window();
    frequencies.push_back(maximize_fourier_integral());
    //frequencies.push_back(fft_frequency);
    return frequencies;
  }
  
  double_vec get_f(double_vec &init_data_x, double_vec &init_data_xp) {
    f_counter=0;
    while (f_counter<3) {
      get_f1(init_data_x, init_data_xp);
      subtract_contribution();  
      f_counter++;   
    }
  return frequencies;
  }

};

