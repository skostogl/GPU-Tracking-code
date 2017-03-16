#pragma once

#include <fftw3.h>
#include <boost/math/tools/minima.hpp>

#include "3Dvec.h"

typedef double Tfloat;
typedef std::vector<double> double_vec;
typedef std::vector<std::complex<double>> complex_vec;
const std::complex<double> z (0,1);

class NAFF {
  private:
    
  size_t fft_size, f_counter;
  double fft_frequency, fft_phase, max_index, RMS; 
  double_vec frequencies;
  fftw_plan fftw_plan_;
  Signal signal;
  std::vector<ComponentVector> norm_vectors;
  bool f_found = true; 
  WindowFunc window;
  std::string merit_func;

  //////// Initialize signal and window
  void input(double_vec &init_data_x, double_vec &init_data_xp) {
    for (size_t i=0; i<init_data_x.size(); i++) { 
      signal.data.emplace_back(std::complex<double>(init_data_x[i], init_data_xp[i]));
    }
    window.compute(signal.size());
  }

  //////// Fast Fourier Transform
  void FFTw () {
    std::vector<fftw_complex> fftw_(signal.size());
    fftw_plan_ = fftw_plan_dft_1d(signal.size(), reinterpret_cast<fftw_complex*>(&signal.data[0]), fftw_.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fftw_plan_);
    fft_size = fftw_.size();
    max_fft_frequency(fftw_);
  }
  
  //////// First estimation of the peak frequency from FFT 
  void max_fft_frequency (std::vector<fftw_complex> &fftw_) {
    double_vec amps;
    double max_amplitude = sqrt(fftw_[0][0]*fftw_[0][0]+fftw_[0][1]*fftw_[0][1]);
    amps.push_back(max_amplitude);
    max_index = 0.0;
    for (size_t i=1; 2*i<fft_size; i++ ) {
      const double current_amplitude = sqrt(fftw_[i][0]*fftw_[i][0]+fftw_[i][1]*fftw_[i][1]);
      amps.push_back(current_amplitude);
      if (current_amplitude > max_amplitude) {
        max_amplitude = current_amplitude;
	max_index = i;
      }      
    }  
    fft_frequency = ((max_index*1.0)/ (fft_size*1.0));
    fft_phase = atan2(fftw_[max_index][1], fftw_[max_index][0]);
    double threshold = 10.0;
    RMS = cmp_RMS(amps);
    if (max_amplitude < threshold * RMS) 
      f_found = false;     
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
    double d = med_value+resphi*(max_value-med_value);
    if (f(d) < f(med_value))
          return cmp_min(f, med_value, d, max_value,precision);
    else
          return cmp_min(f, d, med_value, min_value, precision);
   }

  ///////////////////// Second method: BOOST Brent
  double maximize_fourier_integral () {
    auto y = [this](double f) { 
      Component c(f,signal.size());
      return (-abs(inner_product(signal, c, window))) ;};

    double step = 1.0/signal.size();
    double min = fft_frequency - 1.0*step;
    double max = fft_frequency + 1.0*step;
 
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = boost::math::tools::brent_find_minima(y, min, max, bits);
    /// Save merit function in file
    /*if (f_counter == 0) {
      write_file_merit("hann.dat", y,0.001 ,0.016 , 1e-6);
    }*/
    ///
    std::cout<<"Merit function: Maximize fourier integral"<<std::endl;
    return r.first;   
  }
  
  ///////////////////// Subtraction of frequency components from signal 
  void subtract_frequency(Signal& signal, double& frequency) {
    Component v_i(frequency, signal.size());
    ComponentVector u_i(v_i);
    if (f_counter!= 0) {
      for (size_t i=0; i<f_counter; i++) {
	u_i -= projection (v_i, norm_vectors[i], window);
      } 
    }
    signal_projection (signal, u_i, window);
    signal -= u_i;
    norm_vectors.push_back(u_i);
  }

  ////////////////////// Keep frequency which results in the minimum RMS in time domain
  double minimize_RMS_time () {
    auto y = [this](double current_frequency) {
      Signal signal_copy = signal;    
      Component v_i(current_frequency, signal_copy.size());
      subtract_frequency(signal_copy, current_frequency);
      double step = 1.0/signal.size();
      double min = fft_frequency - 1.0*step;
      double max = fft_frequency + 1.0*step;
      double stepp = (max-min)/100.0;
      double sum = 0;
      for (double i = min; i <= max; i+=stepp) {
        auto curr = abs(signal_copy[i]);
        sum += pow((curr),2);
       }
      return sum;
    };
    
    double step = 1.0/signal.size();
    double min = fft_frequency - 1.0*step;
    double max = fft_frequency + 1.0*step;

    /*if (f_counter == 0) {
      //write_file_merit("time_zoom3.dat", y, 0.001, 0.005, 1e-5);
      write_file_merit("time.dat", y,0.002999992 ,0.0030000095 , 1e-11);
    } */
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = boost::math::tools::brent_find_minima(y, min, max, bits);
    std::cout<<"Merit function: Time domain power"<<std::endl;
    return r.first;
  }



  ////////////////////// Keep frequency which results in the minimum RMS in frequency domain
 double minimize_RMS_frequency () {
    auto y = [this](double current_frequency) {
      Signal signal_copy = signal;    
      Component v_i(current_frequency, signal_copy.size());
      subtract_frequency(signal_copy, current_frequency);
      double step = 1.0/signal.size();
      double min = fft_frequency - 1.0*step;
      double max = fft_frequency + 1.0*step;
      double stepp = (max-min)/100.0;

      auto fourier_integral = [&signal_copy, this](double f) { 
        Component c(f,signal_copy.size());
	return (-abs(inner_product(signal_copy, c,window))); };
      
      auto area= [this, &stepp, &max, &min, &fourier_integral, &current_frequency] () {
        double sum = 0;
        for (double i = min+stepp; i <= max; i+=stepp) {
	  auto curr = fourier_integral(i);
	  sum += pow((curr),2);
        }
        return sum;
      };
      return area();
    };
      
    /*if (f_counter ==1) {
      //write_file_merit("area_zoom_hann3.dat", y, 0.001, 0.004, 1e-5);
      write_file_merit("area.dat", y,0.002999992 ,0.0030000095 , 1e-11);
      }*/
     
    double step = 1.0/signal.size();
    double min = fft_frequency - 1.*step;
    double max = fft_frequency + 1.*step;
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = boost::math::tools::brent_find_minima(y, min, max, bits);
    std::cout<<"Merit function: Fourier integral after subtraction"<<std::endl;
    return r.first;
  }

  public:
  size_t fmax = 4;

  ~NAFF() { 
    fftw_destroy_plan(fftw_plan_);
  } 
  
  void set_window_parameter(const double p, const char tp) {
    window.parameter = p;
    window.type = tp;
  }

  double get_window_parameter() const {
    return window.parameter;
  }

  void set_merit_function(const std::string m) {
    merit_func = m;    
  }
  
  double_vec get_f1 (double_vec &init_data_x,double_vec &init_data_xp) {
    if (frequencies.size() == 0) {
      input(init_data_x, init_data_xp);
    }
    FFTw();
    if (f_found == true) { 
      if (merit_func == "minimize_RMS_frequency") {
        frequencies.push_back(minimize_RMS_frequency());
      }
      else if (merit_func == "minimize_RMS_time") {
        frequencies.push_back(minimize_RMS_time());
      }
      ////////Default: maximize fourier integral
      else 
        frequencies.push_back(maximize_fourier_integral());
    }
    return frequencies;
  }
  
  double_vec get_f(double_vec &init_data_x, double_vec &init_data_xp) {
    f_counter = 0;
    while ((f_counter<fmax) && f_found == true) {
      std::cout<<"Frequency: "<<f_counter+1<<std::endl;
      get_f1(init_data_x, init_data_xp);
      if (f_found == true) {
        subtract_frequency(signal, frequencies.back());
      }
    f_counter++;    
    }    
    std::cout<< "Total number of frequencies found: " <<f_counter<<std::endl;
    return frequencies;
  }
};

