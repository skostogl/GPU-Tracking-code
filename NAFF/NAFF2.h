#include <fftw3.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <utility>
#include <fstream>
#include <memory>
#include <algorithm>
#include <tuple>
#include <math.h>

const double PI = boost::math::constants::pi<double>();
typedef double Tfloat;

namespace NAFF {

  // RMS computation of FFT
  int counter=0;
  double RMS_init = 0;
  bool cmp_RMS (std::vector<double> amps, double threshold) {
    bool f_found = false; 
    counter ++;
    double mean = accumulate (amps.begin(), amps.end(),0.0)/ amps.size();
    std::vector<double> diff (amps.size());
    // sqrt(Sigma (x-mean)^2 /N)
    std::transform(amps.begin(), amps.end(), diff.begin(), [mean](double x) {return x-mean;});
    double sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double RMS= sqrt(sum/amps.size());
    if (counter == 1) {
      RMS_init = RMS;
    }
    double max_amp = *std::max_element(amps.begin(), amps.end());
//    if ((max_amp >= threshold *RMS) && (RMS>0.05*RMS_init) ) ////To be optimized
      f_found=true;
    return f_found;
  }
  
  // First estimate of frequency and phase from FFT
  std::tuple<double, double, bool> find_max_freq(const std::vector<fftw_complex> &fft, double threshold)  {
    size_t maxIndex = 0;
    const size_t N = fft.size();
    std::tuple< std::vector<double>, std::vector<double>, std::vector<double> > amp_freq_phase;
    for (size_t i=0; 2*i<N; i++) {
      // abs (amplitude)
      std::get<0>(amp_freq_phase).push_back(sqrt(fft[i][0]*fft[i][0]+ fft[i][1]*fft[i][1]));
      // frequencies
      std::get<1>(amp_freq_phase).push_back(i/double(N));//-1.0));
      // phase
     //std::get<2>(amp_freq_phase).push_back(atan2(fft[i][1],-fft[i][0])+PI);
      //std::get<2>(amp_freq_phase).push_back(atan2(fft[i][1],fft[i][0])+PI*0.5);
      //std::get<2>(amp_freq_phase).push_back(atan2(fft[i][1],fft[i][0]));
      std::get<2>(amp_freq_phase).push_back(6.28*0.1);
    }
    bool f_found = cmp_RMS (std::get<0>(amp_freq_phase),threshold);
    if (f_found == true) {
      double maxVal  = std::get<0>(amp_freq_phase)[0];
      const size_t M = std::get<0>(amp_freq_phase).size();
      for (size_t j = 1; j < M; j++) {
        const double v = std::get<0>(amp_freq_phase)[j];
        if( v > maxVal ) {
          maxVal   = v;
          maxIndex = j;
        }
      }
      std::tuple<double, double, bool> freq_phase (std::get<1>(amp_freq_phase)[maxIndex], std::get<2>(amp_freq_phase)[maxIndex], f_found);
      std::cout<<fft[maxIndex][0] <<" "<<fft[maxIndex][1]<<" frequency "<<std::get<1>(amp_freq_phase)[maxIndex]<<" phase "<<std::get<2>(amp_freq_phase)[maxIndex]<<" freq,phase "<<std::endl;
      return freq_phase ;
    }
    else { 
      std::tuple<double, double, bool> no_f_found (0.,0., false);
      return no_f_found ;
    }
  }
 
  // Absolute value of fourier integral
  double abs_fourier_integral(const double f, const std::vector<std::complex<Tfloat>> &data, double phase) {
    std::vector<std::complex<double>> sum_vec;
    sum_vec.reserve(data.size());
    const std::complex<double> z (0,-1);
    for (size_t t=0; t<data.size(); t++) {
      sum_vec.push_back((std::conj(data[t])*std::exp((2.0*PI*f*t+phase)*z)));//phase//////////////////
    }  
    return -abs(std::accumulate(sum_vec.begin(),sum_vec.end(), std::complex<double>(0.0)))/data.size();
  }
  // Exponential factor exp(2*pi*f*t), equal to delta(f-f0) in frequency domain
  std::vector<std::complex<double>> positive_exp (double f, int N,double phase) {
    std::vector<std::complex<double>> exp_vec;
    exp_vec.reserve(N);
    const std::complex<double> z (0,1);
    for (int i=0;i<N;i++){
      exp_vec.push_back(std::exp((2.0*PI*f*i+phase)*z));//phase//////////////////
    }
    return exp_vec;
  } 
 
  // Exponential factor exp(-2*pi*f*t), equal to delta(f+f0) in frequency domain
  std::vector<std::complex<double>> negative_exp (double f, int N, double phase){
    std::vector<std::complex<double>> exp_vec;
    exp_vec.reserve(N);
    const std::complex<double> z (0,-1);
    for (int i=0;i<N;i++){
      exp_vec.push_back(std::exp((2.0*PI*f*i+phase)*z));//phase/////////////
    }
    return exp_vec;
  } 

  // Absolute value of fourier integral as a function of f
  auto fourier_integral_function(const std::vector<std::complex<Tfloat> > &data, double phase) {
    auto y= [&data,phase](double f){ return NAFF:: abs_fourier_integral(f, data,phase); };
    double a=0.03+1e-2;
    //double a=0.03+1e-7;
    double b=0.03-1e-2;
    //double b=0.03-1e-7;
    std::ofstream myfile;
    myfile.open ("/home/skostogl/cuTrack/fft.dat");
    int i=0;
    while (b<=a) { 
      myfile<<b<<" "<<y(b)<<std::endl;
      //b=b+1e-7*i*1e-7;
      b+=1e-6;
      i++;
    }
    return y;
  }

  // Fourier integral
  std::complex<Tfloat> fourier_integral(const std::vector<std::complex<Tfloat>> &a,const std::vector<std::complex<Tfloat>> &b){
    std::vector<std::complex<double>> integral (a.size());
    for (size_t i=0;i<a.size();i++){
      integral.push_back(std::conj(a[i])*b[i]);
    }
    return std::accumulate(integral.begin(), integral.end(), std::complex<double>(0.0))/(1.0*a.size());
  }



  // Hann window
  std::vector<std::complex<Tfloat>> apply_hann_window(const std::vector<Tfloat> &init_data_x, const std::vector<Tfloat> &init_data_xp) {
    const size_t windowSize = init_data_x.size();
    std::vector<std::complex<Tfloat>> data;
    data.reserve(init_data_x.size());
    for (size_t i=0; i<windowSize; i++) {
      //const double multiplier = 0.5*(1-cos(2*PI*i/(windowSize-1)));
      const double multiplier = 1;
      //data.emplace_back(init_data_x[i]* multiplier*2.0, init_data_xp[i]*multiplier*2.0);
      data.emplace_back(init_data_x[i]* multiplier, init_data_xp[i]*multiplier);
    }
    return data;
  }

  // Maximize Abs fourier integral function
  double maximize_fourier_integral(auto &y, size_t N, const std::tuple<double,double,bool> &max_freq) {
    double small_step = 1./N;
    double a = std::get<0>(max_freq) - small_step;
    double b = std::get<0>(max_freq) + small_step;
    int bits = std::numeric_limits<double>::digits;
    std::cout<< a<<" "<<b<<std::endl;
    std::pair<double, double> r = boost::math::tools::brent_find_minima(y, a, b, bits);
    return r.first;
  }
  
  // Bundle two real vectors into one compelex vector
  std::vector<std::complex<Tfloat>> convert_vectors_to_complex(const std::vector<Tfloat> &init_data_x, const std::vector<Tfloat> &init_data_xp) {
    std::vector<std::complex<Tfloat>> data;
    data.reserve(init_data_x.size());
    for (size_t i=0; i<init_data_x.size(); i++) {
      data.emplace_back(init_data_x[i], init_data_xp[i]);
    }
    return data;
  }
  
  // Split complex vector into a real and imag pair of vectors
  std::pair<std::vector<Tfloat>,std::vector<Tfloat>> split_vector(const std::vector<std::complex<Tfloat>> &data) {
    std::pair<std::vector<Tfloat>,std::vector<Tfloat>> split_data ;
    for (size_t t=0;t<data.size();t++){
      split_data.first.push_back(data[t].real());
      split_data.second.push_back(data[t].imag());
    } 
    return split_data;
  }

  // Subtract signal with f component and move on to the next frequency
/*  void subtract_amp (std::vector<Tfloat> & data_x,std::vector<Tfloat> &data_xp, double frequency, double phase,std::vector<std::complex<double>> &subtr) {
    //std::vector<std::complex<double>> data = convert_vectors_to_complex(data_x,data_xp);
    //std::vector<std::complex<double>> data = apply_hann_window(data_x,data_xp);
    std::vector<std::complex<double>> data = convert_vectors_to_complex(data_x,data_xp);
    std::pair<std::vector<std::complex<double>>,std::vector<std::complex<double>>> exp_factor;
    std::complex<double> amplitude, amplitude_neg;
    exp_factor.first  = ( positive_exp(frequency, data_x.size(),phase) );
    exp_factor.second = ( negative_exp(frequency, data_x.size(),phase) );
    // amplitude = Sigma (data* exp(i*2*pi*f))
    amplitude  = fourier_integral(exp_factor.first,data);
    amplitude_neg.real(amplitude.real());
    amplitude_neg.imag(-amplitude.imag());
    for (size_t i=0; i<data_x.size(); i++) {
      // data= data - amplitude * Sigma (exp(i*2*pi*f))
      std::complex<double> a = (amplitude*exp_factor.first[i]) + (amplitude_neg*exp_factor.second[i]) ;
      const std::complex<double> z (0,1);
      subtr.push_back(a);:x

      data[i] = data[i]-a;
    }
    // return new data
    std::pair<std::vector<Tfloat>,std::vector<Tfloat>> new_data = split_vector(data);
    data_x  = new_data.first;
    data_xp = new_data.second;
  }*/

  
} //namespace NAFF

std::vector<fftw_complex> FFT(const std::vector<Tfloat> & re, const std::vector<Tfloat> & im) {
  const size_t N = re.size();
  //alloc memory and prepare data
  std::vector<fftw_complex> out(N);
  std::vector<fftw_complex> in(N); 
  for (size_t i = 0; i < N; i++) {
    in[i][0] = re[i];
    in[i][1] = im[i];
  }
  //make plan and call fftw
  fftw_plan p;
  p = fftw_plan_dft_1d(N, in.data(), out.data(), FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  return out;
}


double NAFF_f1( const std::vector<Tfloat> & init_data_x,const std::vector<Tfloat> & init_data_xp)  {
  // FFT
  auto fft = FFT(init_data_x, init_data_xp);
  // Hann window
  std::vector<std::complex<Tfloat>> wdata = NAFF::apply_hann_window(init_data_x,init_data_xp); 
  // first estimate from FFT
  double threshold=3.0;
  std::tuple<double,double,bool> max_freq_phase=NAFF::find_max_freq(fft,threshold);
  if (std::get<2>(max_freq_phase)==true) {
    // optimize frequency determination
    auto y= NAFF::fourier_integral_function (wdata, std::get<1>(max_freq_phase));
    return NAFF::maximize_fourier_integral(y,init_data_x.size(), max_freq_phase);
  }
  else {
    std::cout<<" No frequency component found " <<std::endl;
    return 0;
  }
}


/*std::vector<double> NAFF_f( std::vector<Tfloat> init_data_x,const std::vector<Tfloat> init_data_xp)  {
  int components = 0;
  bool f_found   = true;
  ///////////////// To be removed 
  std::vector<Tfloat> data_x  = init_data_x;
  std::vector<Tfloat> data_xp = init_data_xp;
  /////////////////////////
  std::vector<double> frequencies;
  while (f_found == true) {
    // FFT
    auto fft = FFT(data_x, data_xp);
    // Hann window
    std::vector<std::complex<Tfloat>> wdata = NAFF::apply_hann_window(data_x, data_xp); 
    // first estimate from FFT
    double threshold = 3.0;
    std::tuple<double,double,bool> max_freq_phase = NAFF::find_max_freq(fft, threshold);
    if (std::get<2>(max_freq_phase) == true) {
      // optimize frequency determination
      auto y= NAFF::fourier_integral_function (wdata,std::get<1>(max_freq_phase));
      ////////////// To be removed
      std::vector<Tfloat> initial_data_x; 
      std::vector<Tfloat> initial_data_xp; 
      if (components==0) {
        for (auto i: wdata) {
          initial_data_x.push_back(i.real());
          initial_data_xp.push_back(i.imag());
        }
      }
      ///////////////////////
      frequencies.push_back(maximize_fourier_integral(y, data_x.size(), max_freq_phase));
      // remove frequency component from signal and repeat the process(subtr parameter to be removed)
      std::vector<std::complex<double>> subtr;
      NAFF::subtract_amp(data_x, data_xp, frequencies.back(),std::get<1>(max_freq_phase),subtr);
      /////////////// To be removed
      if (components == 0) {
        std::ofstream myfile;
        myfile.open ("/home/skostogl/cuTrack/text.dat");
        for (int i=0;i<init_data_x.size();i++)
          myfile << initial_data_x[i]<<" "<<initial_data_xp[i]<<" "<<subtr[i].real()<<" "<<subtr[i].imag()<<" "<<data_x[i]<<" "<<data_xp[i]<<std::endl;
        myfile.close();
      }
      ////////////////////////////
      components ++;
      //////////////// To be removed
      if (components == 5) 
        f_found = false;
    }
     //////////////////////
    else {
       if (components == 0)    
         std::cout<<" No frequency component found " <<std::endl;
       else 
	 std::cout<<components<<" frequencies found " <<std::endl;
       f_found=false;
    }  
  }
 return frequencies;
}*/



