#include <fftw3.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <utility>
#include <fstream>

const double PI = boost::math::constants::pi<double>();
typedef double Tfloat;
int counter=0;

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
    return sum_of_elems;
  }
///
//std::ofstream myfile;
//
//
//
  std::vector<std::complex<Tfloat>> apply_hann_window(const std::vector<Tfloat> &init_data_x, const std::vector<Tfloat> &init_data_xp) {
    const size_t windowSize = init_data_x.size();
    //if (counter==1)
      //myfile.open("/home/skostogl/cuTrack/dat_files/hann_window.dat");
    std::vector<std::complex<Tfloat>> data;
    for (int i=0; i<windowSize; i++) {

      const double multiplier = 0.5*(1-cos(2*PI*i/(windowSize-1)));
      //const double multiplier = 1.0;
      data.emplace_back(init_data_x[i]* multiplier, init_data_xp[i]*multiplier);
      //
      /*if (counter==1){
        myfile<<i<<" "<<multiplier<<std::endl;
      }*/
    
   }
    //
    /*if (counter==1){
      std::cout<<"hi"<<std::endl;
      std::ofstream myfile2;
      
      myfile2.open("/home/skostogl/cuTrack/dat_files/hann.dat");
      for (size_t i=0;i<windowSize;i++){  
        myfile2<<init_data_x[i]<<" "<<init_data_xp[i]<<std::endl;
      }
      myfile2.close();
    }*/
    //
    /*if (counter==1)
      myfile.close();
    // */
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
  for (size_t i=0; 2*i<N; i++) {
    amps.push_back(sqrt(out[i][0]*out[i][0]+ out[i][1]*out[i][1]));
    freqs.push_back(i/(N-1.0));
  } 
 /* if (counter==1){
	  std::cout<<"hi"<<std::endl;
    std::ofstream myfile;
    myfile.open("/home/skostogl/cuTrack/dat_files/fft_O3.dat");
    for (size_t i=0;2*i<N;i++){  
      myfile<<freqs[i]<<" "<<amps[i]<<std::endl;
    }
    myfile.close();
  }*/

  fftw_destroy_plan(p);
  return std::make_pair(amps,freqs);
}

double NAFF_f1( const std::vector<Tfloat> & init_data_x,const std::vector<Tfloat> & init_data_xp)  {
  counter=counter+1;
  std::vector<std::complex<Tfloat>> data; 
  for (int i=0;i<init_data_x.size();i++) 
    data.emplace_back(init_data_x[i], init_data_xp[i]);
  auto fft = FFT(data);
  std::vector<std::complex<Tfloat>> wdata( NAFF::apply_hann_window(init_data_x,init_data_xp) ); 
  std::vector<double> & amps = fft.first;
  std::vector<double> & freqs = fft.second;
  double peak_frequency=NAFF::find_peak(amps,freqs);
  auto y = [&wdata](double f){ return NAFF::stat_exp_value(f, wdata); };
  double small_step = 1./init_data_x.size();
  double a = peak_frequency-small_step;
  double b = peak_frequency+small_step;
  /*if (counter==1){
      std::cout<<"hi"<<std::endl;
      std::ofstream myfile2;    
      myfile2.open("/home/skostogl/cuTrack/dat_files/without_hann_o3.dat");
      double i=a-0.01;
      while (i<=b+0.01){  
        myfile2<<i<<" "<<-y(i)<<std::endl;
	i=i+0.0001;
      }
      myfile2.close();
    }*/ 
  int bits = std::numeric_limits<double>::digits;
  std::pair<double, double> r = boost::math::tools::brent_find_minima(y, a, b, bits);
  std::cout << std::setprecision(14) << " NAFF: " << r.first<<"\n";
  return r.first;
}


   
  
  
  
  


