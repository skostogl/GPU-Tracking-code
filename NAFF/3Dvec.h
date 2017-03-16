#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
#include <cmath>

#include "windows.h"
#include "consts.h"
#include "spline_interpolation.h"

#include <boost/numeric/odeint.hpp>

class Signal;
class Component;
class ComponentVector;

///////// Component: exp(i*2*pi*f*t) for one frequency
class Component {
  public:
    double freq;
    size_t signal_size; 
    std::complex<double> ampl;
    
    Component (): ampl(1.,0.) {}    
    Component(double freq, size_t signal_size): freq(freq), signal_size(signal_size), ampl(1.,0.) {}
    ~Component () {}
    
    std::complex<double> operator[] (size_t t) const {
      return (exp(2*pi*freq*t*imag));    
    }
    
    size_t size() const {
      return signal_size;
    }
};

///////// ComponentVector: Components with multiple frequencies
class ComponentVector { 
  public:
    typedef std::vector<std::complex<double>> data_t;
    std::complex<double> ampl;
    std::vector<Component> cs;

    ComponentVector():ampl(1,0) {}
    ComponentVector(const Component &other): ampl(1.,0.) { cs.emplace_back(other);}
    ~ComponentVector() {}

    size_t size() const{
      return cs[0].size();
    }

    std::complex<double> operator[] (size_t t) const {
      std::complex<double> res;
      for (const auto & c: cs) {
        res += ampl*c.ampl*c[t];
      }
      return res;
    }

    ComponentVector operator-=(const ComponentVector& other) {
      if (size() != other.size() ) throw std::runtime_error("Inner product sizes not equal!");
      for (size_t i = 0; i<size(); ++i) {
        (*this)[i] -= other[i];
      }
      return *this;
    }

};

///////// Signal with data in time domain
class Signal {
  public:
    ComponentVector::data_t data; 
    Signal() {}
    
    Signal(const std::vector<std::complex<double>> & v){
      for (const auto i:v) data.emplace_back(i);
    }
    Signal(const Signal &) = default;
    Signal& operator=(const Signal &) = default;

    ~Signal() {}
    
    std::complex<double> operator[] (size_t t) const {
      return data[t];
    }    

    std::complex<double> operator[] (double t) const {
      ///////////// Linear interpolation of data	    
      /*size_t a = (int)(t);
      size_t b = (int)(t)+1;
      return data[a]+(data[b]-data[a])*(t-a)/(1.0*(b-a));*/
      ///////////// Spline interpolation of data	    
      return spline(t, data);
     }
           
    size_t size() const{
      return data.size();
    }

    Signal operator-=(const ComponentVector& other) {
      if (size() != other.size() ) throw std::runtime_error("Inner product sizes not equal!");
      for (size_t i = 0; i<size(); ++i) {
        data[i] -= other[i];
        data[i] -= std::conj(other[i]);
      }
      return *this;
    }    
};

//////// Inner product
template <typename T1, typename T2>
std::complex<double> inner_product(const T1& t1, const T2& t2, const WindowFunc & window_func){
  if ( t1.size() != t2.size()) throw std::runtime_error("Inner product sizes not equal!");
  std::complex<double> result;
  for (size_t i =0; i<t1.size(); i++) {
  ////////// Double i for interpolation of signal	  
  //for (double i =0.0; i<t1.size(); i+=0.01) {
    result += t1[i]*std::conj(t2[i])*window_func(i,t1.size());
  }
  return result/(t1.size()*1.0);
}

//////// Projection of v on u
ComponentVector projection (const Component& v, const ComponentVector& u, const WindowFunc & window_func ){
  std::complex<double> num = inner_product(v, u, window_func);
  std::complex<double> den = inner_product(u, u ,window_func);
  ComponentVector s = v;
  s.ampl *= num/den;
  return s;
}

////////Projection of v on u 
void signal_projection (const Signal& v, ComponentVector& u, const WindowFunc & window_func){
  std::complex<double> num = inner_product(v, u, window_func);
  std::complex<double> den = inner_product(u, u, window_func);
  u.ampl = num/den;
}

//////// RMS for FFT
template <typename T>
double cmp_RMS (const T& data) {
  double res = abs(std::inner_product(data.begin(), data.end(), data.begin(),std::complex<double>(0,0))/(data.size()*1.0));
  return sqrt(res);  
}

////////Print component
std::ostream& operator <<(std::ostream& os, Component c) {
  os << "Frequency: " << c.freq << " Amplitude: " << c.ampl<< " ABS amplitude: "<< abs(c.ampl) << std::endl;
  return os;
}

////////Print signal
std::ostream& operator <<(std::ostream& os, ComponentVector s) {
  for (size_t i=0;i<s.size();i++) {
    os << s[i] << std::endl;
  }
  return os;
}
////////Print signal
std::ostream& operator <<(std::ostream& os, Signal s) {
  for (size_t i=0;i<s.size();i++) {
    os << s[i] << std::endl;
  }
  return os;
}
template <typename T1, typename T2>
void write_file (const std::string& file_name, const std::vector<std::complex<double>>& var1 ,const T1& min, const T2& max, const T1& step ) {
  std::ofstream myfile;
  myfile.open(file_name);
  for (T2 i=min; i<max; i+=step) {
    myfile<<std::setprecision(15)<<var1[i].real()<<" "<<std::setprecision(15)<<var1[i].imag()<<std::endl;
  }
}

template <typename T1, typename T2>
void write_file_merit (const std::string& file_name, const T1& y,const T2& min, const T2& max, const T2& step ) {
  std::ofstream myfile;
  myfile.open(file_name);
  for (T2 i=min; i<max; i+=step) {
    myfile<<std::setprecision(14)<<i<<" "<<-y(i)<<std::endl;
  }
}
