#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
#include <cmath>

const double pi = 3.141592653589793238462643383279;
const std::complex<double> imag(0.,1.);
class Signal;
class Component;

class Component {
  public:
    double freq;
    std::vector<double> frequencies;
    size_t signal_size; 
    std::complex<double> ampl;
    Component () {}    
    Component(std::vector<double> frequencies,size_t signal_size): freq(frequencies.back()), frequencies(frequencies), signal_size(signal_size){}
    Component(double freq,size_t signal_size): freq(freq), signal_size(signal_size){}
    ~Component() {}

    std::complex<double> operator[](size_t t) const{
      return (exp(2*pi*freq*t*imag));    
    }
    
    size_t size() const{
      return signal_size;
    }

    std::vector<std::complex<double>> exp_ () {
      std::vector<std::complex<double>> res;
      for (size_t i=0;i<size();i++) {
        res.push_back((*this)[i]);
      }
      return res;
    }

};

class Signal { 
  public:
    typedef std::vector<std::complex<double>> data_t;
    data_t data;
    Signal() {}
    Signal(const data_t & data):data(data){}
    ~Signal() {}
    std::complex<double> ampl;

    size_t size() const{
      return data.size();
    }

    std::complex<double> operator[] (size_t t) const {
      return data[t];
    }

    void hann_filter () {
      for (size_t i=0; i<size(); i++) {
        //const double multiplier = 0.5 * (1-cos(2*pi*i/(size()-1.0)));
        const double multiplier =1 + cos(pi*i/(size()*1.0));
        data[i]*=multiplier;
      }
    }
    std::vector<std::complex<double>> data_ () {
      std::vector<std::complex<double>> res;
      for (size_t i=0;i<size();i++) {
        res.push_back((*this)[i]);
      }
      return res;
    }

};


//////// Inner product
template <typename T1, typename T2>
std::complex<double> operator% (const T1& t1, const T2& t2){
  if ( t1.size() != t2.size()) throw std::runtime_error("Inner product sizes not equal!");
  std::complex<double>result;
  for (size_t i =0; i<t1.size(); i++) {
    result += t1[i]*std::conj(t2[i]);
  }
  return result/(t1.size()*1.0);
}

Signal projection (const Component& v, Signal& u){
  std::complex<double> num = v%u;
  std::complex<double> den = u%u;
  std::cout<<u.ampl<<" Initially " <<std::endl;
  u.ampl *= num/den;
  for (size_t i=0;i<v.size();i++){
    u.data[i] = v[i]; 
  }
  std::cout<<u.ampl<<" Then " <<std::endl;
  return u;
}

Signal projection (const Signal& v, Signal& u){
  std::complex<double> num = v%u;
  std::complex<double> den = u%u;
  std::cout<<u.ampl<<" Initially " <<std::endl;
  if (u.ampl.real()!=0) {
    u.ampl *= num/den;
  }
  else {
    u.ampl = num/den;
  }
  std::cout<<u.ampl<<" Then " <<std::endl;
  return u;
}
//Projection of v on u
//template <typename T>
//Component projection (const T& v, Component u){
template <typename T1, typename T2>
Component projection (const T1& v,  T2& u){
  std::complex<double> num = v%u;
  std::complex<double> den = u%u;
  std::cout<<u.ampl<<" Initially " <<std::endl;
  if (u.ampl.real()!=0) {
    u.ampl *= num/den;
  }
  else {
    u.ampl = num/den;
  }
  std::cout<<u.ampl<<" Then " <<std::endl;
  return u;
}

template <typename T>
Component normalized (T& t1) {
  for (size_t i=0;i<t1.size();i++ ) {
    t1[i]/=std::norm(t1[i]);
  }
  return t1;
}

template <typename T>
Signal operator-(const Signal& t1, const T& t2) {
  if ( t1.size() != t2.size() ) throw std::runtime_error("Inner product sizes not equal!");
  Signal::data_t data;
  for (size_t i = 0; i<t1.size(); ++i) {
    data.emplace_back(t1[i]-t2[i]);
  }
  return Signal(data);
}

/*template <typename T>
Signal operator-(const Component& t1, const T& t2) {
  if ( t1.size() != t2.size() ) throw std::runtime_error("Inner product sizes not equal!");
  Signal::data_t data;
  for (size_t i = 0; i<t1.size(); ++i) {
    data.emplace_back(t1[i]-t2[i]);
  }
  return Signal(data);
}*/


/*template <typename T>
Signal operator-(const Signal& t1, const T& t2) {
  if ( t1.size() != t2.size() ) throw std::runtime_error("Inner product sizes not equal!");
  Signal::data_t data;
  for (size_t i = 0; i<t1.size(); ++i) {
    data.emplace_back(t1[i]-t2[i]);
  }
  return Signal(data);
}*/

std::vector<std::complex<double>> operator* (const std::complex<double> &ck, const std::vector<std::complex<double>> &exp_){
  std::vector<std::complex<double>>result;
  for (size_t i =0; i<exp_.size(); i++) {
    result.push_back(ck*exp_[i]);
  }
  return result;
}



/*/////Subtract previous frequency
std::vector <std::complex<double>> & operator -(std::vector<std::complex<double>> &s, const Component &proj) {
  std::ofstream myfile2;    
  std::string name = "/home/skostogl/cuTrack/signal.dat";
  myfile2.open(name);
  for (size_t i = 0; i<s.size(); i++) {
    ///
    myfile2<<s[i].real()<<" "<<s[i].imag()<<" ";
    ///
    s[i] -= proj.ampl*proj[i];
    s[i] -= proj.ampl*conj(proj[i]);
    ///
    myfile2<<(proj.ampl*proj[i]).real()<<" "<<(proj.ampl*proj[i]).imag()<<" "<<s[i].real()<<" "<<s[i].imag()<<std::endl;
  }
  myfile2.close();
  return s;
}*/

////////Print component
std::ostream& operator <<(std::ostream& os, Component c) {
  os<<"Frequency: "<<c.freq <<" Amplitude: " <<c.ampl<<" ABS amplitude: "<<abs(c.ampl) <<std::endl;
  return os;
}


//////// Print signal
std::ostream& operator <<(std::ostream& os, Signal s) {
  for (auto i:s.data){
    os<<i<<std::endl;
  }
  return os;
}
