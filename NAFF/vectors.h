#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#pragma once

class Vectors {
  public:
  std::vector<std::complex<double>> vec;
  std::vector<Vectors> basis_func;
  double f_counter=0; 
  double wmega;
  size_t signal_size;
  ~Vectors() {}
  Vectors()  {}
  Vectors(std::vector<std::complex<double>> vec, double wmega): vec(vec), wmega(wmega) {cmp_basis_func(wmega,vec.size());}
  Vectors(double wmega,size_t signal_size): wmega(wmega), signal_size(signal_size) {cmp_basis_func(wmega, signal_size);}
 
  size_t size() const {
    return this->vec.size();
  } 
  
  Vectors operator+(const Vectors& vec) {
    Vectors vector;
    for (size_t i=0;i<vec.size();i++) {
      vector.vec.push_back(this->vec[i]+ vec.vec[i]);
    }
    return vector;
  }
  
  Vectors operator-(const Vectors& vec) {
    Vectors vector;
    for (size_t i=0;i<vec.size();i++) {
      vector.vec.push_back(this->vec[i]- vec.vec[i]);
    }
    return vector;
  }
  
  Vectors operator*(const Vectors &vec) {
    Vectors vector;
    for (size_t i=0;i<vec.size();i++) {
      vector.vec.push_back(this->vec[i]*vec.vec[i]);
    }
    return vector;
  }
  
   
  Vectors conjugate(const Vectors &vec) {
    Vectors vector;
    for (size_t i=0;i<vec.size();i++) {
      vector.vec.push_back(conj(this->vec[i]));
    }
    return vector;
  }

  std::complex<double> inproduct(const Vectors &vec) {
    Vectors vector;
    vector=(conjugate(*this))*vec;
    return std::accumulate(vector.vec.begin(), vector.vec.end(), std::complex<double>(0.0));
  }
  
  std::complex<double> amplitude(const Vectors &vec) {
    return (inproduct(vec))/(inproduct(*this));
  }

  Vectors projection(const Vectors &vec) {
    Vectors vector;
    for (size_t i=0;i<vec.size();i++){
      vector.vec.push_back(amplitude(vec)*vec.vec[i]);
    }
    return vector;
  }
  
  void cmp_basis_func (double frequency, size_t n) {
    std::complex<double> j(0,1);
    for (size_t i=0;i<n;i++){
      basis_func[f_counter].vec.push_back(exp(i*frequency*j));
  }
    f_counter++;
 }
  
 /* void orthog () {
    size_t j=2;
    while (j<f_counter){
      basis_func[f_counter]=basis_func[f_counter]-projection(basis_func[j-1])
    }

  }*/

  std::ostream & operator<<( std::ostream & out) {
   for (size_t j=0;j<this->size();j++){
      out << this->vec[j]<<std::endl;
    }
    return out;
  }
};


