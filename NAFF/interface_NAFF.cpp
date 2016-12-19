#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/constants/constants.hpp>
#include <fftw3.h>

#include "NAFF.h"

BOOST_PYTHON_MODULE(NAFF)
{
  using namespace boost::python;

  class_<std::vector<Tfloat>>("Vec_cpp")
	  .def(vector_indexing_suite<std::vector<Tfloat>>())
   ;	 

  def("NAFF_f1",NAFF_f1);
  def("NAFF_f",NAFF_f);
  def("FFT",FFT);

}

