#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/constants/constants.hpp>

#include "NAFF.h"

BOOST_PYTHON_MODULE(NAFF)
{
  using namespace boost::python;

  class_<std::vector<Tfloat>>("Vec_cpp")
	.def(vector_indexing_suite<std::vector<Tfloat>>())
  ;	 

  class_<NAFF>("NAFF")
	.def("get_f1",&NAFF::get_f1)
 	.def("get_f",&NAFF::get_f)
  ;


}


