#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "lattice.h"


BOOST_PYTHON_MODULE(tracker)
{
  using namespace boost::python;
  
  class_< std::vector<Tfloat> >("vec_Tfloat")
    .def(vector_indexing_suite< std::vector<Tfloat> >())
  ;

  class_< std::vector<HostBunch> >("vec_HostBunch")
    .def(vector_indexing_suite< std::vector<HostBunch> >())
    .def("x",  &get_x_from_bunches ) 
    .def("xp", &get_xp_from_bunches)
    .def("y",  &get_y_from_bunches )
    .def("yp", &get_yp_from_bunches)
    .def("d",  &get_d_from_bunches )
    .def("z",  &get_z_from_bunches )
  ;

  void (Lattice::*lattice_track_H)(HostBunch&)   = &Lattice::track;
  void (Lattice::*lattice_track_D)(DeviceBunch&) = &Lattice::track;

  class_<Lattice, boost::noncopyable>("Lattice")
    .def(init<const std::string &>())
    .def("read_twiss_table", &Lattice::read_twiss_table)
    .def("optimise", &Lattice::optimise)
    .def("add", &Lattice::add)
    .def("n_elements", &Lattice::get_n_elements)
    .def("compile", &Lattice::compile)
    .def("write_ptx", &Lattice::write_ptx)
    .def("read_ptx", &Lattice::read_ptx)
    .def("make_matched_bunch", &Lattice::make_matched_bunch)
    .def("track", lattice_track_H)
    .def("track", lattice_track_D)
    .def("sigma_x", &Lattice::sigma_x)
    .def("sigma_y", &Lattice::sigma_y)
    .def("rel_gamma", &Lattice::rel_gamma)
    .def_readwrite("turns", &Lattice::turn_by_turn_data)
    .def_readwrite("n_turns", &Lattice::n_turns)
    .def_readwrite("collect_tbt_data", &Lattice::collect_tbt_data)
    .def_readwrite("BETX", &Lattice::BETX)
    .def_readwrite("BETY", &Lattice::BETY)
    .def_readwrite("ALFX", &Lattice::ALFX)
    .def_readwrite("ALFY", &Lattice::ALFY)
    .def_readwrite("DX"  , &Lattice::DX)
    .def_readwrite("DY"  , &Lattice::DY)
    .def_readwrite("DPX" , &Lattice::DPX)
    .def_readwrite("DPY" , &Lattice::DPY)
    .def_readwrite("X"   , &Lattice::X) //closed orbit
    .def_readwrite("Y"   , &Lattice::Y)
    .def_readwrite("PX"  , &Lattice::PX)
    .def_readwrite("PY"  , &Lattice::PY)
    .def_readwrite("norm_emit_x", &Lattice::norm_emit_x)
    .def_readwrite("norm_emit_y", &Lattice::norm_emit_y)
    .def_readwrite("bunch_length", &Lattice::bunch_length)
    .def_readwrite("bunch_energy_spread", &Lattice::bunch_energy_spread)
    .def_readwrite("energy", &Lattice::energy)
    .def_readwrite("mass", &Lattice::mass)
  ;


  class_<HostBunch>("HostBunch", init<const size_t>())
    .def("size", &HostBunch::size)
    .def("fromDeviceBunch", &HostBunch::copyFromDeviceBunch)
    .def("toDeviceBunch", &HostBunch::copyToDeviceBunch)
    .def_readwrite("x" , &HostBunch::x )
    .def_readwrite("xp", &HostBunch::xp)
    .def_readwrite("y" , &HostBunch::y )
    .def_readwrite("yp", &HostBunch::yp)
    .def_readwrite("z" , &HostBunch::z )
    .def_readwrite("d" , &HostBunch::d )
  ;
 


  class_<DeviceBunch>("DeviceBunch",init<const HostBunch &>())
  ;
  

}

