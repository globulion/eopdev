// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "../libutil/cphf.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_cphf(py::module &m) {

    // Function pointer types
    typedef std::shared_ptr<psi::Matrix> (oepdev::CPHF::*mat_void)(void) const;
    typedef std::shared_ptr<psi::Matrix> (oepdev::CPHF::*mat_int )(int ) const;

    /* Class oepdev::CPHF */
    py::class_<oepdev::CPHF, std::shared_ptr<oepdev::CPHF>>(m, "CPHF", "docstring") 
	.def(py::init<std::shared_ptr<psi::Wavefunction>, psi::Options&>())
	.def("compute", &oepdev::CPHF::compute, "")
	.def("nocc", &oepdev::CPHF::nocc, "", py::return_value_policy::copy)
	.def("print", &oepdev::CPHF::print, "")
	.def("polarizability", mat_void(&oepdev::CPHF::polarizability), "")
	.def("polarizability", mat_int (&oepdev::CPHF::polarizability), "")
	    ;
}


} // EndNameSpace psi
