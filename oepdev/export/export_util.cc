// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "../libutil/util.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_util(py::module &m) {

    m.def("calculate_JK", &oepdev::calculate_JK, "Calculatr J and K matrices in MO basis", py::return_value_policy::take_ownership);
    m.def("calculate_JK_r", &oepdev::calculate_JK_r, "Calculatr J and K matrices in MO basis", py::return_value_policy::take_ownership);

}

} // EndNameSpace psi
