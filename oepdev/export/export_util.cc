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

    m.def("calculate_Kij", &oepdev::calculate_Kij, "Calculatr Kij exchange matrix", py::return_value_policy::take_ownership);

}

} // EndNameSpace psi
