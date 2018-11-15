// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "../libsolver/solver.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_solver(py::module &m) {

    // Function pointer types

    /* Class oepdev::OEPDevSolver */
    py::class_<oepdev::OEPDevSolver, std::shared_ptr<oepdev::OEPDevSolver>> Solver(m, "OEPDevSolver", "Implements solver routines of OEPDev.");
    Solver
	.def_static("build", &oepdev::OEPDevSolver::build, "Build a chosen OEPDev solver", py::return_value_policy::take_ownership)
	.def("compute_oep_based",   &oepdev::OEPDevSolver::compute_oep_based  , "Run the OEP-based method")
	.def("compute_benchmark",   &oepdev::OEPDevSolver::compute_benchmark  , "Run the benchmark method")
        ;

    /* Class oepdev::WavefunctionUnion */
    py::class_<oepdev::WavefunctionUnion, std::shared_ptr<oepdev::WavefunctionUnion>> WFN(m, "WavefunctionUnion", "The Hartree-Product of multiple wavefunctions");
    WFN
         .def(py::init<std::shared_ptr<psi::Wavefunction>, psi::Options&>())
    	 .def("transform_integrals",  &oepdev::WavefunctionUnion::transform_integrals, "Transform integrals to MO basis")
        ;
}


} // EndNameSpace psi
