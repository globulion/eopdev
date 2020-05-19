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
	 .def(py::init<
		       std::shared_ptr<psi::Molecule>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
       		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
       		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::BasisSet>, 
		       std::shared_ptr<psi::Wavefunction>,
		       std::shared_ptr<psi::Wavefunction>, 
		       psi::Options&>())
    	 .def("transform_integrals",  &oepdev::WavefunctionUnion::transform_integrals, "Transform integrals to MO basis")
    	 .def("localize_orbitals",  &oepdev::WavefunctionUnion::localize_orbitals, "Localize occupied molecular orbitals")
	 .def("clear_dpd", &oepdev::WavefunctionUnion::clear_dpd, "Clear the DPD instance (necessary at the end of work)")
         .def("l_nmo", &oepdev::WavefunctionUnion::l_nmo, "")
         .def("l_nso", &oepdev::WavefunctionUnion::l_nso, "")
         .def("l_ndocc", &oepdev::WavefunctionUnion::l_ndocc, "")
         .def("l_nvir", &oepdev::WavefunctionUnion::l_nvir, "")
         .def("l_nalpha", &oepdev::WavefunctionUnion::l_nalpha, "")
         .def("l_nbeta", &oepdev::WavefunctionUnion::l_nbeta, "")
         .def("l_nbf", &oepdev::WavefunctionUnion::l_nbf, "")
         .def("l_noffs_ao", &oepdev::WavefunctionUnion::l_noffs_ao, "")
         .def("l_energy", &oepdev::WavefunctionUnion::l_energy, "")
         .def("l_molecule", &oepdev::WavefunctionUnion::l_molecule, "")
         .def("l_primary", &oepdev::WavefunctionUnion::l_primary, "")
         .def("l_auxiliary", &oepdev::WavefunctionUnion::l_auxiliary, "")
         .def("l_intermediate", &oepdev::WavefunctionUnion::l_intermediate, "")
         .def("l_guess", &oepdev::WavefunctionUnion::l_guess, "")
         .def("l_wfn", &oepdev::WavefunctionUnion::l_wfn, "")
         .def("l_localizer", &oepdev::WavefunctionUnion::l_localizer, "")
         .def("integrals", &oepdev::WavefunctionUnion::integrals, "")
         .def("has_localized_orbitals", &oepdev::WavefunctionUnion::has_localized_orbitals, "")
         .def("primary", &oepdev::WavefunctionUnion::primary, "")
         .def("Ca_subset", &oepdev::WavefunctionUnion::Ca_subset, "")
         .def("Cb_subset", &oepdev::WavefunctionUnion::Cb_subset, "")
         .def("C_subset_helper", &oepdev::WavefunctionUnion::C_subset_helper, "")
         .def("epsilon_subset_helper", &oepdev::WavefunctionUnion::epsilon_subset_helper, "")
         .def("print_header", &oepdev::WavefunctionUnion::print_header, "")
         .def("print_mo_integrals", &oepdev::WavefunctionUnion::print_mo_integrals, "")
        ;
}


} // EndNameSpace psi
