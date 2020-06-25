// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "../libgefp/gefp.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_gefp(py::module &m) {

    // Function pointer types
    typedef psi::SharedMatrix (oepdev::GenEffPar::*susceptibility_4)(int, int, int, int) const;
    typedef std::vector<psi::SharedMatrix> (oepdev::GenEffPar::*susceptibility_3)(int, int, int) const;
    typedef std::vector<std::vector<psi::SharedMatrix>> (oepdev::GenEffPar::*susceptibility_2)(int, int) const;
    typedef std::vector<std::vector<psi::SharedMatrix>> (oepdev::GenEffPar::*dpol_0)() const; 
    typedef std::vector<psi::SharedMatrix> (oepdev::GenEffPar::*dpol_1)(int) const; 
    typedef psi::SharedMatrix (oepdev::GenEffPar::*dpol_2)(int, int) const; 
    typedef psi::SharedMatrix (oepdev::GenEffPar::*cdm_f)(psi::SharedVector);
    typedef psi::SharedMatrix (oepdev::GenEffPar::*cdm_3)(double, double, double);
    typedef psi::SharedMatrix (oepdev::GenEffPar::*cdm_F)(std::vector<psi::SharedVector>);
    typedef psi::SharedMatrix (oepdev::GenEffPar::*cdm_FG)(std::vector<psi::SharedVector>, std::vector<psi::SharedMatrix>);
    typedef std::shared_ptr<oepdev::GenEffParFactory> (*build_1)(const std::string&, psi::SharedWavefunction, psi::Options&);
    typedef std::shared_ptr<oepdev::GenEffParFactory> (*build_2)(const std::string&, psi::SharedWavefunction, psi::Options&, psi::SharedBasisSet, psi::SharedBasisSet);

    /* Class oepdev::FragmentedSystem */
    py::class_<oepdev::FragmentedSystem, std::shared_ptr<oepdev::FragmentedSystem>> System(m, "FragmentedSystem", "GEFP System");
    System
	.def_static("build", &oepdev::FragmentedSystem::build, "Build a fragmented system", 
                               py::return_value_policy::take_ownership)
        .def("set_geometry", &oepdev::FragmentedSystem::set_geometry, "Set geometry of all fragments")
        .def("set_primary", &oepdev::FragmentedSystem::set_primary, "Set primary basis sets")
        .def("set_auxiliary", &oepdev::FragmentedSystem::set_auxiliary, "Set auxiliary basis sets")
        .def("superimpose", &oepdev::FragmentedSystem::superimpose, "Superimpose")
        .def("compute_energy", &oepdev::FragmentedSystem::compute_energy, "Compute energy")
        .def("compute_energy_term", &oepdev::FragmentedSystem::compute_energy_term, "Compute one of energy terms")
        ;

    /* Class oepdev::GenEffFrag */
    py::class_<oepdev::GenEffFrag, std::shared_ptr<oepdev::GenEffFrag>> Fragment(m, "GenEffFrag", "GEFP Fragment");
    Fragment
	.def_static("build", &oepdev::GenEffFrag::build, "Build a fragment", 
                               py::return_value_policy::take_ownership)
        .def("clone", &oepdev::GenEffFrag::clone, "Copy the fragment", py::return_value_policy::take_ownership)
        .def("set_nbf", &oepdev::GenEffFrag::set_nbf, "Set the number of primary basis functions")
        .def("set_ndocc", &oepdev::GenEffFrag::set_ndocc, "")
        .def("set_molecule", &oepdev::GenEffFrag::set_molecule, "")
        .def("set_basisset", &oepdev::GenEffFrag::set_basisset, "")
        .def("set_parameters", &oepdev::GenEffFrag::set_parameters, "")
        .def_readwrite("parameters", &oepdev::GenEffFrag::parameters, "")
        ;


    /* Class oepdev::GenEffParFactory */
    py::class_<oepdev::GenEffParFactory, std::shared_ptr<oepdev::GenEffParFactory>> Solver(m, "GenEffParFactory", "GEFP factory of OEPDev.");
    Solver
	.def_static("build", build_1(&oepdev::GenEffParFactory::build), "Build a chosen GEFP solver", py::return_value_policy::take_ownership)
	.def_static("build", build_2(&oepdev::GenEffParFactory::build), "Build a chosen GEFP solver for OEP-based purposes", py::return_value_policy::take_ownership)
	.def("compute",   &oepdev::GenEffParFactory::compute  , "Run the GEFP calculation", py::return_value_policy::take_ownership)
        .def("wfn", &oepdev::GenEffParFactory::wfn, "Retrieve wfn")
        .def("options", &oepdev::GenEffParFactory::options, "Retrieve options")
        .def("cphf_solver", &oepdev::GenEffParFactory::cphf_solver, "Retrieve CPHF object")
        ;

    /* Class oepdev::GenEffPar */
    py::class_<oepdev::GenEffPar, std::shared_ptr<oepdev::GenEffPar>> Parameters(m, "GenEffPar", "GEFP parameters of OEPDev.");
    Parameters
        .def("hasDensityMatrixDipolePolarizability", &oepdev::GenEffPar::hasDensityMatrixDipolePolarizability, "")
        .def("hasDensityMatrixDipoleDipoleHyperpolarizability", &oepdev::GenEffPar::hasDensityMatrixDipoleDipoleHyperpolarizability, "")
        .def("hasDensityMatrixQuadrupolePolarizability", &oepdev::GenEffPar::hasDensityMatrixQuadrupolePolarizability, "")
        .def("centres",  &oepdev::GenEffPar::centres, "")
        .def("centre",  &oepdev::GenEffPar::centre, "")
        .def("clone",  &oepdev::GenEffPar::clone, "Make a deep copy", py::return_value_policy::take_ownership)
        .def("susceptibility", susceptibility_2(&oepdev::GenEffPar::susceptibility), "", py::return_value_policy::take_ownership)
        .def("susceptibility", susceptibility_3(&oepdev::GenEffPar::susceptibility), "", py::return_value_policy::take_ownership)
        .def("susceptibility", susceptibility_4(&oepdev::GenEffPar::susceptibility), "", py::return_value_policy::take_ownership)
        .def("dipole_polarizability", dpol_0(&oepdev::GenEffPar::dipole_polarizability), "", py::return_value_policy::take_ownership)
        .def("dipole_polarizability", dpol_1(&oepdev::GenEffPar::dipole_polarizability), "", py::return_value_policy::take_ownership)
        .def("dipole_polarizability", dpol_2(&oepdev::GenEffPar::dipole_polarizability), "", py::return_value_policy::take_ownership)
        .def("dipole_dipole_hyperpolarizability", dpol_0(&oepdev::GenEffPar::dipole_dipole_hyperpolarizability), "", py::return_value_policy::take_ownership)
        .def("dipole_dipole_hyperpolarizability", dpol_1(&oepdev::GenEffPar::dipole_dipole_hyperpolarizability), "", py::return_value_policy::take_ownership)
        .def("dipole_dipole_hyperpolarizability", dpol_2(&oepdev::GenEffPar::dipole_dipole_hyperpolarizability), "", py::return_value_policy::take_ownership)
        .def("quadrupole_polarizability", dpol_0(&oepdev::GenEffPar::quadrupole_polarizability), "", py::return_value_policy::take_ownership)
        .def("quadrupole_polarizability", dpol_1(&oepdev::GenEffPar::quadrupole_polarizability), "", py::return_value_policy::take_ownership)
        .def("quadrupole_polarizability", dpol_2(&oepdev::GenEffPar::quadrupole_polarizability), "", py::return_value_policy::take_ownership)
        .def("compute_density_matrix", cdm_f(&oepdev::GenEffPar::compute_density_matrix), "", py::return_value_policy::take_ownership)
        .def("compute_density_matrix", cdm_3(&oepdev::GenEffPar::compute_density_matrix), "", py::return_value_policy::take_ownership)
        .def("compute_density_matrix", cdm_F(&oepdev::GenEffPar::compute_density_matrix), "", py::return_value_policy::take_ownership)
        .def("compute_density_matrix", cdm_FG(&oepdev::GenEffPar::compute_density_matrix), "", py::return_value_policy::take_ownership)
        ;




}
} // EndNameSpace psi
