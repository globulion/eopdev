// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "../lib3d/dmtp.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

std::shared_ptr<oepdev::DMTPole> get_camm(std::shared_ptr<psi::Wavefunction> wfn) {
   std::shared_ptr<oepdev::DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn);
   dmtp->compute();
   return dmtp;
}

void export_dmtp(py::module &m) {

    // Function pointer types
    typedef std::shared_ptr<oepdev::DMTPole> (*build_dmtp)(const std::string&, psi::SharedWavefunction, int);
    typedef void (oepdev::DMTPole::*compute_default)(void);
    typedef void (oepdev::DMTPole::*compute_dmatrix)(std::vector<psi::SharedMatrix>, std::vector<bool>);
    typedef psi::SharedMatrix (oepdev::DMTPole::*get_matrix_int)(int) const;
    typedef std::vector<psi::SharedMatrix> (oepdev::DMTPole::*get_vecmat)(void) const;
    typedef void (oepdev::DMTPole::*set_matrix_int)(psi::SharedMatrix, int);
    typedef void (oepdev::DMTPole::*set_vecmat)(std::vector<psi::SharedMatrix>);
    typedef void (oepdev::DMTPole::*recenter_all)(psi::SharedMatrix);

    /* Function get_camm */
    m.def("get_camm", &get_camm, "Compute CAMM", py::return_value_policy::take_ownership);

    /* Class oepdev::MultipoleConvergence */
    py::class_<oepdev::MultipoleConvergence, std::shared_ptr<oepdev::MultipoleConvergence>> MC(m, "MultipoleConvergence", "Handles the convergence of the distributed multipole expansions up to hexadecapole.");
    MC
	.def(py::init<std::shared_ptr<oepdev::DMTPole>, std::shared_ptr<oepdev::DMTPole>, oepdev::MultipoleConvergence::ConvergenceLevel>())
        .def("compute", &oepdev::MultipoleConvergence::compute, "Compute generalized property based on distributed multipole expansion sets")
	.def("level", &oepdev::MultipoleConvergence::level, "", py::return_value_policy::take_ownership)
    ;

    py::enum_<oepdev::MultipoleConvergence::ConvergenceLevel>(MC, "ConvergenceLevel", "Level of convergence - Rn series")
	.value("R1", oepdev::MultipoleConvergence::ConvergenceLevel::R1)
	.value("R2", oepdev::MultipoleConvergence::ConvergenceLevel::R2)
	.value("R3", oepdev::MultipoleConvergence::ConvergenceLevel::R3)
	.value("R4", oepdev::MultipoleConvergence::ConvergenceLevel::R4)
	.value("R5", oepdev::MultipoleConvergence::ConvergenceLevel::R5)
	.export_values();

    py::enum_<oepdev::MultipoleConvergence::Property>(MC, "Property", "Property to evalueate (energy or potential)")
	.value("Energy"   , oepdev::MultipoleConvergence::Property::Energy)
	.value("Potential", oepdev::MultipoleConvergence::Property::Potential)
	.export_values();


    /* Class oepdev::DMTPole */
    py::class_<oepdev::DMTPole, std::shared_ptr<oepdev::DMTPole>>(m, "DMTPole", "Set of Distributed Multipole Expansions")
	.def_static("build"            , build_dmtp     (&oepdev::DMTPole::build )                    , "Build DMTP set of distributions. Use one of the compute methods to complete the build.", py::return_value_policy::take_ownership)
	.def       ("compute"          , compute_default(&oepdev::DMTPole::compute)                   , "Compute CAMM's")
	.def       ("compute"          , compute_dmatrix(&oepdev::DMTPole::compute)                   , "Compute generalized CAMM's from user-provided density matrix")
	.def       ("centres"          ,                 &oepdev::DMTPole::centres                    , "Get DMTP centres")
	.def       ("origins"          ,                 &oepdev::DMTPole::origins                    , "Get DMTP origins")
	.def       ("centre"           ,                 &oepdev::DMTPole::centre                     , "Get DMTP centre")
	.def       ("origin"           ,                 &oepdev::DMTPole::origin                     , "Get DMTP origin")
	.def       ("charges"          , get_matrix_int (&oepdev::DMTPole::charges)                   , "Get DMTP charges")
	.def       ("dipoles"          , get_matrix_int (&oepdev::DMTPole::dipoles)                   , "Get DMTP dipoles")
	.def       ("quadrupoles"      , get_matrix_int (&oepdev::DMTPole::quadrupoles)               , "Get DMTP quadrupoles")
	.def       ("octupoles"        , get_matrix_int (&oepdev::DMTPole::octupoles)                 , "Get DMTP octupoles")
	.def       ("hexadecapoles"    , get_matrix_int (&oepdev::DMTPole::hexadecapoles)             , "Get DMTP hexadecapoles")
	.def       ("charges"          , get_vecmat     (&oepdev::DMTPole::charges)                   , "Get set of DMTP charges")
	.def       ("dipoles"          , get_vecmat     (&oepdev::DMTPole::dipoles)                   , "Get set of DMTP dipoles")
	.def       ("quadrupoles"      , get_vecmat     (&oepdev::DMTPole::quadrupoles)               , "Get set of DMTP quadrupoles")
	.def       ("octupoles"        , get_vecmat     (&oepdev::DMTPole::octupoles)                 , "Get set of DMTP octupoles")
	.def       ("hexadecapoles"    , get_vecmat     (&oepdev::DMTPole::hexadecapoles)             , "Get set of DMTP hexadecapoles")
	.def       ("has_charges"      ,                 &oepdev::DMTPole::has_charges                , "", py::return_value_policy::copy)       
	.def       ("has_dipoles"      ,                 &oepdev::DMTPole::has_dipoles                , "", py::return_value_policy::copy)
	.def       ("has_quadrupoles"  ,                 &oepdev::DMTPole::has_quadrupoles            , "", py::return_value_policy::copy)
	.def       ("has_octupoles"    ,                 &oepdev::DMTPole::has_octupoles              , "", py::return_value_policy::copy)
	.def       ("has_hexadecapoles",                 &oepdev::DMTPole::has_hexadecapoles          , "", py::return_value_policy::copy)
	.def       ("set_charges"          , get_matrix_int (&oepdev::DMTPole::charges)                   , "Get DMTP charges", py::return_value_policy::copy)
	.def       ("dipoles"          , get_matrix_int (&oepdev::DMTPole::dipoles)                   , "Get DMTP dipoles", py::return_value_policy::copy)
	.def       ("set_quadrupoles"      , set_matrix_int (&oepdev::DMTPole::set_quadrupoles)               , "Set DMTP quadrupoles")
	.def       ("set_octupoles"        , set_matrix_int (&oepdev::DMTPole::set_octupoles)                 , "Set DMTP octupoles")
	.def       ("set_hexadecapoles"    , set_matrix_int (&oepdev::DMTPole::set_hexadecapoles)             , "Set DMTP hexadecapoles")
	.def       ("set_charges"          , set_vecmat     (&oepdev::DMTPole::set_charges)                   , "Set set of DMTP charges")
	.def       ("set_dipoles"          , set_vecmat     (&oepdev::DMTPole::set_dipoles)                   , "Set set of DMTP dipoles")
	.def       ("set_quadrupoles"      , set_vecmat     (&oepdev::DMTPole::set_quadrupoles)               , "Set set of DMTP quadrupoles")
	.def       ("set_octupoles"        , set_vecmat     (&oepdev::DMTPole::set_octupoles)                 , "Set set of DMTP octupoles")
	.def       ("set_hexadecapoles"    , set_vecmat     (&oepdev::DMTPole::set_hexadecapoles)             , "Set set of DMTP hexadecapoles")

  	.def       ("n_sites"          ,                 &oepdev::DMTPole::n_sites                    , "Number of DMTP sites", py::return_value_policy::copy)
  	.def       ("n_dmtp"           ,                 &oepdev::DMTPole::n_dmtp                     , "Number of DMTP sets", py::return_value_policy::copy)
  	.def       ("recenter"         , recenter_all(   &oepdev::DMTPole::recenter)                  , "Change origins of the distributed multipole moments of all sets")
	.def       ("translate"        ,                 &oepdev::DMTPole::translate                  , "Translate the DMTP sets")
	.def       ("rotate"           ,                 &oepdev::DMTPole::rotate                     , "Rotate the DMTP sets")
	.def       ("superimpose"      ,                 &oepdev::DMTPole::superimpose                , "Superimpose the DMTP sets")
	.def       ("energy"           ,                 &oepdev::DMTPole::energy                     , "Evaluate the generalized energy", py::return_value_policy::take_ownership)
	.def       ("potential"        ,                 &oepdev::DMTPole::potential                  , "Evaluate the generalized potential", py::return_value_policy::take_ownership)
	.def       ("print_header"     ,                 &oepdev::DMTPole::print_header               , "Print header info to output file")
	.def       ("print"            ,                 &oepdev::DMTPole::print                      , "Print the multipoles to output file")
        ;
}


} // EndNameSpace psi
