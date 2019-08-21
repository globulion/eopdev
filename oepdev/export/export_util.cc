// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "../libutil/util.h"
#include "../libutil/cis.h"
#include "../lib3d/dmtp.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_util(py::module &m) {

    typedef std::shared_ptr<psi::Matrix> (*dfi_jk)(std::shared_ptr<psi::IntegralFactory>, std::shared_ptr<psi::IntegralFactory>, std::shared_ptr<psi::Matrix>);
    typedef std::shared_ptr<psi::Matrix> (*dfi_j)(std::shared_ptr<psi::IntegralFactory>, std::shared_ptr<psi::Matrix>);

    m.def("calculate_JK", &oepdev::calculate_JK, "Calculate J and K matrices in MO basis", py::return_value_policy::take_ownership);
    m.def("calculate_JK_r", &oepdev::calculate_JK_r, "Calculate J and K matrices in MO basis", py::return_value_policy::take_ownership);
    m.def("calculate_der_D", &oepdev::calculate_der_D, "Calculate derivatives of E_XC wrt D in MO basis", py::return_value_policy::take_ownership);
    m.def("calculate_e_xc", &oepdev::calculate_e_xc, "Calculate exchange-correlation energy");
    m.def("matrix_power_derivative", &oepdev::matrix_power_derivative, "Calculate derivative of matrix power", py::return_value_policy::take_ownership);
    m.def("calculate_DFI_Vel_JK", &oepdev::calculate_DFI_Vel_JK, "Compute DFI electronic interaction matrix with J and K contributions", py::return_value_policy::take_ownership);
    m.def("calculate_DFI_Vel_J", &oepdev::calculate_DFI_Vel_J, "Compute DFI electronic interaction matrix with J contribution only", py::return_value_policy::take_ownership);

    /* Class oepdev::CISComputer */
    typedef psi::SharedMatrix (oepdev::CISComputer::*m_i)(int) const;
    typedef psi::SharedMatrix (oepdev::CISComputer::*m_ij)(int, int) const;
    typedef std::shared_ptr<oepdev::DMTPole> (oepdev::CISComputer::*trcamm_i)(int, bool) const;
    typedef std::shared_ptr<oepdev::DMTPole> (oepdev::CISComputer::*trcamm_ij)(int, int, bool) const;
    typedef psi::SharedVector (oepdev::CISComputer::*d_i)(int) const;
    typedef psi::SharedVector (oepdev::CISComputer::*d_ij)(int, int) const;
    typedef double (oepdev::CISComputer::*f_i)(int) const;
    typedef double (oepdev::CISComputer::*f_ij)(int, int) const;

    py::class_<oepdev::CISComputer, std::shared_ptr<oepdev::CISComputer>> CISComputer(m, "CISComputer", "Implements the CIS method.");
    CISComputer
	.def_static("build", &oepdev::CISComputer::build, "Build a chosen CIS solver", py::return_value_policy::take_ownership)
	.def("compute", &oepdev::CISComputer::compute, "Run the CIS calculations")
        .def("clear_dpd", &oepdev::CISComputer::clear_dpd, "Clear DPD instance by destroying the IntegralTransform object")
        .def("nstates", &oepdev::CISComputer::nstates, "Extract the number of excited states found")
	.def("eigenvalues", &oepdev::CISComputer::eigenvalues, "Extract CIS excitation energies wrt HF ground state", py::return_value_policy::take_ownership)
	.def("eigenvectors", &oepdev::CISComputer::eigenvectors, "Extract CIS amplitudes", py::return_value_policy::take_ownership)
	.def("E", &oepdev::CISComputer::E, "Extract CIS excitation energies wrt HF ground state", py::return_value_policy::take_ownership)
	.def("U", &oepdev::CISComputer::U, "Extract CIS amplitudes", py::return_value_policy::take_ownership)
	.def("U_homo_lumo", &oepdev::CISComputer::U_homo_lumo, "Extract CIS amplitudes for HOMO->LUMO excitation")
        .def("Da_mo", &oepdev::CISComputer::Da_mo, "", py::return_value_policy::take_ownership)
        .def("Db_mo", &oepdev::CISComputer::Db_mo, "", py::return_value_policy::take_ownership)
        .def("Da_ao", &oepdev::CISComputer::Da_ao, "", py::return_value_policy::take_ownership)
        .def("Db_ao", &oepdev::CISComputer::Db_ao, "", py::return_value_policy::take_ownership)
        .def("camm", &oepdev::CISComputer::camm, "", py::return_value_policy::take_ownership)
        .def("Ta_ao", m_i( &oepdev::CISComputer::Ta_ao), "", py::return_value_policy::take_ownership)
        .def("Ta_ao", m_ij(&oepdev::CISComputer::Ta_ao), "", py::return_value_policy::take_ownership)
        .def("Tb_ao", m_i( &oepdev::CISComputer::Tb_ao), "", py::return_value_policy::take_ownership)
        .def("Tb_ao", m_ij(&oepdev::CISComputer::Tb_ao), "", py::return_value_policy::take_ownership)
        .def("trcamm", trcamm_i( &oepdev::CISComputer::trcamm), "", py::return_value_policy::take_ownership)
        .def("trcamm", trcamm_ij(&oepdev::CISComputer::trcamm), "", py::return_value_policy::take_ownership)
        .def("transition_dipole", d_i( &oepdev::CISComputer::transition_dipole), "", py::return_value_policy::take_ownership)
        .def("transition_dipole", d_ij(&oepdev::CISComputer::transition_dipole), "", py::return_value_policy::take_ownership)
        .def("oscillator_strength", f_i( &oepdev::CISComputer::oscillator_strength), "", py::return_value_policy::take_ownership)
        .def("oscillator_strength", f_ij(&oepdev::CISComputer::oscillator_strength), "", py::return_value_policy::take_ownership)
        ;

}

} // EndNameSpace psi
