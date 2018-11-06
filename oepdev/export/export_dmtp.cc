// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "../lib3d/dmtp.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;
namespace psi{

std::shared_ptr<oepdev::DMTPole> get_camm(std::shared_ptr<psi::Wavefunction> wfn) {
   std::shared_ptr<oepdev::DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn);
   dmtp->compute();
   return dmtp;
}

void export_dmtp(py::module &m) {
    m.def("get_camm", &get_camm, "Compute CAMM");

    typedef std::shared_ptr<oepdev::DMTPole> (*build_dmtp)(const std::string&, std::shared_ptr<psi::Wavefunction>, int);
    typedef void (oepdev::DMTPole::*compute_default)(void);
    typedef void (oepdev::DMTPole::*compute_dmatrix)(std::shared_ptr<psi::Matrix>, bool, int);
    typedef std::shared_ptr<psi::Matrix> (oepdev::DMTPole::*get_charges_mat)(int) const;

    py::class_<oepdev::DMTPole, std::shared_ptr<oepdev::DMTPole>>(m, "DMTPole", "docstring")
	.def_static("build"  , build_dmtp     (&oepdev::DMTPole::build ) , "docstring")
	.def       ("compute", compute_default(&oepdev::DMTPole::compute), "docstring")
	.def       ("compute", compute_dmatrix(&oepdev::DMTPole::compute), "docstring")
	.def       ("centres",                 &oepdev::DMTPole::centres , "docstring")
	.def       ("origins",                 &oepdev::DMTPole::origins , "docstring")
	.def       ("charges", get_charges_mat(&oepdev::DMTPole::charges), "docstring");

}


} // EndNameSpace psi
