// PyBind11
#include "psi4/libmints/wavefunction.h"
#include "../liboep/oep_gdf.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py = pybind11;
namespace psi{

void export_oep(py::module &m) {

    // Function pointer types
    using SharedBasisSet = std::shared_ptr<psi::BasisSet>;
    using SharedMatrix   = std::shared_ptr<psi::Matrix>;
    typedef std::shared_ptr<oepdev::GeneralizedDensityFit> (*build_gdf_single)(                SharedBasisSet, SharedMatrix);
    typedef std::shared_ptr<oepdev::GeneralizedDensityFit> (*build_gdf_double)(SharedBasisSet, SharedBasisSet, SharedMatrix);

    /* Class oepdev::GeneralizedDensityFit */
    py::class_<oepdev::GeneralizedDensityFit, std::shared_ptr<oepdev::GeneralizedDensityFit>> GDF(m, "GeneralizedDensityFit", "Generalized Density Fitting Method");
    GDF
	.def_static("build_single", build_gdf_single(&oepdev::GeneralizedDensityFit::build), "Build GDF calculator. Single GDF scheme.", py::return_value_policy::take_ownership)
	.def_static("build_double", build_gdf_double(&oepdev::GeneralizedDensityFit::build), "Build GDF calculator. Double GDF scheme.", py::return_value_policy::take_ownership)
        .def("compute", &oepdev::GeneralizedDensityFit::compute, "Compute GDF matrix")
        .def("G", &oepdev::GeneralizedDensityFit::G, "Return GDF matrix", py::return_value_policy::take_ownership)
    ;

}


} // EndNameSpace psi
