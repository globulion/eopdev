#include "oep.h"
namespace oepdev{

OEPotential::OEPotential(SharedWavefunction wfn, Options& options) 
   : wfn_(wfn),
     options_(options),
     is_density_fitted(false),
     is_esp_based(true)
{
    common_init();
}

OEPotential::OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options) 
   : wfn_(wfn),
     auxiliary_(auxiliary),
     options_(options),
     is_density_fitted(true),
     is_esp_based(false)
{
    common_init();
}

std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "REPULSION ENERGY")  oep = std::make_shared<RepulsionEnergyOEPotential>(wfn, options);
   else if (category == "EET COUPLING"    )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}

std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, 
                                                                           SharedBasisSet auxiliary, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "REPULSION ENERGY")  oep = std::make_shared<RepulsionEnergyOEPotential>(wfn, auxiliary, options);
   else if (category == "EET COUPLING"    )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, auxiliary, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}



OEPotential::~OEPotential() {}

void OEPotential::common_init(void) 
{
   name_ = "default";
}

void OEPotential::compute(const std::string& oepType) {}
void OEPotential::compute(void) {
   for (int i=0; i<oepTypes_.size(); ++i) {
        this->compute(oepTypes_[i]);
   }
}
void OEPotential::compute_3D(const std::string& oepType, const std::string& fileName) {}
void OEPotential::rotate(const Matrix& rotmat) {}
void OEPotential::translate(const Vector& trans) {}
void OEPotential::superimpose(const Matrix& refGeometry,
                              const std::vector<int>& supList,
                              const std::vector<int>& reordList) {}
void OEPotential::print_header(void) const {}



// <============== Repulsion Energy ==============> //

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, 
                                                       SharedBasisSet auxiliary, Options& options) 
 : OEPotential(wfn, auxiliary, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::~RepulsionEnergyOEPotential() {}
void RepulsionEnergyOEPotential::common_init() 
{
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) {}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const std::string& fileName) {}
void RepulsionEnergyOEPotential::print_header(void) const {}

// <============== EET Coupling ==============> //

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, 
                                               SharedBasisSet auxiliary, Options& options) 
 : OEPotential(wfn, auxiliary, options)
{ 
   common_init();
}

EETCouplingOEPotential::~EETCouplingOEPotential() {}
void EETCouplingOEPotential::common_init() 
{
}

void EETCouplingOEPotential::compute(const std::string& oepType) {}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const std::string& fileName) {}
void EETCouplingOEPotential::print_header(void) const {}

//OEPotential::compute_3D(std::string oepType, std::string fileName) {
//  std::shared_ptr<CubicScalarGrid> grid(new CubicScalarGrid(auxiliary_, options_));
//  grid->write_cube_file(v_oep, oepType);
//}

} // EndNameSpace oepdev
