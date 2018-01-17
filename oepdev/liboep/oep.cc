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

   if      (category == "ELECTROSTATIC ENERGY")  oep = std::make_shared<ElectrostaticEnergyOEPotential>(wfn, options);
   else if (category == "REPULSION ENERGY"    )  oep = std::make_shared<    RepulsionEnergyOEPotential>(wfn, options);
   else if (category == "EET COUPLING"        )  oep = std::make_shared<        EETCouplingOEPotential>(wfn, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}

std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, 
                                                                           SharedBasisSet auxiliary, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY")  throw PSIEXCEPTION("OEPDEV: OEPotential build. DF-based OEP not available for Electrostatic Energy!");
   else if (category == "REPULSION ENERGY"    )  oep = std::make_shared<RepulsionEnergyOEPotential>(wfn, auxiliary, options);
   else if (category == "EET COUPLING"        )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, auxiliary, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}



OEPotential::~OEPotential() {}

void OEPotential::common_init(void) 
{
   name_        = "default";
   primary_     = wfn_->basisset();
   intsFactory_ = std::make_shared<psi::IntegralFactory>(primary_, primary_);
   potMat_      = std::make_shared<psi::Matrix>("Potential Integrals", primary_->nbf(), primary_->nbf());
   potInt_      = std::make_shared<PotentialInt>(intsFactory_->spherical_transform(), primary_, primary_, 0);
}

void OEPotential::compute(const std::string& oepType) {}
void OEPotential::compute(void) { for ( const std::string& type : oepTypes_ ) this->compute(type); }
void OEPotential::write_cube(const std::string& oepType, const std::string& fileName) 
{
   OEPotential3D<OEPotential> cube(60, 60, 60, 10.0, 10.0, 10.0, shared_from_this(), oepType, options_);
   cube.write_cube_file(fileName);
}
void OEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
void OEPotential::rotate(const Matrix& rotmat) {}
void OEPotential::translate(const Vector& trans) {}
void OEPotential::superimpose(const Matrix& refGeometry,
                              const std::vector<int>& supList,
                              const std::vector<int>& reordList) {}
void OEPotential::print_header(void) const {}


// <============== Electrostatic Energy (demo class) ==============> //

ElectrostaticEnergyOEPotential::ElectrostaticEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

ElectrostaticEnergyOEPotential::~ElectrostaticEnergyOEPotential() {}
void ElectrostaticEnergyOEPotential::common_init() 
{
   oepTypes_.push_back("V");
}

void ElectrostaticEnergyOEPotential::compute(const std::string& oepType) {}
void ElectrostaticEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) 
{
  double val = 0.0;
  if (oepType == "V" || oepType == "TOTAL") {
      // ===> Nuclear contribution <=== //
      for (int i=0; i<wfn_->molecule()->natom(); ++i) {
           val+= (double)wfn_->molecule()->Z(i) /
                     sqrt( pow(wfn_->molecule()->x(i)-x, 2.0) 
                         + pow(wfn_->molecule()->y(i)-y, 2.0) 
                         + pow(wfn_->molecule()->z(i)-z, 2.0) );
      }

      // ===> Electronic contribution <=== //
      double v;
      potInt_->set_charge_field(x, y, z);
      OEInt_ = potInt_;
      OEInt_->compute(potMat_);
      for (int i=0; i<primary_->nbf(); ++i) {
           for (int j=0; j<=i; ++j) {
                v= (wfn_->Da()->get(i,j) + wfn_->Db()->get(i,j)) * potMat_->get(i,j);
                val += 2.0*v;
                if (i==j) val -= v;
           }
      }
      potMat_->zero();

   } else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
   }

   // Assign final value
   v = val;
}
void ElectrostaticEnergyOEPotential::print_header(void) const {}


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
   oepTypes_.push_back("S1");
   oepTypes_.push_back("S2");
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) {}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
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
    oepTypes_.push_back("ET1");
    oepTypes_.push_back("ET2");
    oepTypes_.push_back("HT1");
    oepTypes_.push_back("HT2");
    oepTypes_.push_back("CT1");
    oepTypes_.push_back("CT2");
}

void EETCouplingOEPotential::compute(const std::string& oepType) {}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
void EETCouplingOEPotential::print_header(void) const {}


//OEPotential::compute_3D(const std::string& oepType, const std::string& fileName) {
//  std::shared_ptr<CubicScalarGrid> grid = std::make_shared<CubicScalarGrid>(primary_, options_);
//  grid->write_cube_file(v_oep, fileName);
//}

} // EndNameSpace oepdev
