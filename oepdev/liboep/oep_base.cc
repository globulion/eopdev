#include "oep.h"
#include "../lib3d/space3d.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

OEPType::OEPType(std::string name, bool is_density_fitted, int n, SharedMatrix matrix, SharedDMTPole dmtp, SharedCISData cis_data) {
   this->name = name;
   this->is_density_fitted = is_density_fitted;
   this->n = n;
   this->matrix = matrix;
   this->dmtp = dmtp;
   this->cis_data = cis_data;
}
OEPType::OEPType(const OEPType* f) {
   this->name = f->name;
   this->is_density_fitted = f->is_density_fitted;
   this->n = f->n;
   if (f->matrix) this->matrix = f->matrix->clone();
   if (f->dmtp) this->dmtp = f->dmtp->clone();
   if (f->cis_data) this->cis_data = std::make_shared<CISData>(f->cis_data.get());
}

OEPotential::OEPotential(SharedWavefunction wfn, Options& options) 
   : wfn_(wfn),
     options_(options)
{
    common_init();
}
OEPotential::OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options) 
   : wfn_(wfn),
     auxiliary_(auxiliary),
     intermediate_(intermediate),
     options_(options)
{
    common_init();
}
OEPotential::OEPotential(const OEPotential* f) {
  // Do shallow-copy on those:
  options_ = f->options_;
  wfn_ = f->wfn_;
  primary_ = f->primary_;
  auxiliary_ = f->auxiliary_;
  intermediate_ = f->intermediate_;
  localizer_ = f->localizer_;
  // Do deep-copy on those:
  name_ = f->name_;
  use_localized_orbitals = f->use_localized_orbitals;
  if (f->cOcc_) cOcc_ = f->cOcc_->clone();
  if (f->cVir_) cVir_ = f->cVir_->clone();
  if (f->potMat_) {potMat_= f->potMat_->clone();} else {potMat_=nullptr;}
  if (f->lOcc_) {lOcc_ = f->lOcc_->clone();} else {lOcc_=nullptr;}
  if (f->T_) {T_ = f->T_->clone();} else {T_=nullptr;}
  // Set to null
  intsFactory_ = nullptr;
  OEInt_= nullptr;
  potInt_ = nullptr;
  // Copy the rest (lmoc_ and oepTypes_)
  this->copy_from(f);
}
void OEPotential::copy_from(const OEPotential* f) {
  // copy lmoc_
  this->lmoc_ = {nullptr, nullptr, nullptr};
  if (f->lmoc_[0]) {
      this->lmoc_[0] = std::make_shared<psi::Vector>(*(f->lmoc_[0])); 
      this->lmoc_[1] = std::make_shared<psi::Vector>(*(f->lmoc_[1])); 
      this->lmoc_[2] = std::make_shared<psi::Vector>(*(f->lmoc_[2])); 
  }
  // copy oepTypes_
  this->oepTypes_.clear(); 
  for (auto const& x : f->oepTypes_) {
     std::string key = x.first;
     OEPType o = OEPType(x.second);
     this->oepTypes_[key] = o;
  }
}
void OEPotential::common_init(void) 
{
   name_        = "default";
   primary_     = wfn_->basisset();
   intsFactory_ = std::make_shared<psi::IntegralFactory>(primary_);
   potMat_      = std::make_shared<psi::Matrix>("Potential Integrals", primary_->nbf(), primary_->nbf());
   potInt_      = std::make_shared<oepdev::PotentialInt>(intsFactory_->spherical_transform(), primary_, primary_, 0);
   cOcc_        = wfn_->Ca_subset("AO","OCC");
   cVir_        = wfn_->Ca_subset("AO","VIR");
   lOcc_        = nullptr; //std::make_shared<psi::Matrix>();
   localizer_   = nullptr; //psi::Localizer::build(options_.get_str("SOLVER_CT_LOCALIZER"), primary_, cOcc_, options_);
   T_           = nullptr;
   lmoc_        = {nullptr, nullptr, nullptr};
   use_localized_orbitals = false;
}
std::vector<psi::SharedVector> OEPotential::mo_centroids(psi::SharedMatrix C){
   std::vector<std::shared_ptr<psi::Matrix>> Rao;                                                                  
   std::vector<std::shared_ptr<psi::Vector>> Rmo;
   for (int z=0; z<3; ++z) Rao.push_back(std::make_shared<psi::Matrix>("Rao", primary_->nbf(), primary_->nbf()));
   for (int z=0; z<3; ++z) Rmo.push_back(std::make_shared<psi::Vector>("Rmo", C->ncol()));
   psi::IntegralFactory fact(primary_);
   std::shared_ptr<psi::OneBodyAOInt> dipInt(fact.ao_dipole());
   dipInt->compute(Rao);
   for (int z=0; z<3; ++z) {
        Rao[z]->scale(-1.0);
        std::shared_ptr<psi::Matrix> CRC = psi::Matrix::triplet(C, Rao[z], C, true, false, false);
        for (int a=0; a<C->ncol(); ++a) {
             Rmo[z]->set(a, CRC->get(a,a));
        }
        Rao[z].reset(); 
   }
   return Rmo;
}
void OEPotential::localize(void) 
{
   std::string o_loc = options_.get_str("SOLVER_CT_LOCALIZER");
   localizer_ = psi::Localizer::build(o_loc, primary_, cOcc_, options_);
   localizer_->localize();
   // LMO's
   lOcc_ = localizer_->L()->clone();
   // LMO centroids
   lmoc_ = this->mo_centroids(lOcc_);
   // Transformation
   T_ = localizer_->U()->clone();
}
std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY"  )  oep = std::make_shared< ElectrostaticEnergyOEPotential>(wfn, options);
   else if (category == "REPULSION ENERGY"      )  oep = std::make_shared<     RepulsionEnergyOEPotential>(wfn, options);
   else if (category == "CHARGE TRANSFER ENERGY")  oep = std::make_shared<ChargeTransferEnergyOEPotential>(wfn, options);
   else if (category == "EET COUPLING CONSTANT" )  oep = std::make_shared<         EETCouplingOEPotential>(wfn, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}
std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, 
                                                                           SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY"  )  throw PSIEXCEPTION("OEPDEV: OEPotential build. DF-based OEP not available for Electrostatic Energy!");
   else if (category == "REPULSION ENERGY"      )  oep = std::make_shared< RepulsionEnergyOEPotential>(wfn, auxiliary, intermediate, options);
   else if (category == "CHARGE TRANSFER ENERGY")  oep = std::make_shared< ChargeTransferEnergyOEPotential>(wfn, auxiliary, intermediate, options);
   else if (category == "EET COUPLING CONSTANT" )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, auxiliary, intermediate, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}
OEPotential::~OEPotential() {}
void OEPotential::compute(const std::string& oepType) {}
void OEPotential::compute(void) { for ( auto const& oepType : oepTypes_ ) this->compute(oepType.second.name); }
void OEPotential::write_cube(const std::string& oepType, const std::string& fileName) 
{
   OEPotential3D<OEPotential> oeps3d(oepTypes_[oepType].n, 60, 60, 60, 10.0, 10.0, 10.0, shared_from_this(), oepType, options_);
   oeps3d.write_cube_file(fileName);
}
void OEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
std::shared_ptr<OEPotential3D<OEPotential>> OEPotential::make_oeps3d(const std::string& oepType)
{
      std::shared_ptr<OEPotential3D<OEPotential>> oeps3d(new OEPotential3D<OEPotential>(oepTypes_[oepType].n, 
                                        options_.get_int   ("ESP_NPOINTS_PER_ATOM") * wfn_->molecule()->natom(), 
                                        options_.get_double("ESP_PAD_SPHERE"      ),
                                        shared_from_this(), oepType));
      return oeps3d;
}
void OEPotential::rotate(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) {}
void OEPotential::rotate_basic(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) {
  // Rotate orbitals
  psi::SharedMatrix Ri = R_prim->clone(); Ri->invert(); Ri->transpose_this();
  psi::SharedMatrix new_cOcc = psi::Matrix::doublet(Ri, cOcc_, true, false);
  psi::SharedMatrix new_cVir = psi::Matrix::doublet(Ri, cVir_, true, false);
  this->cOcc_->copy(new_cOcc);
  this->cVir_->copy(new_cVir);
  if (lOcc_) {
  psi::SharedMatrix new_lOcc = psi::Matrix::doublet(Ri, lOcc_, true, false);
  this->lOcc_->copy(new_lOcc);
  }

  // Rotate lmoc_
  //TODO - now not necessary!
  // 

  // Rotate potMat_
  if (potMat_) {
  psi::SharedMatrix new_potMat = psi::Matrix::triplet(R_prim, potMat_, R_prim, true, false, false);
  this->potMat_->copy(new_potMat);
  }
}
void OEPotential::translate(psi::SharedVector t) {}
void OEPotential::translate_basic(psi::SharedVector t) {
  // Translate lmoc_
  //TODO - now not necessary!
  // 
}
void OEPotential::superimpose(const Matrix& refGeometry,
                              const std::vector<int>& supList,
                              const std::vector<int>& reordList) {}
void OEPotential::print_header(void) const {}
void OEPotential::print(void) const {
  psi::outfile->Printf("\n\n ===> OEP %10s data <===\n\n", name_.c_str());

  for ( auto const& oepType : oepTypes_ ) {
   psi::outfile->Printf("\n --> Type: %s <--\n\n", oepType.second.name.c_str());

   // ==> Print GDF parameters <== //
   if (oepType.second.is_density_fitted){
   psi::outfile->Printf(  "     Density-Fitted\n\n");
   oepType.second.matrix->print();
   } else {
   psi::outfile->Printf(  "     Multipole-Based\n\n");
   oepType.second.matrix->print();
   if (oepType.second.dmtp) oepType.second.dmtp->print();
   } 
  }
}
