#include "gefp.h"
#include <iostream>
#include <fstream>

using namespace std;

//-- OEP_EFP2_GEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::OEP_EFP2_GEFactory::OEP_EFP2_GEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::EFP2_GEFactory(wfn, opt), 
 auxiliary_(nullptr), intermediate_(nullptr), oep_rep_(nullptr), oep_ct_(nullptr) {}
oepdev::OEP_EFP2_GEFactory::OEP_EFP2_GEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt,
                                               psi::SharedBasisSet auxiliary, psi::SharedBasisSet intermediate) :
 oepdev::OEP_EFP2_GEFactory(wfn, opt) { 
  auxiliary_ = auxiliary;
  intermediate_ = intermediate;
}
oepdev::OEP_EFP2_GEFactory::~OEP_EFP2_GEFactory() {}
std::shared_ptr<oepdev::GenEffPar> oepdev::OEP_EFP2_GEFactory::compute()
{
   // The same as in EFP2
   std::shared_ptr<oepdev::DMTPole> camm = this->compute_dmtp();
   std::shared_ptr<oepdev::CPHF> cphf = this->compute_cphf();

   dmtp_ = camm;
   cphfSolver_ = cphf;

   // Compute additional OEP parameters
   oep_rep_ = oepdev::OEPotential::build(      "REPULSION ENERGY", wfn_, auxiliary_, intermediate_, wfn_->options());
   oep_ct_  = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY", wfn_, auxiliary_, intermediate_, wfn_->options());

   oep_rep_->compute("Murrell-etal.S1");
   oep_rep_->compute("Otto-Ladik.S2.CAMM.a");
   oep_rep_->compute("Otto-Ladik.S2.CAMM.A");
   oep_ct_->compute("Otto-Ladik.V1.GDF");
   oep_ct_->localize();
   oep_ct_->compute("Otto-Ladik.V3.CAMM-nj");

   // Assemble all
   this->assemble_efp2_parameters();
   this->assemble_oep_efp2_parameters();

   return this->EFP2Parameters_;
}
void oepdev::OEP_EFP2_GEFactory::assemble_oep_efp2_parameters() {
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling OEP-REP data...\n");
  this->EFP2Parameters_->set_oep("rep", oep_rep_);
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling OEP-CT data...\n");
  this->EFP2Parameters_->set_oep("ct", oep_ct_);
}
