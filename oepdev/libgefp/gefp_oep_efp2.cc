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
   // The same as in EFP2 but omit QUAMBOs since they are computed by OEP objects
   dmtpSolver_ = this->compute_dmtp();
   cphfSolver_ = this->compute_cphf();

   // Compute additional OEP parameters
   oep_rep_ = oepdev::OEPotential::build(      "REPULSION ENERGY", wfn_, auxiliary_, intermediate_, wfn_->options());
   oep_ct_  = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY", wfn_, auxiliary_, intermediate_, wfn_->options());

   oep_rep_->use_quambo_orbitals = false;
   oep_rep_->use_localized_orbitals = true;
   oep_rep_->compute("Murrell-etal.S1");
   oep_rep_->compute("Otto-Ladik.S2.CAMM.a");
   oep_rep_->compute("Otto-Ladik.S2.CAMM.A");

   oep_ct_->use_quambo_orbitals = options_.get_bool("EFP2_WITH_VVO");
   oep_ct_->use_localized_orbitals = false;
   oep_ct_->compute("Otto-Ladik.V1.GDF");
   oep_ct_->set_occupied_canonical_orbitals(oep_rep_); // ensures occupied canonical orbitals in rep and ct are exactly the same
   oep_ct_->set_localized_orbitals(oep_rep_);          // ensures occupied localized orbitals in rep and ct are exactly the same
   oep_ct_->use_localized_orbitals = true;
   oep_ct_->compute("Otto-Ladik.V3.CAMM-nj");

   // Assemble all parameters
   this->assemble_efp2_parameters();
   this->assemble_oep_efp2_parameters();

   return this->EFP2Parameters_;
}
void oepdev::OEP_EFP2_GEFactory::assemble_oep_efp2_parameters() {
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling OEP-LMO data...\n");
  this->assemble_oep_lmo_centroids();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling OEP-REP data...\n");
  this->EFP2Parameters_->set_oep("rep", oep_rep_);
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling OEP-CT data...\n");
  this->EFP2Parameters_->set_oep("ct", oep_ct_);
}
void oepdev::OEP_EFP2_GEFactory::assemble_oep_lmo_centroids() {
  /* Note: 
   * For OEP calculations we use different set of LMOs than in EFP2.
   * this->oep_rep_ is always computed during parameter calculations
   * with use_localized_orbitals set to True. It is used as a source
   * of localized orbitals
   */
  psi::SharedMatrix C = oep_rep_->lOcc()->clone();
  this->EFP2Parameters_->set_matrix("lmoo-oep", C);

  psi::SharedMatrix T = oep_rep_->T()->clone();
  this->EFP2Parameters_->set_matrix("cmoo2lmoo", T);

  psi::SharedMatrix lmoc = std::make_shared<psi::Matrix>("LMO-OEP Centroids", cphfSolver_->nocc(), 3);
  std::vector<psi::SharedVector> lmoc_as_vecs_3_n = oep_rep_->mo_centroids(C);
  for (int z=0; z<3; ++z) {
       lmoc->set_column(0, z, lmoc_as_vecs_3_n[z]);
  }
  this->EFP2Parameters_->set_matrix("lmoc-oep", lmoc);
}
void oepdev::OEP_EFP2_GEFactory::assemble_canonical_orbitals() {
  /* 
     Set original HF or QUAMBO-based orbitals depending on the oep_ct_ settings
   */
  psi::SharedMatrix Ca_occ_canonical = oep_ct_->cOcc()->clone();
  psi::SharedMatrix Ca_vir_canonical = oep_ct_->cVir()->clone();
  psi::SharedVector eps_a_occ_canonical = std::make_shared<psi::Vector>(*(oep_ct_->epsOcc()));
  psi::SharedVector eps_a_vir_canonical = std::make_shared<psi::Vector>(*(oep_ct_->epsVir()));

  this->EFP2Parameters_->set_matrix("cmoo", Ca_occ_canonical);
  this->EFP2Parameters_->set_matrix("cmov", Ca_vir_canonical);

  this->EFP2Parameters_->set_vector("epso", eps_a_occ_canonical);
  this->EFP2Parameters_->set_vector("epsv", eps_a_vir_canonical);
}

