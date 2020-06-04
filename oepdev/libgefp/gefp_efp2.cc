#include "gefp.h"
#include <iostream>
#include <fstream>

using namespace std;

//-- EFP2_GEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::EFP2_GEFactory::EFP2_GEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::GenEffParFactory(wfn, opt),
 EFP2Parameters_(std::make_shared<oepdev::GenEffPar>("EFP2"))
{ 
  // Assert if orbital localization for CPHF solver is switched on
  if (options_.get_bool("CPHF_LOCALIZE") == false) 
      throw psi::PSIEXCEPTION(" Error! EFP2 requires orbital localization. Switch CPHF_LOCALIZE to true.");
}
oepdev::EFP2_GEFactory::~EFP2_GEFactory() { }

std::shared_ptr<oepdev::GenEffPar> oepdev::EFP2_GEFactory::compute()
{
   std::shared_ptr<oepdev::DMTPole> camm = this->compute_dmtp();
   std::shared_ptr<oepdev::CPHF> cphf = this->compute_cphf();

   dmtp_ = camm;
   cphfSolver_ = cphf;

   this->assemble_parameters();

   return this->EFP2Parameters_;
}

std::shared_ptr<oepdev::DMTPole> oepdev::EFP2_GEFactory::compute_dmtp() {
  psi::outfile->Printf(" @EFP2_GEFactory: Calculating CAMM...\n");
  std::shared_ptr<oepdev::DMTPole> camm = oepdev::DMTPole::build("CAMM", wfn_);
  camm->compute();
  psi::outfile->Printf(" @EFP2_GEFactory: CAMM Done.\n");
  return camm;
}

void oepdev::EFP2_GEFactory::compute_lmoc() {}
std::shared_ptr<oepdev::CPHF> oepdev::EFP2_GEFactory::compute_cphf() {
  psi::outfile->Printf(" @EFP2_GEFactory: Solving CPHF Equations...\n");
  std::shared_ptr<oepdev::CPHF> cphf = std::make_shared<oepdev::CPHF>(wfn_, options_);
  cphf->compute();
  psi::outfile->Printf(" @EFP2_GEFactory: CPHF Done.\n");
  return cphf;
}
void oepdev::EFP2_GEFactory::assemble_parameters() {
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Geometry data...\n");
  this->assemble_geometry_data();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling DMTP data...\n");
  this->assemble_dmtp_data();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling LMO centroids...\n");
  this->assemble_lmo_centroids();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Fock Matrix...\n");
  this->assemble_fock_matrix();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Distributed Polarizabilities...\n");
  this->assemble_distributed_polarizabilities();
}

void oepdev::EFP2_GEFactory::assemble_geometry_data() {
  psi::SharedMatrix geom = std::make_shared<psi::Matrix>(wfn_->molecule()->geometry());
  this->EFP2Parameters_->set_matrix("POS", geom);
}
void oepdev::EFP2_GEFactory::assemble_dmtp_data() {
  this->EFP2Parameters_->set_dmtp("CAMM", dmtp_);
}
void oepdev::EFP2_GEFactory::assemble_lmo_centroids() {
  psi::SharedMatrix C = cphfSolver_->localizer()->L();

  psi::SharedMatrix lmoc = std::make_shared<psi::Matrix>("LMO Centroids", cphfSolver_->nocc(), 3);
  for (int i=0; i<cphfSolver_->nocc(); ++i) {
       lmoc->set_row(0, i, cphfSolver_->lmo_centroid(i));
  }
  this->EFP2Parameters_->set_matrix("LMOC", lmoc);
}

void oepdev::EFP2_GEFactory::assemble_fock_matrix() {
  psi::SharedMatrix C = cphfSolver_->localizer()->L();
  psi::SharedMatrix Fa_ao = wfn_->Fa();
  psi::SharedMatrix Fa_mo = psi::Matrix::triplet(C, Fa_ao, C, true, false, false);

  this->EFP2Parameters_->set_matrix("FOCK", Fa_mo);

  psi::SharedMatrix Ca_occ_canonical = cphfSolver_->Cocc();
  psi::SharedMatrix Ca_vir_canonical = cphfSolver_->Cvir();

  this->EFP2Parameters_->set_matrix("CMOO", Ca_occ_canonical);
  this->EFP2Parameters_->set_matrix("CMOV", Ca_vir_canonical);
}

void oepdev::EFP2_GEFactory::assemble_distributed_polarizabilities() {
  std::vector<psi::SharedMatrix> dpol;
  for (int i=0; i<cphfSolver_->nocc(); ++i) {
       dpol.push_back(cphfSolver_->polarizability(i));
   }

  this->EFP2Parameters_->set_dpol("DPOL-0", dpol);
}
