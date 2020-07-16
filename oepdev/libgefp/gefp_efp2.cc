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
   dmtpSolver_   = this->compute_dmtp();
   cphfSolver_   = this->compute_cphf(); if (options_.get_bool("EFP2_WITH_VVO"))
   quamboSolver_ = this->compute_quambo(); 

   this->assemble_efp2_parameters();

   return this->EFP2Parameters_;
}

oepdev::SharedDMTPole oepdev::EFP2_GEFactory::compute_dmtp() {
  psi::outfile->Printf(" @EFP2_GEFactory: Calculating CAMM...\n");
  oepdev::SharedDMTPole camm = oepdev::DMTPole::build("CAMM", wfn_);
  camm->compute();
  psi::outfile->Printf(" @EFP2_GEFactory: CAMM Done.\n");
  return camm;
}

void oepdev::EFP2_GEFactory::compute_lmoc() {}

oepdev::SharedCPHF oepdev::EFP2_GEFactory::compute_cphf() {
  psi::outfile->Printf(" @EFP2_GEFactory: Solving CPHF Equations...\n");
  oepdev::SharedCPHF cphf = std::make_shared<oepdev::CPHF>(wfn_, options_);
  cphf->compute();
  psi::outfile->Printf(" @EFP2_GEFactory: CPHF Done.\n");
  return cphf;
}

oepdev::SharedQUAMBO oepdev::EFP2_GEFactory::compute_quambo() {
  psi::outfile->Printf(" @EFP2_GEFactory: Computing QUAMBO and VVO...\n");
  oepdev::SharedQUAMBO solver = std::make_shared<QUAMBO>(this->wfn_, true);
  solver->compute();
  psi::outfile->Printf(" @EFP2_GEFactory: QUAMBO and VVO Done.\n");
  return solver;
}

void oepdev::EFP2_GEFactory::assemble_efp2_parameters() {
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Geometry Data...\n");
  this->assemble_geometry_data();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling DMTP Data...\n");
  this->assemble_dmtp_data();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling LMO Dentroids...\n");
  this->assemble_lmo_centroids();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Fock Matrix...\n");
  this->assemble_fock_matrix();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling HF Canonical Orbitals...\n");
  this->assemble_canonical_orbitals();
  psi::outfile->Printf(" @EFP2_GEFactory: Assembling Distributed Polarizabilities...\n");
  this->assemble_distributed_polarizabilities();
}

void oepdev::EFP2_GEFactory::assemble_geometry_data() {
  psi::SharedMatrix geom = std::make_shared<psi::Matrix>(wfn_->molecule()->geometry());
  this->EFP2Parameters_->set_matrix("pos", geom);

  psi::SharedVector Z = std::make_shared<psi::Vector>("Atomic numbers", wfn_->molecule()->natom());
  for (int i=0; i<wfn_->molecule()->natom(); ++i) {
       Z->set(0, (double)wfn_->molecule()->Z(i));
  }
  this->EFP2Parameters_->set_vector("atno", Z);
}
void oepdev::EFP2_GEFactory::assemble_dmtp_data() {
  this->EFP2Parameters_->set_dmtp("camm", dmtpSolver_);
}
void oepdev::EFP2_GEFactory::assemble_lmo_centroids() {
  psi::SharedMatrix C = cphfSolver_->localizer()->L()->clone();
  this->EFP2Parameters_->set_matrix("lmoo", C);

  psi::SharedMatrix lmoc = std::make_shared<psi::Matrix>("LMO Centroids", cphfSolver_->nocc(), 3);
  for (int i=0; i<cphfSolver_->nocc(); ++i) {
       lmoc->set_row(0, i, cphfSolver_->lmo_centroid(i));
  }
  this->EFP2Parameters_->set_matrix("lmoc", lmoc);
}
void oepdev::EFP2_GEFactory::assemble_fock_matrix() {
  psi::SharedMatrix C = cphfSolver_->localizer()->L();
  psi::SharedMatrix Fa_ao = wfn_->Fa()->clone();
  psi::SharedMatrix Fa_mo = psi::Matrix::triplet(C, Fa_ao, C, true, false, false);

  this->EFP2Parameters_->set_matrix("fock_lmo", Fa_mo);
  this->EFP2Parameters_->set_matrix("fock_ao", Fa_ao);
}
void oepdev::EFP2_GEFactory::assemble_canonical_orbitals() {
  psi::SharedMatrix Ca_occ_canonical;
  psi::SharedMatrix Ca_vir_canonical;
  psi::SharedVector eps_a_occ_canonical;
  psi::SharedVector eps_a_vir_canonical;
  if (!this->options_.get_bool("EFP2_WITH_VVO")) {
       // Original HF canonical orbitals
       Ca_occ_canonical = cphfSolver_->Cocc();      
       Ca_vir_canonical = cphfSolver_->Cvir();
       eps_a_occ_canonical = cphfSolver_->epsocc();
       eps_a_vir_canonical = cphfSolver_->epsvir();
  } else {
       // Original HF canonical Occupied Orbitals + VVOs from QUAMBO
       Ca_occ_canonical = quamboSolver_->Ca_subset("AO","OCC");
       Ca_vir_canonical = quamboSolver_->Ca_subset("AO","VIR");
       eps_a_occ_canonical = quamboSolver_->epsilon_a_subset("MO","OCC");
       eps_a_vir_canonical = quamboSolver_->epsilon_a_subset("MO","VIR");
  }

  this->EFP2Parameters_->set_matrix("cmoo", Ca_occ_canonical);
  this->EFP2Parameters_->set_matrix("cmov", Ca_vir_canonical);

  this->EFP2Parameters_->set_vector("epso", eps_a_occ_canonical);
  this->EFP2Parameters_->set_vector("epsv", eps_a_vir_canonical);
}
void oepdev::EFP2_GEFactory::assemble_distributed_polarizabilities() {
  std::vector<psi::SharedMatrix> dpol;
  for (int i=0; i<cphfSolver_->nocc(); ++i) {
       dpol.push_back(cphfSolver_->polarizability(i));
   }

  this->EFP2Parameters_->set_dpol("0", dpol);
}
