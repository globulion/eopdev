#include "gefp.h"
#include <iostream>
#include <fstream>

using namespace std;

//-- EFP2_GEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::EFP2_GEFactory::EFP2_GEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::GenEffParFactory(wfn, opt) { }
oepdev::EFP2_GEFactory::~EFP2_GEFactory() { }

std::shared_ptr<oepdev::GenEffPar> oepdev::EFP2_GEFactory::compute()
{
   // TODO
   std::shared_ptr<oepdev::DMTPole> camm = this->compute_dmtp();
   //this->compute_lmoc();
   std::shared_ptr<oepdev::CPHF> cphf = this->compute_cphf();
   this->assemble_parameters();
}

std::shared_ptr<oepdev::DMTPole> oepdev::EFP2_GEFactory::compute_dmtp() {
  psi::outfile->Printf(" @EFP2_GEFactory: Calculating CAMM...");
  std::shared_ptr<oepdev::DMTPole> camm = oepdev::DMTPole::build("CAMM", wfn_);
  camm->compute();
  return camm;
}

void oepdev::EFP2_GEFactory::compute_lmoc() {}
std::shared_ptr<oepdev::CPHF> oepdev::EFP2_GEFactory::compute_cphf() {
  std::shared_ptr<oepdev::CPHF> cphf = std::make_shared<oepdev::CPHF>(wfn_, options_);
  cphf->compute();
  return cphf;
}
void oepdev::EFP2_GEFactory::assemble_parameters() {}
