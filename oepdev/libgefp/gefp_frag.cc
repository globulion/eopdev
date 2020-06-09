#include <iostream>
#include <random>
#include "gefp.h"

using namespace std;

oepdev::GenEffFrag::GenEffFrag(std::string name) : 
  name_(name), frag_(nullptr),
  densityMatrixSusceptibilityGEF_(std::make_shared<oepdev::GenEffPar>("POLARIZATION")) // deprecate this!!!/TODO
//electrostaticEnergyGEF_(nullptr),
//repulsionEnergyGEF_(nullptr),
//chargeTransferEnergyGEF_(nullptr),
//EETCouplingConstantGEF_(nullptr)
{
// Deprecate all below:
//parameters["POLARIZATION"   ] = densityMatrixSusceptibilityGEF_;
//parameters["COULOMBIC   "   ] = electrostaticEnergyGEF_;
//parameters["REPULSION"      ] = repulsionEnergyGEF_;
//parameters["CHARGE_TRANSFER"] = chargeTransferEnergyGEF_;
//parameters["EET_COUPLING"   ] = EETCouplingConstantGEF_;
}
oepdev::GenEffFrag::GenEffFrag(const GenEffFrag* f) {
 name_ = f->name_;
 this->set_molecule(f->frag_);
 densityMatrixSusceptibilityGEF_ = f->densityMatrixSusceptibilityGEF_;//TODO ->copy it, not pass the pointer
}
oepdev::GenEffFrag::~GenEffFrag() {}
void oepdev::GenEffFrag::rotate(std::shared_ptr<psi::Matrix> R)
{
  for (auto const& x : this->parameters) {
       x.second->rotate(R);
  }
}
void oepdev::GenEffFrag::translate(std::shared_ptr<psi::Vector> T)
{
  for (auto const& x : this->parameters) {
       x.second->translate(T);
  }
}
void oepdev::GenEffFrag::superimpose(std::shared_ptr<psi::Matrix> targetXYZ, std::vector<int> supList)
{
  for (auto const& x : this->parameters) {
       outfile->Printf(" Superimposing now %s in Fragment %s\n", x.first.c_str(), this->name_.c_str());
       x.second->superimpose(targetXYZ, supList);
  }
}
void oepdev::GenEffFrag::superimpose(psi::SharedMolecule mol, std::vector<int> supList)
{
 psi::SharedMatrix xyz = std::make_shared<psi::Matrix>(mol->geometry());
 this->superimpose(xyz, supList);
}
void oepdev::GenEffFrag::superimpose(void)
{
 psi::SharedMatrix xyz = std::make_shared<psi::Matrix>(this->frag_->geometry());
 this->superimpose(xyz, {});
}
// 
double oepdev::GenEffFrag::energy(std::string theory, std::shared_ptr<GenEffFrag> other) {
 double e_tot = this->compute_pairwise_energy(theory, other);
 return e_tot;
}

double oepdev::GenEffFrag::compute_energy(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments, bool manybody)
{
  double e_tot = 0.0;
  if (manybody) {e_tot += oepdev::GenEffFrag::compute_many_body_energy(theory, fragments);}
  else {

     const int n = fragments.size();                                            
     for (int i=0; i<n; ++i) {
     for (int j=0; j<i; ++j) {
          e_tot += fragments[i]->compute_pairwise_energy(theory, fragments[j]);
     }}
  }
  return e_tot;
}
double oepdev::GenEffFrag::compute_many_body_energy(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments)
{
  double e;

  if (theory == "EFP2:IND" || theory == "OEPa-EFP2:IND" || theory == "OEPb-EFP2:IND") {
      e = 0.0;
      //TODO
  } else {
     throw psi::PSIEXCEPTION("Wrong manybody theory chosen.");
  }
  return e;
}
//double oepdev::GenEffFrag::compute_pairwise_total_energy(std::vector<std::shared_ptr<GenEffFrag>> fragments,
//          double (oepdev::GenEffFrag::*pairwise_energy_computer)(void )) {
//  const int n = fragments.size();
//  double e_tot = 0.0;
//  for (int i=0; i<n; ++i) {
//  for (int j=0; j<i; ++j) {
//       e_tot += fragments[i]->
//       e_tot += pairwise_energy_computer(fragments[i], fragments[j]);
//  }
//  }
//  return e_tot;
//}
double oepdev::GenEffFrag::compute_pairwise_energy(std::string theory, std::shared_ptr<GenEffFrag> other) {
  double e;
 if (theory == "EFP2:COUL" || theory == "OEPa-EFP2:COUL" || theory == "OEPb-EFP2:COUL") 
     {e = this->compute_pairwise_energy_efp2_coul(other);}
 else if (theory == "EFP2:IND" || theory == "OEPa-EFP2:IND" || theory == "OEPb-EFP2:IND") 
     {e = this->compute_pairwise_energy_efp2_ind(other);}
 else if (theory == "EFP2:EXREP" || theory == "OEPa-EFP2:EXREP") 
     {e = this->compute_pairwise_energy_efp2_exrep(other);}
 else if (theory == "EFP2:CT")
     {e = this->compute_pairwise_energy_efp2_ct(other);}
 else if (theory == "EFP2:DISP" || theory == "OEPa-EFP2:DISP" || theory == "OEPb-EFP2:DISP")
     {e = this->compute_pairwise_energy_efp2_disp(other);}
 else if (theory == "OEPb-EFP2:EXREP") 
     {e = this->compute_pairwise_energy_oep_efp2_exrep(other);}
 else if (theory == "OEPa-EFP2:CT" || theory == "OEPb-EFP2:CT") 
     {e = this->compute_pairwise_energy_oep_efp2_ct(other);}
 else {
     throw psi::PSIEXCEPTION("Wrong pairwise theory chosen.");
 }
 return e;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_coul(std::shared_ptr<GenEffFrag> other) {
 oepdev::MultipoleConvergence::ConvergenceLevel clevel = oepdev::DMTPole::determine_dmtp_convergence_level("DMTP_CONVER");
 oepdev::SharedDMTPole camm_1 = this->parameters["efp2"]->dmtp("camm");
 oepdev::SharedDMTPole camm_2 =other->parameters["efp2"]->dmtp("camm");
 double e = camm_1->energy(camm_2)->level(clevel)->get(0,0);
 return e;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 outfile->Printf(" Computing EFP2-ExRep interaction energy\n");
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ind(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_disp(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
