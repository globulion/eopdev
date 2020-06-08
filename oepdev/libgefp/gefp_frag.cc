#include <iostream>
#include <random>
#include "gefp.h"

using namespace std;

oepdev::GenEffFrag::GenEffFrag(std::string name) : 
  name_(name), frag_(nullptr),
  densityMatrixSusceptibilityGEF_(std::make_shared<oepdev::GenEffPar>("POLARIZATION"))
//electrostaticEnergyGEF_(nullptr),
//repulsionEnergyGEF_(nullptr),
//chargeTransferEnergyGEF_(nullptr),
//EETCouplingConstantGEF_(nullptr)
{
  parameters["POLARIZATION"   ] = densityMatrixSusceptibilityGEF_;
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
double oepdev::GenEffFrag::energy(std::string theory, std::shared_ptr<GenEffFrag> other) {
 double e_tot;
 if (theory == "EFP2") {
     double e_coul = this->compute_energy_efp2_coul(other);
     double e_exrep= this->compute_energy_efp2_exrep(other);
     double e_ind  = this->compute_energy_efp2_ind(other);
     double e_ct   = this->compute_energy_efp2_ct(other);
     double e_disp = this->compute_energy_efp2_disp(other);
     e_tot = e_coul + e_exrep + e_ind + e_ct + e_disp;
 } else if (theory == "OEP_EFP2-a") {
     double e_coul = this->compute_energy_efp2_coul(other);
     double e_exrep= this->compute_energy_oep_efp2_exrep(other);
     double e_ind  = this->compute_energy_efp2_ind(other);
     double e_ct   = this->compute_energy_oep_efp2_ct(other);
     double e_disp = this->compute_energy_efp2_disp(other);
     e_tot = e_coul + e_exrep + e_ind + e_ct + e_disp;
 } else if (theory == "OEP_EFP2-b") {
     double e_coul = this->compute_energy_efp2_coul(other);
     double e_exrep= this->compute_energy_efp2_exrep(other);
     double e_ind  = this->compute_energy_efp2_ind(other);
     double e_ct   = this->compute_energy_oep_efp2_ct(other);
     double e_disp = this->compute_energy_efp2_disp(other);
     e_tot = e_coul + e_exrep + e_ind + e_ct + e_disp;
 } else {
     throw psi::PSIEXCEPTION("Wrong theory chosen.");
 }
 return e_tot;
}
double oepdev::GenEffFrag::compute_energy_efp2_coul(std::shared_ptr<GenEffFrag> other) {
 //TODO
 oepdev::MultipoleConvergence::ConvergenceLevel clevel = oepdev::DMTPole::determine_dmtp_convergence_level("DMTP_CONVER");
 oepdev::SharedDMTPole camm_1 = this->parameters["efp2"]->dmtp("camm");
 oepdev::SharedDMTPole camm_2 =other->parameters["efp2"]->dmtp("camm");
 outfile->Printf(" Computing CAMM interaction energy\n");
 double e = camm_1->energy(camm_2)->level(clevel)->get(0,0);
 return e;
}
double oepdev::GenEffFrag::compute_energy_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_energy_efp2_ind(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_energy_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_energy_efp2_disp(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_energy_oep_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_energy_oep_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
