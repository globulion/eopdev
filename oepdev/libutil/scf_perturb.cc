#include "psi4/libmints/integral.h"
#include "../libpsi/potential.h"
#include "scf_perturb.h"

namespace oepdev {

RHFPerturbed::RHFPerturbed(std::shared_ptr<psi::Wavefunction> ref_wfn, std::shared_ptr<psi::SuperFunctional> functional)
 : psi::scf::RHF(ref_wfn, functional), 
   perturbField_(std::make_shared<psi::Vector>("Perturbing Electric Field",3)),
   perturbCharges_(std::make_shared<PerturbCharges>())
{
}
RHFPerturbed::RHFPerturbed(std::shared_ptr<psi::Wavefunction> ref_wfn, std::shared_ptr<psi::SuperFunctional> functional,
    psi::Options& options, std::shared_ptr<psi::PSIO> psio) 
 : psi::scf::RHF(ref_wfn, functional, options, psio),
   perturbField_(std::make_shared<psi::Vector>("Perturbing Electric Field",3)),
   perturbCharges_(std::make_shared<PerturbCharges>())
{
}
RHFPerturbed::~RHFPerturbed()
{
}
double RHFPerturbed::compute_energy()
{
  initialize();
  perturb_Hcore();
  iterations();
  return finalize_E();
}
void RHFPerturbed::set_perturbation(std::shared_ptr<psi::Vector> field)
{
  perturbField_->copy(*field);
}
void RHFPerturbed::set_perturbation(const double& fx, const double& fy, const double& fz)
{
  std::shared_ptr<psi::Vector> f = std::make_shared<psi::Vector>("", 3);
  f->set(0, fx);
  f->set(1, fy);
  f->set(2, fz);
  set_perturbation(f);
}
void RHFPerturbed::set_perturbation(std::shared_ptr<psi::Vector> position, const double& charge)
{
  perturbCharges_->charges  .push_back(charge);
  perturbCharges_->positions.push_back(std::shared_ptr<psi::Vector>(position->clone()));
}
void RHFPerturbed::set_perturbation(const double& rx, const double& ry, const double& rz, const double& charge)
{
  std::shared_ptr<psi::Vector> r = std::make_shared<psi::Vector>("", 3);
  r->set(0, rx);
  r->set(1, ry);
  r->set(2, rz);
  set_perturbation(r, charge);
}
void RHFPerturbed::perturb_Hcore()
{
  // Initialize the perturbation Hcore matrix
  std::shared_ptr<psi::Matrix> Hadd = std::make_shared<psi::Matrix>("Hcore perturbation", basisset_->nbf(), basisset_->nbf());

  // Build-up perturbations due to electric field
  std::vector<std::shared_ptr<psi::Matrix>> Dip;
  for (int z=0; z<3; ++z) Dip.push_back(std::make_shared<psi::Matrix>("DipInt",basisset_->nbf(), basisset_->nbf()));
  std::shared_ptr<psi::OneBodyAOInt> dipInt(integral_->ao_dipole());
  dipInt->compute(Dip);
  for (int z=0; z<3; ++z) Hadd->axpy(-perturbField_->get(z), Dip[z]);

  // Build-up perturbations due to point charges
  std::shared_ptr<oepdev::PotentialInt> potInt = std::make_shared<oepdev::PotentialInt>(integral_->spherical_transform(), basisset_, basisset_, 0);
  std::shared_ptr<psi::OneBodyAOInt> oneInt;
  std::shared_ptr<psi::Matrix> V = std::make_shared<psi::Matrix>("V", basisset_->nbf(), basisset_->nbf());

  for (int n=0; n<perturbCharges_->charges.size(); ++n){
       potInt->set_charge_field(perturbCharges_->positions[n]->get(0),
                                perturbCharges_->positions[n]->get(1),
                                perturbCharges_->positions[n]->get(2));
       oneInt = potInt;
       oneInt->compute(V);
       Hadd->axpy(perturbCharges_->charges[n], V);
       V->zero();
  }

  // Add perturbation to Hcore matrix
  H_->add(Hadd);

  // Re-calculate guess
  psi::timer_on("HF: Guess-Perturbed");
  guess();
  psi::timer_off("HF: Guess-Perturbed");

  // Add the contribution from nuclei
  for (int z=0; z<3; ++z) nuclearrep_ -= perturbField_->get(z) * molecule_->nuclear_dipole().get(z);
  for (int I=0; I<molecule_->natom(); ++I) {
       for (int n=0; n<perturbCharges_->charges.size(); ++n){
            double x = molecule_->x(I) - perturbCharges_->positions[n]->get(0);
            double y = molecule_->y(I) - perturbCharges_->positions[n]->get(1);
            double z = molecule_->z(I) - perturbCharges_->positions[n]->get(2);
            double r_In = sqrt(x*x + y*y + z*z);
            nuclearrep_ += perturbCharges_->charges[n] * (double)(molecule_->Z(I)) / r_In;
       }
  }
}

} // EndNameSpace oepdev
