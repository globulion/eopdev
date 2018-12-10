#include <iostream>
#include "test.h"
#include "../libutil/util.h"
#include "../libutil/scf_perturb.h"

using namespace std;


double oepdev::test::Test::test_scf_perturb()
{
  /* Test includes the contribution from nuclear charges */
  double result = 0.0; 
  double const AtoBohr = 1.889725989;
  const double energy_field_ref =-74.9710884012145584; // From Psi4 ref/scf_perturb/field.inp
  const double energy_charge_ref=-74.9698931724711173; // From Psi4 ref/scf_perturb/charges.inp
  const double Fx = 0.024, Rx1 = 1.4000*AtoBohr, Rx2 =-0.9000*AtoBohr;
  const double Fy =-0.019, Ry1 = 0.0939*AtoBohr, Ry2 =-0.3939*AtoBohr;
  const double Fz = 0.009, Rz1 = 3.0030*AtoBohr, Rz2 =-1.0480*AtoBohr;
  const double q1 = 0.9100, q2 =-0.6200;

  std::shared_ptr<psi::SuperFunctional> func = oepdev::create_superfunctional("HF", options_);

  // Solve SCF in external electric field
  psi::outfile->Printf("\n ==> Computing for uniform electric field <==\n\n");
  std::shared_ptr<oepdev::RHFPerturbed> scf_field = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
  scf_field->set_perturbation(Fx, Fy, Fz);
  scf_field->compute_energy();
  const double energy_field = scf_field->reference_energy();

  // Solve SCF with one point charge
  psi::outfile->Printf("\n ==> Computing for 2 point charges <==\n\n");
  std::shared_ptr<oepdev::RHFPerturbed> scf_charge = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
  scf_charge->set_perturbation(Rx1, Ry1, Rz1, q1);
  scf_charge->set_perturbation(Rx2, Ry2, Rz2, q2);
  scf_charge->compute_energy();
  const double energy_charge = scf_charge->reference_energy();

  // Accumulate errors
  result += pow(energy_field - energy_field_ref , 2.0);
  result += pow(energy_charge- energy_charge_ref, 2.0);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
