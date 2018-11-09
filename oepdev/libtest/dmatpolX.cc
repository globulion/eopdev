#include <iostream>
#include "test.h"
#include "../libutil/cphf.h"
#include "../libutil/util.h"
#include "../libgefp/gefp.h"

using namespace std;


double oepdev::test::Test::test_dmatPolX(void)
{
  /* In field perturbation nuclear dipole - field interaction was not added to the total energy */
  //std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn_, options_);
  //solver->compute();

  psi::timer_on (" Test: Computation of Dmat Susc-X");
  std::shared_ptr<oepdev::GenEffParFactory> factory = oepdev::GenEffParFactory::build("POLARIZATION", wfn_, options_);
  std::shared_ptr<oepdev::GenEffPar> par = factory->compute();
  psi::timer_off(" Test: Computation of Dmat Susc-X");

  // Accumulate errors
  double r_sum = 0.0;

  // Compute SCF in the external field
  double Fx = options_.get_double("DMATPOL_TEST_FIELD_X");
  double Fy = options_.get_double("DMATPOL_TEST_FIELD_Y");
  double Fz = options_.get_double("DMATPOL_TEST_FIELD_Z");

  std::shared_ptr<psi::SuperFunctional> func = oepdev::create_superfunctional("HF", options_);
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
  scf->set_perturbation(Fx, Fy, Fz);
  scf->compute_energy();

  // Compute dipole integrals
  std::vector<std::shared_ptr<psi::Matrix>> Mao;
  int nbf = wfn_->basisset()->nbf();
  for (int z=0;z<3;++z) Mao.push_back(std::make_shared<psi::Matrix>("",nbf,nbf));
  psi::IntegralFactory ints(wfn_->basisset());
  std::shared_ptr<psi::OneBodyAOInt> dipInt(ints.ao_dipole());
  dipInt->compute(Mao);
  std::shared_ptr<psi::Matrix> MF = std::make_shared<psi::Matrix>("",nbf,nbf);
  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            MF->set(i,j, Fx * Mao[0]->get(i,j) + Fy * Mao[1]->get(i,j) + Fz * Mao[2]->get(i,j));
       }
  }
  MF->scale(-1.0);
  std::shared_ptr<psi::Matrix> S = wfn_->S();

  std::shared_ptr<psi::Matrix> D0 = wfn_->Da();
  std::shared_ptr<psi::Matrix> D  = scf->Da();
  std::shared_ptr<psi::Matrix> dD = D->clone();
  dD->subtract(D0);
  dD->set_name("Difference Density Matrix: Exact");
  
  std::shared_ptr<psi::Matrix> dD_m = dD->clone();
  dD_m->copy(par->compute_density_matrix(Fx, Fy, Fz));
  dD_m->set_name("Difference Density Matrix: Model");

  dD->print();
  dD_m->print();

  cout << " Traces with Hcore" << endl;
  cout << dD  ->vector_dot(wfn_->H()) << endl;
  cout << dD_m->vector_dot(wfn_->H()) << endl;
  cout << " Traces with Da" << endl;
  cout << dD  ->vector_dot(wfn_->Da()) << endl;
  cout << dD_m->vector_dot(wfn_->Da()) << endl;
  cout << " Traces with Fa" << endl;
  cout << dD  ->vector_dot(wfn_->Fa()) << endl;
  cout << dD_m->vector_dot(wfn_->Fa()) << endl;
  cout << " Traces with S" << endl;
  cout << dD  ->vector_dot(wfn_->S ()) << endl;
  cout << dD_m->vector_dot(wfn_->S ()) << endl;
  cout << " Traces with MF" << endl;
  cout << dD  ->vector_dot(MF)        << endl;
  cout << dD_m->vector_dot(MF)        << endl;

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  dD->subtract(dD_m);
  std::cout << " Difference RMS: " << dD->rms() << std::endl;

  // Return
  return r_sum;
}

