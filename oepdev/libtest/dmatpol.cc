#include <iostream>
#include "test.h"
#include "../libutil/cphf.h"
#include "../libutil/util.h"
#include "../libgefp/gefp.h"


using namespace std;

double oepdev::test::Test::test_dmatPol(void)
{
  /* In field perturbation nuclear dipole - field interaction was not added to the total energy */
  //std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn_, options_);
  //solver->compute();

  psi::timer_on (" Test: Computation of Dmat Susc");
  std::shared_ptr<oepdev::GenEffParFactory> factory = std::make_shared<oepdev::AbInitioPolarGEFactory>(wfn_, options_);
  std::shared_ptr<oepdev::GenEffPar> par = factory->compute();
  psi::timer_off(" Test: Computation of Dmat Susc");

  par->dipole_polarizability(0,0)->print();
  std::shared_ptr<oepdev::CPHF> solver = factory->cphf_solver();

  // Accumulate errors
  double r_sum = 0.0;

  // Compute SCF in the external field
  double Fx = 0.002;
  double Fy = 0.04;
  double Fz =-0.01;
  double F[3] = {Fx, Fy, Fz};

  std::shared_ptr<psi::Wavefunction> scf_base = solve_scf(wfn_->molecule(),wfn_->basisset(),wfn_->get_basisset("DF_BASIS_SCF"),
                      oepdev::create_superfunctional("HF", options_),
                      options_, wfn_->psio());
  std::shared_ptr<psi::scf::RHF> scf = std::make_shared<psi::scf::RHF>(scf_base, 
                      oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());

  scf->initialize();

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

  double** H = scf->H()->pointer();
  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            H[i][j] += MF->get(i,j);
       }
  }
  scf->iterations();
  scf->finalize_E();

  std::shared_ptr<psi::Matrix> D0 = wfn_->Da();
  std::shared_ptr<psi::Matrix> D  = scf->Da();
  std::shared_ptr<psi::Matrix> dD = D->clone();
  dD->subtract(D0);
  dD->set_name("Difference Density Matrix: Exact");
  
  std::shared_ptr<psi::Matrix> dD_m = dD->clone();
  dD_m->zero();
  dD_m->set_name("Difference Density Matrix: Model");

  std::vector<std::shared_ptr<psi::Vector>> ind_dipoles;
  std::shared_ptr<psi::Vector> ind_dipole = std::make_shared<psi::Vector>("INDUCED DIPOLE REFERENCE",3);

  for (int o=0; o<solver->nocc(); ++o) {
       ind_dipoles.push_back(std::make_shared<psi::Vector>(oepdev::string_sprintf("IndDip -%d-", o+1),3));
       ind_dipoles[o]->set(0, solver->polarizability(o)->get(0,0) * Fx
                             +solver->polarizability(o)->get(0,1) * Fy
                             +solver->polarizability(o)->get(0,2) * Fz);
       ind_dipoles[o]->set(1, solver->polarizability(o)->get(1,0) * Fx
                             +solver->polarizability(o)->get(1,1) * Fy
                             +solver->polarizability(o)->get(1,2) * Fz);
       ind_dipoles[o]->set(2, solver->polarizability(o)->get(2,0) * Fx
                             +solver->polarizability(o)->get(2,1) * Fy
                             +solver->polarizability(o)->get(2,2) * Fz);

       ind_dipole->add(ind_dipoles[o]);

       ind_dipoles[o]->print();

       std::shared_ptr<psi::Matrix> temp = std::make_shared<psi::Matrix>("",nbf,nbf);
       temp->copy(par->dipole_polarizability(o, 0)); 
       temp->scale(Fx);
       dD_m->add(temp);
       temp->zero();

       temp->copy(par->dipole_polarizability(o, 1)); 
       temp->scale(Fy);
       dD_m->add(temp);
       temp->zero();

       temp->copy(par->dipole_polarizability(o, 2)); 
       temp->scale(Fz);
       dD_m->add(temp);
       temp->zero();
  }
  ind_dipole->print();
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

  //double e_coul_0 = 4.0 * wfn_->Da()->vector_dot(MF);
  ////for (int i=0; i<wfn_->molecule()->natom(); ++i) {
  ////     e_coul_0 -= (double)wfn_->molecule()->Z(i) * (
  ////                         wfn_->molecule()->x(i) * Fx +
  ////                         wfn_->molecule()->y(i) * Fy +
  ////                         wfn_->molecule()->z(i) * Fz
  ////                 );
  ////}
  //cout << e_coul_0 << endl;

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}

