#include <iostream>
#include "test.h"
#include "psi4/libmints/matrix.h"
#include "../libutil/cphf.h"
#include "../libutil/util.h"
#include "../libutil/unitary_optimizer.h"
#include "../libgefp/gefp.h"
#include "../lib3d/dmtp.h"

using namespace std;

oepdev::test::Test::Test(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options) :
 wfn_(wfn), options_(options)
{
}
oepdev::test::Test::~Test() 
{
}
double oepdev::test::Test::run(void)
{
  double result;
  if      (options_.get_str("OEPDEV_TEST_NAME")=="BASIC"  ) result = test_basic  ();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CPHF"   ) result = test_cphf   ();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL") result = test_dmatPol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL_X") result = test_dmatPolX();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_1_1") result = test_eri_1_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_3_1") result = test_eri_3_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER") result = test_unitaryOptimizer();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER_4_2") result = test_unitaryOptimizer_4_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="SCF_PERTURB") result = test_scf_perturb();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CAMM") result = test_camm();
  else throw psi::PSIEXCEPTION("Incorrect test name specified!");
  return result;
}
double oepdev::test::Test::test_basic(void)
{
  // Read reference data
  const double ref_energy = -74.965002974633;
  const double ref_D[49]  = {     2.10706290,
      -0.45100674,     0.08310605,     0.06647334,    -0.02146005,    -0.02536451,
      -0.02539307,    -0.45100674,     1.98400112,    -0.46570102,    -0.37261261,
       0.12028659,    -0.04161382,    -0.04174240,     0.08310605,    -0.46570102,
       1.01664358,     0.23894684,    -0.01948050,     0.71844075,     0.02334713,
       0.06647334,    -0.37261261,     0.23894684,     0.99697621,     0.25666998,
      -0.10823398,     0.70131230,    -0.02146005,     0.12028659,    -0.01948050,
       0.25666998,     1.93244525,    -0.00345464,    -0.18802698,    -0.02536451,
      -0.04161382,     0.71844075,    -0.10823398,    -0.00345464,     0.61777967,
      -0.18193560,    -0.02539307,    -0.04174240,     0.02334713,     0.70131230,
      -0.18802698,    -0.18193560,     0.61752184};
  
  const double ref_S[49]  = {     1.00000002,
       0.23670394,     0.00000000,    -0.00000000,    -0.00000000,     0.04988317,
       0.04995228,     0.23670394,     1.00000003,     0.00000000,     0.00000000,
      -0.00000000,     0.45329197,     0.45366441,     0.00000000,     0.00000000,
       1.00000002,     0.00000000,    -0.00000000,     0.37260686,    -0.01487029,
      -0.00000000,     0.00000000,     0.00000000,     1.00000002,     0.00000000,
      -0.08261482,     0.36889805,    -0.00000000,    -0.00000000,    -0.00000000,
       0.00000000,     1.00000002,     0.00526845,    -0.09768319,     0.04988317,
       0.45329197,     0.37260686,    -0.08261482,     0.00526845,     0.99999999,
       0.23348072,     0.04995228,     0.45366441,    -0.01487029,     0.36889805,
      -0.09768319,     0.23348072,     0.99999999};
  
  const double ref_F[49]  = {   -20.23721180,
      -5.16293310,    -0.02178888,    -0.01744041,     0.00562972,    -1.10550392,
      -1.10702115,    -5.16293310,    -2.43580923,    -0.08638505,    -0.06927099,
       0.02235340,    -0.96807077,    -0.96899727,    -0.02178888,    -0.08638505,
      -0.31709770,    -0.01787751,     0.00151548,    -0.52223226,    -0.02449398,
      -0.01744041,    -0.06927099,    -0.01787751,    -0.31581677,    -0.01880895,
       0.07126681,    -0.50871736,     0.00562972,     0.02235340,     0.00151548,
      -0.01880895,    -0.38443448,     0.00448625,     0.13673079,    -1.10550392,
      -0.96807077,    -0.52223226,     0.07126681,     0.00448625,    -0.53571599,
      -0.35502764,    -1.10702115,    -0.96899727,    -0.02449398,    -0.50871736,
       0.13673079,    -0.35502764,    -0.53644448};
  
  const double ref_H[49]  = {   -32.68421101,
      -7.60411655,    -0.01359058,    -0.01088508,     0.00351329,    -1.61526323,
      -1.61749754,    -7.60411655,    -9.30126043,    -0.16169184,    -0.12947179,
       0.04179037,    -3.53236969,    -3.53577701,    -0.01359058,    -0.16169184,
      -7.54061793,     0.02880600,    -0.00271303,    -2.45623567,     0.04537431,
      -0.01088508,    -0.12947179,     0.02880600,    -7.54370603,     0.02897248,
       0.49288256,    -2.42256523,     0.00351329,     0.04179037,    -0.00271303,
       0.02897248,    -7.43777932,    -0.02094154,     0.64383980,    -1.61526323,
      -3.53236969,    -2.45623567,     0.49288256,    -0.02094154,    -4.94164328,
      -1.46364788,    -1.61749754,    -3.53577701,     0.04537431,    -2.42256523,
       0.64383980,    -1.46364788,    -4.94376112};

  // Error accumulator
  double r_sum = 0.0;

  int nbf = wfn_->basisset()->nbf();
  // Matrices
  std::shared_ptr<psi::Matrix> C = wfn_->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> D = wfn_->Da();
  std::shared_ptr<psi::Matrix> F = wfn_->Fa();
  D->scale(2.0); // Total density matrix

  // 1-Electron Integrals
  psi::IntegralFactory fact(wfn_->basisset());

  std::shared_ptr<psi::Matrix> S = std::make_shared<psi::Matrix>("S", nbf, nbf);
  std::shared_ptr<psi::Matrix> T = std::make_shared<psi::Matrix>("T", nbf, nbf);
  std::shared_ptr<psi::Matrix> V = std::make_shared<psi::Matrix>("V", nbf, nbf);

  std::shared_ptr<psi::OneBodyAOInt> sI(fact.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> tI(fact.ao_kinetic());
  std::shared_ptr<psi::OneBodyAOInt> vI(fact.ao_potential());

  sI->compute(S);
  tI->compute(T);
  vI->compute(V);
  T->add(V); // Construct the core Hamiltonian

  // Accumulate the error
  r_sum = sqrt(pow(wfn_->reference_energy() - ref_energy, 2.0));
  double** pD = D->pointer();
  double** pF = F->pointer();
  double** pS = S->pointer();
  double** pH = T->pointer();
  const double* pref_D = ref_D;
  const double* pref_F = ref_F;
  const double* pref_S = ref_S;
  const double* pref_H = ref_H;

  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            r_sum += sqrt(pow(pD[i][j] - *pref_D, 2.0));
            r_sum += sqrt(pow(pF[i][j] - *pref_F, 2.0));
            r_sum += sqrt(pow(pS[i][j] - *pref_S, 2.0));
            r_sum += sqrt(pow(pH[i][j] - *pref_H, 2.0));
            pref_D++;
            pref_F++;
            pref_S++;
            pref_H++;
       }
  }

  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test SCF energy= " << wfn_->reference_energy() << std::endl;
  std::cout << " Ref  SCF energy= " <<             ref_energy   << std::endl;
 
  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}
double oepdev::test::Test::test_cphf(void)
{
  // Reference data for MeNH2 at RHF/6-311++G** (6D)
  const double pol_ref[9] = { 22.375420,  -0.545125,   0.630563,
                              -0.545355,  20.939757,  -0.263902,  
                               0.630892,  -0.263922,  20.901605};

  std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn_, options_);
  solver->compute();
  for (int i=0; i<solver->nocc(); i++) {
       solver->lmo_centroid(i)->print();
       solver->polarizability(i)->print();
  }
  solver->polarizability()->print();

  // Accumulate errors
  double r_sum = 0.0;
  double** p = solver->polarizability()->pointer();
  const double*  r = pol_ref;
  for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
            r_sum += sqrt(pow(p[i][j] - *r, 2.0));
            r++;
       }
  }

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}
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

  std::shared_ptr<psi::Wavefunction> scf_base = solve_scf(wfn_->molecule(),wfn_->basisset(),
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

double oepdev::test::Test::test_eri_2_2(void)
{
  // Sizing
  int nbf = wfn_->basisset()->nbf();
  size_t size = nbf*nbf*nbf*nbf;

  // Residual error buffer
  double* errors = new double[size];
  memset(errors, 0, sizeof(double)*size);

  // OEPDev implementation of ERI's
  psi::timer_on(" Test: Computation of OepDev ERI_2_2");
  oepdev::IntegralFactory fact_oepdev(wfn_->basisset());
  std::shared_ptr<psi::TwoBodyAOInt> eri_2_2(fact_oepdev.eri_2_2());

  oepdev::SharedShellsIterator shellIter = oepdev::ShellCombinationsIterator::build(fact_oepdev, "ALL");
  int i, j, k, l; 
  double integral;
  const double * buffer_2_2 = eri_2_2->buffer();

  int icount = 0;
  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(eri_2_2);
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            int i = intsIter->i();int j = intsIter->j();int k = intsIter->k();int l = intsIter->l();
            double integral = buffer_2_2[intsIter->index()];
            errors[icount] = integral;
            ++icount;
       }
  }
  psi::timer_off(" Test: Computation of OepDev ERI_2_2");

  // Psi4 implementation of ERI's
  psi::timer_on(" Test: Computation of PSI4   ERI_2_2");
  oepdev::IntegralFactory fact_psi(wfn_->basisset());
  std::shared_ptr<psi::TwoBodyAOInt> eri(fact_psi.eri());
  const double* buffer = eri->buffer();

  icount = 0;
  shellIter = oepdev::ShellCombinationsIterator::build(fact_psi, "ALL");
  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(eri);
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            int i = intsIter->i();int j = intsIter->j();int k = intsIter->k();int l = intsIter->l();
            double integral = buffer[intsIter->index()];
            errors[icount]-=integral;
            ++icount;
       }
  }
  psi::timer_off(" Test: Computation of PSI4   ERI_2_2");

  // Compute the residual error sum
  double r_sum = 0.0;
  for (int i=0; i<size; ++i) {
       r_sum += errors[i]*errors[i];
  }

  // Clean up
  delete[] errors;

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;
  // Return
  return r_sum;
}
double oepdev::test::Test::test_eri_1_1(void)
{
  // Sizing
  int nbf = wfn_->basisset()->nbf();
  size_t size = nbf*nbf;

  // Store ERI's and square errors
  psi::Matrix eri_1_1_oepdev("ERI_1_1: OepDev", nbf, nbf);
  psi::Matrix eri_1_1_psi4  ("ERI_1_1: Psi4  ", nbf, nbf);
  psi::Matrix error         ("Squared Errors ", nbf, nbf);

  // OEPDev implementation of ERI's
  psi::timer_on(" Test: Computation of OepDev ERI_1_1");
  oepdev::IntegralFactory fact_oepdev(wfn_->basisset());
  std::shared_ptr<oepdev::TwoBodyAOInt> eri_1_1(fact_oepdev.eri_1_1());

  oepdev::SharedShellsIterator shellIter; 
  eri_1_1->compute(eri_1_1_oepdev);
  //oepdev::SharedShellsIterator shellIter = oepdev::ShellCombinationsIterator::build(fact_oepdev, "ALL", 2);
  //const double * buffer_1_1 = eri_1_1->buffer();

  //for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  //{
  //     shellIter->compute_shell(eri_1_1);
  //     oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
  //     for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
  //     {
  //          int i = intsIter->i();int j = intsIter->j();
  //          double integral = buffer_1_1[intsIter->index()];
  //          eri_1_1_oepdev.set(i,j,integral);
  //     }
  //}
  psi::timer_off(" Test: Computation of OepDev ERI_1_1");

  // Psi4 implementation of ERI's
  psi::timer_on(" Test: Computation of Psi4   ERI_1_1");
  oepdev::IntegralFactory fact_psi4(wfn_->basisset(),psi::BasisSet::zero_ao_basis_set());

  std::shared_ptr<psi::TwoBodyAOInt> eri_2_2(fact_psi4.eri());

  shellIter = oepdev::ShellCombinationsIterator::build(fact_psi4, "ALL", 4);
  const double * buffer_2_2 = eri_2_2->buffer();

  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(eri_2_2);
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            int i = intsIter->i();int k = intsIter->k(); //int j = intsIter->j(); int l = intsIter->l();
            double integral = buffer_2_2[intsIter->index()]; 
            //psi::outfile->Printf("(%d%d|%d%d) = %13.6f\n", i, j, k, l, integral);
            eri_1_1_psi4.set(i,k,integral);
       }
  }
  psi::timer_off(" Test: Computation of Psi4   ERI_1_1");

  // Compute squared errors
  error.copy(eri_1_1_oepdev.clone());
  error.subtract(&eri_1_1_psi4); 
  double result = error.sum_of_squares();

  // Print ERI's and errors
  eri_1_1_oepdev.print();
  eri_1_1_psi4.print();
  error.print();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  // Finish
  return result;
}
double oepdev::test::Test::test_eri_3_1(void)
{
  // Sizing
  int nbf = wfn_->basisset()->nbf();
  size_t size = nbf*nbf;

  // Store ERI's and square errors
  psi::Matrix eri_3_1_oepdev("ERI_3_1: OepDev", nbf, nbf);
  psi::Matrix eri_3_1_psi4  ("ERI_3_1: Psi4  ", nbf, nbf);
  psi::Matrix error         ("Squared Errors ", nbf, nbf);

  // OEPDev implementation of ERI's
  psi::timer_on(" Test: Computation of OepDev ERI_3_1");
  oepdev::IntegralFactory fact_oepdev(wfn_->basisset(),psi::BasisSet::zero_ao_basis_set(),
                                      psi::BasisSet::zero_ao_basis_set(),wfn_->basisset());
  std::shared_ptr<oepdev::TwoBodyAOInt> eri_3_1(fact_oepdev.eri_3_1());

  oepdev::SharedShellsIterator shellIter = oepdev::ShellCombinationsIterator::build(fact_oepdev, "ALL", 4);
  const double * buffer_3_1 = eri_3_1->buffer();

  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(eri_3_1);
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            int i = intsIter->i();int j = intsIter->j();int k = intsIter->k();int l = intsIter->l();
            double integral = buffer_3_1[intsIter->index()];
            eri_3_1_oepdev.set(i,l,integral);
       }
  }
  psi::timer_off(" Test: Computation of OepDev ERI_3_1");

  // Psi4 implementation of ERI's
  psi::timer_on(" Test: Computation of Psi4   ERI_3_1");
  oepdev::IntegralFactory fact_psi4(wfn_->basisset(), psi::BasisSet::zero_ao_basis_set());

  std::shared_ptr<psi::TwoBodyAOInt> eri_2_2(fact_psi4.eri());

  shellIter = oepdev::ShellCombinationsIterator::build(fact_psi4, "ALL", 4);
  const double * buffer_2_2 = eri_2_2->buffer();

  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(eri_2_2);
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            int i = intsIter->i();int k = intsIter->k();
            double integral = buffer_2_2[intsIter->index()]; 
            eri_3_1_psi4.set(i,k,integral);
       }
  }
  psi::timer_off(" Test: Computation of Psi4   ERI_3_1");


  // Compute squared errors
  error.copy(eri_3_1_oepdev.clone());
  error.subtract(&eri_3_1_psi4); 
  double result = error.sum_of_squares();

  // Print ERI's and errors
  eri_3_1_oepdev.print();
  eri_3_1_psi4.print();
  error.print();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  // Finish
  return result;
}
double oepdev::test::Test::test_unitaryOptimizer()
{
  double result = 0.0;
  const double P[5] = {-0.2563, -0.95325, 1.535253, 0.0553, 2.04};
  const double R[25]= {0.242412, -0.3324216, -0.61421416, -0.9023525, 0.072185, 
                       -1.4562215367, -1.52156, -0.58757, 0.305301,
                        0.45, 0.9, 1.1, -0.5, 0.01, 0.4, 0.723, -0.7512, 0.04, -0.1121, 
                        -1.533, 0.553, -0.04566, 0.05654, 0.1432, -1.04};
  const double Xmin_ref[25] = { -8.74267861e-01,  -3.17344968e-01,  -3.67352484e-01,  -6.75418301e-05,
                                 1.14274386e-05,   2.71298173e-01,  -7.21888732e-01,  -2.19297103e-02,    
                                -6.27836730e-01,   1.03024980e-01,  -3.39007266e-01,  -1.42396112e-01,
                                 9.29823205e-01,  -1.48347999e-02,   2.47344382e-03,   2.13741091e-01,
                                -5.89011051e-01,   2.07712921e-06,   7.77760943e-01,   4.96864340e-02,
                                -3.79713378e-02,   1.04681613e-01,  -3.67038294e-05,   2.62477372e-02,
                                 9.93433921e-01};
  const double Xmax_ref[25] = {-0.96081674, -0.17487789,  0.08508092, -0.03368997,  0.19461540, -0.07491281,
                               -0.40602029, -0.82986219, -0.01744834, -0.37491303, -0.18717579,  0.87781524,  
                               -0.43742598, -0.00431968,  0.05518537, -0.05791162,  0.01419595,  0.02906414,
                                0.99123309, -0.11426684, -0.18119255,  0.18385879,  0.33452004, -0.12648636,
                               -0.89747569};

  std::shared_ptr<psi::Matrix> Rm = std::make_shared<psi::Matrix>("R matrix",5,5);
  std::shared_ptr<psi::Vector> Pm = std::make_shared<psi::Vector>("P vector",5);
  for (int i=0; i< 5; ++i) {
       Pm->set(i, P[i]);
       for (int j=0; j<5; ++j) {
            Rm->set(i,j,R[5*i+j]);
       }
  }
  oepdev::UnitaryOptimizer optimizer(Rm, Pm, 1.0e-8, 100, true);

  psi::outfile->Printf(" ==> Unitary Minimization in 5 Dimensions <==\n");
  bool success_min = optimizer.minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer.X();

  psi::outfile->Printf(" ==> Unitary Maximization in 5 Dimensions <==\n");
  bool success_max = optimizer.maximize();
  std::shared_ptr<psi::Matrix> X_max = optimizer.X();

  // Accumulate errors
  double** xmin = X_min->pointer();
  double** xmax = X_max->pointer();
  const double*  xmin_ref = Xmin_ref;
  const double*  xmax_ref = Xmax_ref;
  for (int i=0; i<5; ++i) {
       for (int j=0; j<5; ++j) {
            result += sqrt(pow(std::abs(xmin[i][j]) - std::abs(*xmin_ref), 2.0));
            result += sqrt(pow(std::abs(xmax[i][j]) - std::abs(*xmax_ref), 2.0));
            xmin_ref++; 
            xmax_ref++;
       }
  }
  X_min->print();
  X_max->print();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
double oepdev::test::Test::test_unitaryOptimizer_4_2()
{
  double result = 0.0;

  double R[729] = {  -5.8298E-01,
    -2.7968E-01,  -9.9989E-01,  -6.9767E-01,  -8.5324E-01,  -9.0766E-01,  -8.1374E-01,  -6.5444E-01,
    -6.0323E-01,  -4.6118E-01,  -5.8081E-01,  -3.1478E-01,  -7.9555E-01,  -1.2188E-01,  -9.7261E-01,
    -3.2953E-01,  -5.8270E-01,  -4.4131E-01,  -8.5961E-01,  -8.0190E-01,  -1.9926E-01,  -3.1738E-02,
    -6.8658E-01,  -3.0768E-01,  -1.2361E-01,  -1.0539E-01,  -9.1496E-01,  -9.6095E-01,  -8.3017E-01,
    -1.2186E-01,  -9.0165E-01,  -5.7889E-01,  -4.2110E-02,  -4.6683E-01,  -3.0812E-01,  -6.8448E-01,
    -3.1350E-01,  -1.6537E-01,  -9.8171E-01,  -2.4986E-01,  -1.1139E-02,  -2.5183E-01,  -7.1956E-01,
    -2.1072E-01,  -8.9677E-01,  -5.5211E-01,  -9.1404E-02,  -7.0639E-01,  -7.1222E-01,  -8.6997E-01,
    -9.8063E-01,  -3.2116E-01,  -7.8837E-01,  -7.3445E-01,  -5.0843E-01,  -9.4664E-01,  -4.2588E-01,
    -8.5327E-01,  -4.1069E-01,  -3.0024E-01,  -8.9767E-01,  -5.8594E-01,  -3.0560E-01,  -5.8582E-01,
    -9.5005E-01,  -4.6410E-01,  -3.3621E-01,  -4.8511E-01,  -5.5405E-02,  -4.1344E-01,  -9.6598E-02,
    -8.6253E-01,  -8.6072E-01,  -1.9261E-01,  -6.0232E-01,  -8.3465E-01,  -7.2491E-02,  -6.5223E-01,
    -2.4919E-01,  -2.7400E-01,  -1.1669E-01,  -3.7633E-01,  -2.4906E-01,  -6.5110E-01,  -7.3007E-01,
    -1.0411E-01,  -5.7191E-01,  -3.5160E-02,  -3.3656E-01,  -3.7830E-01,  -8.8525E-01,  -5.0511E-02,
    -5.5009E-01,  -4.2161E-01,  -5.9186E-01,  -7.6297E-01,  -9.6620E-02,  -4.2632E-01,  -9.9713E-01,
    -3.8286E-01,  -6.7336E-01,  -4.7294E-01,  -1.1406E-01,  -6.4273E-01,  -9.1465E-02,  -3.7664E-01,
    -9.8418E-01,  -7.0563E-02,  -3.0910E-01,  -2.6771E-03,  -8.2766E-01,  -8.6286E-01,  -6.7405E-02,
    -3.0318E-01,  -9.3400E-01,  -2.4454E-01,  -2.4612E-01,  -7.6975E-02,  -2.8848E-01,  -8.7573E-01,
    -9.8012E-01,  -9.7379E-01,  -9.7169E-01,  -7.5379E-01,  -1.3997E-01,  -4.6117E-01,  -4.4718E-01,
    -1.5797E-01,  -8.7583E-01,  -7.2082E-01,  -4.1424E-01,  -3.0404E-02,  -4.3897E-01,  -9.8135E-01,
    -1.9937E-01,  -7.6703E-01,  -1.9289E-01,  -6.1214E-01,  -1.3646E-01,  -2.5288E-01,  -4.4376E-01,
    -8.6354E-01,  -9.4008E-01,  -8.7866E-01,  -9.5545E-01,  -8.9251E-01,  -7.7429E-01,  -2.8701E-01,
    -4.4028E-01,  -9.8744E-01,  -9.2803E-01,  -3.2724E-02,  -4.3190E-01,  -7.9671E-01,  -7.4767E-01,
    -2.5617E-01,  -8.0457E-01,  -4.1864E-01,  -2.9980E-02,  -1.5317E-01,  -7.6015E-01,  -5.0623E-01,
    -3.8004E-01,  -1.7102E-01,  -8.4321E-01,  -9.8142E-01,  -9.2998E-01,  -5.1365E-01,  -3.9367E-01,
    -4.3115E-01,  -6.8264E-01,  -1.1384E-02,  -4.2025E-01,  -6.1986E-01,  -4.4905E-01,  -2.5467E-01,
    -3.3077E-01,  -7.3508E-01,  -9.3367E-01,  -6.2992E-01,  -3.7028E-01,  -7.8983E-01,  -2.4724E-01,
    -9.3346E-01,  -7.3968E-01,  -1.9525E-01,  -8.0657E-01,  -3.6054E-01,  -4.7533E-01,  -7.5192E-02,
    -7.3670E-01,  -9.3404E-01,  -2.6493E-01,  -2.2782E-01,  -9.2184E-02,  -6.8028E-02,  -9.8605E-01,
    -7.6564E-01,  -3.8322E-01,  -5.0984E-02,  -4.9824E-02,  -4.4335E-01,  -8.4394E-02,  -3.5843E-01,
    -6.0999E-01,  -5.1401E-01,  -3.9569E-01,  -4.5045E-01,  -7.3819E-02,  -8.1267E-02,  -6.0512E-01,
    -3.6737E-02,  -8.2604E-01,  -8.7367E-01,  -8.6492E-01,  -4.9434E-01,  -9.7848E-01,  -5.2030E-02,
    -1.7288E-01,  -9.8498E-01,  -8.2380E-01,  -6.6794E-01,  -8.6900E-01,  -1.9051E-01,  -6.5526E-01,
    -5.9893E-02,  -4.1799E-01,  -1.2117E-01,  -1.5527E-01,  -9.4608E-02,  -5.4012E-01,  -4.5365E-01,
    -2.0140E-01,  -7.1428E-01,  -5.0975E-01,  -4.0089E-01,  -9.8447E-01,  -4.0652E-01,  -5.6632E-01,
    -1.9264E-01,  -6.8476E-01,  -1.0711E-01,  -4.2214E-01,  -8.1599E-01,  -2.1207E-01,  -3.8797E-01,
    -9.4609E-01,  -5.7981E-01,  -3.2093E-01,  -8.1398E-02,  -9.9960E-01,  -2.3241E-02,  -6.2342E-01,
    -2.6216E-02,  -3.9528E-01,  -1.7115E-01,  -4.2529E-01,  -3.7192E-01,  -7.1442E-01,  -4.1317E-01,
    -2.4998E-01,  -1.4169E-01,  -2.4492E-01,  -3.0194E-01,  -1.3552E-01,  -6.7732E-01,  -3.2921E-01,
    -5.4913E-01,  -6.1790E-01,  -5.8919E-01,  -5.9852E-01,  -6.8262E-01,  -3.7808E-01,  -5.6975E-01,
    -2.6198E-02,  -3.2220E-01,  -8.0143E-01,  -5.7330E-01,  -6.5665E-01,  -2.0236E-01,  -1.2000E-01,
    -9.6158E-02,  -3.3728E-01,  -7.2979E-01,  -7.4763E-01,  -1.4510E-01,  -4.7229E-01,  -1.9784E-01,
    -4.2751E-01,  -2.6686E-01,  -4.8099E-01,  -2.2912E-01,  -4.3114E-01,  -5.3429E-01,  -6.5731E-01,
    -9.3179E-01,  -6.2208E-01,  -9.2037E-01,  -1.7183E-02,  -8.1839E-01,  -1.8814E-01,  -1.2504E-01,
    -3.1159E-01,  -4.3051E-01,  -8.3903E-01,  -5.3312E-01,  -6.5483E-01,  -7.7496E-01,  -4.0749E-01,
    -6.8773E-01,  -8.3694E-02,  -9.0364E-02,  -7.4288E-01,  -8.8911E-01,  -8.0704E-01,  -5.0042E-01,
    -2.7141E-01,  -7.9181E-01,  -7.5197E-01,  -1.4833E-01,  -5.8415E-01,  -3.8331E-01,  -7.6633E-01,
    -8.9803E-01,  -4.8414E-01,  -5.2286E-01,  -8.4733E-01,  -3.7819E-01,  -4.5599E-01,  -3.4586E-01,
    -8.5545E-01,  -2.4847E-01,  -7.7795E-01,  -4.8065E-01,  -2.1470E-01,  -9.7767E-01,  -6.7564E-01,
    -1.2708E-01,  -1.5529E-01,  -4.6156E-01,  -1.3339E-01,  -5.0194E-02,  -1.7359E-01,  -1.4588E-01,
    -9.0126E-01,  -3.4870E-01,  -2.9648E-01,  -3.8976E-01,  -2.0038E-01,  -9.6543E-01,  -2.2976E-01,
    -2.6827E-01,  -7.4030E-01,  -7.4293E-01,  -3.6770E-01,  -6.5470E-01,  -2.0341E-01,  -5.5385E-01,
    -2.1725E-01,  -9.5282E-03,  -6.9975E-01,  -8.5699E-01,  -9.8692E-02,  -4.5844E-01,  -2.5260E-02,
    -3.6340E-01,  -6.0870E-03,  -4.5393E-01,  -4.7357E-01,  -8.6457E-01,  -6.4429E-01,  -9.7378E-01,
    -8.3960E-01,  -2.5436E-01,  -9.6960E-01,  -6.3346E-01,  -1.3765E-01,  -3.0732E-01,  -3.0906E-01,
    -8.1136E-01,  -5.5810E-01,  -4.1842E-01,  -1.0248E-02,  -7.9609E-01,  -7.5227E-01,  -7.3783E-01,
    -2.4983E-01,  -5.4302E-01,  -9.4307E-01,  -4.9148E-01,  -7.8804E-01,  -2.0140E-01,  -7.0267E-01,
    -9.7239E-01,  -4.0657E-01,  -1.5616E-01,  -6.1898E-01,  -2.5014E-01,  -4.8886E-01,  -4.5905E-01,
    -4.0566E-02,  -1.9604E-01,  -9.6768E-01,  -2.9061E-01,  -5.3500E-01,  -5.2451E-02,  -7.7857E-01,
    -7.3293E-01,  -9.1853E-01,  -5.7138E-01,  -8.9098E-01,  -3.6621E-01,  -1.9704E-01,  -3.0320E-01,
    -2.3379E-01,  -6.5755E-01,  -1.5415E-01,  -5.7123E-01,  -1.7599E-01,  -3.7350E-01,  -8.5658E-01,
    -9.2161E-01,  -9.8167E-01,  -9.3328E-01,  -5.4142E-01,  -8.8666E-01,  -9.7222E-01,  -2.4514E-01,
    -6.0515E-01,  -2.5306E-01,  -5.4760E-01,  -5.4991E-01,  -5.2193E-01,  -5.2600E-01,  -1.9684E-01,
    -5.9761E-01,  -9.5314E-02,  -9.6294E-01,  -2.2613E-01,  -8.7436E-01,  -3.8149E-01,  -9.8964E-01,
    -4.6137E-01,  -9.9698E-01,  -4.8806E-02,  -9.4598E-02,  -2.0403E-01,  -8.4726E-02,  -8.5444E-01,
    -8.4227E-01,  -8.1237E-01,  -3.7750E-01,  -9.4191E-02,  -1.0045E-02,  -2.8888E-01,  -2.6820E-01,
    -9.0707E-02,  -5.9913E-01,  -7.5015E-01,  -8.2657E-01,  -8.8054E-01,  -1.8739E-01,  -8.5321E-01,
    -7.3570E-01,  -1.8091E-01,  -6.8941E-01,  -1.7583E-02,  -7.3336E-01,  -4.6635E-01,  -6.8553E-01,
    -8.9227E-02,  -6.3344E-01,  -5.6641E-01,  -4.8771E-01,  -6.1114E-02,  -9.6905E-01,  -2.8312E-01,
    -1.0898E-01,  -9.7271E-01,  -4.7795E-01,  -6.7401E-01,  -1.4051E-01,  -4.4148E-01,  -3.0977E-01,
    -5.4715E-01,  -3.7169E-01,  -7.0990E-01,  -9.9065E-01,  -4.2324E-01,  -6.8856E-01,  -4.8273E-01,
    -8.3594E-02,  -5.7353E-01,  -7.5260E-01,  -6.2871E-01,  -6.8139E-02,  -6.3132E-02,  -1.5567E-01,
    -7.9793E-02,  -7.7210E-01,  -9.1252E-01,  -7.7269E-01,  -6.8562E-01,  -8.2523E-01,  -3.9291E-01,
    -5.8641E-01,  -1.8365E-01,  -8.1487E-01,  -2.9812E-01,  -7.5964E-01,  -4.2578E-01,  -6.5101E-01,
    -9.4304E-01,  -7.7119E-01,  -3.3590E-01,  -5.0275E-01,  -4.8098E-01,  -8.2528E-01,  -4.2928E-01,
    -3.2466E-03,  -1.8316E-01,  -4.0563E-01,  -2.4011E-02,  -9.8437E-02,  -4.0439E-01,  -9.6757E-01,
    -9.0642E-01,  -9.3463E-01,  -5.4827E-01,  -6.2457E-01,  -2.4650E-02,  -8.3202E-01,  -2.7212E-02,
    -2.3253E-01,  -1.7576E-01,  -3.6738E-01,  -3.3127E-01,  -5.2312E-01,  -9.8686E-01,  -6.4699E-01,
    -5.0793E-01,  -2.6991E-01,  -5.3137E-01,  -5.4260E-01,  -8.6234E-01,  -9.8911E-01,  -2.4172E-01,
    -6.8005E-01,  -1.5617E-02,  -7.7977E-01,  -6.6129E-01,  -4.7610E-01,  -2.4511E-01,  -5.3614E-01,
    -8.7518E-01,  -6.8750E-01,  -4.9548E-01,  -3.2615E-01,  -2.2985E-01,  -8.6966E-01,  -9.7708E-01,
    -4.8092E-01,  -1.9001E-01,  -9.8740E-01,  -3.2753E-01,  -3.1319E-01,  -5.5075E-01,  -8.5211E-02,
    -3.5564E-01,  -9.9476E-01,  -5.1557E-01,  -1.4068E-01,  -1.6960E-01,  -3.5085E-01,  -3.2630E-01,
    -4.2150E-01,  -7.2588E-01,  -4.3947E-01,  -3.2827E-01,  -6.4757E-01,  -1.4417E-01,  -8.0496E-01,
    -2.5268E-01,  -7.1040E-01,  -2.2620E-01,  -5.7226E-01,  -1.9230E-01,  -6.4647E-01,  -7.8631E-01,
    -2.3272E-01,  -6.9136E-01,  -2.6675E-01,  -2.5553E-01,  -7.7860E-01,  -7.8589E-01,  -8.0105E-01,
    -8.5748E-01,  -6.2292E-01,  -9.7337E-01,  -8.8908E-01,  -3.2544E-01,  -2.0022E-01,  -9.1947E-01,
    -7.6830E-01,  -7.9237E-01,  -8.2666E-02,  -2.8869E-01,  -4.4612E-01,  -6.9548E-01,  -1.6515E-01,
    -5.6469E-01,  -7.6544E-02,  -2.9395E-01,  -5.2197E-01,  -8.7379E-01,  -2.3956E-02,  -8.4017E-01,
    -7.9740E-01,  -5.6882E-01,  -5.9580E-01,  -8.5325E-01,  -2.7068E-01,  -8.1125E-01,  -3.5610E-01,
    -2.4569E-01,  -7.8927E-01,  -3.9905E-01,  -2.5107E-01,  -3.6178E-01,  -4.0287E-01,  -7.0452E-01,
    -2.6839E-01,  -5.4692E-02,  -5.7444E-01,  -2.1782E-01,  -9.4386E-01,  -1.6473E-01,  -8.0775E-01,
    -6.0490E-01,  -6.9992E-01,  -9.1990E-01,  -9.5369E-02,  -6.2985E-01,  -4.6930E-01,  -5.0588E-01,
    -8.6784E-01,  -7.9355E-01,  -9.2381E-01,  -4.9208E-01,  -7.3845E-01,  -6.4294E-01,  -8.9193E-01,
    -2.1245E-01,  -8.9342E-01,  -1.4291E-02,  -8.2284E-01,  -4.2759E-01,  -9.5515E-01,  -2.1288E-01,
    -8.1039E-01,  -4.7210E-01,  -2.5992E-01,  -8.5007E-01,  -4.4891E-01,  -7.8338E-01,  -2.4080E-01,
    -2.7708E-01,  -8.2345E-01,  -1.3803E-01,  -9.8022E-01,  -1.3976E-01,  -4.4110E-01,  -5.9678E-01,
    -2.4125E-01,  -2.8307E-01,  -1.2674E-02,  -7.2191E-01,  -9.9621E-01,  -6.6097E-02,  -1.4210E-01,
    -2.7115E-01,  -4.8331E-01,  -2.9304E-01,  -2.1947E-01,  -6.2512E-01,  -2.2968E-01,  -2.4938E-01,
    -3.8679E-01,  -5.9813E-01,  -3.0269E-01,  -9.9689E-01,  -2.2510E-01,  -1.0358E-01,  -7.6068E-01,
    -8.7923E-01,  -7.7972E-01,  -6.9790E-01,  -1.1697E-01,  -4.5683E-01,  -7.1329E-01,  -8.6165E-01,
    -7.0986E-01,  -3.8613E-01,  -6.7586E-01,  -5.4264E-01,  -5.5588E-01,  -1.7186E-01,  -5.7365E-01,
    -6.5430E-01,  -3.2503E-01,  -7.7852E-01,  -5.3275E-01,  -6.8523E-01,  -3.7314E-01,  -1.2264E-01,
    -5.5231E-01,  -2.1554E-01,  -5.4303E-01,  -3.4377E-01,  -8.6816E-01,  -5.6702E-01,  -9.0688E-02,
    -3.9452E-01,  -2.3323E-01,  -4.9530E-01,  -5.0194E-01,  -1.5710E-01,  -9.3219E-01,  -4.2673E-01 };
  
  double P[27] = {  -5.7237E-02,
    -4.8214E-01,  -8.0553E-01,  -1.5206E-01,  -7.4836E-01,  -2.9927E-01,  -4.5974E-01,  -5.1164E-02,
    -3.7566E-01,  -1.6202E-01,  -9.9207E-01,  -1.0660E-02,  -9.2229E-01,  -6.7787E-01,  -5.3848E-02,
    -9.9106E-01,  -1.7727E-01,  -1.3879E-01,  -5.6017E-01,  -7.4425E-01,  -1.9731E-01,  -5.2214E-01,
    -8.6566E-01,  -7.2151E-02,  -1.0403E-01,  -5.0845E-01,  -1.4330E-01 };
  
  const double Xmin_ref[9] = {   4.4346E-01,
    -6.1498E-01,   6.5203E-01,  -7.4667E-01,  -6.5590E-01,  -1.1081E-01,   4.9581E-01,  -4.3771E-01,
    -7.5006E-01 };
  
  const double Xmax_ref[9] = {   6.9716E-01,
     2.5450E-01,  -6.7022E-01,   3.4562E-01,  -9.3837E-01,   3.1831E-03,  -6.2811E-01,  -2.3386E-01,
    -7.4215E-01 };

  oepdev::UnitaryOptimizer_4_2 optimizer(R, P, 3, 1.0e-8, 100, true);

  psi::outfile->Printf(" ==> Unitary Minimization 4_2 in 3 Dimensions <==\n");
  bool success_min = optimizer.minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer.X();

  psi::outfile->Printf(" ==> Unitary Maximization 4_2 in 3 Dimensions <==\n");
  bool success_max = optimizer.maximize();
  std::shared_ptr<psi::Matrix> X_max = optimizer.X();

  // Accumulate errors
  double** xmin = X_min->pointer();
  double** xmax = X_max->pointer();
  const double*  xmin_ref = Xmin_ref;
  const double*  xmax_ref = Xmax_ref;
  for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
            result += sqrt(pow(std::abs(xmin[i][j]) - std::abs(*xmin_ref), 2.0));
            result += sqrt(pow(std::abs(xmax[i][j]) - std::abs(*xmax_ref), 2.0));
            xmin_ref++; 
            xmax_ref++;
       }
  }
  X_min->print();
  X_max->print();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
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
  result += sqrt(pow(energy_field - energy_field_ref , 2.0));
  result += sqrt(pow(energy_charge- energy_charge_ref, 2.0));

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
double oepdev::test::Test::test_camm(void) {
  // This test is for H2O at HF/6-31* molecule
  double result = 0.0;

  // Reference CAMM values
                          /* */
  const double c_ref[3] = {-3.366373E-01,  
                            1.682078E-01,  
                            1.684295E-01};
                          /* X              Y              Z */
  const double m_ref[9] = { 1.751974E-01,  1.403032E-01, -4.528554E-02,
                            1.331952E-02,  2.019633E-02, -5.983574E-03,
                            2.305333E-02,  8.621490E-03, -3.335370E-03};
                          /* XX             XY             XZ             YY             YZ             ZZ */
  const double q_ref[18]= {-3.461364E+00, -4.200223E-02,  5.218767E-03, -3.455108E+00, -3.522603E-02, -3.585162E+00,
                           -2.246864E-01, -8.106905E-02,  1.017712E-02, -4.150175E-01, -1.340822E-02, -4.758673E-01,
                           -4.343364E-01, -4.277004E-02,  9.403507E-03, -2.190033E-01, -6.733931E-02, -4.612635E-01};
                          /* XXX            XXY            XXZ            XYY            XYZ            XZZ            YYY */
                          /* YYZ            YZZ            ZZZ */
  const double o_ref[30]= {-4.118915E-01,  1.322068E-01, -1.378869E-02,  2.361463E-02, -2.364673E-03,  3.349728E-02,  -4.066315E-01,
                            1.213214E-01, -8.290932E-03, -1.620214E-02,
                           -7.924737E-01,  2.163772E-01, -2.529971E-02, -1.338179E-01,  1.787315E-02, -4.012411E-02,   1.979517E-01,
                           -2.819298E-02,  3.948758E-02, -2.403612E-02, 
                            1.433252E-01, -1.065108E-01,  2.453663E-02,  1.064634E-01, -1.776971E-02,  3.281251E-02,  -8.072710E-01,
                            1.886268E-01, -8.783347E-02,  4.116774E-02};
  // Compute CAMM
  std::shared_ptr<DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn_);
  dmtp->compute();
  dmtp->charges(0)      ->print();
  dmtp->dipoles(0)      ->print();
  dmtp->quadrupoles(0)  ->print();
  dmtp->octupoles(0)    ->print();
  dmtp->hexadecapoles(0)->print();
  std::shared_ptr<psi::Matrix> c = dmtp->charges      (0);
  std::shared_ptr<psi::Matrix> m = dmtp->dipoles      (0);
  std::shared_ptr<psi::Matrix> q = dmtp->quadrupoles  (0);
  std::shared_ptr<psi::Matrix> o = dmtp->octupoles    (0);
  std::shared_ptr<psi::Matrix> h = dmtp->hexadecapoles(0);
 
  // Accumulate errors
  for (int n=0; n<dmtp->n_sites(); ++n) {
       result += sqrt(pow(c->get(n, 0) - c_ref[n] , 2.0));
       //
       result += sqrt(pow(m->get(n, 0) - m_ref[3*n+0] , 2.0));
       result += sqrt(pow(m->get(n, 1) - m_ref[3*n+1] , 2.0));
       result += sqrt(pow(m->get(n, 2) - m_ref[3*n+2] , 2.0));
       //
       result += sqrt(pow(q->get(n, 0) - q_ref[6*n+0] , 2.0));
       result += sqrt(pow(q->get(n, 1) - q_ref[6*n+1] , 2.0));
       result += sqrt(pow(q->get(n, 2) - q_ref[6*n+2] , 2.0));
       result += sqrt(pow(q->get(n, 3) - q_ref[6*n+3] , 2.0));
       result += sqrt(pow(q->get(n, 4) - q_ref[6*n+4] , 2.0));
       result += sqrt(pow(q->get(n, 5) - q_ref[6*n+5] , 2.0));
       //
       result += sqrt(pow(o->get(n, 0) - o_ref[10*n+0] , 2.0));
       result += sqrt(pow(o->get(n, 1) - o_ref[10*n+1] , 2.0));
       result += sqrt(pow(o->get(n, 2) - o_ref[10*n+2] , 2.0));
       result += sqrt(pow(o->get(n, 3) - o_ref[10*n+3] , 2.0));
       result += sqrt(pow(o->get(n, 4) - o_ref[10*n+4] , 2.0));
       result += sqrt(pow(o->get(n, 5) - o_ref[10*n+5] , 2.0));
       result += sqrt(pow(o->get(n, 6) - o_ref[10*n+6] , 2.0));
       result += sqrt(pow(o->get(n, 7) - o_ref[10*n+7] , 2.0));
       result += sqrt(pow(o->get(n, 8) - o_ref[10*n+8] , 2.0));
       result += sqrt(pow(o->get(n, 9) - o_ref[10*n+9] , 2.0));
       //
  }

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}

