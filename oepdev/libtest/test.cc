#include <iostream>
#include "test.h"
#include "../libutil/cphf.h"
#include "psi4/libmints/matrix.h"
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
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_1_1") result = test_eri_1_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_3_1") result = test_eri_3_1();
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
