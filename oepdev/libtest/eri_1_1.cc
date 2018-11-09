#include <iostream>
#include "test.h"

using namespace std;


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
