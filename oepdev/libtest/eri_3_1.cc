#include <iostream>
#include "test.h"

using namespace std;


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
