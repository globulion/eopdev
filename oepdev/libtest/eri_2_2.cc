#include <iostream>
#include "test.h"

using namespace std;


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

