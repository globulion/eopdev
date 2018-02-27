#include "test.h"
#include <iostream>
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
  if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
  else throw psi::PSIEXCEPTION("Incorrect test name specified!");
  return result;
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
        psi::timer_on(" Test: Computation of OEPDEV ERI");
        oepdev::IntegralFactory fact_oepdev(wfn_->basisset());
        std::shared_ptr<psi::TwoBodyAOInt> eri_2_2(fact_oepdev.eri_2_2());

        oepdev::AllAOShellCombinationsIterator shellIter(fact_oepdev);
        int i, j, k, l; 
        double integral;
        const double * buffer_2_2 = eri_2_2->buffer();

        int icount = 0;
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
        {
             shellIter.compute_shell(eri_2_2);
             oepdev::AllAOIntegralsIterator intsIter(shellIter);
             for (intsIter.first(); intsIter.is_done() == false; intsIter.next())
             {
                  int i = intsIter.i();int j = intsIter.j();int k = intsIter.k();int l = intsIter.l();
                  double integral = buffer_2_2[intsIter.index()];
                  errors[icount] = integral;
                  ++icount;
             }
        }
        psi::timer_off(" Test: Computation of OEPDEV ERI");

        // Psi4 implementation of ERI's
        psi::timer_on(" Test: Computation of PSI4 ERI");
        oepdev::IntegralFactory fact_psi(wfn_->basisset());
        std::shared_ptr<psi::TwoBodyAOInt> eri(fact_psi.eri());
        const double* buffer = eri->buffer();

        icount = 0;
        oepdev::AllAOShellCombinationsIterator sIter(fact_psi);
        for (sIter.first(); sIter.is_done() == false; sIter.next())
        {
             sIter.compute_shell(eri);
             oepdev::AllAOIntegralsIterator intsIter(sIter);
             for (intsIter.first(); intsIter.is_done() == false; intsIter.next())
             {
                  int i = intsIter.i();int j = intsIter.j();int k = intsIter.k();int l = intsIter.l();
                  double integral = buffer[intsIter.index()];
                  errors[icount]-=integral;
                  ++icount;
             }
        }
        psi::timer_off(" Test: Computation of PSI4 ERI");

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
}
