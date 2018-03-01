#include <iostream>
#include "test.h"
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
  if      (options_.get_str("OEPDEV_TEST_NAME")=="ERI_1_1") result = test_eri_1_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
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

        oepdev::SharedShellsIterator shellIter = oepdev::ShellCombinationsIterator::build(fact_oepdev, "ALL", 2);
        const double * buffer_1_1 = eri_1_1->buffer();

        for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
        {
             shellIter->compute_shell(eri_1_1);
             oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");
             for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
             {
                  int i = intsIter->i();int j = intsIter->j();
                  double integral = buffer_1_1[intsIter->index()];
                  eri_1_1_oepdev.set(i,j,integral);
             }
        }
        psi::timer_off(" Test: Computation of OepDev ERI_1_1");

        // Psi4 implementation of ERI's
        psi::timer_on(" Test: Computation of Psi4   ERI_1_1");
        oepdev::IntegralFactory fact_psi4(wfn_->basisset(),psi::BasisSet::zero_ao_basis_set());

        std::shared_ptr<oepdev::TwoBodyAOInt> eri_2_2(fact_psi4.eri_2_2());

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


        // Print ERI's and errors
        eri_1_1_oepdev.print();
        eri_1_1_psi4.print();
        error.print();

        // Compute squared errors
        error.copy(eri_1_1_oepdev.clone());
        error.subtract(&eri_1_1_psi4); 
        double result = error.sum_of_squares();

        // Print result
        std::cout << std::fixed;
        std::cout.precision(8);
        std::cout << " Test result= " << result << std::endl;

        // Finish
        return result;
}
