#include "psi4/libpsi4util/exception.h"
#include "integral.h"
#include "../libints/eri.h"

using namespace std;

oepdev::TwoBodyAOInt::TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv)
 : psi::TwoBodyAOInt(intsfactory, deriv)
{
}
oepdev::TwoBodyAOInt::TwoBodyAOInt(const TwoBodyAOInt & rhs)
 : psi::TwoBodyAOInt(rhs)
{
}
oepdev::TwoBodyAOInt::~TwoBodyAOInt()
{
}
void oepdev::TwoBodyAOInt::compute(std::shared_ptr<psi::Matrix>& result, int ibs1, int ibs2) 
{
    compute(*result, ibs1, ibs2);
}
void oepdev::TwoBodyAOInt::compute(psi::Matrix& result, int ibs1, int ibs2) 
{
  std::shared_ptr<psi::BasisSet> bs1, bs2;
  if      (ibs1==0) bs1 = bs1_;
  else if (ibs1==1) bs1 = bs2_;
  else if (ibs1==2) bs1 = bs3_;
  else if (ibs1==3) bs1 = bs4_;
  else throw psi::PSIEXCEPTION(" Incorrect first basis set axis ID chosen!");
  if      (ibs2==0) bs2 = bs1_;
  else if (ibs2==1) bs2 = bs2_;
  else if (ibs2==2) bs2 = bs3_;
  else if (ibs2==3) bs2 = bs4_;
  else throw psi::PSIEXCEPTION(" Incorrect second basis set axis ID chosen!");
    // Do not worry about zeroing out result
    int ns1 = bs1->nshell();
    int ns2 = bs2->nshell();

    int i_offset = 0;
    double *location;

    // Leave as this full double for loop. We could be computing nonsymmetric integrals
    for (int i = 0; i < ns1; ++i) {
        int ni = bs1->shell(i).ncartesian();
        int j_offset = 0;
        for (int j = 0; j < ns2; ++j) {
            int nj = bs2->shell(j).ncartesian();

            // Compute the shell (automatically transforms to pure am in needed)
            compute_shell(i, j);

            // For each integral that we got put in its contribution
            location = target_full_;
            for (int p = 0; p < ni; ++p) {
                for (int q = 0; q < nj; ++q) {
                    result.add(0, i_offset + p, j_offset + q, *location);
                    location++;
                }
            }
            j_offset += nj;
        }
        i_offset += ni;
    }
}
size_t oepdev::TwoBodyAOInt::compute_shell(int, int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell(int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell(int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv1(int, int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv2(int, int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv1(int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv2(int, int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv1(int, int) {}
size_t oepdev::TwoBodyAOInt::compute_shell_deriv2(int, int) {}
oepdev::IntegralFactory::IntegralFactory(std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2,
                                         std::shared_ptr<psi::BasisSet> bs3, std::shared_ptr<psi::BasisSet> bs4) :
 psi::IntegralFactory(bs1, bs2, bs3, bs4) 
{
}
oepdev::IntegralFactory::IntegralFactory(std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2) :
 psi::IntegralFactory(bs1, bs2) 
{
}
oepdev::IntegralFactory::IntegralFactory(std::shared_ptr<psi::BasisSet> bs1) :
 psi::IntegralFactory(bs1)
{
}
oepdev::IntegralFactory::~IntegralFactory() 
{
}
oepdev::TwoBodyAOInt* oepdev::IntegralFactory::eri_1_1(int deriv, bool use_shell_pairs)
{
  return new oepdev::ERI_1_1(this, deriv, use_shell_pairs);
}
oepdev::TwoBodyAOInt* oepdev::IntegralFactory::eri_2_1(int deriv, bool use_shell_pairs)
{
  throw psi::PSIEXCEPTION("ERI 2_1 is not implemented yet!");
}
oepdev::TwoBodyAOInt* oepdev::IntegralFactory::eri_2_2(int deriv, bool use_shell_pairs)
{
  return new oepdev::ERI_2_2(this, deriv, use_shell_pairs);
}
oepdev::TwoBodyAOInt* oepdev::IntegralFactory::eri_3_1(int deriv, bool use_shell_pairs)
{
  return new oepdev::ERI_3_1(this, deriv, use_shell_pairs);
}
