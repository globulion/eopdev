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
