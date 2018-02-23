#include "psi4/libpsi4util/exception.h"
#include "integral.h"
#include "../libints/eri.h"

using namespace std;

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
psi::TwoBodyAOInt* oepdev::IntegralFactory::eri_2_2(int deriv, bool use_shell_pairs)
{
  return new oepdev::ERI_2_2(this, deriv, use_shell_pairs);
}
psi::TwoBodyAOInt* oepdev::IntegralFactory::eri_3_1(int deriv, bool use_shell_pairs)
{
  return new oepdev::ERI_3_1(this, deriv, use_shell_pairs);
}
psi::TwoBodyAOInt* oepdev::IntegralFactory::eri_2_1(int deriv, bool use_shell_pairs)
{
  throw psi::PSIEXCEPTION("ERI 2_1 is not implemented yet!");
}
psi::TwoBodyAOInt* oepdev::IntegralFactory::eri_1_1(int deriv, bool use_shell_pairs)
{
  throw psi::PSIEXCEPTION("ERI 1_1 is not implemented yet!");
}
