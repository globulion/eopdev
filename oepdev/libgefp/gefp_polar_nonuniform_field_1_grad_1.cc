#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // Three blocks: Electric field parameters (1 block), Electric field gradient parameters (2 blocks)
  nBlocks_ = 3;
  nParameters_ = 3*(nSites_ * 3);
  for (int z=0; z<nBlocks_; ++z) nParametersBlock_.push_back(nSites_ * 3);
}
oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::LinearGradientNonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::LinearGradientNonUniformEFieldPolarGEFactory::~LinearGradientNonUniformEFieldPolarGEFactory()
{
}
// Implementations of abstract methods from base class
void oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute(void)
//{
//}
