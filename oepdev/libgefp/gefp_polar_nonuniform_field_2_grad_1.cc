#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // Five blocks: Electric field parameters (1 block), Squared electric field parameters (2 blocks), Electric field gradient parameters (2 blocks)
  nBlocks_ = 5;
  nParameters_ = 5*(nSites_ * 3);
  for (int z=0; z<nBlocks_; ++z) nParametersBlock_.push_back(nSites_ * 3);
  // Allocate memory
  allocate();
}
oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::~QuadraticGradientNonUniformEFieldPolarGEFactory()
{
}
// Implementations of abstract methods from base class
void oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
  // TODO
}
void oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
  // TODO
}
