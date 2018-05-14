#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticUniformEFieldPolarGEFactory::QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::UniformEFieldPolarGEFactory(cphf, opt)
{
  // Three blocks: Electric field parameters (1 block), Squared electric field parameters (2 blocks)
  nBlocks_ = 3;
  nParameters_ = 9;
  for (int z=0; z<nBlocks_; ++z) nParametersBlock_.push_back(3);
  // Allocate memory
  allocate();
}
oepdev::QuadraticUniformEFieldPolarGEFactory::QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::QuadraticUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::QuadraticUniformEFieldPolarGEFactory::~QuadraticUniformEFieldPolarGEFactory()
{
}
// Implementations of abstract methods from base class
void oepdev::QuadraticUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
  // TODO
}
void oepdev::QuadraticUniformEFieldPolarGEFactory::compute_hessian(void)
{
  // TODO
}
