#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // Three blocks: Electric field parameters (1 block), Electric field gradient parameters (2 blocks)
  hasDipolePolarizability_ = true;
  hasQuadrupolePolarizability_ = true;
  // Allocate memory
  allocate();
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
  // TODO
}
void oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
  // TODO
}
