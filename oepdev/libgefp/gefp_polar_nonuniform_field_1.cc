#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearNonUniformEFieldPolarGEFactory::LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // One block: Electric field parameters
  hasDipolePolarizability_ = true;
  // Allocate memory
  allocate();
}
oepdev::LinearNonUniformEFieldPolarGEFactory::LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::LinearNonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::LinearNonUniformEFieldPolarGEFactory::~LinearNonUniformEFieldPolarGEFactory()
{
}
// Implementations of abstract methods from base class
void oepdev::LinearNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
  // TODO
}
void oepdev::LinearNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
  // TODO
}
