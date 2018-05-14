#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticNonUniformEFieldPolarGEFactory::QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // Two blocks: Electric field parameters (1 block), Squared electric field parameters (1 block)
  hasDipolePolarizability_ = true;
  hasDipoleDipoleHyperpolarizability_ = true;
  // Allocate memory
  allocate();
}
oepdev::QuadraticNonUniformEFieldPolarGEFactory::QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::QuadraticNonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::QuadraticNonUniformEFieldPolarGEFactory::~QuadraticNonUniformEFieldPolarGEFactory()
{
}
// Implementations of abstract methods from base class
void oepdev::QuadraticNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
  // TODO
}
void oepdev::QuadraticNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
  // TODO
}
