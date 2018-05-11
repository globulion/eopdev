#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearNonUniformEFieldPolarGEFactory::LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
  // One block: Electric field parameters
  nBlocks_ = 1;
  nParameters_ = nSites_ * 3;
  nParametersBlock_.push_back(nSites_ * 3);
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
}
void oepdev::LinearNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::LinearNonUniformEFieldPolarGEFactory::compute(void)
//{
//}
