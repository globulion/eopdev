#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(cphf, opt)
{
   // Atoms are assumed to be distributed centres
   nSites_ = wfn_->molecule()->natom();
}
oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::NonUniformEFieldPolarGEFactory::~NonUniformEFieldPolarGEFactory()
{
}
// implementations of abstract methods from base
void oepdev::NonUniformEFieldPolarGEFactory::compute_samples(void)
{
  // TODO
}
// abstract methods
void oepdev::NonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::NonUniformEFieldPolarGEFactory::compute_hessian(void)
{
}
