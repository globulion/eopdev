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
   for (int s=0; s<nSites_; ++s) {
   for (int z=0; z<3; ++z) {
        int sz = 3*s + z;
        double g = 0.0;
        for (int n=0; n<nSamples_; ++n) {                            
             double dij = -referenceStatisticalSet_.DensityMatrixSet[n]->get(i, j);
             g += dij * electricFieldSet_[n][s]->get(z);
        }
        Gradient_->set(sz, 0, g);
   }}
   Gradient_->scale(2.0);

}
void oepdev::LinearNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
   double** H = Hessian_->pointer();
   for (int s1=0; s1<nSites_; ++s1) {
   for (int z1=0; z1<3; ++z1) {
        int s1z1 = 3*s1 + z1;
        for (int s2=0; s2<nSites_; ++s2) {
        for (int z2=0; z2<3; ++z2) {
             int s2z2 = 3*s2 + z2;
             double v = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  v += electricFieldSet_[n][s1]->get(z1) * electricFieldSet_[n][s2]->get(z2);
             }
             H[s1z1][s2z2] = v;
        }}
   }}
   Hessian_->scale(2.0);
}
