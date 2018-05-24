#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearUniformEFieldPolarGEFactory::LinearUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::UniformEFieldPolarGEFactory(wfn, opt)
{
  // One block: Electric field parameters
  hasDipolePolarizability_ = true;

  // Allocate memory
  allocate();
}
oepdev::LinearUniformEFieldPolarGEFactory::~LinearUniformEFieldPolarGEFactory()
{

}
// implementations of abstract methods from base
void oepdev::LinearUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
   double g_x = 0.0;
   double g_y = 0.0;
   double g_z = 0.0;
   for (int n=0; n<nSamples_; ++n) {
        double dij = -referenceStatisticalSet_.DensityMatrixSet[n]->get(i, j);
        g_x += dij * electricFieldSet_[n][0]->get(0);
        g_y += dij * electricFieldSet_[n][0]->get(1);
        g_z += dij * electricFieldSet_[n][0]->get(2);
   }
   Gradient_->set(0, 0, g_x);
   Gradient_->set(1, 0, g_y);
   Gradient_->set(2, 0, g_z);
   Gradient_->scale(2.0);
}
void oepdev::LinearUniformEFieldPolarGEFactory::compute_hessian(void)
{
   double** H = Hessian_->pointer();
   for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             double v = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  v += electricFieldSet_[n][0]->get(z1) * electricFieldSet_[n][0]->get(z2);
             }
             H[z1][z2] = v;
        }
   }
   Hessian_->scale(2.0);
}
