#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearUniformEFieldPolarGEFactory::LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::UniformEFieldPolarGEFactory(cphf, opt)
{
  hasDipolePolarizability_ = true;
  // One block: Electric field parameters
  nBlocks_ = 1;
  nParameters_ = 3;
  nParametersBlock_.push_back(3);
  // Allocate memory
  allocate();//-> add later
  //Gradient_   = std::make_shared<psi::Matrix>("Gradient"  , 3, 1);
  //Hessian_    = std::make_shared<psi::Matrix>("Hessian"   , 3, 3);
  //Parameters_ = std::make_shared<psi::Matrix>("Parameters", 3, 1);
}
oepdev::LinearUniformEFieldPolarGEFactory::LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::LinearUniformEFieldPolarGEFactory(cphf, cphf->options())
{
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
        double dij = guessDensityMatrixSet_[n]->get(i, j) - referenceDensityMatrixSet_[n]->get(i, j);
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
//std::shared_ptr<oepdev::GenEffPar> oepdev::LinearUniformEFieldPolarGEFactory::compute(void)
//{
//}
