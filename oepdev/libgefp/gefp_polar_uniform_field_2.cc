#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticUniformEFieldPolarGEFactory::QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::UniformEFieldPolarGEFactory(cphf, opt)
{
  // Two blocks: Electric field parameters (1 block), Squared electric field parameters (1 block)
  hasDipolePolarizability_ = true;
  hasDipoleDipoleHyperpolarizability_ = true;
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
   double g1_x = 0.0;
   double g1_y = 0.0;
   double g1_z = 0.0;
   double g2_x = 0.0;
   double g2_y = 0.0;
   double g2_z = 0.0;
   for (int n=0; n<nSamples_; ++n) {
        double fx = electricFieldSet_[n][0]->get(0);
        double fy = electricFieldSet_[n][0]->get(1);
        double fz = electricFieldSet_[n][0]->get(2);
        double fs = electricFieldSumSet_[n][0] * 2.0 * mField_;
        double dij=-referenceDensityMatrixSet_[n]->get(i, j);
        g1_x += dij * fx;
        g1_y += dij * fy;
        g1_z += dij * fz;
        g2_x += dij * fx * fs;
        g2_y += dij * fy * fs;
        g2_z += dij * fz * fs;
   }
   Gradient_->set(0, 0, g1_x);
   Gradient_->set(1, 0, g1_y);
   Gradient_->set(2, 0, g1_z);
   Gradient_->set(3, 0, g2_x);
   Gradient_->set(4, 0, g2_y);
   Gradient_->set(5, 0, g2_z);
   Gradient_->scale(2.0);
}
void oepdev::QuadraticUniformEFieldPolarGEFactory::compute_hessian(void)
{
   double** H = Hessian_->pointer();
   for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             double v_AA = 0.0;
             double v_BB = 0.0;
             double v_AB = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  double fz1 = electricFieldSet_[n][0]->get(z1);
                  double fz2 = electricFieldSet_[n][0]->get(z2);
                  double fs  = electricFieldSumSet_[n][0] * 2.0 * mField_;
                  v_AA += fz1 * fz2;
                  v_AB += fz1 * fz2 * fs;
                  v_BB += fz1 * fz2 * fs * fs;
             }
             H[z1  ][z2  ] = v_AA;
             H[z1+3][z2+3] = v_BB;
             H[z1  ][z2+3] = v_AB;
             H[z1+3][z2  ] = v_AB;
        }
   }
   Hessian_->scale(2.0);
}
