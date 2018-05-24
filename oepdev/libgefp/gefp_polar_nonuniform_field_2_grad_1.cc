#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(wfn, opt)
{
  // Three blocks: 
  //  - Electric field parameters
  //  - Squared electric field parameters
  //  - Electric field gradient parameters
  hasDipolePolarizability_ = true;
  hasDipoleDipoleHyperpolarizability_ = true;
  hasQuadrupolePolarizability_ = true;

  // Allocate memory
  allocate();
}
oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::~QuadraticGradientNonUniformEFieldPolarGEFactory()
{

}
// Implementations of abstract methods from base class
void oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
  const int d = nParametersBlock_[0];

  for (int s=0; s<nSites_; ++s) {
  for (int z=0; z<3; ++z) {
       int sz = 3*s + z;
       double g1 = 0.0;
       double g2 = 0.0;
       double g3 = 0.0;
       for (int n=0; n<nSamples_; ++n) {                                                                  
            double dij = -referenceStatisticalSet_.DensityMatrixSet[n]->get(i, j);
            double fz = electricFieldSet_[n][s]->get(z);
            double fs = electricFieldSumSet_[n][s] * 2.0 * mField_;
            double gs = electricFieldGradientSumSet_[n][s]->get(z) * 2.0;
            g1 += dij * fz;
            g2 += dij * fz * fs;
            g3 += dij * gs;
       }
       Gradient_->set(sz    , 0, g1);
       Gradient_->set(sz+d  , 0, g2);
       Gradient_->set(sz+2*d, 0, g3);
  }}
  Gradient_->scale(2.0);
}
void oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
   double** H = Hessian_->pointer();
   const int d = nParametersBlock_[0];

   for (int s1=0; s1<nSites_; ++s1) {
   for (int z1=0; z1<3; ++z1) {
        int s1z1 = 3*s1 + z1;
        for (int s2=0; s2<nSites_; ++s2) {
        for (int z2=0; z2<3; ++z2) {
             int s2z2 = 3*s2 + z2;
             double v_AA = 0.0;
             double v_BB = 0.0;
             double v_CC = 0.0;
             double v_AB = 0.0;
             double v_BA = 0.0;
             double v_AC = 0.0;
             double v_CA = 0.0;
             double v_BC = 0.0;
             double v_CB = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  double fz1 = electricFieldSet_[n][s1]->get(z1);
                  double fz2 = electricFieldSet_[n][s2]->get(z2);
                  double fs1 = electricFieldSumSet_[n][s1] * 2.0 * mField_;
                  double fs2 = electricFieldSumSet_[n][s2] * 2.0 * mField_;
                  double gs1 = electricFieldGradientSumSet_[n][s1]->get(z1) * 2.0;
                  double gs2 = electricFieldGradientSumSet_[n][s2]->get(z2) * 2.0;

                  v_AA += fz1 * fz2;
                  v_BB += fz1 * fz2 * fs1 * fs2;
                  v_CC += gs1 * gs2;
                  v_AB += fz1 * fz2 * fs2;
                  v_BA += fz1 * fz2 * fs1;
                  v_AC += fz1 * gs2;
                  v_CA += gs1 * fz2;
                  v_BC += fz1 * fs1 * gs2;
                  v_CB += gs1 * fz2 * fs2;
             }
             H[s1z1    ][s2z2    ] = v_AA;
             H[s1z1+d  ][s2z2+d  ] = v_BB;
             H[s1z1+d*2][s2z2+d*2] = v_CC;
             H[s1z1    ][s2z2+d  ] = v_AB;
             H[s1z1+d  ][s2z2    ] = v_BA;
             H[s1z1    ][s2z2+d*2] = v_AC;
             H[s1z1+d*2][s2z2    ] = v_CA;
             H[s1z1+d  ][s2z2+d*2] = v_BC;
             H[s1z1+d*2][s2z2+d  ] = v_CB;
        }}
   }}
   Hessian_->scale(2.0);
}
