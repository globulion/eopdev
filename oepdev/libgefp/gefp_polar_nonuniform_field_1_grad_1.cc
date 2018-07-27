#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(wfn, opt)
{
  // Two blocks: 
  //  - Electric field parameters
  //  - Electric field gradient parameters
  hasDipolePolarizability_ = true;
  hasQuadrupolePolarizability_ = true;

  // Allocate memory
  allocate();
}
oepdev::LinearGradientNonUniformEFieldPolarGEFactory::~LinearGradientNonUniformEFieldPolarGEFactory()
{

}
// Implementations of abstract methods from base class
void oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
   const int d = nParametersBlock_[0];

   for (int s=0; s<nSites_; ++s) {
        double g1_x = 0.0;                                                        
        double g1_y = 0.0;
        double g1_z = 0.0;
        double g2_xx= 0.0;
        double g2_yy= 0.0;
        double g2_zz= 0.0;
        double g2_xy= 0.0;
        double g2_xz= 0.0;
        double g2_yz= 0.0;
                                                                                  
        for (int n=0; n<nSamples_; ++n) {
             double fx = electricFieldSet_[n][s]->get(0);
             double fy = electricFieldSet_[n][s]->get(1);
             double fz = electricFieldSet_[n][s]->get(2);
             double fxx= electricFieldGradientSet_[n][s]->get(0, 0);
             double fxy= electricFieldGradientSet_[n][s]->get(0, 1);
             double fxz= electricFieldGradientSet_[n][s]->get(0, 2);
             double fyy= electricFieldGradientSet_[n][s]->get(1, 1);
             double fyz= electricFieldGradientSet_[n][s]->get(1, 2);
             double fzz= electricFieldGradientSet_[n][s]->get(2, 2);
             double dij=-referenceStatisticalSet_.DensityMatrixSet[n]->get(i, j);
             g1_x += dij * fx;
             g1_y += dij * fy;
             g1_z += dij * fz;
             g2_xx += dij * fxx;
             g2_yy += dij * fyy;
             g2_zz += dij * fzz;
             g2_xy += dij * fxy * 2.0;
             g2_xz += dij * fxz * 2.0;
             g2_yz += dij * fyz * 2.0;
        }
        Gradient_->set(3*s+0, 0, g1_x);
        Gradient_->set(3*s+1, 0, g1_y);
        Gradient_->set(3*s+2, 0, g1_z);
        Gradient_->set(d+6*s+0, 0, g2_xx);
        Gradient_->set(d+6*s+1, 0, g2_xy);
        Gradient_->set(d+6*s+2, 0, g2_xz);
        Gradient_->set(d+6*s+3, 0, g2_yy);
        Gradient_->set(d+6*s+4, 0, g2_yz);
        Gradient_->set(d+6*s+5, 0, g2_zz);
   }
   Gradient_->scale(2.0);
   Gradient_->print();

   //for (int s=0; s<nSites_; ++s) {
   //for (int z=0; z<3; ++z) {
   //     int sz = 3*s + z;
   //     double g1 = 0.0;
   //     double g2 = 0.0;
   //     for (int n=0; n<nSamples_; ++n) {                                                                  
   //          double dij = -referenceStatisticalSet_.DensityMatrixSet[n]->get(i, j);
   //          double fz = electricFieldSet_[n][s]->get(z);
   //          double gs = electricFieldGradientSumSet_[n][s]->get(z) * 2.0;
   //          g1 += dij * fz;
   //          g2 += dij * gs;
   //     }
   //     Gradient_->set(sz  , 0, g1);
   //     Gradient_->set(sz+d, 0, g2);
   //}}
   //Gradient_->scale(2.0);
}
void oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute_hessian(void)
{
   double** H = Hessian_->pointer();
   const int d = nParametersBlock_[0];

   // Block AA
   for (int s1=0; s1<nSites_; ++s1) {
   for (int s2=0; s2<nSites_; ++s2) {

   for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             double v_AA = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  double fz1 = electricFieldSet_[n][s1]->get(z1);
                  double fz2 = electricFieldSet_[n][s2]->get(z2);
                  v_AA += fz1 * fz2;
             }
             H[3*s1 + z1  ][3*s2 + z2  ] = v_AA;
        }
   }
   // Block AB
   for (int z1=0; z1<3; ++z1) {
       int i = 0;
       for (int z2=0; z2<3; ++z2) {
       for (int z3=z2; z3<3; ++z3) {
            double v_AB = 0.0;
            for (int n=0; n<nSamples_; ++n) {
                 double fz1 = electricFieldSet_[n][s1]->get(z1);
                 double fz2z3 = electricFieldGradientSet_[n][s2]->get(z2,z3);
                 double r = symmetryNumber_[i];
                 v_AB += fz1 * fz2z3 * r;
            }
            H[3*s1 + z1][d+6*s2 + i] = v_AB;
            i+=1;
       }}
   }
   // Block BA
   {int i = 0;
   for (int z1=0 ; z1<3; ++z1) {
   for (int z2=z1; z2<3; ++z2) {
        for (int z3=0; z3<3; ++z3) {
             double v_BA = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  double fz1z2 = electricFieldGradientSet_[n][s1]->get(z1,z2);
                  double fz3 = electricFieldSet_[n][s2]->get(z3);
                  double r = symmetryNumber_[i];
                  v_BA += fz1z2 * fz3 * r;
             }
             H[d+6*s1 + i][3*s2 + z3] = v_BA;
        }
        i+=1;
   }}
   }
   // Block BB
   {int i = 0;
   for (int z1=0 ; z1<3; ++z1) {
   for (int z2=z1; z2<3; ++z2) {
        int j = 0;
        for (int z3=0 ;z3<3; ++z3) {
        for (int z4=z3;z4<3; ++z4) {
             double v_BB = 0.0;
             for (int n=0; n<nSamples_; ++n) {
                  double fz1z2 = electricFieldGradientSet_[n][s1]->get(z1,z2);
                  double fz3z4 = electricFieldGradientSet_[n][s2]->get(z3,z4);
                  double r = symmetryNumber_[i];
                  double s = symmetryNumber_[j];
                  v_BB += fz1z2 * fz3z4 * r * s;
             }
             H[d+6*s1 + i][d+6*s2 + j] = v_BB;
             j += 1;
        }} 
        i += 1;
   }}
   }

   }}

   //for (int s1=0; s1<nSites_; ++s1) {
   //for (int z1=0; z1<3; ++z1) {
   //     int s1z1 = 3*s1 + z1;
   //     for (int s2=0; s2<nSites_; ++s2) {
   //     for (int z2=0; z2<3; ++z2) {
   //          int s2z2 = 3*s2 + z2;
   //          double v_AA = 0.0;
   //          double v_BB = 0.0;
   //          double v_AB = 0.0;
   //          double v_BA = 0.0;
   //          for (int n=0; n<nSamples_; ++n) {
   //               double fz1 = electricFieldSet_[n][s1]->get(z1);
   //               double fz2 = electricFieldSet_[n][s2]->get(z2);
   //               double gs1 = electricFieldGradientSumSet_[n][s1]->get(z1) * 2.0;
   //               double gs2 = electricFieldGradientSumSet_[n][s2]->get(z2) * 2.0;

   //               v_AA += fz1 * fz2;
   //               v_AB += fz1 * gs2;
   //               v_BA += gs1 * fz2;
   //               v_BB += gs1 * gs2;
   //          }
   //          H[s1z1  ][s2z2  ] = v_AA;
   //          H[s1z1+d][s2z2+d] = v_BB;
   //          H[s1z1  ][s2z2+d] = v_AB;
   //          H[s1z1+d][s2z2  ] = v_BA;
   //     }}
   //}}
   Hessian_->scale(2.0);
}
