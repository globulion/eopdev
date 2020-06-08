#include "dmtp.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libpsi4util/process.h"
#include "../libutil/kabsch_superimposer.h"
#include <cassert>
#include <iostream>


namespace oepdev{

using namespace std;

// MultipoleConvergence

MultipoleConvergence::MultipoleConvergence(std::shared_ptr<DMTPole> dmtp1, std::shared_ptr<DMTPole> dmtp2,
                                           MultipoleConvergence::ConvergenceLevel max_clevel) 
 : dmtp_1_(dmtp1),
   dmtp_2_(dmtp2),
   max_clevel_(max_clevel),
   convergenceList_{}
{
   convergenceList_["qq"] = std::make_shared<psi::Matrix>("q-q term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["qD"] = std::make_shared<psi::Matrix>("q-D term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["DD"] = std::make_shared<psi::Matrix>("D-D term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["qQ"] = std::make_shared<psi::Matrix>("q-Q term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["DQ"] = std::make_shared<psi::Matrix>("D-Q term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["qO"] = std::make_shared<psi::Matrix>("q-O term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["QQ"] = std::make_shared<psi::Matrix>("Q-Q term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["DO"] = std::make_shared<psi::Matrix>("D-O term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
   convergenceList_["qH"] = std::make_shared<psi::Matrix>("q-H term",dmtp_1_->nDMTPs_,dmtp_2_->nDMTPs_);
}
MultipoleConvergence::~MultipoleConvergence() 
{
}
void MultipoleConvergence::compute(MultipoleConvergence::Property property)
{
  if      (property == MultipoleConvergence::Property::Energy   ) {this->compute_energy   ();}
  else if (property == MultipoleConvergence::Property::Potential) {this->compute_potential();}
}
std::shared_ptr<psi::Matrix> MultipoleConvergence::level(MultipoleConvergence::ConvergenceLevel clevel)
{
  std::shared_ptr<psi::Matrix> result = std::make_shared<psi::Matrix>("Result", dmtp_1_->nDMTPs_, dmtp_2_->nDMTPs_);
  // R-1
  result->add(convergenceList_["qq"]);
  // R-2
  if ((clevel > MultipoleConvergence::ConvergenceLevel::R1) && (clevel <= max_clevel_)) {
      result->add(convergenceList_["qD"]);
  }
  // R-3
  if ((clevel > MultipoleConvergence::ConvergenceLevel::R2) && (clevel <= max_clevel_)) {
      result->add(convergenceList_["qQ"]);
      result->add(convergenceList_["DD"]);
  }
  // R-4
  if ((clevel > MultipoleConvergence::ConvergenceLevel::R3) && (clevel <= max_clevel_)) {
      result->add(convergenceList_["qO"]);
      result->add(convergenceList_["DQ"]);
  }
  // R-5
  if ((clevel > MultipoleConvergence::ConvergenceLevel::R4) && (clevel <= max_clevel_)) {
      result->add(convergenceList_["qH"]);
      result->add(convergenceList_["DO"]);
      result->add(convergenceList_["QQ"]);
  }

  return result;
}
void MultipoleConvergence::compute_energy()
{
   double** r_A_ = dmtp_1_->centres()->pointer();
   double** r_B_ = dmtp_2_->centres()->pointer();

   for (int N1=0; N1<dmtp_1_->nDMTPs_; ++N1) {
   for (int N2=0; N2<dmtp_2_->nDMTPs_; ++N2) {
        double** q_A_ = dmtp_1_->charges_[N1]->pointer();
        double** q_B_ = dmtp_2_->charges_[N2]->pointer();
        double** D_A_ = dmtp_1_->dipoles_[N1]->pointer();
        double** D_B_ = dmtp_2_->dipoles_[N2]->pointer();
        double** Q_A_ = dmtp_1_->quadrupoles_[N1]->pointer();
        double** Q_B_ = dmtp_2_->quadrupoles_[N2]->pointer();
        double** O_A_ = dmtp_1_->octupoles_[N1]->pointer();
        double** O_B_ = dmtp_2_->octupoles_[N2]->pointer();
        double** H_A_ = dmtp_1_->hexadecapoles_[N1]->pointer();
        double** H_B_ = dmtp_2_->hexadecapoles_[N2]->pointer();

        double qq = 0.0; // R-1
        double qD = 0.0; // R-2
        double DD = 0.0; // R-3
        double qQ = 0.0; // R-3
        double DQ = 0.0; // R-4
        double qO = 0.0; // R-4
        double QQ = 0.0; // R-5
        double DO = 0.0; // R-5
        double qH = 0.0; // R-5
              
        for (int i=0; i<dmtp_1_->nSites_; ++i) {
             double rix = r_A_[i][0];
             double riy = r_A_[i][1];
             double riz = r_A_[i][2];

             double qi = q_A_[i][0];
             double dix, diy, diz;
             double Qixx, Qixy, Qixz, Qiyy, Qiyz, Qizz;
             double Oixxx, Oixxy, Oixxz, Oixyy, Oixyz, Oixzz, Oiyyy, Oiyyz, Oiyzz, Oizzz;
             double Hixxxx, Hixxxy, Hixxxz, Hixxyy, Hixxyz, Hixxzz, Hixyyy, Hixyyz, Hixyzz, 
                    Hixzzz, Hiyyyy, Hiyyyz, Hiyyzz, Hiyzzz, Hizzzz;

             if (dmtp_1_->hasDipoles_) {
                 dix = D_A_[i][0]; 
                 diy = D_A_[i][1];
                 diz = D_A_[i][2];
                 if (dmtp_1_->hasQuadrupoles_) {
                     Qixx = Q_A_[i][0] * 1.5; 
                     Qixy = Q_A_[i][1] * 1.5;
                     Qixz = Q_A_[i][2] * 1.5;
                     Qiyy = Q_A_[i][3] * 1.5;
                     Qiyz = Q_A_[i][4] * 1.5;
                     Qizz = Q_A_[i][5] * 1.5;
                     double t = (Qixx + Qiyy + Qizz) / 3.0;
                     Qixx -= t; Qiyy -= t; Qizz -= t;
                     if (dmtp_1_->hasOctupoles_) {
                         Oixxx = O_A_[i][0] * 2.5;
                         Oixxy = O_A_[i][1] * 2.5; 
                         Oixxz = O_A_[i][2] * 2.5; 
                         Oixyy = O_A_[i][3] * 2.5; 
                         Oixyz = O_A_[i][4] * 2.5; 
                         Oixzz = O_A_[i][5] * 2.5; 
                         Oiyyy = O_A_[i][6] * 2.5; 
                         Oiyyz = O_A_[i][7] * 2.5; 
                         Oiyzz = O_A_[i][8] * 2.5; 
                         Oizzz = O_A_[i][9] * 2.5; 
                         double tx = (Oixxx + Oixyy + Oixzz) / 5.0;
                         double ty = (Oixxy + Oiyyy + Oiyzz) / 5.0;
                         double tz = (Oixxz + Oiyyz + Oizzz) / 5.0;
                         Oixxx -= tx * 3.0;
                         Oixxy -= ty;
                         Oixxz -= tz;
                         Oixyy -= tx;
                         Oixzz -= tx;
                         Oiyyy -= ty * 3.0;
                         Oiyyz -= tz;
                         Oiyzz -= ty;
                         Oizzz -= tz * 3.0;
                         if (dmtp_1_->hasHexadecapoles_) {
                             Hixxxx = H_A_[i][ 0] * 4.375; 
                             Hixxxy = H_A_[i][ 1] * 4.375;
                             Hixxxz = H_A_[i][ 2] * 4.375;
                             Hixxyy = H_A_[i][ 3] * 4.375;
                             Hixxyz = H_A_[i][ 4] * 4.375;
                             Hixxzz = H_A_[i][ 5] * 4.375;
                             Hixyyy = H_A_[i][ 6] * 4.375;
                             Hixyyz = H_A_[i][ 7] * 4.375;
                             Hixyzz = H_A_[i][ 8] * 4.375;
                             Hixzzz = H_A_[i][ 9] * 4.375;
                             Hiyyyy = H_A_[i][10] * 4.375;
                             Hiyyyz = H_A_[i][11] * 4.375;
                             Hiyyzz = H_A_[i][12] * 4.375;
                             Hiyzzz = H_A_[i][13] * 4.375;
                             Hizzzz = H_A_[i][14] * 4.375;
                             double txx = (Hixxxx + Hixxyy + Hixxzz) / 7.0;
                             double txy = (Hixxxy + Hixyyy + Hixyzz) / 7.0;
                             double txz = (Hixxxz + Hixyyz + Hixzzz) / 7.0;
                             double tyy = (Hixxyy + Hiyyyy + Hiyyzz) / 7.0;
                             double tyz = (Hixxyz + Hiyyyz + Hiyzzz) / 7.0;
                             double tzz = (Hixxzz + Hiyyzz + Hizzzz) / 7.0; 
                             double th  = (txx + tyy + tzz) / 5.0;
                             Hixxxx-= 6.0 * txx - 3.0 * th;
                             Hixxxy-= 3.0 * txy;
                             Hixxxz-= 3.0 * txz;
                             Hixxyy-= txx + tyy - th;
                             Hixxyz-= tyz;
                             Hixxzz-= txx + tzz - th;
                             Hixyyy-= 3.0 * txy;
                             Hixyyz-= txz; 
                             Hixyzz-= txy;
                             Hixzzz-= 3.0 * txz;
                             Hiyyyy-= 6.0 * tyy - 3.0 * th;
                             Hiyyyz-= 3.0 * tyz;
                             Hiyyzz-= tyy + tzz - th;
                             Hiyzzz-= 3.0 * tyz;
                             Hizzzz-= 6.0 * tzz - 3.0 * th;
                         }
                     }
                 }
             }

             for (int j=0; j<dmtp_2_->nSites_; ++j) {
                  double rjx = r_B_[j][0]; 
                  double rjy = r_B_[j][1];
                  double rjz = r_B_[j][2];
                                          
                  double qj = q_B_[j][0];

                  double rjix = rjx - rix;
                  double rjiy = rjy - riy;
                  double rjiz = rjz - riz;

                  double rji1 = 1.0 / sqrt(rjix*rjix+rjiy*rjiy+rjiz*rjiz);

                  // R-1 TERMS
                  qq += qi * qj * rji1;

                  // R-2 TERMS
                  if (dmtp_2_->hasDipoles_) {
                    double djx = D_B_[j][0];                                                             
                    double djy = D_B_[j][1];
                    double djz = D_B_[j][2];
                                                                                                         
                    double rji2 = rji1 * rji1;
                    double rji3 = rji1 * rji2;
                                                                                                         
                    double dirji = dix * rjix + diy * rjiy + diz * rjiz;
                    double djrji = djx * rjix + djy * rjiy + djz * rjiz;
                    qD += ( dirji * qj - djrji * qi ) * rji3;
                                                                                                         
                    // R-3 TERMS
                    if (dmtp_2_->hasQuadrupoles_) {
                      double Qjxx = Q_B_[j][0] * 1.5;
                      double Qjxy = Q_B_[j][1] * 1.5;
                      double Qjxz = Q_B_[j][2] * 1.5;
                      double Qjyy = Q_B_[j][3] * 1.5;
                      double Qjyz = Q_B_[j][4] * 1.5;
                      double Qjzz = Q_B_[j][5] * 1.5;
                      double t = (Qjxx + Qjyy + Qjzz) / 3.0;
                      Qjxx -= t; Qjyy -= t; Qjzz -= t;
                                                                                                           
                      double rji5 = rji3 * rji2;
                      double didj = dix * djx + diy * djy + diz * djz;
                      double Qirji2 = Qixx * rjix * rjix + Qiyy * rjiy * rjiy + Qizz * rjiz * rjiz +
                               2.0 * (Qixy * rjix * rjiy + Qixz * rjix * rjiz + Qiyz * rjiy * rjiz);
                      double Qjrji2 = Qjxx * rjix * rjix + Qjyy * rjiy * rjiy + Qjzz * rjiz * rjiz +
                               2.0 * (Qjxy * rjix * rjiy + Qjxz * rjix * rjiz + Qjyz * rjiy * rjiz);
                      DD += -3.0 * dirji * djrji * rji5 + didj * rji3;
                                                                                                           
                      qQ += (qj * Qirji2 + qi * Qjrji2) * rji5;

                      // R-4 TERMS
                      if (dmtp_2_->hasOctupoles_) {
                        double rji7 = rji5 * rji2;                                                            
                        double Ojxxx = O_B_[j][0] * 2.5;
                        double Ojxxy = O_B_[j][1] * 2.5; 
                        double Ojxxz = O_B_[j][2] * 2.5; 
                        double Ojxyy = O_B_[j][3] * 2.5; 
                        double Ojxyz = O_B_[j][4] * 2.5; 
                        double Ojxzz = O_B_[j][5] * 2.5; 
                        double Ojyyy = O_B_[j][6] * 2.5; 
                        double Ojyyz = O_B_[j][7] * 2.5; 
                        double Ojyzz = O_B_[j][8] * 2.5; 
                        double Ojzzz = O_B_[j][9] * 2.5; 
                        double tx = (Ojxxx + Ojxyy + Ojxzz) / 5.0;
                        double ty = (Ojxxy + Ojyyy + Ojyzz) / 5.0;
                        double tz = (Ojxxz + Ojyyz + Ojzzz) / 5.0;
                        Ojxxx -= tx * 3.0;
                        Ojxxy -= ty;
                        Ojxxz -= tz;
                        Ojxyy -= tx;
                        Ojxzz -= tx;
                        Ojyyy -= ty * 3.0;
                        Ojyyz -= tz;
                        Ojyzz -= ty;
                        Ojzzz -= tz * 3.0;

                        double Ojrji3 = Ojxxx * rjix * rjix * rjix       +
                                        Ojxxy * rjix * rjix * rjiy * 3.0 +
                                        Ojxxz * rjix * rjix * rjiz * 3.0 +  
                                        Ojxyy * rjix * rjiy * rjiy * 3.0 + 
                                        Ojxyz * rjix * rjiy * rjiz * 6.0 + 
                                        Ojxzz * rjix * rjiz * rjiz * 3.0 + 
                                        Ojyyy * rjiy * rjiy * rjiy       + 
                                        Ojyyz * rjiy * rjiy * rjiz * 3.0 + 
                                        Ojyzz * rjiy * rjiz * rjiz * 3.0 + 
                                        Ojzzz * rjiz * rjiz * rjiz       ; 
                        double Oirji3 = Oixxx * rjix * rjix * rjix       +
                                        Oixxy * rjix * rjix * rjiy * 3.0 +
                                        Oixxz * rjix * rjix * rjiz * 3.0 +  
                                        Oixyy * rjix * rjiy * rjiy * 3.0 + 
                                        Oixyz * rjix * rjiy * rjiz * 6.0 + 
                                        Oixzz * rjix * rjiz * rjiz * 3.0 + 
                                        Oiyyy * rjiy * rjiy * rjiy       + 
                                        Oiyyz * rjiy * rjiy * rjiz * 3.0 + 
                                        Oiyzz * rjiy * rjiz * rjiz * 3.0 + 
                                        Oizzz * rjiz * rjiz * rjiz       ; 

                        DQ += -2.0 * (Qjxx * dix * rjix + Qjyy * diy * rjiy + Qjzz * diz * rjiz +
                                      Qjxy *(dix * rjiy + diy * rjix) +
                                      Qjxz *(dix * rjiz + diz * rjix) +
                                      Qjyz *(diy * rjiz + diz * rjiy) -
                                      Qixx * djx * rjix - Qiyy * djy * rjiy - Qizz * djz * rjiz -
                                      Qixy *(djx * rjiy + djy * rjix) -
                                      Qixz *(djx * rjiz + djz * rjix) -
                                      Qiyz *(djy * rjiz + djz * rjiy)  ) * rji5
                              +5.0 * ( dirji * Qjrji2 - djrji * Qirji2 ) * rji7;
                        qO += ( qj * Oirji3 - qi * Ojrji3 ) * rji7;

                        // R-5 TERMS
                        if (dmtp_2_->hasHexadecapoles_) {
                          double Hjxxxx = H_B_[j][ 0] * 4.375; 
                          double Hjxxxy = H_B_[j][ 1] * 4.375;
                          double Hjxxxz = H_B_[j][ 2] * 4.375;
                          double Hjxxyy = H_B_[j][ 3] * 4.375;
                          double Hjxxyz = H_B_[j][ 4] * 4.375;
                          double Hjxxzz = H_B_[j][ 5] * 4.375;
                          double Hjxyyy = H_B_[j][ 6] * 4.375;
                          double Hjxyyz = H_B_[j][ 7] * 4.375;
                          double Hjxyzz = H_B_[j][ 8] * 4.375;
                          double Hjxzzz = H_B_[j][ 9] * 4.375;
                          double Hjyyyy = H_B_[j][10] * 4.375;
                          double Hjyyyz = H_B_[j][11] * 4.375;
                          double Hjyyzz = H_B_[j][12] * 4.375;
                          double Hjyzzz = H_B_[j][13] * 4.375;
                          double Hjzzzz = H_B_[j][14] * 4.375;
                          double txx = (Hjxxxx + Hjxxyy + Hjxxzz) / 7.0;
                          double txy = (Hjxxxy + Hjxyyy + Hjxyzz) / 7.0;
                          double txz = (Hjxxxz + Hjxyyz + Hjxzzz) / 7.0;
                          double tyy = (Hjxxyy + Hjyyyy + Hjyyzz) / 7.0;
                          double tyz = (Hjxxyz + Hjyyyz + Hjyzzz) / 7.0;
                          double tzz = (Hjxxzz + Hjyyzz + Hjzzzz) / 7.0; 
                          double th  = (txx + tyy + tzz) / 5.0;
                          Hjxxxx-= 6.0 * txx - 3.0 * th;
                          Hjxxxy-= 3.0 * txy;
                          Hjxxxz-= 3.0 * txz;
                          Hjxxyy-= txx + tyy - th;
                          Hjxxyz-= tyz;
                          Hjxxzz-= txx + tzz - th;
                          Hjxyyy-= 3.0 * txy;
                          Hjxyyz-= txz; 
                          Hjxyzz-= txy;
                          Hjxzzz-= 3.0 * txz;
                          Hjyyyy-= 6.0 * tyy - 3.0 * th;
                          Hjyyyz-= 3.0 * tyz;
                          Hjyyzz-= tyy + tzz - th;
                          Hjyzzz-= 3.0 * tyz;
                          Hjzzzz-= 6.0 * tzz - 3.0 * th;

                          double rji9 = rji7 * rji2;

                          double QiQj = Qixx * Qjxx + Qiyy * Qjyy + Qizz * Qjzz +
                                  2.0 *(Qixy * Qjxy + Qixz * Qjxz + Qiyz * Qjyz);
                          double Qirjix = Qixx * rjix + Qixy * rjiy + Qixz * rjiz;
                          double Qirjiy = Qixy * rjix + Qiyy * rjiy + Qiyz * rjiz;
                          double Qirjiz = Qixz * rjix + Qiyz * rjiy + Qizz * rjiz;
                          double Qjrjix = Qjxx * rjix + Qjxy * rjiy + Qjxz * rjiz;
                          double Qjrjiy = Qjxy * rjix + Qjyy * rjiy + Qjyz * rjiz;
                          double Qjrjiz = Qjxz * rjix + Qjyz * rjiy + Qjzz * rjiz;

                          double Ojrji2di = Ojxxx * rjix * rjix * dix       +   
                                            Ojxxy *(rjix * rjix * diy + 2.0 * rjix * rjiy * dix) +
                                            Ojxxz *(rjix * rjix * diz + 2.0 * rjix * rjiz * dix) +
                                            Ojxyy *(rjiy * rjiy * dix + 2.0 * rjix * rjiy * diy) +
                                            Ojxyz *(rjix * rjiy * diz * 2.0 + rjix * rjiz * diy * 2.0 + rjiy * rjiz * dix * 2.0) +
                                            Ojxzz *(rjiz * rjiz * dix + 2.0 * rjix * rjiz * diz) +
                                            Ojyyy * rjiy * rjiy * diy       + 
                                            Ojyyz *(rjiy * rjiy * diz + 2.0 * rjiy * rjiz * diy) +
                                            Ojyzz *(rjiz * rjiz * diy + 2.0 * rjiy * rjiz * diz) +
                                            Ojzzz * rjiz * rjiz * diz       ;
                          double Oirji2dj = Oixxx * rjix * rjix * djx       +   
                                            Oixxy *(rjix * rjix * djy + 2.0 * rjix * rjiy * djx) +
                                            Oixxz *(rjix * rjix * djz + 2.0 * rjix * rjiz * djx) +
                                            Oixyy *(rjiy * rjiy * djx + 2.0 * rjix * rjiy * djy) +
                                            Oixyz *(rjix * rjiy * djz * 2.0 + rjix * rjiz * djy * 2.0 + rjiy * rjiz * djx * 2.0) +
                                            Oixzz *(rjiz * rjiz * djx + 2.0 * rjix * rjiz * djz) +
                                            Oiyyy * rjiy * rjiy * djy       + 
                                            Oiyyz *(rjiy * rjiy * djz + 2.0 * rjiy * rjiz * djy) +
                                            Oiyzz *(rjiz * rjiz * djy + 2.0 * rjiy * rjiz * djz) +
                                            Oizzz * rjiz * rjiz * djz       ;
                          double Hjrji4 = Hjxxxx * rjix * rjix * rjix * rjix       +     
                                          Hjxxxy * rjix * rjix * rjix * rjiy * 4.0 +     
                                          Hjxxxz * rjix * rjix * rjix * rjiz * 4.0 +     
                                          Hjxxyy * rjix * rjix * rjiy * rjiy * 6.0 +     
                                          Hjxxyz * rjix * rjix * rjiy * rjiz *12.0 +     
                                          Hjxxzz * rjix * rjix * rjiz * rjiz * 6.0 +     
                                          Hjxyyy * rjix * rjiy * rjiy * rjiy * 4.0 +
                                          Hjxyyz * rjix * rjiy * rjiy * rjiz *12.0 +
                                          Hjxyzz * rjix * rjiy * rjiz * rjiz *12.0 +
                                          Hjxzzz * rjix * rjiz * rjiz * rjiz * 4.0 + 
                                          Hjyyyy * rjiy * rjiy * rjiy * rjiy       +
                                          Hjyyyz * rjiy * rjiy * rjiy * rjiz * 4.0 +
                                          Hjyyzz * rjiy * rjiy * rjiz * rjiz * 6.0 +
                                          Hjyzzz * rjiy * rjiz * rjiz * rjiz * 4.0 +
                                          Hjzzzz * rjiz * rjiz * rjiz * rjiz       ;
                          double Hirji4 = Hixxxx * rjix * rjix * rjix * rjix       +     
                                          Hixxxy * rjix * rjix * rjix * rjiy * 4.0 +     
                                          Hixxxz * rjix * rjix * rjix * rjiz * 4.0 +     
                                          Hixxyy * rjix * rjix * rjiy * rjiy * 6.0 +     
                                          Hixxyz * rjix * rjix * rjiy * rjiz *12.0 +     
                                          Hixxzz * rjix * rjix * rjiz * rjiz * 6.0 +     
                                          Hixyyy * rjix * rjiy * rjiy * rjiy * 4.0 +
                                          Hixyyz * rjix * rjiy * rjiy * rjiz *12.0 +
                                          Hixyzz * rjix * rjiy * rjiz * rjiz *12.0 +
                                          Hixzzz * rjix * rjiz * rjiz * rjiz * 4.0 + 
                                          Hiyyyy * rjiy * rjiy * rjiy * rjiy       +
                                          Hiyyyz * rjiy * rjiy * rjiy * rjiz * 4.0 +
                                          Hiyyzz * rjiy * rjiy * rjiz * rjiz * 6.0 +
                                          Hiyzzz * rjiy * rjiz * rjiz * rjiz * 4.0 +
                                          Hizzzz * rjiz * rjiz * rjiz * rjiz       ;

                          QQ += (35.0 / 3.0) * Qirji2 * Qjrji2 * rji9
                               -(20.0 / 3.0) *(Qirjix * Qjrjix + Qirjiy * Qjrjiy + Qirjiz * Qjrjiz) * rji7
                               +( 2.0 / 3.0) * QiQj * rji5;
                          DO += -7.0 * (djrji * Oirji3 + dirji * Ojrji3) * rji9
                                +3.0 * (Oirji2dj + Ojrji2di) * rji7;
                          qH += (qi * Hjrji4 + qj * Hirji4) * rji9;
                        } // EndIfBhasHexadecapoles
                      } // EndIfBhasOctupoles
                    } // EndIfBhasQuadrupoles
                  } // EndIfBhasDipoles

             } // EndForDMTPCentres_B
        } // EndForDMTPCentres_A

        convergenceList_["qq"]->set(N1, N2, qq);
        convergenceList_["qD"]->set(N1, N2, qD);
        convergenceList_["DD"]->set(N1, N2, DD);
        convergenceList_["qQ"]->set(N1, N2, qQ);
        convergenceList_["DQ"]->set(N1, N2, DQ);
        convergenceList_["qO"]->set(N1, N2, qO);
        convergenceList_["QQ"]->set(N1, N2, QQ);
        convergenceList_["DO"]->set(N1, N2, DO);
        convergenceList_["qH"]->set(N1, N2, qH);

   }} // EndForDMTPs
}
void MultipoleConvergence::compute_potential()
{
  throw psi::PSIEXCEPTION("The potential from DMTP's is not implemented yet.");
}


// DMTPole

DMTPole::DMTPole(psi::SharedWavefunction wfn, int n) 
 : mol_(wfn->molecule()), 
   wfn_(wfn), 
   primary_(wfn->basisset()),
   nDMTPs_(n),
   name_("none"),
   order_(0),
   nSites_(0),
   hasCharges_(false),
   hasDipoles_(false),
   hasQuadrupoles_(false),
   hasOctupoles_(false),
   hasHexadecapoles_(false),
   centres_(nullptr),
   origins_(nullptr),
   charges_({}),
   dipoles_({}),
   quadrupoles_({}),
   octupoles_({}),
   hexadecapoles_({}),
   mpInts_({})
{

}
DMTPole::DMTPole(const DMTPole* d) {
 // Do shallow copy on those:
 mol_ = d->mol_;
 wfn_ = d->wfn_;
 primary_ = d->primary_;
 // Do deep copy on those:
 nDMTPs_ = d->nDMTPs_;
 name_ = d->name_;
 order_ = d->order_;
 nSites_ = d->nSites_;
 hasCharges_ = d->hasCharges_;
 hasDipoles_ = d->hasDipoles_;
 hasQuadrupoles_ = d->hasQuadrupoles_;
 hasOctupoles_ = d->hasOctupoles_;
 hasHexadecapoles_ = d->hasHexadecapoles_;
 centres_ = std::make_shared<psi::Matrix>(d->centres_);
 origins_ = std::make_shared<psi::Matrix>(d->origins_);
 this->copy_from(d);
}

void DMTPole::copy_from(const DMTPole* d) {
 //
 charges_.clear();
 for (unsigned int i=0; i<d->nDMTPs_; ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->charges_[i]);
      charges_.push_back(m);
 }
 //
 dipoles_.clear();
 for (unsigned int i=0; i<d->nDMTPs_; ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->dipoles_[i]);
      dipoles_.push_back(m);
 }
 //
 quadrupoles_.clear();
 for (unsigned int i=0; i<d->nDMTPs_; ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->quadrupoles_[i]);
      quadrupoles_.push_back(m);
 }
 //
 octupoles_.clear();
 for (unsigned int i=0; i<d->nDMTPs_; ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->octupoles_[i]);
      octupoles_.push_back(m);
 }
 //
 hexadecapoles_.clear();
 for (unsigned int i=0; i<d->nDMTPs_; ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->hexadecapoles_[i]);
      hexadecapoles_.push_back(m);
 }
 //
 mpInts_.clear();
 for (unsigned int i=0; i<d->mpInts_.size(); ++i) {
      psi::SharedMatrix m = std::make_shared<psi::Matrix>(d->mpInts_[i]);
      mpInts_.push_back(m);
 }
}
DMTPole::~DMTPole()
{

}
std::shared_ptr<DMTPole> DMTPole::build(const std::string& type, 
                                        std::shared_ptr<psi::Wavefunction> wfn, 
                                        int n)
{
  std::shared_ptr<DMTPole> dmtp;
  if (type == "CAMM") dmtp = std::make_shared<oepdev::CAMM>(wfn, n);
  else throw psi::PSIEXCEPTION("Invalid DMTP type requested.");
  return dmtp;
}
psi::SharedVector DMTPole::centre(int x) const {
 psi::SharedVector r = std::make_shared<psi::Vector>("", 3);
 for (int z=0; z<3; ++z) r->set(z, centres_->get(x, z));
 return r;
}
psi::SharedVector DMTPole::origin(int x) const {
 psi::SharedVector r = std::make_shared<psi::Vector>("", 3);
 for (int z=0; z<3; ++z) r->set(z, origins_->get(x, z));
 return r;
}
MultipoleConvergence::ConvergenceLevel DMTPole::determine_dmtp_convergence_level(const std::string& option)
{
  MultipoleConvergence::ConvergenceLevel clevel;
  if      (psi::Process::environment.options.get_str(option) == "R1") clevel = MultipoleConvergence::ConvergenceLevel::R1;
  else if (psi::Process::environment.options.get_str(option) == "R2") clevel = MultipoleConvergence::ConvergenceLevel::R2;
  else if (psi::Process::environment.options.get_str(option) == "R3") clevel = MultipoleConvergence::ConvergenceLevel::R3;
  else if (psi::Process::environment.options.get_str(option) == "R4") clevel = MultipoleConvergence::ConvergenceLevel::R4;
  else if (psi::Process::environment.options.get_str(option) == "R5") clevel = MultipoleConvergence::ConvergenceLevel::R5;
  else {throw psi::PSIEXCEPTION("Incorrect convergence level specified!");}
  return clevel;
}
void DMTPole::allocate()
{
  for (int i=0; i<nDMTPs_; ++i) {
       if (hasCharges_      )  {charges_      .push_back( std::make_shared<psi::Matrix>("DMTP 0-th order tensor", nSites_  , 1 ) );}
       else {charges_      .push_back( std::make_shared<psi::Matrix>() );}
       if (hasDipoles_      )  {dipoles_      .push_back( std::make_shared<psi::Matrix>("DMTP 1-th order tensor", nSites_  , 3 ) );}
       else {dipoles_      .push_back( std::make_shared<psi::Matrix>() );}
       if (hasQuadrupoles_  )  {quadrupoles_  .push_back( std::make_shared<psi::Matrix>("DMTP 2-th order tensor", nSites_  , 6 ) );}
       else {quadrupoles_  .push_back( std::make_shared<psi::Matrix>() );}
       if (hasOctupoles_    )  {octupoles_    .push_back( std::make_shared<psi::Matrix>("DMTP 3-th order tensor", nSites_  , 10) );}
       else {octupoles_    .push_back( std::make_shared<psi::Matrix>() );}
       if (hasHexadecapoles_)  {hexadecapoles_.push_back( std::make_shared<psi::Matrix>("DMTP 4-th order tensor", nSites_  , 15) );}
       else {hexadecapoles_.push_back( std::make_shared<psi::Matrix>() );}
  }
  centres_ = std::make_shared<psi::Matrix>("DMTP Centres", nSites_, 3);
  origins_ = std::make_shared<psi::Matrix>("DMTP Origins", nSites_, 3);
}
void DMTPole::compute(std::vector<psi::SharedMatrix> D, std::vector<bool> transition) {
 if (D.size() != nDMTPs_) throw psi::PSIEXCEPTION("The number of OED's does not match the allocated size of DMTP object!");
 for (int i=0; i<nDMTPs_; ++i) this->compute(D.at(i), transition.at(i), i);
}
void DMTPole::compute(void) {
  psi::SharedMatrix D = wfn_->Da()->clone(); D->add(wfn_->Db());
  this->compute(D, false, 0);
}
void DMTPole::compute_order(void) {
 order_ = (int)hasCharges_ + (int)hasDipoles_ + (int)hasQuadrupoles_ + (int)hasOctupoles_ + (int)hasHexadecapoles_ - 1;
}
void DMTPole::compute_integrals(void) {
 this->compute_order(); // computes order_
 std::shared_ptr<psi::IntegralFactory> integral = std::make_shared<psi::IntegralFactory>(wfn_->basisset());
 psi::MultipoleSymmetry mpsymm(order_, mol_, integral, wfn_->matrix_factory());
 mpInts_ = mpsymm.create_matrices("Multipole Integrals", true);
 std::shared_ptr<psi::OneBodyAOInt> aompOBI(integral->ao_multipoles(order_));
 aompOBI->set_origin(psi::Vector3());
 aompOBI->compute(mpInts_);
}
void DMTPole::recenter(psi::SharedMatrix new_origins, int i)
{
 psi::SharedMatrix dipoles_new = std::make_shared<psi::Matrix>(dipoles_[i]); 
 psi::SharedMatrix qdpoles_new = std::make_shared<psi::Matrix>(quadrupoles_[i]); 
 psi::SharedMatrix ocpoles_new = std::make_shared<psi::Matrix>(octupoles_[i]); 
 psi::SharedMatrix hdpoles_new = std::make_shared<psi::Matrix>(hexadecapoles_[i]); 
 double** cp = charges_[i]->pointer();
 double** mp = dipoles_[i]->pointer();
 double** mp_= dipoles_new->pointer();
 double** qp = quadrupoles_[i]->pointer();
 double** qp_= qdpoles_new->pointer();
 double** op = octupoles_[i]->pointer();
 double** op_= ocpoles_new->pointer();
 double** hp = hexadecapoles_[i]->pointer();
 double** hp_= hdpoles_new->pointer();

 // Recenter
 for (int ic=0; ic<nSites_; ++ic) {
      if (hasDipoles_) {
      double rx_o = origins_->get(ic, 0);
      double ry_o = origins_->get(ic, 1);
      double rz_o = origins_->get(ic, 2);
      double rx_n = new_origins->get(ic, 0);
      double ry_n = new_origins->get(ic, 1);
      double rz_n = new_origins->get(ic, 2);

      double d1x = rx_n - rx_o;
      double d1y = ry_n - ry_o;
      double d1z = rz_n - rz_o;

      double c = cp[ic][0];

      double mx_o = mp[ic][0];
      double my_o = mp[ic][1];
      double mz_o = mp[ic][2];

      double mx_n = mx_o - c * d1x;
      double my_n = my_o - c * d1y;
      double mz_n = mz_o - c * d1z;

      // Collect
      mp_[ic][0] = mx_n;
      mp_[ic][1] = my_n;
      mp_[ic][2] = mz_n;

      if (hasQuadrupoles_) {
      double rxx_n = rx_n * rx_n;  double rxx_o = rx_o * rx_o;
      double rxy_n = rx_n * ry_n;  double rxy_o = rx_o * ry_o;
      double rxz_n = rx_n * rz_n;  double rxz_o = rx_o * rz_o;
      double ryy_n = ry_n * ry_n;  double ryy_o = ry_o * ry_o;
      double ryz_n = ry_n * rz_n;  double ryz_o = ry_o * rz_o;
      double rzz_n = rz_n * rz_n;  double rzz_o = rz_o * rz_o;
      
      double d2xx = rxx_n - rxx_o;
      double d2xy = rxy_n - rxy_o;
      double d2xz = rxz_n - rxz_o;
      double d2yy = ryy_n - ryy_o;
      double d2yz = ryz_n - ryz_o;
      double d2zz = rzz_n - rzz_o;

      double qxx_o = qp[ic][0];
      double qxy_o = qp[ic][1];
      double qxz_o = qp[ic][2];
      double qyy_o = qp[ic][3];
      double qyz_o = qp[ic][4];
      double qzz_o = qp[ic][5];

      double qxx_n = qxx_o + c * d2xx - 2.0 * mx_o * d1x              - 2.0 * c * d1x * rx_o;
      double qxy_n = qxy_o + c * d2xy -       mx_o * d1y - my_o * d1x - c * (d1x * ry_o + d1y * rx_o);
      double qxz_n = qxz_o + c * d2xz -       mx_o * d1z - mz_o * d1x - c * (d1x * rz_o + d1z * rx_o);
      double qyy_n = qyy_o + c * d2yy - 2.0 * my_o * d1y              - 2.0 * c * d1y * ry_o;
      double qyz_n = qyz_o + c * d2yz -       my_o * d1z - mz_o * d1y - c * (d1y * rz_o + d1z * ry_o); 
      double qzz_n = qzz_o + c * d2zz - 2.0 * mz_o * d1z              - 2.0 * c * d1z * rz_o;

      // Collect
      qp_[ic][0] = qxx_n;
      qp_[ic][1] = qxy_n;
      qp_[ic][2] = qxz_n;
      qp_[ic][3] = qyy_n;
      qp_[ic][4] = qyz_n;
      qp_[ic][5] = qzz_n;

      if (hasOctupoles_) {
      double rxxx_n = rxx_n * rx_n;  double rxxx_o = rxx_o * rx_o;  
      double rxxy_n = rxx_n * ry_n;  double rxxy_o = rxx_o * ry_o;
      double rxxz_n = rxx_n * rz_n;  double rxxz_o = rxx_o * rz_o;
      double rxyy_n = rxy_n * ry_n;  double rxyy_o = rxy_o * ry_o;
      double rxyz_n = rxy_n * rz_n;  double rxyz_o = rxy_o * rz_o;
      double rxzz_n = rxz_n * rz_n;  double rxzz_o = rxz_o * rz_o;
      double ryyy_n = ryy_n * ry_n;  double ryyy_o = ryy_o * ry_o;
      double ryyz_n = ryy_n * rz_n;  double ryyz_o = ryy_o * rz_o;
      double ryzz_n = ryz_n * rz_n;  double ryzz_o = ryz_o * rz_o;
      double rzzz_n = rzz_n * rz_n;  double rzzz_o = rzz_o * rz_o;

      double d3xxx = rxxx_n - rxxx_o; 
      double d3xxy = rxxy_n - rxxy_o;
      double d3xxz = rxxz_n - rxxz_o;
      double d3xyy = rxyy_n - rxyy_o;
      double d3xyz = rxyz_n - rxyz_o;
      double d3xzz = rxzz_n - rxzz_o;
      double d3yyy = ryyy_n - ryyy_o;
      double d3yyz = ryyz_n - ryyz_o;
      double d3yzz = ryzz_n - ryzz_o;
      double d3zzz = rzzz_n - rzzz_o;

      double oxxx_o = op[ic][0];
      double oxxy_o = op[ic][1];
      double oxxz_o = op[ic][2];
      double oxyy_o = op[ic][3];
      double oxyz_o = op[ic][4];
      double oxzz_o = op[ic][5];
      double oyyy_o = op[ic][6];
      double oyyz_o = op[ic][7];
      double oyzz_o = op[ic][8];
      double ozzz_o = op[ic][9];

      double oxxx_n = oxxx_o  - c * d3xxx  + mx_o * d2xx + mx_o * d2xx + mx_o * d2xx  - qxx_o * d1x - qxx_o * d1x - qxx_o * d1x;   
      double oxxy_n = oxxy_o  - c * d3xxy  + mx_o * d2xy + my_o * d2xx + mx_o * d2xy  - qxy_o * d1x - qxx_o * d1y - qxy_o * d1x;
      double oxxz_n = oxxz_o  - c * d3xxz  + mx_o * d2xz + mz_o * d2xx + mx_o * d2xz  - qxz_o * d1x - qxx_o * d1z - qxz_o * d1x;
      double oxyy_n = oxyy_o  - c * d3xyy  + mx_o * d2yy + my_o * d2xy + my_o * d2xy  - qxy_o * d1y - qxy_o * d1y - qyy_o * d1x;
      double oxyz_n = oxyz_o  - c * d3xyz  + mx_o * d2yz + mz_o * d2xy + my_o * d2xz  - qxz_o * d1y - qxy_o * d1z - qyz_o * d1x;
      double oxzz_n = oxzz_o  - c * d3xzz  + mx_o * d2zz + mz_o * d2xz + mz_o * d2xz  - qxz_o * d1z - qxz_o * d1z - qzz_o * d1x;
      double oyyy_n = oyyy_o  - c * d3yyy  + my_o * d2yy + my_o * d2yy + my_o * d2yy  - qyy_o * d1y - qyy_o * d1y - qyy_o * d1y;
      double oyyz_n = oyyz_o  - c * d3yyz  + my_o * d2yz + mz_o * d2yy + my_o * d2yz  - qyz_o * d1y - qyy_o * d1z - qyz_o * d1y;
      double oyzz_n = oyzz_o  - c * d3yzz  + my_o * d2zz + mz_o * d2yz + mz_o * d2yz  - qyz_o * d1z - qyz_o * d1z - qzz_o * d1y;
      double ozzz_n = ozzz_o  - c * d3zzz  + mz_o * d2zz + mz_o * d2zz + mz_o * d2zz  - qzz_o * d1z - qzz_o * d1z - qzz_o * d1z;

      oxxx_n += 3.0 * c * rx_o * d2xx - 3.0 * d1x * rx_o * (2.0 * mx_o + c * rx_o);
      oxxy_n += c * (2.0 * rx_o * d2xy + ry_o * d2xx) - 2.0 * d1x * (c * rx_o * ry_o + rx_o * my_o + ry_o * mx_o) - d1y * rx_o * (2.0 * mx_o + c * rx_o);
      oxxz_n += c * (2.0 * rx_o * d2xz + rz_o * d2xx) - 2.0 * d1x * (c * rx_o * rz_o + rx_o * mz_o + rz_o * mx_o) - d1z * rx_o * (2.0 * mx_o + c * rx_o);
      oxyy_n += c * (2.0 * ry_o * d2xy + rx_o * d2yy) - 2.0 * d1y * (c * ry_o * rx_o + ry_o * mx_o + rx_o * my_o) - d1x * ry_o * (2.0 * my_o + c * ry_o);
      oxyz_n += c * (d2yz * rx_o + d2xy * rz_o + d2xz * ry_o) - d1x * (c * ry_o * rz_o + ry_o * mz_o + rz_o * my_o)
	                                                      - d1y * (c * rx_o * rz_o + rx_o * mz_o + rz_o * mx_o)
							      - d1z * (c * rx_o * ry_o + rx_o * my_o + ry_o * mx_o);
      oxzz_n += c * (2.0 * rz_o * d2xz + rx_o * d2zz) - 2.0 * d1z * (c * rz_o * rx_o + rz_o * mx_o + rx_o * mz_o) - d1x * rz_o * (2.0 * mz_o + c * rz_o);
      oyyy_n += 3.0 * c * ry_o * d2yy - 3.0 * d1y * ry_o * (2.0 * my_o + c * ry_o);
      oyyz_n += c * (2.0 * ry_o * d2yz + rz_o * d2yy) - 2.0 * d1y * (c * ry_o * rz_o + ry_o * mz_o + rz_o * my_o) - d1z * ry_o * (2.0 * my_o + c * ry_o);
      oyzz_n += c * (2.0 * rz_o * d2yz + ry_o * d2zz) - 2.0 * d1z * (c * rz_o * ry_o + rz_o * my_o + ry_o * mz_o) - d1y * rz_o * (2.0 * mz_o + c * rz_o);
      ozzz_n += 3.0 * c * rz_o * d2zz - 3.0 * d1z * rz_o * (2.0 * mz_o + c * rz_o);       

      // Collect
      op_[ic][0] = oxxx_n;
      op_[ic][1] = oxxy_n;
      op_[ic][2] = oxxz_n;
      op_[ic][3] = oxyy_n;
      op_[ic][4] = oxyz_n;
      op_[ic][5] = oxzz_n;
      op_[ic][6] = oyyy_n;
      op_[ic][7] = oyyz_n;
      op_[ic][8] = oyzz_n;
      op_[ic][9] = ozzz_n;

      if (hasHexadecapoles_) {
      double rxxxx_n = rxxx_n * rx_n;  double rxxxx_o = rxxx_o * rx_o;
      double rxxxy_n = rxxx_n * ry_n;  double rxxxy_o = rxxx_o * ry_o;
      double rxxxz_n = rxxx_n * rz_n;  double rxxxz_o = rxxx_o * rz_o;
      double rxxyy_n = rxxy_n * ry_n;  double rxxyy_o = rxxy_o * ry_o;
      double rxxyz_n = rxxy_n * rz_n;  double rxxyz_o = rxxy_o * rz_o;
      double rxxzz_n = rxxz_n * rz_n;  double rxxzz_o = rxxz_o * rz_o;
      double rxyyy_n = rxyy_n * ry_n;  double rxyyy_o = rxyy_o * ry_o;
      double rxyyz_n = rxyy_n * rz_n;  double rxyyz_o = rxyy_o * rz_o;
      double rxyzz_n = rxyz_n * rz_n;  double rxyzz_o = rxyz_o * rz_o;
      double rxzzz_n = rxzz_n * rz_n;  double rxzzz_o = rxzz_o * rz_o;
      double ryyyy_n = ryyy_n * ry_n;  double ryyyy_o = ryyy_o * ry_o;
      double ryyyz_n = ryyy_n * rz_n;  double ryyyz_o = ryyy_o * rz_o;
      double ryyzz_n = ryyz_n * rz_n;  double ryyzz_o = ryyz_o * rz_o;
      double ryzzz_n = ryzz_n * rz_n;  double ryzzz_o = ryzz_o * rz_o;
      double rzzzz_n = rzzz_n * rz_n;  double rzzzz_o = rzzz_o * rz_o;

      double d4xxxx = rxxxx_n - rxxxx_o;
      double d4xxxy = rxxxy_n - rxxxy_o;
      double d4xxxz = rxxxz_n - rxxxz_o;
      double d4xxyy = rxxyy_n - rxxyy_o;
      double d4xxyz = rxxyz_n - rxxyz_o;
      double d4xxzz = rxxzz_n - rxxzz_o;
      double d4xyyy = rxyyy_n - rxyyy_o;
      double d4xyyz = rxyyz_n - rxyyz_o;
      double d4xyzz = rxyzz_n - rxyzz_o;
      double d4xzzz = rxzzz_n - rxzzz_o;
      double d4yyyy = ryyyy_n - ryyyy_o;
      double d4yyyz = ryyyz_n - ryyyz_o;
      double d4yyzz = ryyzz_n - ryyzz_o;
      double d4yzzz = ryzzz_n - ryzzz_o;
      double d4zzzz = rzzzz_n - rzzzz_o;

      double hxxxx_o = hp[ic][ 0];
      double hxxxy_o = hp[ic][ 1];
      double hxxxz_o = hp[ic][ 2];
      double hxxyy_o = hp[ic][ 3];
      double hxxyz_o = hp[ic][ 4];
      double hxxzz_o = hp[ic][ 5];
      double hxyyy_o = hp[ic][ 6];
      double hxyyz_o = hp[ic][ 7];
      double hxyzz_o = hp[ic][ 8];
      double hxzzz_o = hp[ic][ 9];
      double hyyyy_o = hp[ic][10];
      double hyyyz_o = hp[ic][11];
      double hyyzz_o = hp[ic][12];
      double hyzzz_o = hp[ic][13];
      double hzzzz_o = hp[ic][14];

      double hxxxx_n = hxxxx_o + c * d4xxxx - mx_o * d3xxx - mx_o * d3xxx - mx_o * d3xxx - mx_o * d3xxx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx - oxxx_o * d1x - oxxx_o * d1x - oxxx_o * d1x - oxxx_o * d1x; 
      double hxxxy_n = hxxxy_o + c * d4xxxy - mx_o * d3xxy - mx_o * d3xxy - mx_o * d3xxy - my_o * d3xxx + qxx_o * d2xy + qxx_o * d2xy + qxy_o * d2xx + qxy_o * d2xx + qxy_o * d2xx + qxx_o * d2xy - oxxx_o * d1y - oxxy_o * d1x - oxxy_o * d1x - oxxy_o * d1x;
      double hxxxz_n = hxxxz_o + c * d4xxxz - mx_o * d3xxz - mx_o * d3xxz - mx_o * d3xxz - mz_o * d3xxx + qxx_o * d2xz + qxx_o * d2xz + qxz_o * d2xx + qxz_o * d2xx + qxz_o * d2xx + qxx_o * d2xz - oxxx_o * d1z - oxxz_o * d1x - oxxz_o * d1x - oxxz_o * d1x;
      double hxxyy_n = hxxyy_o + c * d4xxyy - mx_o * d3xyy - mx_o * d3xyy - my_o * d3xxy - my_o * d3xxy + qxx_o * d2yy + qxy_o * d2xy + qxy_o * d2xy + qyy_o * d2xx + qxy_o * d2xy + qxy_o * d2xy - oxxy_o * d1y - oxxy_o * d1y - oxyy_o * d1x - oxyy_o * d1x;
      double hxxyz_n = hxxyz_o + c * d4xxyz - mx_o * d3xyz - mx_o * d3xyz - my_o * d3xxz - mz_o * d3xxy + qxx_o * d2yz + qxy_o * d2xz + qxz_o * d2xy + qyz_o * d2xx + qxz_o * d2xy + qxy_o * d2xz - oxxy_o * d1z - oxxz_o * d1y - oxyz_o * d1x - oxyz_o * d1x;
      double hxxzz_n = hxxzz_o + c * d4xxzz - mx_o * d3xzz - mx_o * d3xzz - mz_o * d3xxz - mz_o * d3xxz + qxx_o * d2zz + qxz_o * d2xz + qxz_o * d2xz + qzz_o * d2xx + qxz_o * d2xz + qxz_o * d2xz - oxxz_o * d1z - oxxz_o * d1z - oxzz_o * d1x - oxzz_o * d1x;
      double hxyyy_n = hxyyy_o + c * d4xyyy - mx_o * d3yyy - my_o * d3xyy - my_o * d3xyy - my_o * d3xyy + qxy_o * d2yy + qxy_o * d2yy + qyy_o * d2xy + qyy_o * d2xy + qxy_o * d2yy + qyy_o * d2xy - oxyy_o * d1y - oxyy_o * d1y - oxyy_o * d1y - oyyy_o * d1x;
      double hxyyz_n = hxyyz_o + c * d4xyyz - mx_o * d3yyz - my_o * d3xyz - my_o * d3xyz - mz_o * d3xyy + qxy_o * d2yz + qxy_o * d2yz + qyz_o * d2xy + qyz_o * d2xy + qxz_o * d2yy + qyy_o * d2xz - oxyy_o * d1z - oxyz_o * d1y - oxyz_o * d1y - oyyz_o * d1x;
      double hxyzz_n = hxyzz_o + c * d4xyzz - mx_o * d3yzz - my_o * d3xzz - mz_o * d3xyz - mz_o * d3xyz + qxy_o * d2zz + qxz_o * d2yz + qyz_o * d2xz + qzz_o * d2xy + qxz_o * d2yz + qyz_o * d2xz - oxyz_o * d1z - oxyz_o * d1z - oxzz_o * d1y - oyzz_o * d1x;
      double hxzzz_n = hxzzz_o + c * d4xzzz - mx_o * d3zzz - mz_o * d3xzz - mz_o * d3xzz - mz_o * d3xzz + qxz_o * d2zz + qxz_o * d2zz + qzz_o * d2xz + qzz_o * d2xz + qxz_o * d2zz + qzz_o * d2xz - oxzz_o * d1z - oxzz_o * d1z - oxzz_o * d1z - ozzz_o * d1x;
      double hyyyy_n = hyyyy_o + c * d4yyyy - my_o * d3yyy - my_o * d3yyy - my_o * d3yyy - my_o * d3yyy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy - oyyy_o * d1y - oyyy_o * d1y - oyyy_o * d1y - oyyy_o * d1y;
      double hyyyz_n = hyyyz_o + c * d4yyyz - my_o * d3yyz - my_o * d3yyz - my_o * d3yyz - mz_o * d3yyy + qyy_o * d2yz + qyy_o * d2yz + qyz_o * d2yy + qyz_o * d2yy + qyz_o * d2yy + qyy_o * d2yz - oyyy_o * d1z - oyyz_o * d1y - oyyz_o * d1y - oyyz_o * d1y;
      double hyyzz_n = hyyzz_o + c * d4yyzz - my_o * d3yzz - my_o * d3yzz - mz_o * d3yyz - mz_o * d3yyz + qyy_o * d2zz + qyz_o * d2yz + qyz_o * d2yz + qzz_o * d2yy + qyz_o * d2yz + qyz_o * d2yz - oyyz_o * d1z - oyyz_o * d1z - oyzz_o * d1y - oyzz_o * d1y;
      double hyzzz_n = hyzzz_o + c * d4yzzz - my_o * d3zzz - mz_o * d3yzz - mz_o * d3yzz - mz_o * d3yzz + qyz_o * d2zz + qyz_o * d2zz + qzz_o * d2yz + qzz_o * d2yz + qyz_o * d2zz + qzz_o * d2yz - oyzz_o * d1z - oyzz_o * d1z - oyzz_o * d1z - ozzz_o * d1y;
      double hzzzz_n = hzzzz_o + c * d4zzzz - mz_o * d3zzz - mz_o * d3zzz - mz_o * d3zzz - mz_o * d3zzz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz - ozzz_o * d1z - ozzz_o * d1z - ozzz_o * d1z - ozzz_o * d1z;

      // TODO!
      hxxxx_o += -4.0 * d3xxx * c * rx_o + 6.0 * d2xx * rx_o * (c * rx_o + 2.0 * mx_o) - 4.0 * d1x * rx_o * (c * rx_o * rx_o + 3.0 * mx_o * rx_o + 3.0 * qxx_o);
      // other correction terms here!

      // Collect
      hp_[ic][ 0] = hxxxx_n; 
      hp_[ic][ 1] = hxxxy_n;
      hp_[ic][ 2] = hxxxz_n;
      hp_[ic][ 3] = hxxyy_n;
      hp_[ic][ 4] = hxxyz_n;
      hp_[ic][ 5] = hxxzz_n;
      hp_[ic][ 6] = hxyyy_n;
      hp_[ic][ 7] = hxyyz_n;
      hp_[ic][ 8] = hxyzz_n;
      hp_[ic][ 9] = hxzzz_n;
      hp_[ic][10] = hyyyy_n;
      hp_[ic][11] = hyyyz_n;
      hp_[ic][12] = hyyzz_n;
      hp_[ic][13] = hyzzz_n;
      hp_[ic][14] = hzzzz_n;

      }}}} // EndIFHasMultipoles
 }

 // Save
 if (hasDipoles_      ) dipoles_      [i]->copy(dipoles_new);
 if (hasQuadrupoles_  ) quadrupoles_  [i]->copy(qdpoles_new);
 if (hasOctupoles_    ) octupoles_    [i]->copy(ocpoles_new);
 if (hasHexadecapoles_) hexadecapoles_[i]->copy(hdpoles_new);
}
void DMTPole::recenter(psi::SharedMatrix new_origins)
{
 for (int i=0; i<nDMTPs_; ++i) this->recenter(new_origins, i);
 origins_->copy(new_origins);
}
std::shared_ptr<MultipoleConvergence> DMTPole::energy(std::shared_ptr<DMTPole> other, MultipoleConvergence::ConvergenceLevel max_clevel)
{
  std::shared_ptr<MultipoleConvergence> convergence = std::make_shared<MultipoleConvergence>(shared_from_this(), other, max_clevel);
  convergence->compute(MultipoleConvergence::Energy);
  return convergence;
}
std::shared_ptr<MultipoleConvergence> DMTPole::potential(std::shared_ptr<DMTPole> other, MultipoleConvergence::ConvergenceLevel max_clevel)
{
  std::shared_ptr<MultipoleConvergence> convergence = std::make_shared<MultipoleConvergence>(shared_from_this(), other, max_clevel);
  convergence->compute(MultipoleConvergence::Potential);
  //std::vector<double> potentials;
  //for (int i=0; i<nDMTPs_; ++i) potentials.push_back(0.0);
  //TODO
  // ... if (type == "R-5") ...
  // Return
  return convergence;
}
/// Translate the DMTP sets
void DMTPole::translate(psi::SharedVector transl) 
{
   double** c = centres_->pointer();
   double** o = origins_->pointer();
   for (int i=0; i<nSites_; ++i) {
	c[i][0] += transl->get(0); o[i][0] += transl->get(0); 
	c[i][1] += transl->get(1); o[i][1] += transl->get(1);
	c[i][2] += transl->get(2); o[i][2] += transl->get(2);
   }
}
/// Rotate the DMTP sets
void DMTPole::rotate(psi::SharedMatrix rotmat) 
{
  psi::SharedMatrix new_centres = psi::Matrix::doublet(this->centres_, rotmat, false, false);
  psi::SharedMatrix new_origins = psi::Matrix::doublet(this->origins_, rotmat, false, false);
  this->centres_->copy(new_centres);
  this->origins_->copy(new_origins);

  double** rr = rotmat->pointer();
  const double RXX = rr[0][0];
  const double RYX = rr[1][0];
  const double RZX = rr[2][0];
  const double RXY = rr[0][1];
  const double RYY = rr[1][1];
  const double RZY = rr[2][1];
  const double RXZ = rr[0][2];
  const double RYZ = rr[1][2];
  const double RZZ = rr[2][2];

  for (int n=0; n<this->nDMTPs_; ++n) {
      
       // Dipoles 
       psi::SharedMatrix new_dipoles = psi::Matrix::doublet(this->dipoles_[n], rotmat, false, false);
       this->dipoles_[n]->copy(new_dipoles);

       for (int i=0; i<this->nSites_; ++i) {
            // Quadrupoles                                                          
            double** Q = this->quadrupoles_[n]->pointer();
            double QXX = Q[i][0]; 
            double QXY = Q[i][1];
            double QXZ = Q[i][2];
            double QYY = Q[i][3];
            double QYZ = Q[i][4];
            double QZZ = Q[i][5];
                                                                                    
            double rQXX = RXX * RXX * QXX +  
                          RXX * RYX * QXY +
                          RXX * RZX * QXZ +
                          RYX * RXX * QXY +
                          RYX * RYX * QYY +
                          RYX * RZX * QYZ +
                          RZX * RXX * QXZ +
                          RZX * RYX * QYZ +
                          RZX * RZX * QZZ;
            double rQYY = RXY * RXY * QXX +
                          RXY * RYY * QXY +
                          RXY * RZY * QXZ +
                          RYY * RXY * QXY +
                          RYY * RYY * QYY +
                          RYY * RZY * QYZ +
                          RZY * RXY * QXZ +
                          RZY * RYY * QYZ +
                          RZY * RZY * QZZ;
            double rQZZ = RXZ * RXZ * QXX +
                          RXZ * RYZ * QXY +
                          RXZ * RZZ * QXZ +
                          RYZ * RXZ * QXY +
                          RYZ * RYZ * QYY +
                          RYZ * RZZ * QYZ +
                          RZZ * RXZ * QXZ +
                          RZZ * RYZ * QYZ +
                          RZZ * RZZ * QZZ;
            double rQXY = RXX * RXY * QXX +
                          RXX * RYY * QXY +
                          RXX * RZY * QXZ +
                          RYX * RXY * QXY +
                          RYX * RYY * QYY +
                          RYX * RZY * QYZ +
                          RZX * RXY * QXZ +
                          RZX * RYY * QYZ +
                          RZX * RZY * QZZ;
            double rQXZ = RXX * RXZ * QXX +
                          RXX * RYZ * QXY +
                          RXX * RZZ * QXZ +
                          RYX * RXZ * QXY +
                          RYX * RYZ * QYY +
                          RYX * RZZ * QYZ +
                          RZX * RXZ * QXZ +
                          RZX * RYZ * QYZ +
                          RZX * RZZ * QZZ;
            double rQYZ = RXY * RXZ * QXX + 
                          RXY * RYZ * QXY + 
                          RXY * RZZ * QXZ +
                          RYY * RXZ * QXY +
                          RYY * RYZ * QYY +
                          RYY * RZZ * QYZ +
                          RZY * RXZ * QXZ +
                          RZY * RYZ * QYZ +
                          RZY * RZZ * QZZ;

            Q[i][0] = rQXX; 
            Q[i][1] = rQXY;
            Q[i][2] = rQXZ;
            Q[i][3] = rQYY;
            Q[i][4] = rQYZ;
            Q[i][5] = rQZZ;

            // Octupoles
            double** O = this->octupoles_[n]->pointer();
            double OXXX = O[i][0];
            double OXXY = O[i][1]; 
            double OXXZ = O[i][2]; 
            double OXYY = O[i][3]; 
            double OXYZ = O[i][4]; 
            double OXZZ = O[i][5]; 
            double OYYY = O[i][6]; 
            double OYYZ = O[i][7]; 
            double OYZZ = O[i][8]; 
            double OZZZ = O[i][9]; 

            double rOXXX = RXX * RXX * RXX * OXXX + 
                           RXX * RXX * RYX * OXXY +
                           RXX * RXX * RZX * OXXZ +
                           RXX * RYX * RXX * OXXY +
                           RXX * RYX * RYX * OXYY +
                           RXX * RYX * RZX * OXYZ +
                           RXX * RZX * RXX * OXXZ +
                           RXX * RZX * RYX * OXYZ +
                           RXX * RZX * RZX * OXZZ +
                           RYX * RXX * RXX * OXXY +
                           RYX * RXX * RYX * OXYY +
                           RYX * RXX * RZX * OXYZ +
                           RYX * RYX * RXX * OXYY +
                           RYX * RYX * RYX * OYYY +
                           RYX * RYX * RZX * OYYZ +
                           RYX * RZX * RXX * OXYZ +
                           RYX * RZX * RYX * OYYZ +
                           RYX * RZX * RZX * OYZZ +
                           RZX * RXX * RXX * OXXZ +
                           RZX * RXX * RYX * OXYZ +
                           RZX * RXX * RZX * OXZZ +
                           RZX * RYX * RXX * OXYZ +
                           RZX * RYX * RYX * OYYZ +
                           RZX * RYX * RZX * OYZZ +
                           RZX * RZX * RXX * OXZZ +
                           RZX * RZX * RYX * OYZZ +
                           RZX * RZX * RZX * OZZZ;
            double rOYYY = RXY * RXY * RXY * OXXX +
                           RXY * RXY * RYY * OXXY +
                           RXY * RXY * RZY * OXXZ +
                           RXY * RYY * RXY * OXXY +
                           RXY * RYY * RYY * OXYY +
                           RXY * RYY * RZY * OXYZ +
                           RXY * RZY * RXY * OXXZ +
                           RXY * RZY * RYY * OXYZ +
                           RXY * RZY * RZY * OXZZ +
                           RYY * RXY * RXY * OXXY +
                           RYY * RXY * RYY * OXYY +
                           RYY * RXY * RZY * OXYZ +
                           RYY * RYY * RXY * OXYY +
                           RYY * RYY * RYY * OYYY +
                           RYY * RYY * RZY * OYYZ +
                           RYY * RZY * RXY * OXYZ +
                           RYY * RZY * RYY * OYYZ +
                           RYY * RZY * RZY * OYZZ +
                           RZY * RXY * RXY * OXXZ +
                           RZY * RXY * RYY * OXYZ +
                           RZY * RXY * RZY * OXZZ +
                           RZY * RYY * RXY * OXYZ +
                           RZY * RYY * RYY * OYYZ +
                           RZY * RYY * RZY * OYZZ +
                           RZY * RZY * RXY * OXZZ +
                           RZY * RZY * RYY * OYZZ +
                           RZY * RZY * RZY * OZZZ;
            double rOZZZ = RXZ * RXZ * RXZ * OXXX +
                           RXZ * RXZ * RYZ * OXXY +
                           RXZ * RXZ * RZZ * OXXZ +
                           RXZ * RYZ * RXZ * OXXY +
                           RXZ * RYZ * RYZ * OXYY +
                           RXZ * RYZ * RZZ * OXYZ +
                           RXZ * RZZ * RXZ * OXXZ +
                           RXZ * RZZ * RYZ * OXYZ +
                           RXZ * RZZ * RZZ * OXZZ +
                           RYZ * RXZ * RXZ * OXXY +
                           RYZ * RXZ * RYZ * OXYY +
                           RYZ * RXZ * RZZ * OXYZ +
                           RYZ * RYZ * RXZ * OXYY +
                           RYZ * RYZ * RYZ * OYYY +
                           RYZ * RYZ * RZZ * OYYZ +
                           RYZ * RZZ * RXZ * OXYZ +
                           RYZ * RZZ * RYZ * OYYZ +
                           RYZ * RZZ * RZZ * OYZZ +
                           RZZ * RXZ * RXZ * OXXZ +
                           RZZ * RXZ * RYZ * OXYZ +
                           RZZ * RXZ * RZZ * OXZZ +
                           RZZ * RYZ * RXZ * OXYZ +
                           RZZ * RYZ * RYZ * OYYZ +
                           RZZ * RYZ * RZZ * OYZZ +
                           RZZ * RZZ * RXZ * OXZZ +
                           RZZ * RZZ * RYZ * OYZZ +
                           RZZ * RZZ * RZZ * OZZZ;
            double rOXXY = RXX * RXX * RXY * OXXX +
                           RXX * RXX * RYY * OXXY +
                           RXX * RXX * RZY * OXXZ +
                           RXX * RYX * RXY * OXXY +
                           RXX * RYX * RYY * OXYY +
                           RXX * RYX * RZY * OXYZ +
                           RXX * RZX * RXY * OXXZ +
                           RXX * RZX * RYY * OXYZ +
                           RXX * RZX * RZY * OXZZ +
                           RYX * RXX * RXY * OXXY +
                           RYX * RXX * RYY * OXYY +
                           RYX * RXX * RZY * OXYZ +
                           RYX * RYX * RXY * OXYY +
                           RYX * RYX * RYY * OYYY +
                           RYX * RYX * RZY * OYYZ +
                           RYX * RZX * RXY * OXYZ +
                           RYX * RZX * RYY * OYYZ +
                           RYX * RZX * RZY * OYZZ +
                           RZX * RXX * RXY * OXXZ +
                           RZX * RXX * RYY * OXYZ +
                           RZX * RXX * RZY * OXZZ +
                           RZX * RYX * RXY * OXYZ +
                           RZX * RYX * RYY * OYYZ +
                           RZX * RYX * RZY * OYZZ +
                           RZX * RZX * RXY * OXZZ +
                           RZX * RZX * RYY * OYZZ +
                           RZX * RZX * RZY * OZZZ;
            double rOXXZ = RXX * RXX * RXZ * OXXX +
                           RXX * RXX * RYZ * OXXY +
                           RXX * RXX * RZZ * OXXZ +
                           RXX * RYX * RXZ * OXXY +
                           RXX * RYX * RYZ * OXYY +
                           RXX * RYX * RZZ * OXYZ +
                           RXX * RZX * RXZ * OXXZ +
                           RXX * RZX * RYZ * OXYZ +
                           RXX * RZX * RZZ * OXZZ +
                           RYX * RXX * RXZ * OXXY +
                           RYX * RXX * RYZ * OXYY +
                           RYX * RXX * RZZ * OXYZ +
                           RYX * RYX * RXZ * OXYY +
                           RYX * RYX * RYZ * OYYY +
                           RYX * RYX * RZZ * OYYZ +
                           RYX * RZX * RXZ * OXYZ +
                           RYX * RZX * RYZ * OYYZ +
                           RYX * RZX * RZZ * OYZZ +
                           RZX * RXX * RXZ * OXXZ +
                           RZX * RXX * RYZ * OXYZ +
                           RZX * RXX * RZZ * OXZZ +
                           RZX * RYX * RXZ * OXYZ +
                           RZX * RYX * RYZ * OYYZ +
                           RZX * RYX * RZZ * OYZZ +
                           RZX * RZX * RXZ * OXZZ +
                           RZX * RZX * RYZ * OYZZ +
                           RZX * RZX * RZZ * OZZZ;
            double rOXYY = RXX * RXY * RXY * OXXX +
                           RXX * RXY * RYY * OXXY +
                           RXX * RXY * RZY * OXXZ +
                           RXX * RYY * RXY * OXXY +
                           RXX * RYY * RYY * OXYY +
                           RXX * RYY * RZY * OXYZ +
                           RXX * RZY * RXY * OXXZ +
                           RXX * RZY * RYY * OXYZ +
                           RXX * RZY * RZY * OXZZ +
                           RYX * RXY * RXY * OXXY +
                           RYX * RXY * RYY * OXYY +
                           RYX * RXY * RZY * OXYZ +
                           RYX * RYY * RXY * OXYY +
                           RYX * RYY * RYY * OYYY +
                           RYX * RYY * RZY * OYYZ +
                           RYX * RZY * RXY * OXYZ +
                           RYX * RZY * RYY * OYYZ +
                           RYX * RZY * RZY * OYZZ +
                           RZX * RXY * RXY * OXXZ +
                           RZX * RXY * RYY * OXYZ +
                           RZX * RXY * RZY * OXZZ +
                           RZX * RYY * RXY * OXYZ +
                           RZX * RYY * RYY * OYYZ +
                           RZX * RYY * RZY * OYZZ +
                           RZX * RZY * RXY * OXZZ +
                           RZX * RZY * RYY * OYZZ +
                           RZX * RZY * RZY * OZZZ;
            double rOYYZ = RXY * RXY * RXZ * OXXX +
                           RXY * RXY * RYZ * OXXY +
                           RXY * RXY * RZZ * OXXZ +
                           RXY * RYY * RXZ * OXXY +
                           RXY * RYY * RYZ * OXYY +
                           RXY * RYY * RZZ * OXYZ +
                           RXY * RZY * RXZ * OXXZ +
                           RXY * RZY * RYZ * OXYZ +
                           RXY * RZY * RZZ * OXZZ +
                           RYY * RXY * RXZ * OXXY +
                           RYY * RXY * RYZ * OXYY +
                           RYY * RXY * RZZ * OXYZ +
                           RYY * RYY * RXZ * OXYY +
                           RYY * RYY * RYZ * OYYY +
                           RYY * RYY * RZZ * OYYZ +
                           RYY * RZY * RXZ * OXYZ +
                           RYY * RZY * RYZ * OYYZ +
                           RYY * RZY * RZZ * OYZZ +
                           RZY * RXY * RXZ * OXXZ +
                           RZY * RXY * RYZ * OXYZ +
                           RZY * RXY * RZZ * OXZZ +
                           RZY * RYY * RXZ * OXYZ +
                           RZY * RYY * RYZ * OYYZ +
                           RZY * RYY * RZZ * OYZZ +
                           RZY * RZY * RXZ * OXZZ +
                           RZY * RZY * RYZ * OYZZ +
                           RZY * RZY * RZZ * OZZZ;
            double rOXZZ = RXX * RXZ * RXZ * OXXX +
                           RXX * RXZ * RYZ * OXXY +
                           RXX * RXZ * RZZ * OXXZ +
                           RXX * RYZ * RXZ * OXXY +
                           RXX * RYZ * RYZ * OXYY +
                           RXX * RYZ * RZZ * OXYZ +
                           RXX * RZZ * RXZ * OXXZ +
                           RXX * RZZ * RYZ * OXYZ +
                           RXX * RZZ * RZZ * OXZZ +
                           RYX * RXZ * RXZ * OXXY +
                           RYX * RXZ * RYZ * OXYY +
                           RYX * RXZ * RZZ * OXYZ +
                           RYX * RYZ * RXZ * OXYY +
                           RYX * RYZ * RYZ * OYYY +
                           RYX * RYZ * RZZ * OYYZ +
                           RYX * RZZ * RXZ * OXYZ +
                           RYX * RZZ * RYZ * OYYZ +
                           RYX * RZZ * RZZ * OYZZ +
                           RZX * RXZ * RXZ * OXXZ +
                           RZX * RXZ * RYZ * OXYZ +
                           RZX * RXZ * RZZ * OXZZ +
                           RZX * RYZ * RXZ * OXYZ +
                           RZX * RYZ * RYZ * OYYZ +
                           RZX * RYZ * RZZ * OYZZ +
                           RZX * RZZ * RXZ * OXZZ +
                           RZX * RZZ * RYZ * OYZZ +
                           RZX * RZZ * RZZ * OZZZ;
            double rOYZZ = RXY * RXZ * RXZ * OXXX +
                           RXY * RXZ * RYZ * OXXY +
                           RXY * RXZ * RZZ * OXXZ +
                           RXY * RYZ * RXZ * OXXY +
                           RXY * RYZ * RYZ * OXYY +
                           RXY * RYZ * RZZ * OXYZ +
                           RXY * RZZ * RXZ * OXXZ +
                           RXY * RZZ * RYZ * OXYZ +
                           RXY * RZZ * RZZ * OXZZ +
                           RYY * RXZ * RXZ * OXXY +
                           RYY * RXZ * RYZ * OXYY +
                           RYY * RXZ * RZZ * OXYZ +
                           RYY * RYZ * RXZ * OXYY +
                           RYY * RYZ * RYZ * OYYY +
                           RYY * RYZ * RZZ * OYYZ +
                           RYY * RZZ * RXZ * OXYZ +
                           RYY * RZZ * RYZ * OYYZ +
                           RYY * RZZ * RZZ * OYZZ +
                           RZY * RXZ * RXZ * OXXZ +
                           RZY * RXZ * RYZ * OXYZ +
                           RZY * RXZ * RZZ * OXZZ +
                           RZY * RYZ * RXZ * OXYZ +
                           RZY * RYZ * RYZ * OYYZ +
                           RZY * RYZ * RZZ * OYZZ +
                           RZY * RZZ * RXZ * OXZZ +
                           RZY * RZZ * RYZ * OYZZ +
                           RZY * RZZ * RZZ * OZZZ;
            double rOXYZ = RXX * RXY * RXZ * OXXX +
                           RXX * RXY * RYZ * OXXY +
                           RXX * RXY * RZZ * OXXZ +
                           RXX * RYY * RXZ * OXXY +
                           RXX * RYY * RYZ * OXYY +
                           RXX * RYY * RZZ * OXYZ +
                           RXX * RZY * RXZ * OXXZ +
                           RXX * RZY * RYZ * OXYZ +
                           RXX * RZY * RZZ * OXZZ +
                           RYX * RXY * RXZ * OXXY +
                           RYX * RXY * RYZ * OXYY +
                           RYX * RXY * RZZ * OXYZ +
                           RYX * RYY * RXZ * OXYY +
                           RYX * RYY * RYZ * OYYY +
                           RYX * RYY * RZZ * OYYZ +
                           RYX * RZY * RXZ * OXYZ +
                           RYX * RZY * RYZ * OYYZ +
                           RYX * RZY * RZZ * OYZZ +
                           RZX * RXY * RXZ * OXXZ +
                           RZX * RXY * RYZ * OXYZ +
                           RZX * RXY * RZZ * OXZZ +
                           RZX * RYY * RXZ * OXYZ +
                           RZX * RYY * RYZ * OYYZ +
                           RZX * RYY * RZZ * OYZZ +
                           RZX * RZY * RXZ * OXZZ +
                           RZX * RZY * RYZ * OYZZ +
                           RZX * RZY * RZZ * OZZZ;
        
            O[i][0] = rOXXX;
            O[i][1] = rOXXY; 
            O[i][2] = rOXXZ; 
            O[i][3] = rOXYY; 
            O[i][4] = rOXYZ; 
            O[i][5] = rOXZZ; 
            O[i][6] = rOYYY; 
            O[i][7] = rOYYZ; 
            O[i][8] = rOYZZ; 
            O[i][9] = rOZZZ; 

            // Hexadecapoles
            double** H = this->hexadecapoles_[n]->pointer();

            double HXXXX = H[i][ 0]; 
            double HXXXY = H[i][ 1];
            double HXXXZ = H[i][ 2];
            double HXXYY = H[i][ 3];
            double HXXYZ = H[i][ 4];
            double HXXZZ = H[i][ 5];
            double HXYYY = H[i][ 6];
            double HXYYZ = H[i][ 7];
            double HXYZZ = H[i][ 8];
            double HXZZZ = H[i][ 9];
            double HYYYY = H[i][10];
            double HYYYZ = H[i][11];
            double HYYZZ = H[i][12];
            double HYZZZ = H[i][13];
            double HZZZZ = H[i][14];

            double rHXXXX = 
                            RXX * RXX * RXX * RXX * HXXXX +
                            RXX * RXX * RXX * RYX * HXXXY +
                            RXX * RXX * RXX * RZX * HXXXZ +
                            RXX * RXX * RYX * RXX * HXXXY +
                            RXX * RXX * RYX * RYX * HXXYY +
                            RXX * RXX * RYX * RZX * HXXYZ +
                            RXX * RXX * RZX * RXX * HXXXZ +
                            RXX * RXX * RZX * RYX * HXXYZ +
                            RXX * RXX * RZX * RZX * HXXZZ +
                            RXX * RYX * RXX * RXX * HXXXY +
                            RXX * RYX * RXX * RYX * HXXYY +
                            RXX * RYX * RXX * RZX * HXXYZ +
                            RXX * RYX * RYX * RXX * HXXYY +
                            RXX * RYX * RYX * RYX * HXYYY +
                            RXX * RYX * RYX * RZX * HXYYZ +
                            RXX * RYX * RZX * RXX * HXXYZ +
                            RXX * RYX * RZX * RYX * HXYYZ +
                            RXX * RYX * RZX * RZX * HXYZZ +
                            RXX * RZX * RXX * RXX * HXXXZ +
                            RXX * RZX * RXX * RYX * HXXYZ +
                            RXX * RZX * RXX * RZX * HXXZZ +
                            RXX * RZX * RYX * RXX * HXXYZ +
                            RXX * RZX * RYX * RYX * HXYYZ +
                            RXX * RZX * RYX * RZX * HXYZZ +
                            RXX * RZX * RZX * RXX * HXXZZ +
                            RXX * RZX * RZX * RYX * HXYZZ +
                            RXX * RZX * RZX * RZX * HXZZZ +
                            RYX * RXX * RXX * RXX * HXXXY +
                            RYX * RXX * RXX * RYX * HXXYY +
                            RYX * RXX * RXX * RZX * HXXYZ +
                            RYX * RXX * RYX * RXX * HXXYY +
                            RYX * RXX * RYX * RYX * HXYYY +
                            RYX * RXX * RYX * RZX * HXYYZ +
                            RYX * RXX * RZX * RXX * HXXYZ +
                            RYX * RXX * RZX * RYX * HXYYZ +
                            RYX * RXX * RZX * RZX * HXYZZ +
                            RYX * RYX * RXX * RXX * HXXYY +
                            RYX * RYX * RXX * RYX * HXYYY +
                            RYX * RYX * RXX * RZX * HXYYZ +
                            RYX * RYX * RYX * RXX * HXYYY +
                            RYX * RYX * RYX * RYX * HYYYY +
                            RYX * RYX * RYX * RZX * HYYYZ +
                            RYX * RYX * RZX * RXX * HXYYZ +
                            RYX * RYX * RZX * RYX * HYYYZ +
                            RYX * RYX * RZX * RZX * HYYZZ +
                            RYX * RZX * RXX * RXX * HXXYZ +
                            RYX * RZX * RXX * RYX * HXYYZ +
                            RYX * RZX * RXX * RZX * HXYZZ +
                            RYX * RZX * RYX * RXX * HXYYZ +
                            RYX * RZX * RYX * RYX * HYYYZ +
                            RYX * RZX * RYX * RZX * HYYZZ +
                            RYX * RZX * RZX * RXX * HXYZZ +
                            RYX * RZX * RZX * RYX * HYYZZ +
                            RYX * RZX * RZX * RZX * HYZZZ +
                            RZX * RXX * RXX * RXX * HXXXZ +
                            RZX * RXX * RXX * RYX * HXXYZ +
                            RZX * RXX * RXX * RZX * HXXZZ +
                            RZX * RXX * RYX * RXX * HXXYZ +
                            RZX * RXX * RYX * RYX * HXYYZ +
                            RZX * RXX * RYX * RZX * HXYZZ +
                            RZX * RXX * RZX * RXX * HXXZZ +
                            RZX * RXX * RZX * RYX * HXYZZ +
                            RZX * RXX * RZX * RZX * HXZZZ +
                            RZX * RYX * RXX * RXX * HXXYZ +
                            RZX * RYX * RXX * RYX * HXYYZ +
                            RZX * RYX * RXX * RZX * HXYZZ +
                            RZX * RYX * RYX * RXX * HXYYZ +
                            RZX * RYX * RYX * RYX * HYYYZ +
                            RZX * RYX * RYX * RZX * HYYZZ +
                            RZX * RYX * RZX * RXX * HXYZZ +
                            RZX * RYX * RZX * RYX * HYYZZ +
                            RZX * RYX * RZX * RZX * HYZZZ +
                            RZX * RZX * RXX * RXX * HXXZZ +
                            RZX * RZX * RXX * RYX * HXYZZ +
                            RZX * RZX * RXX * RZX * HXZZZ +
                            RZX * RZX * RYX * RXX * HXYZZ +
                            RZX * RZX * RYX * RYX * HYYZZ +
                            RZX * RZX * RYX * RZX * HYZZZ +
                            RZX * RZX * RZX * RXX * HXZZZ +
                            RZX * RZX * RZX * RYX * HYZZZ +
                            RZX * RZX * RZX * RZX * HZZZZ ;
            double rHXXXY = 
                            RXX * RXX * RXX * RXY * HXXXX +
                            RXX * RXX * RXX * RYY * HXXXY +
                            RXX * RXX * RXX * RZY * HXXXZ +
                            RXX * RXX * RYX * RXY * HXXXY +
                            RXX * RXX * RYX * RYY * HXXYY +
                            RXX * RXX * RYX * RZY * HXXYZ +
                            RXX * RXX * RZX * RXY * HXXXZ +
                            RXX * RXX * RZX * RYY * HXXYZ +
                            RXX * RXX * RZX * RZY * HXXZZ +
                            RXX * RYX * RXX * RXY * HXXXY +
                            RXX * RYX * RXX * RYY * HXXYY +
                            RXX * RYX * RXX * RZY * HXXYZ +
                            RXX * RYX * RYX * RXY * HXXYY +
                            RXX * RYX * RYX * RYY * HXYYY +
                            RXX * RYX * RYX * RZY * HXYYZ +
                            RXX * RYX * RZX * RXY * HXXYZ +
                            RXX * RYX * RZX * RYY * HXYYZ +
                            RXX * RYX * RZX * RZY * HXYZZ +
                            RXX * RZX * RXX * RXY * HXXXZ +
                            RXX * RZX * RXX * RYY * HXXYZ +
                            RXX * RZX * RXX * RZY * HXXZZ +
                            RXX * RZX * RYX * RXY * HXXYZ +
                            RXX * RZX * RYX * RYY * HXYYZ +
                            RXX * RZX * RYX * RZY * HXYZZ +
                            RXX * RZX * RZX * RXY * HXXZZ +
                            RXX * RZX * RZX * RYY * HXYZZ +
                            RXX * RZX * RZX * RZY * HXZZZ +
                            RYX * RXX * RXX * RXY * HXXXY +
                            RYX * RXX * RXX * RYY * HXXYY +
                            RYX * RXX * RXX * RZY * HXXYZ +
                            RYX * RXX * RYX * RXY * HXXYY +
                            RYX * RXX * RYX * RYY * HXYYY +
                            RYX * RXX * RYX * RZY * HXYYZ +
                            RYX * RXX * RZX * RXY * HXXYZ +
                            RYX * RXX * RZX * RYY * HXYYZ +
                            RYX * RXX * RZX * RZY * HXYZZ +
                            RYX * RYX * RXX * RXY * HXXYY +
                            RYX * RYX * RXX * RYY * HXYYY +
                            RYX * RYX * RXX * RZY * HXYYZ +
                            RYX * RYX * RYX * RXY * HXYYY +
                            RYX * RYX * RYX * RYY * HYYYY +
                            RYX * RYX * RYX * RZY * HYYYZ +
                            RYX * RYX * RZX * RXY * HXYYZ +
                            RYX * RYX * RZX * RYY * HYYYZ +
                            RYX * RYX * RZX * RZY * HYYZZ +
                            RYX * RZX * RXX * RXY * HXXYZ +
                            RYX * RZX * RXX * RYY * HXYYZ +
                            RYX * RZX * RXX * RZY * HXYZZ +
                            RYX * RZX * RYX * RXY * HXYYZ +
                            RYX * RZX * RYX * RYY * HYYYZ +
                            RYX * RZX * RYX * RZY * HYYZZ +
                            RYX * RZX * RZX * RXY * HXYZZ +
                            RYX * RZX * RZX * RYY * HYYZZ +
                            RYX * RZX * RZX * RZY * HYZZZ +
                            RZX * RXX * RXX * RXY * HXXXZ +
                            RZX * RXX * RXX * RYY * HXXYZ +
                            RZX * RXX * RXX * RZY * HXXZZ +
                            RZX * RXX * RYX * RXY * HXXYZ +
                            RZX * RXX * RYX * RYY * HXYYZ +
                            RZX * RXX * RYX * RZY * HXYZZ +
                            RZX * RXX * RZX * RXY * HXXZZ +
                            RZX * RXX * RZX * RYY * HXYZZ +
                            RZX * RXX * RZX * RZY * HXZZZ +
                            RZX * RYX * RXX * RXY * HXXYZ +
                            RZX * RYX * RXX * RYY * HXYYZ +
                            RZX * RYX * RXX * RZY * HXYZZ +
                            RZX * RYX * RYX * RXY * HXYYZ +
                            RZX * RYX * RYX * RYY * HYYYZ +
                            RZX * RYX * RYX * RZY * HYYZZ +
                            RZX * RYX * RZX * RXY * HXYZZ +
                            RZX * RYX * RZX * RYY * HYYZZ +
                            RZX * RYX * RZX * RZY * HYZZZ +
                            RZX * RZX * RXX * RXY * HXXZZ +
                            RZX * RZX * RXX * RYY * HXYZZ +
                            RZX * RZX * RXX * RZY * HXZZZ +
                            RZX * RZX * RYX * RXY * HXYZZ +
                            RZX * RZX * RYX * RYY * HYYZZ +
                            RZX * RZX * RYX * RZY * HYZZZ +
                            RZX * RZX * RZX * RXY * HXZZZ +
                            RZX * RZX * RZX * RYY * HYZZZ +
                            RZX * RZX * RZX * RZY * HZZZZ ;
            double rHXXXZ = 
                            RXX * RXX * RXX * RXZ * HXXXX +
                            RXX * RXX * RXX * RYZ * HXXXY +
                            RXX * RXX * RXX * RZZ * HXXXZ +
                            RXX * RXX * RYX * RXZ * HXXXY +
                            RXX * RXX * RYX * RYZ * HXXYY +
                            RXX * RXX * RYX * RZZ * HXXYZ +
                            RXX * RXX * RZX * RXZ * HXXXZ +
                            RXX * RXX * RZX * RYZ * HXXYZ +
                            RXX * RXX * RZX * RZZ * HXXZZ +
                            RXX * RYX * RXX * RXZ * HXXXY +
                            RXX * RYX * RXX * RYZ * HXXYY +
                            RXX * RYX * RXX * RZZ * HXXYZ +
                            RXX * RYX * RYX * RXZ * HXXYY +
                            RXX * RYX * RYX * RYZ * HXYYY +
                            RXX * RYX * RYX * RZZ * HXYYZ +
                            RXX * RYX * RZX * RXZ * HXXYZ +
                            RXX * RYX * RZX * RYZ * HXYYZ +
                            RXX * RYX * RZX * RZZ * HXYZZ +
                            RXX * RZX * RXX * RXZ * HXXXZ +
                            RXX * RZX * RXX * RYZ * HXXYZ +
                            RXX * RZX * RXX * RZZ * HXXZZ +
                            RXX * RZX * RYX * RXZ * HXXYZ +
                            RXX * RZX * RYX * RYZ * HXYYZ +
                            RXX * RZX * RYX * RZZ * HXYZZ +
                            RXX * RZX * RZX * RXZ * HXXZZ +
                            RXX * RZX * RZX * RYZ * HXYZZ +
                            RXX * RZX * RZX * RZZ * HXZZZ +
                            RYX * RXX * RXX * RXZ * HXXXY +
                            RYX * RXX * RXX * RYZ * HXXYY +
                            RYX * RXX * RXX * RZZ * HXXYZ +
                            RYX * RXX * RYX * RXZ * HXXYY +
                            RYX * RXX * RYX * RYZ * HXYYY +
                            RYX * RXX * RYX * RZZ * HXYYZ +
                            RYX * RXX * RZX * RXZ * HXXYZ +
                            RYX * RXX * RZX * RYZ * HXYYZ +
                            RYX * RXX * RZX * RZZ * HXYZZ +
                            RYX * RYX * RXX * RXZ * HXXYY +
                            RYX * RYX * RXX * RYZ * HXYYY +
                            RYX * RYX * RXX * RZZ * HXYYZ +
                            RYX * RYX * RYX * RXZ * HXYYY +
                            RYX * RYX * RYX * RYZ * HYYYY +
                            RYX * RYX * RYX * RZZ * HYYYZ +
                            RYX * RYX * RZX * RXZ * HXYYZ +
                            RYX * RYX * RZX * RYZ * HYYYZ +
                            RYX * RYX * RZX * RZZ * HYYZZ +
                            RYX * RZX * RXX * RXZ * HXXYZ +
                            RYX * RZX * RXX * RYZ * HXYYZ +
                            RYX * RZX * RXX * RZZ * HXYZZ +
                            RYX * RZX * RYX * RXZ * HXYYZ +
                            RYX * RZX * RYX * RYZ * HYYYZ +
                            RYX * RZX * RYX * RZZ * HYYZZ +
                            RYX * RZX * RZX * RXZ * HXYZZ +
                            RYX * RZX * RZX * RYZ * HYYZZ +
                            RYX * RZX * RZX * RZZ * HYZZZ +
                            RZX * RXX * RXX * RXZ * HXXXZ +
                            RZX * RXX * RXX * RYZ * HXXYZ +
                            RZX * RXX * RXX * RZZ * HXXZZ +
                            RZX * RXX * RYX * RXZ * HXXYZ +
                            RZX * RXX * RYX * RYZ * HXYYZ +
                            RZX * RXX * RYX * RZZ * HXYZZ +
                            RZX * RXX * RZX * RXZ * HXXZZ +
                            RZX * RXX * RZX * RYZ * HXYZZ +
                            RZX * RXX * RZX * RZZ * HXZZZ +
                            RZX * RYX * RXX * RXZ * HXXYZ +
                            RZX * RYX * RXX * RYZ * HXYYZ +
                            RZX * RYX * RXX * RZZ * HXYZZ +
                            RZX * RYX * RYX * RXZ * HXYYZ +
                            RZX * RYX * RYX * RYZ * HYYYZ +
                            RZX * RYX * RYX * RZZ * HYYZZ +
                            RZX * RYX * RZX * RXZ * HXYZZ +
                            RZX * RYX * RZX * RYZ * HYYZZ +
                            RZX * RYX * RZX * RZZ * HYZZZ +
                            RZX * RZX * RXX * RXZ * HXXZZ +
                            RZX * RZX * RXX * RYZ * HXYZZ +
                            RZX * RZX * RXX * RZZ * HXZZZ +
                            RZX * RZX * RYX * RXZ * HXYZZ +
                            RZX * RZX * RYX * RYZ * HYYZZ +
                            RZX * RZX * RYX * RZZ * HYZZZ +
                            RZX * RZX * RZX * RXZ * HXZZZ +
                            RZX * RZX * RZX * RYZ * HYZZZ +
                            RZX * RZX * RZX * RZZ * HZZZZ ;
            double rHXXYY = 
                            RXX * RXX * RXY * RXY * HXXXX +
                            RXX * RXX * RXY * RYY * HXXXY +
                            RXX * RXX * RXY * RZY * HXXXZ +
                            RXX * RXX * RYY * RXY * HXXXY +
                            RXX * RXX * RYY * RYY * HXXYY +
                            RXX * RXX * RYY * RZY * HXXYZ +
                            RXX * RXX * RZY * RXY * HXXXZ +
                            RXX * RXX * RZY * RYY * HXXYZ +
                            RXX * RXX * RZY * RZY * HXXZZ +
                            RXX * RYX * RXY * RXY * HXXXY +
                            RXX * RYX * RXY * RYY * HXXYY +
                            RXX * RYX * RXY * RZY * HXXYZ +
                            RXX * RYX * RYY * RXY * HXXYY +
                            RXX * RYX * RYY * RYY * HXYYY +
                            RXX * RYX * RYY * RZY * HXYYZ +
                            RXX * RYX * RZY * RXY * HXXYZ +
                            RXX * RYX * RZY * RYY * HXYYZ +
                            RXX * RYX * RZY * RZY * HXYZZ +
                            RXX * RZX * RXY * RXY * HXXXZ +
                            RXX * RZX * RXY * RYY * HXXYZ +
                            RXX * RZX * RXY * RZY * HXXZZ +
                            RXX * RZX * RYY * RXY * HXXYZ +
                            RXX * RZX * RYY * RYY * HXYYZ +
                            RXX * RZX * RYY * RZY * HXYZZ +
                            RXX * RZX * RZY * RXY * HXXZZ +
                            RXX * RZX * RZY * RYY * HXYZZ +
                            RXX * RZX * RZY * RZY * HXZZZ +
                            RYX * RXX * RXY * RXY * HXXXY +
                            RYX * RXX * RXY * RYY * HXXYY +
                            RYX * RXX * RXY * RZY * HXXYZ +
                            RYX * RXX * RYY * RXY * HXXYY +
                            RYX * RXX * RYY * RYY * HXYYY +
                            RYX * RXX * RYY * RZY * HXYYZ +
                            RYX * RXX * RZY * RXY * HXXYZ +
                            RYX * RXX * RZY * RYY * HXYYZ +
                            RYX * RXX * RZY * RZY * HXYZZ +
                            RYX * RYX * RXY * RXY * HXXYY +
                            RYX * RYX * RXY * RYY * HXYYY +
                            RYX * RYX * RXY * RZY * HXYYZ +
                            RYX * RYX * RYY * RXY * HXYYY +
                            RYX * RYX * RYY * RYY * HYYYY +
                            RYX * RYX * RYY * RZY * HYYYZ +
                            RYX * RYX * RZY * RXY * HXYYZ +
                            RYX * RYX * RZY * RYY * HYYYZ +
                            RYX * RYX * RZY * RZY * HYYZZ +
                            RYX * RZX * RXY * RXY * HXXYZ +
                            RYX * RZX * RXY * RYY * HXYYZ +
                            RYX * RZX * RXY * RZY * HXYZZ +
                            RYX * RZX * RYY * RXY * HXYYZ +
                            RYX * RZX * RYY * RYY * HYYYZ +
                            RYX * RZX * RYY * RZY * HYYZZ +
                            RYX * RZX * RZY * RXY * HXYZZ +
                            RYX * RZX * RZY * RYY * HYYZZ +
                            RYX * RZX * RZY * RZY * HYZZZ +
                            RZX * RXX * RXY * RXY * HXXXZ +
                            RZX * RXX * RXY * RYY * HXXYZ +
                            RZX * RXX * RXY * RZY * HXXZZ +
                            RZX * RXX * RYY * RXY * HXXYZ +
                            RZX * RXX * RYY * RYY * HXYYZ +
                            RZX * RXX * RYY * RZY * HXYZZ +
                            RZX * RXX * RZY * RXY * HXXZZ +
                            RZX * RXX * RZY * RYY * HXYZZ +
                            RZX * RXX * RZY * RZY * HXZZZ +
                            RZX * RYX * RXY * RXY * HXXYZ +
                            RZX * RYX * RXY * RYY * HXYYZ +
                            RZX * RYX * RXY * RZY * HXYZZ +
                            RZX * RYX * RYY * RXY * HXYYZ +
                            RZX * RYX * RYY * RYY * HYYYZ +
                            RZX * RYX * RYY * RZY * HYYZZ +
                            RZX * RYX * RZY * RXY * HXYZZ +
                            RZX * RYX * RZY * RYY * HYYZZ +
                            RZX * RYX * RZY * RZY * HYZZZ +
                            RZX * RZX * RXY * RXY * HXXZZ +
                            RZX * RZX * RXY * RYY * HXYZZ +
                            RZX * RZX * RXY * RZY * HXZZZ +
                            RZX * RZX * RYY * RXY * HXYZZ +
                            RZX * RZX * RYY * RYY * HYYZZ +
                            RZX * RZX * RYY * RZY * HYZZZ +
                            RZX * RZX * RZY * RXY * HXZZZ +
                            RZX * RZX * RZY * RYY * HYZZZ +
                            RZX * RZX * RZY * RZY * HZZZZ ;
            double rHXXYZ = 
                            RXX * RXX * RXY * RXZ * HXXXX +
                            RXX * RXX * RXY * RYZ * HXXXY +
                            RXX * RXX * RXY * RZZ * HXXXZ +
                            RXX * RXX * RYY * RXZ * HXXXY +
                            RXX * RXX * RYY * RYZ * HXXYY +
                            RXX * RXX * RYY * RZZ * HXXYZ +
                            RXX * RXX * RZY * RXZ * HXXXZ +
                            RXX * RXX * RZY * RYZ * HXXYZ +
                            RXX * RXX * RZY * RZZ * HXXZZ +
                            RXX * RYX * RXY * RXZ * HXXXY +
                            RXX * RYX * RXY * RYZ * HXXYY +
                            RXX * RYX * RXY * RZZ * HXXYZ +
                            RXX * RYX * RYY * RXZ * HXXYY +
                            RXX * RYX * RYY * RYZ * HXYYY +
                            RXX * RYX * RYY * RZZ * HXYYZ +
                            RXX * RYX * RZY * RXZ * HXXYZ +
                            RXX * RYX * RZY * RYZ * HXYYZ +
                            RXX * RYX * RZY * RZZ * HXYZZ +
                            RXX * RZX * RXY * RXZ * HXXXZ +
                            RXX * RZX * RXY * RYZ * HXXYZ +
                            RXX * RZX * RXY * RZZ * HXXZZ +
                            RXX * RZX * RYY * RXZ * HXXYZ +
                            RXX * RZX * RYY * RYZ * HXYYZ +
                            RXX * RZX * RYY * RZZ * HXYZZ +
                            RXX * RZX * RZY * RXZ * HXXZZ +
                            RXX * RZX * RZY * RYZ * HXYZZ +
                            RXX * RZX * RZY * RZZ * HXZZZ +
                            RYX * RXX * RXY * RXZ * HXXXY +
                            RYX * RXX * RXY * RYZ * HXXYY +
                            RYX * RXX * RXY * RZZ * HXXYZ +
                            RYX * RXX * RYY * RXZ * HXXYY +
                            RYX * RXX * RYY * RYZ * HXYYY +
                            RYX * RXX * RYY * RZZ * HXYYZ +
                            RYX * RXX * RZY * RXZ * HXXYZ +
                            RYX * RXX * RZY * RYZ * HXYYZ +
                            RYX * RXX * RZY * RZZ * HXYZZ +
                            RYX * RYX * RXY * RXZ * HXXYY +
                            RYX * RYX * RXY * RYZ * HXYYY +
                            RYX * RYX * RXY * RZZ * HXYYZ +
                            RYX * RYX * RYY * RXZ * HXYYY +
                            RYX * RYX * RYY * RYZ * HYYYY +
                            RYX * RYX * RYY * RZZ * HYYYZ +
                            RYX * RYX * RZY * RXZ * HXYYZ +
                            RYX * RYX * RZY * RYZ * HYYYZ +
                            RYX * RYX * RZY * RZZ * HYYZZ +
                            RYX * RZX * RXY * RXZ * HXXYZ +
                            RYX * RZX * RXY * RYZ * HXYYZ +
                            RYX * RZX * RXY * RZZ * HXYZZ +
                            RYX * RZX * RYY * RXZ * HXYYZ +
                            RYX * RZX * RYY * RYZ * HYYYZ +
                            RYX * RZX * RYY * RZZ * HYYZZ +
                            RYX * RZX * RZY * RXZ * HXYZZ +
                            RYX * RZX * RZY * RYZ * HYYZZ +
                            RYX * RZX * RZY * RZZ * HYZZZ +
                            RZX * RXX * RXY * RXZ * HXXXZ +
                            RZX * RXX * RXY * RYZ * HXXYZ +
                            RZX * RXX * RXY * RZZ * HXXZZ +
                            RZX * RXX * RYY * RXZ * HXXYZ +
                            RZX * RXX * RYY * RYZ * HXYYZ +
                            RZX * RXX * RYY * RZZ * HXYZZ +
                            RZX * RXX * RZY * RXZ * HXXZZ +
                            RZX * RXX * RZY * RYZ * HXYZZ +
                            RZX * RXX * RZY * RZZ * HXZZZ +
                            RZX * RYX * RXY * RXZ * HXXYZ +
                            RZX * RYX * RXY * RYZ * HXYYZ +
                            RZX * RYX * RXY * RZZ * HXYZZ +
                            RZX * RYX * RYY * RXZ * HXYYZ +
                            RZX * RYX * RYY * RYZ * HYYYZ +
                            RZX * RYX * RYY * RZZ * HYYZZ +
                            RZX * RYX * RZY * RXZ * HXYZZ +
                            RZX * RYX * RZY * RYZ * HYYZZ +
                            RZX * RYX * RZY * RZZ * HYZZZ +
                            RZX * RZX * RXY * RXZ * HXXZZ +
                            RZX * RZX * RXY * RYZ * HXYZZ +
                            RZX * RZX * RXY * RZZ * HXZZZ +
                            RZX * RZX * RYY * RXZ * HXYZZ +
                            RZX * RZX * RYY * RYZ * HYYZZ +
                            RZX * RZX * RYY * RZZ * HYZZZ +
                            RZX * RZX * RZY * RXZ * HXZZZ +
                            RZX * RZX * RZY * RYZ * HYZZZ +
                            RZX * RZX * RZY * RZZ * HZZZZ ;
            double rHXXZZ = 
                            RXX * RXX * RXZ * RXZ * HXXXX +
                            RXX * RXX * RXZ * RYZ * HXXXY +
                            RXX * RXX * RXZ * RZZ * HXXXZ +
                            RXX * RXX * RYZ * RXZ * HXXXY +
                            RXX * RXX * RYZ * RYZ * HXXYY +
                            RXX * RXX * RYZ * RZZ * HXXYZ +
                            RXX * RXX * RZZ * RXZ * HXXXZ +
                            RXX * RXX * RZZ * RYZ * HXXYZ +
                            RXX * RXX * RZZ * RZZ * HXXZZ +
                            RXX * RYX * RXZ * RXZ * HXXXY +
                            RXX * RYX * RXZ * RYZ * HXXYY +
                            RXX * RYX * RXZ * RZZ * HXXYZ +
                            RXX * RYX * RYZ * RXZ * HXXYY +
                            RXX * RYX * RYZ * RYZ * HXYYY +
                            RXX * RYX * RYZ * RZZ * HXYYZ +
                            RXX * RYX * RZZ * RXZ * HXXYZ +
                            RXX * RYX * RZZ * RYZ * HXYYZ +
                            RXX * RYX * RZZ * RZZ * HXYZZ +
                            RXX * RZX * RXZ * RXZ * HXXXZ +
                            RXX * RZX * RXZ * RYZ * HXXYZ +
                            RXX * RZX * RXZ * RZZ * HXXZZ +
                            RXX * RZX * RYZ * RXZ * HXXYZ +
                            RXX * RZX * RYZ * RYZ * HXYYZ +
                            RXX * RZX * RYZ * RZZ * HXYZZ +
                            RXX * RZX * RZZ * RXZ * HXXZZ +
                            RXX * RZX * RZZ * RYZ * HXYZZ +
                            RXX * RZX * RZZ * RZZ * HXZZZ +
                            RYX * RXX * RXZ * RXZ * HXXXY +
                            RYX * RXX * RXZ * RYZ * HXXYY +
                            RYX * RXX * RXZ * RZZ * HXXYZ +
                            RYX * RXX * RYZ * RXZ * HXXYY +
                            RYX * RXX * RYZ * RYZ * HXYYY +
                            RYX * RXX * RYZ * RZZ * HXYYZ +
                            RYX * RXX * RZZ * RXZ * HXXYZ +
                            RYX * RXX * RZZ * RYZ * HXYYZ +
                            RYX * RXX * RZZ * RZZ * HXYZZ +
                            RYX * RYX * RXZ * RXZ * HXXYY +
                            RYX * RYX * RXZ * RYZ * HXYYY +
                            RYX * RYX * RXZ * RZZ * HXYYZ +
                            RYX * RYX * RYZ * RXZ * HXYYY +
                            RYX * RYX * RYZ * RYZ * HYYYY +
                            RYX * RYX * RYZ * RZZ * HYYYZ +
                            RYX * RYX * RZZ * RXZ * HXYYZ +
                            RYX * RYX * RZZ * RYZ * HYYYZ +
                            RYX * RYX * RZZ * RZZ * HYYZZ +
                            RYX * RZX * RXZ * RXZ * HXXYZ +
                            RYX * RZX * RXZ * RYZ * HXYYZ +
                            RYX * RZX * RXZ * RZZ * HXYZZ +
                            RYX * RZX * RYZ * RXZ * HXYYZ +
                            RYX * RZX * RYZ * RYZ * HYYYZ +
                            RYX * RZX * RYZ * RZZ * HYYZZ +
                            RYX * RZX * RZZ * RXZ * HXYZZ +
                            RYX * RZX * RZZ * RYZ * HYYZZ +
                            RYX * RZX * RZZ * RZZ * HYZZZ +
                            RZX * RXX * RXZ * RXZ * HXXXZ +
                            RZX * RXX * RXZ * RYZ * HXXYZ +
                            RZX * RXX * RXZ * RZZ * HXXZZ +
                            RZX * RXX * RYZ * RXZ * HXXYZ +
                            RZX * RXX * RYZ * RYZ * HXYYZ +
                            RZX * RXX * RYZ * RZZ * HXYZZ +
                            RZX * RXX * RZZ * RXZ * HXXZZ +
                            RZX * RXX * RZZ * RYZ * HXYZZ +
                            RZX * RXX * RZZ * RZZ * HXZZZ +
                            RZX * RYX * RXZ * RXZ * HXXYZ +
                            RZX * RYX * RXZ * RYZ * HXYYZ +
                            RZX * RYX * RXZ * RZZ * HXYZZ +
                            RZX * RYX * RYZ * RXZ * HXYYZ +
                            RZX * RYX * RYZ * RYZ * HYYYZ +
                            RZX * RYX * RYZ * RZZ * HYYZZ +
                            RZX * RYX * RZZ * RXZ * HXYZZ +
                            RZX * RYX * RZZ * RYZ * HYYZZ +
                            RZX * RYX * RZZ * RZZ * HYZZZ +
                            RZX * RZX * RXZ * RXZ * HXXZZ +
                            RZX * RZX * RXZ * RYZ * HXYZZ +
                            RZX * RZX * RXZ * RZZ * HXZZZ +
                            RZX * RZX * RYZ * RXZ * HXYZZ +
                            RZX * RZX * RYZ * RYZ * HYYZZ +
                            RZX * RZX * RYZ * RZZ * HYZZZ +
                            RZX * RZX * RZZ * RXZ * HXZZZ +
                            RZX * RZX * RZZ * RYZ * HYZZZ +
                            RZX * RZX * RZZ * RZZ * HZZZZ ;
            double rHXYYY = 
                            RXX * RXY * RXY * RXY * HXXXX +
                            RXX * RXY * RXY * RYY * HXXXY +
                            RXX * RXY * RXY * RZY * HXXXZ +
                            RXX * RXY * RYY * RXY * HXXXY +
                            RXX * RXY * RYY * RYY * HXXYY +
                            RXX * RXY * RYY * RZY * HXXYZ +
                            RXX * RXY * RZY * RXY * HXXXZ +
                            RXX * RXY * RZY * RYY * HXXYZ +
                            RXX * RXY * RZY * RZY * HXXZZ +
                            RXX * RYY * RXY * RXY * HXXXY +
                            RXX * RYY * RXY * RYY * HXXYY +
                            RXX * RYY * RXY * RZY * HXXYZ +
                            RXX * RYY * RYY * RXY * HXXYY +
                            RXX * RYY * RYY * RYY * HXYYY +
                            RXX * RYY * RYY * RZY * HXYYZ +
                            RXX * RYY * RZY * RXY * HXXYZ +
                            RXX * RYY * RZY * RYY * HXYYZ +
                            RXX * RYY * RZY * RZY * HXYZZ +
                            RXX * RZY * RXY * RXY * HXXXZ +
                            RXX * RZY * RXY * RYY * HXXYZ +
                            RXX * RZY * RXY * RZY * HXXZZ +
                            RXX * RZY * RYY * RXY * HXXYZ +
                            RXX * RZY * RYY * RYY * HXYYZ +
                            RXX * RZY * RYY * RZY * HXYZZ +
                            RXX * RZY * RZY * RXY * HXXZZ +
                            RXX * RZY * RZY * RYY * HXYZZ +
                            RXX * RZY * RZY * RZY * HXZZZ +
                            RYX * RXY * RXY * RXY * HXXXY +
                            RYX * RXY * RXY * RYY * HXXYY +
                            RYX * RXY * RXY * RZY * HXXYZ +
                            RYX * RXY * RYY * RXY * HXXYY +
                            RYX * RXY * RYY * RYY * HXYYY +
                            RYX * RXY * RYY * RZY * HXYYZ +
                            RYX * RXY * RZY * RXY * HXXYZ +
                            RYX * RXY * RZY * RYY * HXYYZ +
                            RYX * RXY * RZY * RZY * HXYZZ +
                            RYX * RYY * RXY * RXY * HXXYY +
                            RYX * RYY * RXY * RYY * HXYYY +
                            RYX * RYY * RXY * RZY * HXYYZ +
                            RYX * RYY * RYY * RXY * HXYYY +
                            RYX * RYY * RYY * RYY * HYYYY +
                            RYX * RYY * RYY * RZY * HYYYZ +
                            RYX * RYY * RZY * RXY * HXYYZ +
                            RYX * RYY * RZY * RYY * HYYYZ +
                            RYX * RYY * RZY * RZY * HYYZZ +
                            RYX * RZY * RXY * RXY * HXXYZ +
                            RYX * RZY * RXY * RYY * HXYYZ +
                            RYX * RZY * RXY * RZY * HXYZZ +
                            RYX * RZY * RYY * RXY * HXYYZ +
                            RYX * RZY * RYY * RYY * HYYYZ +
                            RYX * RZY * RYY * RZY * HYYZZ +
                            RYX * RZY * RZY * RXY * HXYZZ +
                            RYX * RZY * RZY * RYY * HYYZZ +
                            RYX * RZY * RZY * RZY * HYZZZ +
                            RZX * RXY * RXY * RXY * HXXXZ +
                            RZX * RXY * RXY * RYY * HXXYZ +
                            RZX * RXY * RXY * RZY * HXXZZ +
                            RZX * RXY * RYY * RXY * HXXYZ +
                            RZX * RXY * RYY * RYY * HXYYZ +
                            RZX * RXY * RYY * RZY * HXYZZ +
                            RZX * RXY * RZY * RXY * HXXZZ +
                            RZX * RXY * RZY * RYY * HXYZZ +
                            RZX * RXY * RZY * RZY * HXZZZ +
                            RZX * RYY * RXY * RXY * HXXYZ +
                            RZX * RYY * RXY * RYY * HXYYZ +
                            RZX * RYY * RXY * RZY * HXYZZ +
                            RZX * RYY * RYY * RXY * HXYYZ +
                            RZX * RYY * RYY * RYY * HYYYZ +
                            RZX * RYY * RYY * RZY * HYYZZ +
                            RZX * RYY * RZY * RXY * HXYZZ +
                            RZX * RYY * RZY * RYY * HYYZZ +
                            RZX * RYY * RZY * RZY * HYZZZ +
                            RZX * RZY * RXY * RXY * HXXZZ +
                            RZX * RZY * RXY * RYY * HXYZZ +
                            RZX * RZY * RXY * RZY * HXZZZ +
                            RZX * RZY * RYY * RXY * HXYZZ +
                            RZX * RZY * RYY * RYY * HYYZZ +
                            RZX * RZY * RYY * RZY * HYZZZ +
                            RZX * RZY * RZY * RXY * HXZZZ +
                            RZX * RZY * RZY * RYY * HYZZZ +
                            RZX * RZY * RZY * RZY * HZZZZ ;
            double rHXYYZ = 
                            RXX * RXY * RXY * RXZ * HXXXX +
                            RXX * RXY * RXY * RYZ * HXXXY +
                            RXX * RXY * RXY * RZZ * HXXXZ +
                            RXX * RXY * RYY * RXZ * HXXXY +
                            RXX * RXY * RYY * RYZ * HXXYY +
                            RXX * RXY * RYY * RZZ * HXXYZ +
                            RXX * RXY * RZY * RXZ * HXXXZ +
                            RXX * RXY * RZY * RYZ * HXXYZ +
                            RXX * RXY * RZY * RZZ * HXXZZ +
                            RXX * RYY * RXY * RXZ * HXXXY +
                            RXX * RYY * RXY * RYZ * HXXYY +
                            RXX * RYY * RXY * RZZ * HXXYZ +
                            RXX * RYY * RYY * RXZ * HXXYY +
                            RXX * RYY * RYY * RYZ * HXYYY +
                            RXX * RYY * RYY * RZZ * HXYYZ +
                            RXX * RYY * RZY * RXZ * HXXYZ +
                            RXX * RYY * RZY * RYZ * HXYYZ +
                            RXX * RYY * RZY * RZZ * HXYZZ +
                            RXX * RZY * RXY * RXZ * HXXXZ +
                            RXX * RZY * RXY * RYZ * HXXYZ +
                            RXX * RZY * RXY * RZZ * HXXZZ +
                            RXX * RZY * RYY * RXZ * HXXYZ +
                            RXX * RZY * RYY * RYZ * HXYYZ +
                            RXX * RZY * RYY * RZZ * HXYZZ +
                            RXX * RZY * RZY * RXZ * HXXZZ +
                            RXX * RZY * RZY * RYZ * HXYZZ +
                            RXX * RZY * RZY * RZZ * HXZZZ +
                            RYX * RXY * RXY * RXZ * HXXXY +
                            RYX * RXY * RXY * RYZ * HXXYY +
                            RYX * RXY * RXY * RZZ * HXXYZ +
                            RYX * RXY * RYY * RXZ * HXXYY +
                            RYX * RXY * RYY * RYZ * HXYYY +
                            RYX * RXY * RYY * RZZ * HXYYZ +
                            RYX * RXY * RZY * RXZ * HXXYZ +
                            RYX * RXY * RZY * RYZ * HXYYZ +
                            RYX * RXY * RZY * RZZ * HXYZZ +
                            RYX * RYY * RXY * RXZ * HXXYY +
                            RYX * RYY * RXY * RYZ * HXYYY +
                            RYX * RYY * RXY * RZZ * HXYYZ +
                            RYX * RYY * RYY * RXZ * HXYYY +
                            RYX * RYY * RYY * RYZ * HYYYY +
                            RYX * RYY * RYY * RZZ * HYYYZ +
                            RYX * RYY * RZY * RXZ * HXYYZ +
                            RYX * RYY * RZY * RYZ * HYYYZ +
                            RYX * RYY * RZY * RZZ * HYYZZ +
                            RYX * RZY * RXY * RXZ * HXXYZ +
                            RYX * RZY * RXY * RYZ * HXYYZ +
                            RYX * RZY * RXY * RZZ * HXYZZ +
                            RYX * RZY * RYY * RXZ * HXYYZ +
                            RYX * RZY * RYY * RYZ * HYYYZ +
                            RYX * RZY * RYY * RZZ * HYYZZ +
                            RYX * RZY * RZY * RXZ * HXYZZ +
                            RYX * RZY * RZY * RYZ * HYYZZ +
                            RYX * RZY * RZY * RZZ * HYZZZ +
                            RZX * RXY * RXY * RXZ * HXXXZ +
                            RZX * RXY * RXY * RYZ * HXXYZ +
                            RZX * RXY * RXY * RZZ * HXXZZ +
                            RZX * RXY * RYY * RXZ * HXXYZ +
                            RZX * RXY * RYY * RYZ * HXYYZ +
                            RZX * RXY * RYY * RZZ * HXYZZ +
                            RZX * RXY * RZY * RXZ * HXXZZ +
                            RZX * RXY * RZY * RYZ * HXYZZ +
                            RZX * RXY * RZY * RZZ * HXZZZ +
                            RZX * RYY * RXY * RXZ * HXXYZ +
                            RZX * RYY * RXY * RYZ * HXYYZ +
                            RZX * RYY * RXY * RZZ * HXYZZ +
                            RZX * RYY * RYY * RXZ * HXYYZ +
                            RZX * RYY * RYY * RYZ * HYYYZ +
                            RZX * RYY * RYY * RZZ * HYYZZ +
                            RZX * RYY * RZY * RXZ * HXYZZ +
                            RZX * RYY * RZY * RYZ * HYYZZ +
                            RZX * RYY * RZY * RZZ * HYZZZ +
                            RZX * RZY * RXY * RXZ * HXXZZ +
                            RZX * RZY * RXY * RYZ * HXYZZ +
                            RZX * RZY * RXY * RZZ * HXZZZ +
                            RZX * RZY * RYY * RXZ * HXYZZ +
                            RZX * RZY * RYY * RYZ * HYYZZ +
                            RZX * RZY * RYY * RZZ * HYZZZ +
                            RZX * RZY * RZY * RXZ * HXZZZ +
                            RZX * RZY * RZY * RYZ * HYZZZ +
                            RZX * RZY * RZY * RZZ * HZZZZ ;
            double rHXYZZ = 
                            RXX * RXY * RXZ * RXZ * HXXXX +
                            RXX * RXY * RXZ * RYZ * HXXXY +
                            RXX * RXY * RXZ * RZZ * HXXXZ +
                            RXX * RXY * RYZ * RXZ * HXXXY +
                            RXX * RXY * RYZ * RYZ * HXXYY +
                            RXX * RXY * RYZ * RZZ * HXXYZ +
                            RXX * RXY * RZZ * RXZ * HXXXZ +
                            RXX * RXY * RZZ * RYZ * HXXYZ +
                            RXX * RXY * RZZ * RZZ * HXXZZ +
                            RXX * RYY * RXZ * RXZ * HXXXY +
                            RXX * RYY * RXZ * RYZ * HXXYY +
                            RXX * RYY * RXZ * RZZ * HXXYZ +
                            RXX * RYY * RYZ * RXZ * HXXYY +
                            RXX * RYY * RYZ * RYZ * HXYYY +
                            RXX * RYY * RYZ * RZZ * HXYYZ +
                            RXX * RYY * RZZ * RXZ * HXXYZ +
                            RXX * RYY * RZZ * RYZ * HXYYZ +
                            RXX * RYY * RZZ * RZZ * HXYZZ +
                            RXX * RZY * RXZ * RXZ * HXXXZ +
                            RXX * RZY * RXZ * RYZ * HXXYZ +
                            RXX * RZY * RXZ * RZZ * HXXZZ +
                            RXX * RZY * RYZ * RXZ * HXXYZ +
                            RXX * RZY * RYZ * RYZ * HXYYZ +
                            RXX * RZY * RYZ * RZZ * HXYZZ +
                            RXX * RZY * RZZ * RXZ * HXXZZ +
                            RXX * RZY * RZZ * RYZ * HXYZZ +
                            RXX * RZY * RZZ * RZZ * HXZZZ +
                            RYX * RXY * RXZ * RXZ * HXXXY +
                            RYX * RXY * RXZ * RYZ * HXXYY +
                            RYX * RXY * RXZ * RZZ * HXXYZ +
                            RYX * RXY * RYZ * RXZ * HXXYY +
                            RYX * RXY * RYZ * RYZ * HXYYY +
                            RYX * RXY * RYZ * RZZ * HXYYZ +
                            RYX * RXY * RZZ * RXZ * HXXYZ +
                            RYX * RXY * RZZ * RYZ * HXYYZ +
                            RYX * RXY * RZZ * RZZ * HXYZZ +
                            RYX * RYY * RXZ * RXZ * HXXYY +
                            RYX * RYY * RXZ * RYZ * HXYYY +
                            RYX * RYY * RXZ * RZZ * HXYYZ +
                            RYX * RYY * RYZ * RXZ * HXYYY +
                            RYX * RYY * RYZ * RYZ * HYYYY +
                            RYX * RYY * RYZ * RZZ * HYYYZ +
                            RYX * RYY * RZZ * RXZ * HXYYZ +
                            RYX * RYY * RZZ * RYZ * HYYYZ +
                            RYX * RYY * RZZ * RZZ * HYYZZ +
                            RYX * RZY * RXZ * RXZ * HXXYZ +
                            RYX * RZY * RXZ * RYZ * HXYYZ +
                            RYX * RZY * RXZ * RZZ * HXYZZ +
                            RYX * RZY * RYZ * RXZ * HXYYZ +
                            RYX * RZY * RYZ * RYZ * HYYYZ +
                            RYX * RZY * RYZ * RZZ * HYYZZ +
                            RYX * RZY * RZZ * RXZ * HXYZZ +
                            RYX * RZY * RZZ * RYZ * HYYZZ +
                            RYX * RZY * RZZ * RZZ * HYZZZ +
                            RZX * RXY * RXZ * RXZ * HXXXZ +
                            RZX * RXY * RXZ * RYZ * HXXYZ +
                            RZX * RXY * RXZ * RZZ * HXXZZ +
                            RZX * RXY * RYZ * RXZ * HXXYZ +
                            RZX * RXY * RYZ * RYZ * HXYYZ +
                            RZX * RXY * RYZ * RZZ * HXYZZ +
                            RZX * RXY * RZZ * RXZ * HXXZZ +
                            RZX * RXY * RZZ * RYZ * HXYZZ +
                            RZX * RXY * RZZ * RZZ * HXZZZ +
                            RZX * RYY * RXZ * RXZ * HXXYZ +
                            RZX * RYY * RXZ * RYZ * HXYYZ +
                            RZX * RYY * RXZ * RZZ * HXYZZ +
                            RZX * RYY * RYZ * RXZ * HXYYZ +
                            RZX * RYY * RYZ * RYZ * HYYYZ +
                            RZX * RYY * RYZ * RZZ * HYYZZ +
                            RZX * RYY * RZZ * RXZ * HXYZZ +
                            RZX * RYY * RZZ * RYZ * HYYZZ +
                            RZX * RYY * RZZ * RZZ * HYZZZ +
                            RZX * RZY * RXZ * RXZ * HXXZZ +
                            RZX * RZY * RXZ * RYZ * HXYZZ +
                            RZX * RZY * RXZ * RZZ * HXZZZ +
                            RZX * RZY * RYZ * RXZ * HXYZZ +
                            RZX * RZY * RYZ * RYZ * HYYZZ +
                            RZX * RZY * RYZ * RZZ * HYZZZ +
                            RZX * RZY * RZZ * RXZ * HXZZZ +
                            RZX * RZY * RZZ * RYZ * HYZZZ +
                            RZX * RZY * RZZ * RZZ * HZZZZ ;
            double rHXZZZ = 
                            RXX * RXZ * RXZ * RXZ * HXXXX +
                            RXX * RXZ * RXZ * RYZ * HXXXY +
                            RXX * RXZ * RXZ * RZZ * HXXXZ +
                            RXX * RXZ * RYZ * RXZ * HXXXY +
                            RXX * RXZ * RYZ * RYZ * HXXYY +
                            RXX * RXZ * RYZ * RZZ * HXXYZ +
                            RXX * RXZ * RZZ * RXZ * HXXXZ +
                            RXX * RXZ * RZZ * RYZ * HXXYZ +
                            RXX * RXZ * RZZ * RZZ * HXXZZ +
                            RXX * RYZ * RXZ * RXZ * HXXXY +
                            RXX * RYZ * RXZ * RYZ * HXXYY +
                            RXX * RYZ * RXZ * RZZ * HXXYZ +
                            RXX * RYZ * RYZ * RXZ * HXXYY +
                            RXX * RYZ * RYZ * RYZ * HXYYY +
                            RXX * RYZ * RYZ * RZZ * HXYYZ +
                            RXX * RYZ * RZZ * RXZ * HXXYZ +
                            RXX * RYZ * RZZ * RYZ * HXYYZ +
                            RXX * RYZ * RZZ * RZZ * HXYZZ +
                            RXX * RZZ * RXZ * RXZ * HXXXZ +
                            RXX * RZZ * RXZ * RYZ * HXXYZ +
                            RXX * RZZ * RXZ * RZZ * HXXZZ +
                            RXX * RZZ * RYZ * RXZ * HXXYZ +
                            RXX * RZZ * RYZ * RYZ * HXYYZ +
                            RXX * RZZ * RYZ * RZZ * HXYZZ +
                            RXX * RZZ * RZZ * RXZ * HXXZZ +
                            RXX * RZZ * RZZ * RYZ * HXYZZ +
                            RXX * RZZ * RZZ * RZZ * HXZZZ +
                            RYX * RXZ * RXZ * RXZ * HXXXY +
                            RYX * RXZ * RXZ * RYZ * HXXYY +
                            RYX * RXZ * RXZ * RZZ * HXXYZ +
                            RYX * RXZ * RYZ * RXZ * HXXYY +
                            RYX * RXZ * RYZ * RYZ * HXYYY +
                            RYX * RXZ * RYZ * RZZ * HXYYZ +
                            RYX * RXZ * RZZ * RXZ * HXXYZ +
                            RYX * RXZ * RZZ * RYZ * HXYYZ +
                            RYX * RXZ * RZZ * RZZ * HXYZZ +
                            RYX * RYZ * RXZ * RXZ * HXXYY +
                            RYX * RYZ * RXZ * RYZ * HXYYY +
                            RYX * RYZ * RXZ * RZZ * HXYYZ +
                            RYX * RYZ * RYZ * RXZ * HXYYY +
                            RYX * RYZ * RYZ * RYZ * HYYYY +
                            RYX * RYZ * RYZ * RZZ * HYYYZ +
                            RYX * RYZ * RZZ * RXZ * HXYYZ +
                            RYX * RYZ * RZZ * RYZ * HYYYZ +
                            RYX * RYZ * RZZ * RZZ * HYYZZ +
                            RYX * RZZ * RXZ * RXZ * HXXYZ +
                            RYX * RZZ * RXZ * RYZ * HXYYZ +
                            RYX * RZZ * RXZ * RZZ * HXYZZ +
                            RYX * RZZ * RYZ * RXZ * HXYYZ +
                            RYX * RZZ * RYZ * RYZ * HYYYZ +
                            RYX * RZZ * RYZ * RZZ * HYYZZ +
                            RYX * RZZ * RZZ * RXZ * HXYZZ +
                            RYX * RZZ * RZZ * RYZ * HYYZZ +
                            RYX * RZZ * RZZ * RZZ * HYZZZ +
                            RZX * RXZ * RXZ * RXZ * HXXXZ +
                            RZX * RXZ * RXZ * RYZ * HXXYZ +
                            RZX * RXZ * RXZ * RZZ * HXXZZ +
                            RZX * RXZ * RYZ * RXZ * HXXYZ +
                            RZX * RXZ * RYZ * RYZ * HXYYZ +
                            RZX * RXZ * RYZ * RZZ * HXYZZ +
                            RZX * RXZ * RZZ * RXZ * HXXZZ +
                            RZX * RXZ * RZZ * RYZ * HXYZZ +
                            RZX * RXZ * RZZ * RZZ * HXZZZ +
                            RZX * RYZ * RXZ * RXZ * HXXYZ +
                            RZX * RYZ * RXZ * RYZ * HXYYZ +
                            RZX * RYZ * RXZ * RZZ * HXYZZ +
                            RZX * RYZ * RYZ * RXZ * HXYYZ +
                            RZX * RYZ * RYZ * RYZ * HYYYZ +
                            RZX * RYZ * RYZ * RZZ * HYYZZ +
                            RZX * RYZ * RZZ * RXZ * HXYZZ +
                            RZX * RYZ * RZZ * RYZ * HYYZZ +
                            RZX * RYZ * RZZ * RZZ * HYZZZ +
                            RZX * RZZ * RXZ * RXZ * HXXZZ +
                            RZX * RZZ * RXZ * RYZ * HXYZZ +
                            RZX * RZZ * RXZ * RZZ * HXZZZ +
                            RZX * RZZ * RYZ * RXZ * HXYZZ +
                            RZX * RZZ * RYZ * RYZ * HYYZZ +
                            RZX * RZZ * RYZ * RZZ * HYZZZ +
                            RZX * RZZ * RZZ * RXZ * HXZZZ +
                            RZX * RZZ * RZZ * RYZ * HYZZZ +
                            RZX * RZZ * RZZ * RZZ * HZZZZ ;
            double rHYYYY = 
                            RXY * RXY * RXY * RXY * HXXXX +
                            RXY * RXY * RXY * RYY * HXXXY +
                            RXY * RXY * RXY * RZY * HXXXZ +
                            RXY * RXY * RYY * RXY * HXXXY +
                            RXY * RXY * RYY * RYY * HXXYY +
                            RXY * RXY * RYY * RZY * HXXYZ +
                            RXY * RXY * RZY * RXY * HXXXZ +
                            RXY * RXY * RZY * RYY * HXXYZ +
                            RXY * RXY * RZY * RZY * HXXZZ +
                            RXY * RYY * RXY * RXY * HXXXY +
                            RXY * RYY * RXY * RYY * HXXYY +
                            RXY * RYY * RXY * RZY * HXXYZ +
                            RXY * RYY * RYY * RXY * HXXYY +
                            RXY * RYY * RYY * RYY * HXYYY +
                            RXY * RYY * RYY * RZY * HXYYZ +
                            RXY * RYY * RZY * RXY * HXXYZ +
                            RXY * RYY * RZY * RYY * HXYYZ +
                            RXY * RYY * RZY * RZY * HXYZZ +
                            RXY * RZY * RXY * RXY * HXXXZ +
                            RXY * RZY * RXY * RYY * HXXYZ +
                            RXY * RZY * RXY * RZY * HXXZZ +
                            RXY * RZY * RYY * RXY * HXXYZ +
                            RXY * RZY * RYY * RYY * HXYYZ +
                            RXY * RZY * RYY * RZY * HXYZZ +
                            RXY * RZY * RZY * RXY * HXXZZ +
                            RXY * RZY * RZY * RYY * HXYZZ +
                            RXY * RZY * RZY * RZY * HXZZZ +
                            RYY * RXY * RXY * RXY * HXXXY +
                            RYY * RXY * RXY * RYY * HXXYY +
                            RYY * RXY * RXY * RZY * HXXYZ +
                            RYY * RXY * RYY * RXY * HXXYY +
                            RYY * RXY * RYY * RYY * HXYYY +
                            RYY * RXY * RYY * RZY * HXYYZ +
                            RYY * RXY * RZY * RXY * HXXYZ +
                            RYY * RXY * RZY * RYY * HXYYZ +
                            RYY * RXY * RZY * RZY * HXYZZ +
                            RYY * RYY * RXY * RXY * HXXYY +
                            RYY * RYY * RXY * RYY * HXYYY +
                            RYY * RYY * RXY * RZY * HXYYZ +
                            RYY * RYY * RYY * RXY * HXYYY +
                            RYY * RYY * RYY * RYY * HYYYY +
                            RYY * RYY * RYY * RZY * HYYYZ +
                            RYY * RYY * RZY * RXY * HXYYZ +
                            RYY * RYY * RZY * RYY * HYYYZ +
                            RYY * RYY * RZY * RZY * HYYZZ +
                            RYY * RZY * RXY * RXY * HXXYZ +
                            RYY * RZY * RXY * RYY * HXYYZ +
                            RYY * RZY * RXY * RZY * HXYZZ +
                            RYY * RZY * RYY * RXY * HXYYZ +
                            RYY * RZY * RYY * RYY * HYYYZ +
                            RYY * RZY * RYY * RZY * HYYZZ +
                            RYY * RZY * RZY * RXY * HXYZZ +
                            RYY * RZY * RZY * RYY * HYYZZ +
                            RYY * RZY * RZY * RZY * HYZZZ +
                            RZY * RXY * RXY * RXY * HXXXZ +
                            RZY * RXY * RXY * RYY * HXXYZ +
                            RZY * RXY * RXY * RZY * HXXZZ +
                            RZY * RXY * RYY * RXY * HXXYZ +
                            RZY * RXY * RYY * RYY * HXYYZ +
                            RZY * RXY * RYY * RZY * HXYZZ +
                            RZY * RXY * RZY * RXY * HXXZZ +
                            RZY * RXY * RZY * RYY * HXYZZ +
                            RZY * RXY * RZY * RZY * HXZZZ +
                            RZY * RYY * RXY * RXY * HXXYZ +
                            RZY * RYY * RXY * RYY * HXYYZ +
                            RZY * RYY * RXY * RZY * HXYZZ +
                            RZY * RYY * RYY * RXY * HXYYZ +
                            RZY * RYY * RYY * RYY * HYYYZ +
                            RZY * RYY * RYY * RZY * HYYZZ +
                            RZY * RYY * RZY * RXY * HXYZZ +
                            RZY * RYY * RZY * RYY * HYYZZ +
                            RZY * RYY * RZY * RZY * HYZZZ +
                            RZY * RZY * RXY * RXY * HXXZZ +
                            RZY * RZY * RXY * RYY * HXYZZ +
                            RZY * RZY * RXY * RZY * HXZZZ +
                            RZY * RZY * RYY * RXY * HXYZZ +
                            RZY * RZY * RYY * RYY * HYYZZ +
                            RZY * RZY * RYY * RZY * HYZZZ +
                            RZY * RZY * RZY * RXY * HXZZZ +
                            RZY * RZY * RZY * RYY * HYZZZ +
                            RZY * RZY * RZY * RZY * HZZZZ ;
            double rHYYYZ = 
                            RXY * RXY * RXY * RXZ * HXXXX +
                            RXY * RXY * RXY * RYZ * HXXXY +
                            RXY * RXY * RXY * RZZ * HXXXZ +
                            RXY * RXY * RYY * RXZ * HXXXY +
                            RXY * RXY * RYY * RYZ * HXXYY +
                            RXY * RXY * RYY * RZZ * HXXYZ +
                            RXY * RXY * RZY * RXZ * HXXXZ +
                            RXY * RXY * RZY * RYZ * HXXYZ +
                            RXY * RXY * RZY * RZZ * HXXZZ +
                            RXY * RYY * RXY * RXZ * HXXXY +
                            RXY * RYY * RXY * RYZ * HXXYY +
                            RXY * RYY * RXY * RZZ * HXXYZ +
                            RXY * RYY * RYY * RXZ * HXXYY +
                            RXY * RYY * RYY * RYZ * HXYYY +
                            RXY * RYY * RYY * RZZ * HXYYZ +
                            RXY * RYY * RZY * RXZ * HXXYZ +
                            RXY * RYY * RZY * RYZ * HXYYZ +
                            RXY * RYY * RZY * RZZ * HXYZZ +
                            RXY * RZY * RXY * RXZ * HXXXZ +
                            RXY * RZY * RXY * RYZ * HXXYZ +
                            RXY * RZY * RXY * RZZ * HXXZZ +
                            RXY * RZY * RYY * RXZ * HXXYZ +
                            RXY * RZY * RYY * RYZ * HXYYZ +
                            RXY * RZY * RYY * RZZ * HXYZZ +
                            RXY * RZY * RZY * RXZ * HXXZZ +
                            RXY * RZY * RZY * RYZ * HXYZZ +
                            RXY * RZY * RZY * RZZ * HXZZZ +
                            RYY * RXY * RXY * RXZ * HXXXY +
                            RYY * RXY * RXY * RYZ * HXXYY +
                            RYY * RXY * RXY * RZZ * HXXYZ +
                            RYY * RXY * RYY * RXZ * HXXYY +
                            RYY * RXY * RYY * RYZ * HXYYY +
                            RYY * RXY * RYY * RZZ * HXYYZ +
                            RYY * RXY * RZY * RXZ * HXXYZ +
                            RYY * RXY * RZY * RYZ * HXYYZ +
                            RYY * RXY * RZY * RZZ * HXYZZ +
                            RYY * RYY * RXY * RXZ * HXXYY +
                            RYY * RYY * RXY * RYZ * HXYYY +
                            RYY * RYY * RXY * RZZ * HXYYZ +
                            RYY * RYY * RYY * RXZ * HXYYY +
                            RYY * RYY * RYY * RYZ * HYYYY +
                            RYY * RYY * RYY * RZZ * HYYYZ +
                            RYY * RYY * RZY * RXZ * HXYYZ +
                            RYY * RYY * RZY * RYZ * HYYYZ +
                            RYY * RYY * RZY * RZZ * HYYZZ +
                            RYY * RZY * RXY * RXZ * HXXYZ +
                            RYY * RZY * RXY * RYZ * HXYYZ +
                            RYY * RZY * RXY * RZZ * HXYZZ +
                            RYY * RZY * RYY * RXZ * HXYYZ +
                            RYY * RZY * RYY * RYZ * HYYYZ +
                            RYY * RZY * RYY * RZZ * HYYZZ +
                            RYY * RZY * RZY * RXZ * HXYZZ +
                            RYY * RZY * RZY * RYZ * HYYZZ +
                            RYY * RZY * RZY * RZZ * HYZZZ +
                            RZY * RXY * RXY * RXZ * HXXXZ +
                            RZY * RXY * RXY * RYZ * HXXYZ +
                            RZY * RXY * RXY * RZZ * HXXZZ +
                            RZY * RXY * RYY * RXZ * HXXYZ +
                            RZY * RXY * RYY * RYZ * HXYYZ +
                            RZY * RXY * RYY * RZZ * HXYZZ +
                            RZY * RXY * RZY * RXZ * HXXZZ +
                            RZY * RXY * RZY * RYZ * HXYZZ +
                            RZY * RXY * RZY * RZZ * HXZZZ +
                            RZY * RYY * RXY * RXZ * HXXYZ +
                            RZY * RYY * RXY * RYZ * HXYYZ +
                            RZY * RYY * RXY * RZZ * HXYZZ +
                            RZY * RYY * RYY * RXZ * HXYYZ +
                            RZY * RYY * RYY * RYZ * HYYYZ +
                            RZY * RYY * RYY * RZZ * HYYZZ +
                            RZY * RYY * RZY * RXZ * HXYZZ +
                            RZY * RYY * RZY * RYZ * HYYZZ +
                            RZY * RYY * RZY * RZZ * HYZZZ +
                            RZY * RZY * RXY * RXZ * HXXZZ +
                            RZY * RZY * RXY * RYZ * HXYZZ +
                            RZY * RZY * RXY * RZZ * HXZZZ +
                            RZY * RZY * RYY * RXZ * HXYZZ +
                            RZY * RZY * RYY * RYZ * HYYZZ +
                            RZY * RZY * RYY * RZZ * HYZZZ +
                            RZY * RZY * RZY * RXZ * HXZZZ +
                            RZY * RZY * RZY * RYZ * HYZZZ +
                            RZY * RZY * RZY * RZZ * HZZZZ ;
            double rHYYZZ = 
                            RXY * RXY * RXZ * RXZ * HXXXX +
                            RXY * RXY * RXZ * RYZ * HXXXY +
                            RXY * RXY * RXZ * RZZ * HXXXZ +
                            RXY * RXY * RYZ * RXZ * HXXXY +
                            RXY * RXY * RYZ * RYZ * HXXYY +
                            RXY * RXY * RYZ * RZZ * HXXYZ +
                            RXY * RXY * RZZ * RXZ * HXXXZ +
                            RXY * RXY * RZZ * RYZ * HXXYZ +
                            RXY * RXY * RZZ * RZZ * HXXZZ +
                            RXY * RYY * RXZ * RXZ * HXXXY +
                            RXY * RYY * RXZ * RYZ * HXXYY +
                            RXY * RYY * RXZ * RZZ * HXXYZ +
                            RXY * RYY * RYZ * RXZ * HXXYY +
                            RXY * RYY * RYZ * RYZ * HXYYY +
                            RXY * RYY * RYZ * RZZ * HXYYZ +
                            RXY * RYY * RZZ * RXZ * HXXYZ +
                            RXY * RYY * RZZ * RYZ * HXYYZ +
                            RXY * RYY * RZZ * RZZ * HXYZZ +
                            RXY * RZY * RXZ * RXZ * HXXXZ +
                            RXY * RZY * RXZ * RYZ * HXXYZ +
                            RXY * RZY * RXZ * RZZ * HXXZZ +
                            RXY * RZY * RYZ * RXZ * HXXYZ +
                            RXY * RZY * RYZ * RYZ * HXYYZ +
                            RXY * RZY * RYZ * RZZ * HXYZZ +
                            RXY * RZY * RZZ * RXZ * HXXZZ +
                            RXY * RZY * RZZ * RYZ * HXYZZ +
                            RXY * RZY * RZZ * RZZ * HXZZZ +
                            RYY * RXY * RXZ * RXZ * HXXXY +
                            RYY * RXY * RXZ * RYZ * HXXYY +
                            RYY * RXY * RXZ * RZZ * HXXYZ +
                            RYY * RXY * RYZ * RXZ * HXXYY +
                            RYY * RXY * RYZ * RYZ * HXYYY +
                            RYY * RXY * RYZ * RZZ * HXYYZ +
                            RYY * RXY * RZZ * RXZ * HXXYZ +
                            RYY * RXY * RZZ * RYZ * HXYYZ +
                            RYY * RXY * RZZ * RZZ * HXYZZ +
                            RYY * RYY * RXZ * RXZ * HXXYY +
                            RYY * RYY * RXZ * RYZ * HXYYY +
                            RYY * RYY * RXZ * RZZ * HXYYZ +
                            RYY * RYY * RYZ * RXZ * HXYYY +
                            RYY * RYY * RYZ * RYZ * HYYYY +
                            RYY * RYY * RYZ * RZZ * HYYYZ +
                            RYY * RYY * RZZ * RXZ * HXYYZ +
                            RYY * RYY * RZZ * RYZ * HYYYZ +
                            RYY * RYY * RZZ * RZZ * HYYZZ +
                            RYY * RZY * RXZ * RXZ * HXXYZ +
                            RYY * RZY * RXZ * RYZ * HXYYZ +
                            RYY * RZY * RXZ * RZZ * HXYZZ +
                            RYY * RZY * RYZ * RXZ * HXYYZ +
                            RYY * RZY * RYZ * RYZ * HYYYZ +
                            RYY * RZY * RYZ * RZZ * HYYZZ +
                            RYY * RZY * RZZ * RXZ * HXYZZ +
                            RYY * RZY * RZZ * RYZ * HYYZZ +
                            RYY * RZY * RZZ * RZZ * HYZZZ +
                            RZY * RXY * RXZ * RXZ * HXXXZ +
                            RZY * RXY * RXZ * RYZ * HXXYZ +
                            RZY * RXY * RXZ * RZZ * HXXZZ +
                            RZY * RXY * RYZ * RXZ * HXXYZ +
                            RZY * RXY * RYZ * RYZ * HXYYZ +
                            RZY * RXY * RYZ * RZZ * HXYZZ +
                            RZY * RXY * RZZ * RXZ * HXXZZ +
                            RZY * RXY * RZZ * RYZ * HXYZZ +
                            RZY * RXY * RZZ * RZZ * HXZZZ +
                            RZY * RYY * RXZ * RXZ * HXXYZ +
                            RZY * RYY * RXZ * RYZ * HXYYZ +
                            RZY * RYY * RXZ * RZZ * HXYZZ +
                            RZY * RYY * RYZ * RXZ * HXYYZ +
                            RZY * RYY * RYZ * RYZ * HYYYZ +
                            RZY * RYY * RYZ * RZZ * HYYZZ +
                            RZY * RYY * RZZ * RXZ * HXYZZ +
                            RZY * RYY * RZZ * RYZ * HYYZZ +
                            RZY * RYY * RZZ * RZZ * HYZZZ +
                            RZY * RZY * RXZ * RXZ * HXXZZ +
                            RZY * RZY * RXZ * RYZ * HXYZZ +
                            RZY * RZY * RXZ * RZZ * HXZZZ +
                            RZY * RZY * RYZ * RXZ * HXYZZ +
                            RZY * RZY * RYZ * RYZ * HYYZZ +
                            RZY * RZY * RYZ * RZZ * HYZZZ +
                            RZY * RZY * RZZ * RXZ * HXZZZ +
                            RZY * RZY * RZZ * RYZ * HYZZZ +
                            RZY * RZY * RZZ * RZZ * HZZZZ ;
            double rHYZZZ = 
                            RXY * RXZ * RXZ * RXZ * HXXXX +
                            RXY * RXZ * RXZ * RYZ * HXXXY +
                            RXY * RXZ * RXZ * RZZ * HXXXZ +
                            RXY * RXZ * RYZ * RXZ * HXXXY +
                            RXY * RXZ * RYZ * RYZ * HXXYY +
                            RXY * RXZ * RYZ * RZZ * HXXYZ +
                            RXY * RXZ * RZZ * RXZ * HXXXZ +
                            RXY * RXZ * RZZ * RYZ * HXXYZ +
                            RXY * RXZ * RZZ * RZZ * HXXZZ +
                            RXY * RYZ * RXZ * RXZ * HXXXY +
                            RXY * RYZ * RXZ * RYZ * HXXYY +
                            RXY * RYZ * RXZ * RZZ * HXXYZ +
                            RXY * RYZ * RYZ * RXZ * HXXYY +
                            RXY * RYZ * RYZ * RYZ * HXYYY +
                            RXY * RYZ * RYZ * RZZ * HXYYZ +
                            RXY * RYZ * RZZ * RXZ * HXXYZ +
                            RXY * RYZ * RZZ * RYZ * HXYYZ +
                            RXY * RYZ * RZZ * RZZ * HXYZZ +
                            RXY * RZZ * RXZ * RXZ * HXXXZ +
                            RXY * RZZ * RXZ * RYZ * HXXYZ +
                            RXY * RZZ * RXZ * RZZ * HXXZZ +
                            RXY * RZZ * RYZ * RXZ * HXXYZ +
                            RXY * RZZ * RYZ * RYZ * HXYYZ +
                            RXY * RZZ * RYZ * RZZ * HXYZZ +
                            RXY * RZZ * RZZ * RXZ * HXXZZ +
                            RXY * RZZ * RZZ * RYZ * HXYZZ +
                            RXY * RZZ * RZZ * RZZ * HXZZZ +
                            RYY * RXZ * RXZ * RXZ * HXXXY +
                            RYY * RXZ * RXZ * RYZ * HXXYY +
                            RYY * RXZ * RXZ * RZZ * HXXYZ +
                            RYY * RXZ * RYZ * RXZ * HXXYY +
                            RYY * RXZ * RYZ * RYZ * HXYYY +
                            RYY * RXZ * RYZ * RZZ * HXYYZ +
                            RYY * RXZ * RZZ * RXZ * HXXYZ +
                            RYY * RXZ * RZZ * RYZ * HXYYZ +
                            RYY * RXZ * RZZ * RZZ * HXYZZ +
                            RYY * RYZ * RXZ * RXZ * HXXYY +
                            RYY * RYZ * RXZ * RYZ * HXYYY +
                            RYY * RYZ * RXZ * RZZ * HXYYZ +
                            RYY * RYZ * RYZ * RXZ * HXYYY +
                            RYY * RYZ * RYZ * RYZ * HYYYY +
                            RYY * RYZ * RYZ * RZZ * HYYYZ +
                            RYY * RYZ * RZZ * RXZ * HXYYZ +
                            RYY * RYZ * RZZ * RYZ * HYYYZ +
                            RYY * RYZ * RZZ * RZZ * HYYZZ +
                            RYY * RZZ * RXZ * RXZ * HXXYZ +
                            RYY * RZZ * RXZ * RYZ * HXYYZ +
                            RYY * RZZ * RXZ * RZZ * HXYZZ +
                            RYY * RZZ * RYZ * RXZ * HXYYZ +
                            RYY * RZZ * RYZ * RYZ * HYYYZ +
                            RYY * RZZ * RYZ * RZZ * HYYZZ +
                            RYY * RZZ * RZZ * RXZ * HXYZZ +
                            RYY * RZZ * RZZ * RYZ * HYYZZ +
                            RYY * RZZ * RZZ * RZZ * HYZZZ +
                            RZY * RXZ * RXZ * RXZ * HXXXZ +
                            RZY * RXZ * RXZ * RYZ * HXXYZ +
                            RZY * RXZ * RXZ * RZZ * HXXZZ +
                            RZY * RXZ * RYZ * RXZ * HXXYZ +
                            RZY * RXZ * RYZ * RYZ * HXYYZ +
                            RZY * RXZ * RYZ * RZZ * HXYZZ +
                            RZY * RXZ * RZZ * RXZ * HXXZZ +
                            RZY * RXZ * RZZ * RYZ * HXYZZ +
                            RZY * RXZ * RZZ * RZZ * HXZZZ +
                            RZY * RYZ * RXZ * RXZ * HXXYZ +
                            RZY * RYZ * RXZ * RYZ * HXYYZ +
                            RZY * RYZ * RXZ * RZZ * HXYZZ +
                            RZY * RYZ * RYZ * RXZ * HXYYZ +
                            RZY * RYZ * RYZ * RYZ * HYYYZ +
                            RZY * RYZ * RYZ * RZZ * HYYZZ +
                            RZY * RYZ * RZZ * RXZ * HXYZZ +
                            RZY * RYZ * RZZ * RYZ * HYYZZ +
                            RZY * RYZ * RZZ * RZZ * HYZZZ +
                            RZY * RZZ * RXZ * RXZ * HXXZZ +
                            RZY * RZZ * RXZ * RYZ * HXYZZ +
                            RZY * RZZ * RXZ * RZZ * HXZZZ +
                            RZY * RZZ * RYZ * RXZ * HXYZZ +
                            RZY * RZZ * RYZ * RYZ * HYYZZ +
                            RZY * RZZ * RYZ * RZZ * HYZZZ +
                            RZY * RZZ * RZZ * RXZ * HXZZZ +
                            RZY * RZZ * RZZ * RYZ * HYZZZ +
                            RZY * RZZ * RZZ * RZZ * HZZZZ ;
            double rHZZZZ = 
                            RXZ * RXZ * RXZ * RXZ * HXXXX +
                            RXZ * RXZ * RXZ * RYZ * HXXXY +
                            RXZ * RXZ * RXZ * RZZ * HXXXZ +
                            RXZ * RXZ * RYZ * RXZ * HXXXY +
                            RXZ * RXZ * RYZ * RYZ * HXXYY +
                            RXZ * RXZ * RYZ * RZZ * HXXYZ +
                            RXZ * RXZ * RZZ * RXZ * HXXXZ +
                            RXZ * RXZ * RZZ * RYZ * HXXYZ +
                            RXZ * RXZ * RZZ * RZZ * HXXZZ +
                            RXZ * RYZ * RXZ * RXZ * HXXXY +
                            RXZ * RYZ * RXZ * RYZ * HXXYY +
                            RXZ * RYZ * RXZ * RZZ * HXXYZ +
                            RXZ * RYZ * RYZ * RXZ * HXXYY +
                            RXZ * RYZ * RYZ * RYZ * HXYYY +
                            RXZ * RYZ * RYZ * RZZ * HXYYZ +
                            RXZ * RYZ * RZZ * RXZ * HXXYZ +
                            RXZ * RYZ * RZZ * RYZ * HXYYZ +
                            RXZ * RYZ * RZZ * RZZ * HXYZZ +
                            RXZ * RZZ * RXZ * RXZ * HXXXZ +
                            RXZ * RZZ * RXZ * RYZ * HXXYZ +
                            RXZ * RZZ * RXZ * RZZ * HXXZZ +
                            RXZ * RZZ * RYZ * RXZ * HXXYZ +
                            RXZ * RZZ * RYZ * RYZ * HXYYZ +
                            RXZ * RZZ * RYZ * RZZ * HXYZZ +
                            RXZ * RZZ * RZZ * RXZ * HXXZZ +
                            RXZ * RZZ * RZZ * RYZ * HXYZZ +
                            RXZ * RZZ * RZZ * RZZ * HXZZZ +
                            RYZ * RXZ * RXZ * RXZ * HXXXY +
                            RYZ * RXZ * RXZ * RYZ * HXXYY +
                            RYZ * RXZ * RXZ * RZZ * HXXYZ +
                            RYZ * RXZ * RYZ * RXZ * HXXYY +
                            RYZ * RXZ * RYZ * RYZ * HXYYY +
                            RYZ * RXZ * RYZ * RZZ * HXYYZ +
                            RYZ * RXZ * RZZ * RXZ * HXXYZ +
                            RYZ * RXZ * RZZ * RYZ * HXYYZ +
                            RYZ * RXZ * RZZ * RZZ * HXYZZ +
                            RYZ * RYZ * RXZ * RXZ * HXXYY +
                            RYZ * RYZ * RXZ * RYZ * HXYYY +
                            RYZ * RYZ * RXZ * RZZ * HXYYZ +
                            RYZ * RYZ * RYZ * RXZ * HXYYY +
                            RYZ * RYZ * RYZ * RYZ * HYYYY +
                            RYZ * RYZ * RYZ * RZZ * HYYYZ +
                            RYZ * RYZ * RZZ * RXZ * HXYYZ +
                            RYZ * RYZ * RZZ * RYZ * HYYYZ +
                            RYZ * RYZ * RZZ * RZZ * HYYZZ +
                            RYZ * RZZ * RXZ * RXZ * HXXYZ +
                            RYZ * RZZ * RXZ * RYZ * HXYYZ +
                            RYZ * RZZ * RXZ * RZZ * HXYZZ +
                            RYZ * RZZ * RYZ * RXZ * HXYYZ +
                            RYZ * RZZ * RYZ * RYZ * HYYYZ +
                            RYZ * RZZ * RYZ * RZZ * HYYZZ +
                            RYZ * RZZ * RZZ * RXZ * HXYZZ +
                            RYZ * RZZ * RZZ * RYZ * HYYZZ +
                            RYZ * RZZ * RZZ * RZZ * HYZZZ +
                            RZZ * RXZ * RXZ * RXZ * HXXXZ +
                            RZZ * RXZ * RXZ * RYZ * HXXYZ +
                            RZZ * RXZ * RXZ * RZZ * HXXZZ +
                            RZZ * RXZ * RYZ * RXZ * HXXYZ +
                            RZZ * RXZ * RYZ * RYZ * HXYYZ +
                            RZZ * RXZ * RYZ * RZZ * HXYZZ +
                            RZZ * RXZ * RZZ * RXZ * HXXZZ +
                            RZZ * RXZ * RZZ * RYZ * HXYZZ +
                            RZZ * RXZ * RZZ * RZZ * HXZZZ +
                            RZZ * RYZ * RXZ * RXZ * HXXYZ +
                            RZZ * RYZ * RXZ * RYZ * HXYYZ +
                            RZZ * RYZ * RXZ * RZZ * HXYZZ +
                            RZZ * RYZ * RYZ * RXZ * HXYYZ +
                            RZZ * RYZ * RYZ * RYZ * HYYYZ +
                            RZZ * RYZ * RYZ * RZZ * HYYZZ +
                            RZZ * RYZ * RZZ * RXZ * HXYZZ +
                            RZZ * RYZ * RZZ * RYZ * HYYZZ +
                            RZZ * RYZ * RZZ * RZZ * HYZZZ +
                            RZZ * RZZ * RXZ * RXZ * HXXZZ +
                            RZZ * RZZ * RXZ * RYZ * HXYZZ +
                            RZZ * RZZ * RXZ * RZZ * HXZZZ +
                            RZZ * RZZ * RYZ * RXZ * HXYZZ +
                            RZZ * RZZ * RYZ * RYZ * HYYZZ +
                            RZZ * RZZ * RYZ * RZZ * HYZZZ +
                            RZZ * RZZ * RZZ * RXZ * HXZZZ +
                            RZZ * RZZ * RZZ * RYZ * HYZZZ +
                            RZZ * RZZ * RZZ * RZZ * HZZZZ ;

            H[i][ 0] = rHXXXX; 
            H[i][ 1] = rHXXXY;
            H[i][ 2] = rHXXXZ;
            H[i][ 3] = rHXXYY;
            H[i][ 4] = rHXXYZ;
            H[i][ 5] = rHXXZZ;
            H[i][ 6] = rHXYYY;
            H[i][ 7] = rHXYYZ;
            H[i][ 8] = rHXYZZ;
            H[i][ 9] = rHXZZZ;
            H[i][10] = rHYYYY;
            H[i][11] = rHYYYZ;
            H[i][12] = rHYYZZ;
            H[i][13] = rHYZZZ;
            H[i][14] = rHZZZZ;
       }
  }
}
/// Superimpose the DMTP sets
double DMTPole::superimpose(psi::SharedMatrix ref_xyz, std::vector<int> suplist) 
{
   // Initialize
   KabschSuperimposer sup = KabschSuperimposer();

   // Determine the overlapping structure slices
   psi::SharedMatrix initial_xyz, final_xyz;
   if (suplist.empty()) {
       initial_xyz = this->centres_;
       final_xyz = ref_xyz;
   } else {
       const int n = suplist.size();
       initial_xyz = std::make_shared<psi::Matrix>("", n, 3);
       final_xyz   = std::make_shared<psi::Matrix>("", n, 3);
       for (int i=0; i<n; ++i) {
            initial_xyz->set_row(0, suplist[i], this->centres_->get_row(0, suplist[i]));
            final_xyz  ->set_row(0, suplist[i],        ref_xyz->get_row(0, suplist[i]));
       }
   }

   // Superimpose
   sup.compute(initial_xyz, final_xyz);
   psi::SharedMatrix r = sup.rotation;
   psi::SharedVector t = sup.translation;

   this->rotate(r);
   this->translate(t);

   double rms = sup.rms();
   return rms;
}

void DMTPole::print(void) const
{
	centres_->print();
	origins_->print();
	for (int i=0; i<nDMTPs_; ++i) {
           charges_[i]->print();
	   dipoles_[i]->print();
	}
}
// abstract methods
void DMTPole::compute(psi::SharedMatrix D, bool transition, int i) {/* nothing to implement here */}
void DMTPole::print_header(void) const {/* nothing to implement here */}


} // EndNameSpace oepdev
