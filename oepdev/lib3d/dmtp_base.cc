#include "dmtp.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/multipolesymmetry.h"
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
   throw psi::PSIEXCEPTION("DMTP rotation is not implemented yet.");
}
/// Superimpose the DMTP sets
void DMTPole::superimpose(psi::SharedMatrix ref_xyz, std::vector<int> suplist) 
{
   throw psi::PSIEXCEPTION("DMTP superimposition is not implemented yet.");
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
