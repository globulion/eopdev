#include "basis_rotation.h"


namespace oepdev{

psi::SharedMatrix r6(psi::SharedMatrix r3) {
 psi::SharedMatrix R = std::make_shared<psi::Matrix>("",6,6);
 double** R6 = R->pointer();
 double** r  =r3->pointer();
 const double r00= r[0][0];
 const double r01= r[0][1];
 const double r02= r[0][2];
 const double r10= r[1][0];
 const double r11= r[1][1];
 const double r12= r[1][2];
 const double r20= r[2][0];
 const double r21= r[2][1];
 const double r22= r[2][2];

 // for 0 - XX
 R6[0][0] = r00 * r00             ;        // XX XX
 R6[0][1] = r00 * r10 * 2.0       ;        // XX XY
 R6[0][2] = r00 * r20 * 2.0       ;        // XX XZ
 R6[0][3] = r10 * r10             ;        // XX YY
 R6[0][4] = r10 * r20 * 2.0       ;        // XX YZ
 R6[0][5] = r20 * r20             ;        // XX ZZ
 // for 1 - XY 
 R6[1][0] = r00 * r01             ;        // XY XX
 R6[1][1] =(r00 * r11 + r10 * r01);        // XY XY
 R6[1][2] =(r00 * r21 + r20 * r01);        // XY XZ
 R6[1][3] = r10 * r11             ;        // XY YY
 R6[1][4] =(r10 * r21 + r20 * r11);        // XY YZ
 R6[1][5] = r20 * r21             ;        // XY ZZ
 // for 2 - XZ
 R6[2][0] = r00 * r02             ;        // XZ XX
 R6[2][1] =(r00 * r12 + r10 * r02);        // XZ XY
 R6[2][2] =(r00 * r22 + r20 * r02);        // XZ XZ
 R6[2][3] = r10 * r12             ;        // XZ YY
 R6[2][4] =(r10 * r22 + r20 * r12);        // XZ YZ
 R6[2][5] = r20 * r22             ;        // XZ ZZ
 // for 3 - YY
 R6[3][0] = r01 * r01             ;        // YY XX
 R6[3][1] = r01 * r11 * 2.0       ;        // YY XY
 R6[3][2] = r01 * r21 * 2.0       ;        // YY XZ
 R6[3][3] = r11 * r11             ;        // YY YY
 R6[3][4] = r11 * r21 * 2.0       ;        // YY YZ
 R6[3][5] = r21 * r21             ;        // YY ZZ
 // for 4 - YZ
 R6[4][0] = r01 * r02             ;        // YZ XX
 R6[4][1] =(r01 * r12 + r11 * r02);        // YZ XY
 R6[4][2] =(r01 * r22 + r21 * r02);        // YZ XZ
 R6[4][3] = r11 * r12             ;        // YZ YY
 R6[4][4] =(r11 * r22 + r21 * r12);        // YZ YZ
 R6[4][5] = r21 * r22             ;        // YZ ZZ
 // for 5 - ZZ
 R6[5][0] = r02 * r02             ;        // ZZ XX
 R6[5][1] = r02 * r12 * 2.0       ;        // ZZ XY
 R6[5][2] = r02 * r22 * 2.0       ;        // ZZ XZ
 R6[5][3] = r12 * r12             ;        // ZZ YY
 R6[5][4] = r12 * r22 * 2.0       ;        // ZZ YZ
 R6[5][5] = r22 * r22             ;        // ZZ ZZ

 R->transpose_this();

 return R;
}

psi::SharedMatrix r10(psi::SharedMatrix r3) {
 psi::SharedMatrix R = std::make_shared<psi::Matrix>("",10,10);
 double** R10= R->pointer();
 double** r  =r3->pointer();
 const double r00= r[0][0];
 const double r01= r[0][1];
 const double r02= r[0][2];
 const double r10= r[1][0];
 const double r11= r[1][1];
 const double r12= r[1][2];
 const double r20= r[2][0];
 const double r21= r[2][1];
 const double r22= r[2][2];

 R10[0][0] = r00*r00*r00;
 R10[0][1] = r00*r00*r01;
 R10[0][2] = r00*r00*r02;
 R10[0][3] = r00*r01*r01;
 R10[0][4] = r00*r01*r02;
 R10[0][5] = r00*r02*r02;
 R10[0][6] = r01*r01*r01;
 R10[0][7] = r01*r01*r02;
 R10[0][8] = r01*r02*r02;
 R10[0][9] = r02*r02*r02;
 R10[1][0] = r10*r00*r00 + r00*r10*r00 + r00*r00*r10;
 R10[1][1] = r10*r00*r01 + r00*r10*r01 + r00*r00*r11;
 R10[1][2] = r10*r00*r02 + r00*r10*r02 + r00*r00*r12;
 R10[1][3] = r10*r01*r01 + r00*r11*r01 + r00*r01*r11;
 R10[1][4] = r10*r01*r02 + r00*r11*r02 + r00*r01*r12;
 R10[1][5] = r10*r02*r02 + r00*r12*r02 + r00*r02*r12;
 R10[1][6] = r11*r01*r01 + r01*r11*r01 + r01*r01*r11;
 R10[1][7] = r11*r01*r02 + r01*r11*r02 + r01*r01*r12;
 R10[1][8] = r11*r02*r02 + r01*r12*r02 + r01*r02*r12;
 R10[1][9] = r12*r02*r02 + r02*r12*r02 + r02*r02*r12;
 R10[2][0] = r20*r00*r00 + r00*r20*r00 + r00*r00*r20;
 R10[2][1] = r20*r00*r01 + r00*r20*r01 + r00*r00*r21;
 R10[2][2] = r20*r00*r02 + r00*r20*r02 + r00*r00*r22;
 R10[2][3] = r20*r01*r01 + r00*r21*r01 + r00*r01*r21;
 R10[2][4] = r20*r01*r02 + r00*r21*r02 + r00*r01*r22;
 R10[2][5] = r20*r02*r02 + r00*r22*r02 + r00*r02*r22;
 R10[2][6] = r21*r01*r01 + r01*r21*r01 + r01*r01*r21;
 R10[2][7] = r21*r01*r02 + r01*r21*r02 + r01*r01*r22;
 R10[2][8] = r21*r02*r02 + r01*r22*r02 + r01*r02*r22;
 R10[2][9] = r22*r02*r02 + r02*r22*r02 + r02*r02*r22;
 R10[3][0] = r00*r10*r10 + r10*r00*r10 + r10*r10*r00;
 R10[3][1] = r00*r10*r11 + r10*r00*r11 + r10*r10*r01;
 R10[3][2] = r00*r10*r12 + r10*r00*r12 + r10*r10*r02;
 R10[3][3] = r00*r11*r11 + r10*r01*r11 + r10*r11*r01;
 R10[3][4] = r00*r11*r12 + r10*r01*r12 + r10*r11*r02;
 R10[3][5] = r00*r12*r12 + r10*r02*r12 + r10*r12*r02;
 R10[3][6] = r01*r11*r11 + r11*r01*r11 + r11*r11*r01;
 R10[3][7] = r01*r11*r12 + r11*r01*r12 + r11*r11*r02;
 R10[3][8] = r01*r12*r12 + r11*r02*r12 + r11*r12*r02;
 R10[3][9] = r02*r12*r12 + r12*r02*r12 + r12*r12*r02;
 R10[4][0] = r00*r10*r20 + r00*r20*r10 + r10*r00*r20 + r10*r20*r00 + r20*r00*r10 + r20*r10*r00;
 R10[4][1] = r00*r10*r21 + r00*r20*r11 + r10*r00*r21 + r10*r20*r01 + r20*r00*r11 + r20*r10*r01;
 R10[4][2] = r00*r10*r22 + r00*r20*r12 + r10*r00*r22 + r10*r20*r02 + r20*r00*r12 + r20*r10*r02;
 R10[4][3] = r00*r11*r21 + r00*r21*r11 + r10*r01*r21 + r10*r21*r01 + r20*r01*r11 + r20*r11*r01;
 R10[4][4] = r00*r11*r22 + r00*r21*r12 + r10*r01*r22 + r10*r21*r02 + r20*r01*r12 + r20*r11*r02;
 R10[4][5] = r00*r12*r22 + r00*r22*r12 + r10*r02*r22 + r10*r22*r02 + r20*r02*r12 + r20*r12*r02;
 R10[4][6] = r01*r11*r21 + r01*r21*r11 + r11*r01*r21 + r11*r21*r01 + r21*r01*r11 + r21*r11*r01;
 R10[4][7] = r01*r11*r22 + r01*r21*r12 + r11*r01*r22 + r11*r21*r02 + r21*r01*r12 + r21*r11*r02;
 R10[4][8] = r01*r12*r22 + r01*r22*r12 + r11*r02*r22 + r11*r22*r02 + r21*r02*r12 + r21*r12*r02;
 R10[4][9] = r02*r12*r22 + r02*r22*r12 + r12*r02*r22 + r12*r22*r02 + r22*r02*r12 + r22*r12*r02;
 R10[5][0] = r00*r20*r20 + r20*r00*r20 + r20*r20*r00;
 R10[5][1] = r00*r20*r21 + r20*r00*r21 + r20*r20*r01;
 R10[5][2] = r00*r20*r22 + r20*r00*r22 + r20*r20*r02;
 R10[5][3] = r00*r21*r21 + r20*r01*r21 + r20*r21*r01;
 R10[5][4] = r00*r21*r22 + r20*r01*r22 + r20*r21*r02;
 R10[5][5] = r00*r22*r22 + r20*r02*r22 + r20*r22*r02;
 R10[5][6] = r01*r21*r21 + r21*r01*r21 + r21*r21*r01;
 R10[5][7] = r01*r21*r22 + r21*r01*r22 + r21*r21*r02;
 R10[5][8] = r01*r22*r22 + r21*r02*r22 + r21*r22*r02;
 R10[5][9] = r02*r22*r22 + r22*r02*r22 + r22*r22*r02;
 R10[6][0] = r10*r10*r10;
 R10[6][1] = r10*r10*r11;
 R10[6][2] = r10*r10*r12;
 R10[6][3] = r10*r11*r11;
 R10[6][4] = r10*r11*r12;
 R10[6][5] = r10*r12*r12;
 R10[6][6] = r11*r11*r11;
 R10[6][7] = r11*r11*r12;
 R10[6][8] = r11*r12*r12;
 R10[6][9] = r12*r12*r12;
 R10[7][0] = r20*r10*r10 + r10*r20*r10 + r10*r10*r20;
 R10[7][1] = r20*r10*r11 + r10*r20*r11 + r10*r10*r21;
 R10[7][2] = r20*r10*r12 + r10*r20*r12 + r10*r10*r22;
 R10[7][3] = r20*r11*r11 + r10*r21*r11 + r10*r11*r21;
 R10[7][4] = r20*r11*r12 + r10*r21*r12 + r10*r11*r22;
 R10[7][5] = r20*r12*r12 + r10*r22*r12 + r10*r12*r22;
 R10[7][6] = r21*r11*r11 + r11*r21*r11 + r11*r11*r21;
 R10[7][7] = r21*r11*r12 + r11*r21*r12 + r11*r11*r22;
 R10[7][8] = r21*r12*r12 + r11*r22*r12 + r11*r12*r22;
 R10[7][9] = r22*r12*r12 + r12*r22*r12 + r12*r12*r22;
 R10[8][0] = r10*r20*r20 + r20*r10*r20 + r20*r20*r10;
 R10[8][1] = r10*r20*r21 + r20*r10*r21 + r20*r20*r11;
 R10[8][2] = r10*r20*r22 + r20*r10*r22 + r20*r20*r12;
 R10[8][3] = r10*r21*r21 + r20*r11*r21 + r20*r21*r11;
 R10[8][4] = r10*r21*r22 + r20*r11*r22 + r20*r21*r12;
 R10[8][5] = r10*r22*r22 + r20*r12*r22 + r20*r22*r12;
 R10[8][6] = r11*r21*r21 + r21*r11*r21 + r21*r21*r11;
 R10[8][7] = r11*r21*r22 + r21*r11*r22 + r21*r21*r12;
 R10[8][8] = r11*r22*r22 + r21*r12*r22 + r21*r22*r12;
 R10[8][9] = r12*r22*r22 + r22*r12*r22 + r22*r22*r12;
 R10[9][0] = r20*r20*r20;
 R10[9][1] = r20*r20*r21;
 R10[9][2] = r20*r20*r22;
 R10[9][3] = r20*r21*r21;
 R10[9][4] = r20*r21*r22;
 R10[9][5] = r20*r22*r22;
 R10[9][6] = r21*r21*r21;
 R10[9][7] = r21*r21*r22;
 R10[9][8] = r21*r22*r22;
 R10[9][9] = r22*r22*r22;

 return R;
}

void populate(double** R, double** r, std::vector<int> idx_am, const int& nam) {
  const int n_p_groups = idx_am.size() / nam;
  int g_c = 0;
  for (int group=0; group<n_p_groups; ++group) {
     //int g_n = g_c + nam; -> not needed

       for (int ir=0; ir<nam; ++ir) {
            int i = idx_am[g_c+ir];
            for (int jr=0; jr<nam; ++jr) {
                 int j = idx_am[g_c+jr];
                 R[i][j] = r[ir][jr];
            }
       }
       g_c += nam;
  }
}

psi::SharedMatrix ao_rotation_matrix(psi::SharedMatrix rot, psi::SharedBasisSet bsf) {
 const int nbf = bsf->nbf();
 psi::SharedMatrix R = std::make_shared<psi::Matrix>("AO rotation matrix", nbf, nbf);
 R->identity();
 double** Rp = R->pointer();
 double** rp = rot->pointer();

 std::map<int,int> nam = {};
 if (bsf->has_puream()) {nam[0]=1;nam[1]=3;nam[2]=5;nam[3]= 7;nam[4]= 9;}
 else                   {nam[0]=1;nam[1]=3;nam[2]=6;nam[3]=10;nam[4]=15;}

 const int max_am = bsf->max_am();

 // Sanity checks
 if (max_am > 3) throw psi::PSIEXCEPTION("AO Rotaton: Agnular momentae larger than 3 are not supported!");
 if (bsf->has_puream()) throw psi::PSIEXCEPTION("AO Rotation: Sorry. Only Cartesian basis sets (puream = False) are supported at present.");

 // Indices per angular momentum
 std::vector<std::vector<int>> idx;
 for (int am=0; am<max_am+1; ++am) {
      std::vector<int> l = {};
      idx.push_back(l);
 }
 for (int i=0; i<nbf; ++i) {
      std::vector<int> l;
      int i_shell = bsf->ao_to_shell(i);
      int am      = bsf->shell(i_shell).am();
      idx[am].push_back(i);
 }
 for (int am=0; am<max_am+1; ++am) {
      if (idx[am].size() % nam[am] != 0) throw psi::PSIEXCEPTION("ERROR!!");
 }


 // Populate R matrix
 /* --- s block --- */
 // Nothing

 /* --- p block --- */
 if (max_am > 0) populate(Rp, rp, idx[1], nam[1]);

 /* --- d block --- */
 if (max_am > 1) {
    psi::SharedMatrix R6 = r6(rot);
    populate(Rp, R6->pointer(), idx[2], nam[2]);}

 /* --- f block --- */
 if (max_am > 2) {
    psi::SharedMatrix R10= r10(rot);
    populate(Rp, R10->pointer(), idx[3], nam[3]);}


 return R;
}


} // EndNameSpace oepdev
