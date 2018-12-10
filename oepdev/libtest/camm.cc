#include <iostream>
#include "test.h"
#include "../lib3d/dmtp.h"

using namespace std;


double oepdev::test::Test::test_camm(void) {
  // This test is for H2O at HF/6-31* molecule
  double result = 0.0;

  // Reference CAMM values
                          /* */
  const double c_ref[3] = {-3.366373E-01,  
                            1.682078E-01,  
                            1.684295E-01};
                          /* X              Y              Z */
  const double m_ref[9] = { 1.751974E-01,  1.403032E-01, -4.528554E-02,
                            1.331952E-02,  2.019633E-02, -5.983574E-03,
                            2.305333E-02,  8.621490E-03, -3.335370E-03};
                          /* XX             XY             XZ             YY             YZ             ZZ */
  const double q_ref[18]= {-3.461364E+00, -4.200223E-02,  5.218767E-03, -3.455108E+00, -3.522603E-02, -3.585162E+00,
                           -2.246864E-01, -8.106905E-02,  1.017712E-02, -4.150175E-01, -1.340822E-02, -4.758673E-01,
                           -4.343364E-01, -4.277004E-02,  9.403507E-03, -2.190033E-01, -6.733931E-02, -4.612635E-01};
                          /* XXX            XXY            XXZ            XYY            XYZ            XZZ            YYY 
                             YYZ            YZZ            ZZZ */
  const double o_ref[30]= {-4.118915E-01,  1.322068E-01, -1.378869E-02,  2.361463E-02, -2.364673E-03,  3.349728E-02,  -4.066315E-01,
                            1.213214E-01, -8.290932E-03, -1.620214E-02,
                           -7.924737E-01,  2.163772E-01, -2.529971E-02, -1.338179E-01,  1.787315E-02, -4.012411E-02,   1.979517E-01,
                           -2.819298E-02,  3.948758E-02, -2.403612E-02, 
                            1.433252E-01, -1.065108E-01,  2.453663E-02,  1.064634E-01, -1.776971E-02,  3.281251E-02,  -8.072710E-01,
                            1.886268E-01, -8.783347E-02,  4.116774E-02};
                          /* XXXX           XXXY           XXXZ           XXYY           XXYZ           XXZZ           XYYY    
                             XYYZ           XYZZ           XZZZ           YYYY           YYYZ           YYZZ           YZZZ        
                             ZZZZ  */
  const double h_ref[45]= {-6.827831E+00,  1.096713E-01, -1.807676E-03, -2.096353E+00,  7.727420E-04, -2.072593E+00,  2.711104e-03,
                           -3.225299E-03, -5.645001E-03,  4.940954E-03, -6.783732E+00,  1.527001E-01, -2.111187E+00,  8.412184E-03, 
                           -6.225445E+00,
                            1.363178E+00, -5.627909E-01,  6.022465E-02,  5.578021E-04, -3.748843E-02, -2.031609E-01, -3.411239E-01,
                            3.780946E-02, -7.639469E-02,  3.115663E-02, -7.362549E-01, -6.815813E-02, -3.397376E-01, -3.764279E-02, 
                           -1.202959E+00,
                           -9.077266E-01, -2.140821E-01,  4.834518E-02, -6.139338E-02, -7.566905E-02, -3.484502E-01, -2.644335E-01,
                            4.190326E-02, -4.922713E-02,  2.978663E-02,  1.391233E+00, -5.221911E-01, -9.444388E-02, -1.815657E-01, 
                           -1.131796E+00};
  // Compute CAMM
  psi::timer_on("CAMM   Calculation              ");
  std::shared_ptr<DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn_);
  dmtp->compute();
  psi::timer_off("CAMM   Calculation              ");
  dmtp->charges(0)      ->print();
  dmtp->dipoles(0)      ->print();
  dmtp->quadrupoles(0)  ->print();
  dmtp->octupoles(0)    ->print();
  dmtp->hexadecapoles(0)->print();

  // Recenter to (0, 0, 0) and then back to atomic centres
  //dmtp->recenter(std::make_shared<psi::Matrix>("", dmtp->n_sites(), 3));
  //dmtp->recenter(dmtp->centres());

  std::shared_ptr<psi::Matrix> c = dmtp->charges      (0);
  std::shared_ptr<psi::Matrix> m = dmtp->dipoles      (0);
  std::shared_ptr<psi::Matrix> q = dmtp->quadrupoles  (0);
  std::shared_ptr<psi::Matrix> o = dmtp->octupoles    (0);
  std::shared_ptr<psi::Matrix> h = dmtp->hexadecapoles(0);
 
  // Accumulate errors
  for (int n=0; n<dmtp->n_sites(); ++n) {
       result += pow(c->get(n, 0) - c_ref[n] , 2.0);
       //
       result += pow(m->get(n, 0) - m_ref[3*n+0] , 2.0);
       result += pow(m->get(n, 1) - m_ref[3*n+1] , 2.0);
       result += pow(m->get(n, 2) - m_ref[3*n+2] , 2.0);
       //
       result += pow(q->get(n, 0) - q_ref[6*n+0] , 2.0);
       result += pow(q->get(n, 1) - q_ref[6*n+1] , 2.0);
       result += pow(q->get(n, 2) - q_ref[6*n+2] , 2.0);
       result += pow(q->get(n, 3) - q_ref[6*n+3] , 2.0);
       result += pow(q->get(n, 4) - q_ref[6*n+4] , 2.0);
       result += pow(q->get(n, 5) - q_ref[6*n+5] , 2.0);
       //
       result += pow(o->get(n, 0) - o_ref[10*n+0] , 2.0);
       result += pow(o->get(n, 1) - o_ref[10*n+1] , 2.0);
       result += pow(o->get(n, 2) - o_ref[10*n+2] , 2.0);
       result += pow(o->get(n, 3) - o_ref[10*n+3] , 2.0);
       result += pow(o->get(n, 4) - o_ref[10*n+4] , 2.0);
       result += pow(o->get(n, 5) - o_ref[10*n+5] , 2.0);
       result += pow(o->get(n, 6) - o_ref[10*n+6] , 2.0);
       result += pow(o->get(n, 7) - o_ref[10*n+7] , 2.0);
       result += pow(o->get(n, 8) - o_ref[10*n+8] , 2.0);
       result += pow(o->get(n, 9) - o_ref[10*n+9] , 2.0);
       //
       result += pow(h->get(n, 0) - h_ref[15*n+ 0] , 2.0);
       result += pow(h->get(n, 1) - h_ref[15*n+ 1] , 2.0);
       result += pow(h->get(n, 2) - h_ref[15*n+ 2] , 2.0);
       result += pow(h->get(n, 3) - h_ref[15*n+ 3] , 2.0);
       result += pow(h->get(n, 4) - h_ref[15*n+ 4] , 2.0);
       result += pow(h->get(n, 5) - h_ref[15*n+ 5] , 2.0);
       result += pow(h->get(n, 6) - h_ref[15*n+ 6] , 2.0);
       result += pow(h->get(n, 7) - h_ref[15*n+ 7] , 2.0);
       result += pow(h->get(n, 8) - h_ref[15*n+ 8] , 2.0);
       result += pow(h->get(n, 9) - h_ref[15*n+ 9] , 2.0);
       result += pow(h->get(n,10) - h_ref[15*n+10] , 2.0);
       result += pow(h->get(n,11) - h_ref[15*n+11] , 2.0);
       result += pow(h->get(n,12) - h_ref[15*n+12] , 2.0);
       result += pow(h->get(n,13) - h_ref[15*n+13] , 2.0);
       result += pow(h->get(n,14) - h_ref[15*n+14] , 2.0);
  }
  result = sqrt(result);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}

