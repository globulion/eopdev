#include <iostream>
#include "test.h"
#include "../libutil/cis.h"

using namespace std;


double oepdev::test::Test::test_cis_rhf(void) {
  // This test is for H2O at HF/6-31G* 6D molecule multiplicity singlet
  // The Gaussian 2016 input and output are given in refs
  double result = 0.0;

  // Reference energies
  const double E_ref[11] = {8.212409544202391, 9.216164293339032, 
      10.154430162329664, 10.303295256306203, 10.984749936098458,
      11.689811665332998, 12.011992061283568, 13.466813981526869,
      13.819390225017136, 14.851902189703274, 14.982747984237099};
  const int tr_J[12] {2, 5, 7, 9, 11, 12, 19, 23, 26, 27, 49, 80};
  const double tr_ref[36] = {0.0108,      0.0639,      0.2397,    // 2
                            -0.0000,     -0.0002,     -0.0007,    // 5
                             0.4872,      0.3878,     -0.1253,    // 7
                             0.3411,     -0.3985,      0.0909,    // 9
                            -0.7383,      0.8587,     -0.1957,    // 11
                            -0.6076,     -0.4863,      0.1570,    // 12
                            -0.0183,     -0.1087,     -0.4078,    // 19
                            -0.1138,     -0.0918,      0.0296,    // 23
                             0.1696,      0.1428,     -0.0457,    // 26
                             0.0089,      0.0530,      0.1986,    // 27
                            -0.2243,     -0.1794,      0.0579,    // 49
                            -0.8885,     -0.7081,      0.2287};   // 80

  // Compute CIS(RHF)
  psi::timer_on("CIS RHF Calculation             ");
  std::shared_ptr<oepdev::CISComputer> cis = oepdev::CISComputer::build("RESTRICTED", wfn_, wfn_->options());
  cis->compute();
  psi::timer_off("CIS RHF Calculation             ");

  // Extract excitation energies
  psi::SharedVector E = cis->eigenvalues();

  // Compute all transition dipole moments for 0->j transition
  // for (int i=0; i<E->dim(); ++i) cis->transition_dipole(i)->print_out();

  // Build-up error
  for (int i=0; i<11; ++i) {
       double ei = E->get(i) * 27.21138; // [EV]
       result += pow(E_ref[i] - ei, 2.0);
       psi::outfile->Printf(" Transition 0-%2d energy: %14.5f [EV]\n", i+1, ei      );
       psi::outfile->Printf( "        (g16 reference): %14.5f [EV]\n", i+1, E_ref[i]);
  }
  for (int i=0; i<12; ++i) {
       int j = tr_J[i]-1;
       psi::SharedVector tj = cis->transition_dipole(j);
       double tjx = tj->get(0); double tjx_r = tr_ref[3*i+0];
       double tjy = tj->get(1); double tjy_r = tr_ref[3*i+1];
       double tjz = tj->get(2); double tjz_r = tr_ref[3*i+2];
       result += pow(abs(tjx_r) - abs(tjx), 2.0);
       result += pow(abs(tjy_r) - abs(tjy), 2.0);
       result += pow(abs(tjz_r) - abs(tjz), 2.0);
       psi::outfile->Printf(" Transition 0-%2d  dipole moment: %14.4f %14.4f %14.4f [a.u.]\n", j+1, tjx, tjy, tjz);
       psi::outfile->Printf( "                (g16 reference): %14.4f %14.4f %14.4f [a.u.]\n",      tjx_r, tjy_r, tjz_r);
  }

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
