#include <iostream>
#include "test.h"
#include "../libutil/cis.h"

using namespace std;


double oepdev::test::Test::test_cis_rhf_dl(void) {
  // This test is for H2O at HF/6-31G* 6D molecule multiplicity singlet
  // The Gaussian 2016 input and output are given in refs
  double result = 0.0;

  // Reference energies
  const double E_ref[11] = {8.212409544202391, 9.216164293339032, 
      10.154430162329664, 10.303295256306203, 10.984749936098458,
      11.689811665332998, 12.011992061283568, 13.466813981526869,
      13.819390225017136, 14.851902189703274, 14.982747984237099};
  const int tr_J[5] {2, 5, 7, 9, 11};
  const double tr_ref[15] = {0.0108,      0.0639,      0.2397,    // 2
                            -0.0000,     -0.0002,     -0.0007,    // 5
                             0.4872,      0.3878,     -0.1253,    // 7
                             0.3411,     -0.3985,      0.0909,    // 9
                            -0.7383,      0.8587,     -0.1957};   // 11
  const double f_ref[5] = {0.0139, 0.0000, 0.1187, 0.0960, 0.4848};

  // Compute CIS(RHF)
  psi::timer_on("CIS RHF Calculation             ");
  std::shared_ptr<oepdev::CISComputer> cis = oepdev::CISComputer::build("RESTRICTED", wfn_, wfn_->options(), "RHF");
  cis->compute();
  psi::timer_off("CIS RHF Calculation             ");

  // Extract excitation energies
  psi::SharedVector E = cis->eigenvalues();

  // Compute all transition dipole moments for 0->j transition
  //for (int i=0; i<E->dim(); ++i) cis->transition_dipole(i)->print_out();

  // Build-up error
  for (int i=0; i<11; ++i) {
       double ei = E->get(i) * OEPDEV_AU_EV; // [EV]
       result += pow(E_ref[i] - ei, 2.0);
       psi::outfile->Printf(" Transition 0-%2d energy: %14.5f [EV]\n", i+1, ei      );
       psi::outfile->Printf( "        (g16 reference): %14.5f [EV]\n", i+1, E_ref[i]);
  }
  for (int i=0; i<5; ++i) {
       int j = tr_J[i]-1;
       psi::SharedVector tj = cis->transition_dipole(j);
       double f = cis->oscillator_strength(j);
       double tjx = tj->get(0); double tjx_r = tr_ref[3*i+0];
       double tjy = tj->get(1); double tjy_r = tr_ref[3*i+1];
       double tjz = tj->get(2); double tjz_r = tr_ref[3*i+2];
       result += pow(abs(tjx_r) - abs(tjx), 2.0);
       result += pow(abs(tjy_r) - abs(tjy), 2.0);
       result += pow(abs(tjz_r) - abs(tjz), 2.0);
       result += pow(abs(f)     - abs(f_ref[i]), 2.0);
       psi::outfile->Printf(" Transition 0-%2d  dipole moment: %14.4f %14.4f %14.4f [a.u.] Osc=%14.4f\n", j+1, 
                                                                 tjx, tjy, tjz, f);
       psi::outfile->Printf( "                (g16 reference): %14.4f %14.4f %14.4f [a.u.] Osc=%14.4f\n",      
                                                                 tjx_r, tjy_r, tjz_r, f_ref[i]);
  }

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
