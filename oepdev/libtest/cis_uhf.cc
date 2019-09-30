#include <iostream>
#include "test.h"
#include "../libutil/cis.h"

using namespace std;


double oepdev::test::Test::test_cis_uhf(void) {
  // This test is for H2O at HF/6-31G* 6D molecule multiplicity triplet
  // The Gaussian 2016 input and output are given in refs
  double result = 0.0;

  // Reference energies
  const double E_ref[11] = {2.155183853880288, 2.624280601617274,
       6.181243536823391, 14.442628590338225, 14.541305258786490,
      15.092987981838029, 16.511819294165495, 16.575685580676815,
      17.143246255410293, 20.266527134865147, 20.372973644109411};
  const int tr_J[12] {1, 2, 4, 5,  7,  8,  9, 10, 15, 38, 43,102};
  const double tr_ref[36] = {-1.0709,      1.2477,     -0.2845,    // 1 
                              0.0053,      0.0314,      0.1179,    // 2 
                              0.2434,      0.1913,     -0.0619,    // 4 
                              0.0361,     -0.0657,      0.0159,    // 5
                             -0.1407,     -0.1689,      0.0513,    // 7 
                             -0.4981,      0.5561,     -0.1258,    // 8 
                             -0.4913,      0.5775,     -0.1319,    // 9 
                             -0.4281,      0.5138,     -0.1177,    // 10 
                             -0.5622,     -0.4493,      0.1451,    // 15
                              0.3950,      0.3167,     -0.1022,    // 38
                              0.5229,     -0.6093,      0.1389,    // 43
                             -0.1141,     -0.0910,      0.0294};   // 102
   const double f_ref[12] = {0.1470, 0.0010, 0.0353, 0.0021, 0.0206, 0.2328, 
                             0.2487, 0.2289, 0.3044, 0.2475, 0.6576, 0.0527};


  // Compute CAMM
  psi::timer_on("CIS UHF Calculation             ");
  std::shared_ptr<oepdev::CISComputer> cis = oepdev::CISComputer::build("UNRESTRICTED", wfn_, wfn_->options(), "UHF");
  cis->compute();
  psi::timer_off("CIS UHF Calculation             ");

  // Extract excitation energies
  psi::SharedVector E = cis->eigenvalues();

  // Compute all transition dipole moments for 0->j transition
  // for (int i=0; i<E->dim(); ++i) cis->transition_dipole(i)->print_out();

  // Build-up error
  for (int i=0; i<11; ++i) {
       double ei = E->get(i) * OEPDEV_AU_EV; // [EV]
       result += pow(E_ref[i] - ei, 2.0);
       psi::outfile->Printf(" Transition 0-%2d energy: %14.5f [EV]\n", i+1, ei      );
       psi::outfile->Printf( "        (g16 reference): %14.5f [EV]\n", i+1, E_ref[i]);
  }
  for (int i=0; i<12; ++i) {
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
