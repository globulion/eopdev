#include <iostream>
#include "test.h"
#include "../lib3d/dmtp.h"

using namespace std;


double oepdev::test::Test::test_dmtp_pot_field(void) {
  // This test is for H2O at HF/6-31* molecule
  double result = 0.0;

  // Reference CAMM values
                          /* */
  const double v_ref    = 0.11363663635731694;
  const double f_ref[3] ={-0.00702115,  0.19233044,  0.11597346};

  // Test point
  const double r[3] = {2.0, -6.0, -2.0};

  // Compute CAMM
  psi::timer_on("CAMM   Calculation              ");
  std::shared_ptr<DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn_);
  dmtp->compute();
  psi::timer_off("CAMM   Calculation              ");

  // Compute electrostatic potential at R
  std::shared_ptr<oepdev::MultipoleConvergence> potential = dmtp->potential(r[0], r[1], r[2], oepdev::MultipoleConvergence::ConvergenceLevel::R5);
  double v = potential->level(oepdev::MultipoleConvergence::ConvergenceLevel::R5)->get(0,0);

  // Compute electric field at R
  std::shared_ptr<oepdev::MultipoleConvergence> efield = dmtp->field(r[0], r[1], r[2], oepdev::MultipoleConvergence::ConvergenceLevel::R5);
  double f_efield_x = efield->level(oepdev::MultipoleConvergence::ConvergenceLevel::R5)->get(0,0);
  double f_efield_y = efield->level(oepdev::MultipoleConvergence::ConvergenceLevel::R5)->get(0,1);
  double f_efield_z = efield->level(oepdev::MultipoleConvergence::ConvergenceLevel::R5)->get(0,2);

  psi::outfile->Printf(" Point R: %14.5f %14.5f %14.5f [Bohr]\n\n"  , r[0], r[1], r[2]);

  psi::outfile->Printf(" Calculated Elc. Potential at Point R: %14.5f\n"  , v_ref);
  psi::outfile->Printf(" Reference  Elc. Potential at Point R: %14.5f\n\n", v);

  psi::outfile->Printf(" Calculated Electric Field at Point R: %14.5f %14.5f %14.5f\n", f_efield_x, f_efield_y, f_efield_z);
  psi::outfile->Printf(" Reference  Electric Field at Point R: %14.5f %14.5f %14.5f\n", f_ref[0], f_ref[1], f_ref[2]);

  // Accumulate errors
  result += pow(v_ref    - v         , 2.0);
  result += pow(f_ref[0] - f_efield_x, 2.0);
  result += pow(f_ref[1] - f_efield_y, 2.0);
  result += pow(f_ref[2] - f_efield_z, 2.0);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}

