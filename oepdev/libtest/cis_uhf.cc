#include <iostream>
#include "test.h"
#include "../libutil/cis.h"

using namespace std;


double oepdev::test::Test::test_cis_uhf(void) {
  // This test is for H2O at HF/sto-3g molecule multiplicity triplet
  double result = 0.0;

  // Compute CAMM
  psi::timer_on("CIS UHF Calculation             ");
  std::shared_ptr<oepdev::CISComputer> cis = oepdev::CISComputer::build("UNRESTRICTED", wfn_, wfn_->options());
  cis->compute();
  psi::timer_off("CIS UHF Calculation             ");

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
