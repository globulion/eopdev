#include <iostream>
#include "test.h"
#include "../libutil/cis.h"

using namespace std;


double oepdev::test::Test::test_cis_rhf(void) {
  // This test is for H2O at HF/sto-3g molecule multiplicity singlet
  double result = 0.0;

  // Compute CIS(RHF)
  psi::timer_on("CIS RHF Calculation             ");
  std::shared_ptr<oepdev::CISComputer> cis = oepdev::CISComputer::build("RESTRICTED", wfn_, wfn_->options());
  cis->compute();
  psi::timer_off("CIS RHF Calculation             ");

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
