#include <iostream>
#include "test.h"
#include "../libutil/quambo.h"

using namespace std;


double oepdev::test::Test::test_quambo(void) {
  // This test is for H2O at HF/6-31* molecule
  double result = 0.0;

  // Reference QUAMBO values
  //TODO

  // Compute QUAMBO in ACBS
  psi::timer_on("QUAMBO Calculation              ");
  std::shared_ptr<oepdev::QUAMBO> solver = std::make_shared<oepdev::QUAMBO>(wfn_, true);
  solver->compute();
  psi::timer_off("QUAMBO Calculation              ");
  result = sqrt(result);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}

