#include <iostream>
#include "test.h"
#include "../libutil/quambo.h"

using namespace std;


double oepdev::test::Test::test_quambo(void) {
  // This test is for asymmetrical CH4 at HF/6-31++G**
  double result = 0.0;

  // Reference values: virtual valence orbital (VVO) energies
  const double eps_vir_ref[4] = {0.342667, 0.506218, 0.524898, 0.597949};

  // Compute QUAMBO in ACBS
  psi::timer_on("QUAMBO Calculation              ");
  std::shared_ptr<oepdev::QUAMBO> solver = std::make_shared<oepdev::QUAMBO>(wfn_, true);
  solver->compute();
  psi::timer_off("QUAMBO Calculation              ");

  // Print VVOs and QUAMBOs
  psi::SharedMatrix quambo= solver->quambo("ALPHA", "ORTHOGONAL");
  psi::SharedMatrix cavir = solver->Ca_subset("AO","VIR");
  psi::SharedVector cavir_e = solver->epsilon_a_subset("MO","VIR");
  quambo->print();
  cavir->print();
  cavir_e->print();

  // Accumulate result errors
  //TODO
  if (cavir_e->dim() != 4) throw psi::PSIEXCEPTION("Something wrong in the QUAMBO test!");
  const double* ref = eps_vir_ref;
  for (int i=0; i<4; ++i) {
       result += pow(*ref - cavir_e->get(i), 2.0);
       ref++;
  }
  result = sqrt(result);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}

