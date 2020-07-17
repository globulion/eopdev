#include <iostream>
#include "test.h"
#include "../libutil/quambo.h"

using namespace std;


double oepdev::test::Test::test_quambo(void) {
  // This test is for asymmetrical CH4 at HF/6-31++G**
  double result = 0.0;

  // Reference values: virtual valence orbital (VVO) energies: reference is Python version of the code in gefp.basis.quambo
  const double eps_vir_ref_acbs[4] = {0.342667, 0.506218, 0.524898, 0.597949};
  const double eps_vir_ref_mcbs[4] = {0.347716, 0.509791, 0.528265, 0.600252};

  // Compute QUAMBO in ACBS
  psi::timer_on("QUAMBO Calculation              ");
  std::shared_ptr<oepdev::QUAMBO> solver = std::make_shared<oepdev::QUAMBO>(wfn_, options_.get_bool("QUAMBO_ACBS"));
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
  if (cavir_e->dim() != 4) throw psi::PSIEXCEPTION("Something wrong in the QUAMBO test!");
  const double* ref;
  ref = options_.get_bool("QUAMBO_ACBS") ? eps_vir_ref_acbs : eps_vir_ref_mcbs;
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

