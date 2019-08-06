#include <iostream>
#include "test.h"
#include "../libutil/wavefunction_union.h"
#include "../libsolver/solver.h"

using namespace std;

double oepdev::test::Test::test_ct_energy_benchmark_ol(void)
{
  // Reference data for H2O dimer  
  const double ref_e = -0.00050753;



  psi::timer_on("Test: CT Energy Calculation - Otto-Ladik             ");
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);
  wfn_union->transform_integrals();
  std::shared_ptr<oepdev::ChargeTransferEnergySolver> ct_ene_OL = std::make_shared<oepdev::ChargeTransferEnergySolver>(wfn_union);
  double e = ct_ene_OL->compute_benchmark("OTTO_LADIK"); 
  psi::timer_off("Test: CT Energy Calculation - Otto-Ladik             ");


  // Error
  double r = 0.0;
  r = ref_e - e;

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " CT Energy Otto-Ladik= " << e << std::endl;
  std::cout << " Test result= " << r << std::endl;

  // Return
  return r;


}
