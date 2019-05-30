#include <iostream>
#include "test.h"
#include "../libutil/wavefunction_union.h"
#include "../libsolver/solver.h"

using namespace std;

double oepdev::test::Test::test_ct_energy_benchmark_ol(void)
{
  // Reference data for  
  const double ref_e = -0.029507;



  psi::timer_on("Test: CT Energy Calculation - Otto-Ladik             ");
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union_ = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);
  std::shared_ptr<oepdev::ChargeTransferEnergySolver> ct_ene_OL = std::make_shared<oepdev::ChargeTransferEnergySolver>(wfn_union_);
//  std::shared_ptr<OEPDevSolver> build(const std::string& target, wfn_union_)
  double e = ct_ene_OL->compute_benchmark("OTTO_LADIK"); 
  psi::timer_off("Test: CT Energy Calculation - Otto-Ladik             ");


  // Error
  double r = 0.0;
  r = e - ref_e;

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r << std::endl;

  // Return
  return r;


}
