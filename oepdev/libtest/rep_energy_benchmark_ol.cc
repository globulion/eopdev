#include <iostream>
#include "test.h"
#include "../libutil/wavefunction_union.h"
#include "../libsolver/solver.h"

using namespace std;

double oepdev::test::Test::test_rep_energy_benchmark_ol(void)
{
  // Reference data for H2O dimer 6-311++G** 6D
  const double ref_e = +0.00740481; // Total Exchange-Repulsion Energy [A.U.]

  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);
  wfn_union->localize_orbitals();
  wfn_union->transform_integrals();
  std::shared_ptr<oepdev::RepulsionEnergySolver> solver = std::make_shared<oepdev::RepulsionEnergySolver>(wfn_union);
  double e = solver->compute_benchmark("OTTO_LADIK"); 

  // Error
  double r = 0.0;
  r = std::abs(ref_e - e);

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Energy     = " << e << std::endl;
  std::cout << " Test result= " << r << std::endl;

  // Return
  return r;


}
