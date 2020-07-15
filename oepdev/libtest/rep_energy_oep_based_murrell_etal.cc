#include <iostream>
#include "test.h"
#include "../libutil/wavefunction_union.h"
#include "../libsolver/solver.h"

using namespace std;

double oepdev::test::Test::test_rep_energy_oep_based_murrell_etal(void)
{
  // Reference data for H2O dimer 6-311++G** 6D
  const double ref_e_esp = +0.01123609; // Total Exchange-Repulsion Energy [A.U.] -> before change of orbital management in WfnUnn
  const double ref_e_camm= +0.00910510; // Total Exchange-Repulsion Energy [A.U.] -> before change of orbital management in WfnUnn

  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);
  wfn_union->localize_orbitals();
  wfn_union->transform_integrals();
  std::shared_ptr<oepdev::RepulsionEnergySolver> solver = std::make_shared<oepdev::RepulsionEnergySolver>(wfn_union);
  double e_esp = solver->compute_oep_based("MURRELL_ETAL_GDF_ESP"); 
  double e_camm= solver->compute_oep_based("MURRELL_ETAL_GDF_CAMM");

  // Error
  double r = 0.0;
  r = std::abs(ref_e_esp - e_esp) + std::abs(ref_e_camm - e_camm);

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Energy ESP = " << e_esp << std::endl;
  std::cout << " Energy CAMM= " << e_camm<< std::endl;
  std::cout << " Test result= " << r << std::endl;

  // Return
  return r;


}
