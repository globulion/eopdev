#include <iostream>
#include "test.h"

using namespace std;


oepdev::test::Test::Test(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options) :
 wfn_(wfn), options_(options)
{
}
oepdev::test::Test::~Test() 
{
}
double oepdev::test::Test::run(void)
{
  double result;
  if      (options_.get_str("OEPDEV_TEST_NAME")=="BASIC"  ) result = test_basic  ();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="BASIS_ROTATION") result = test_basis_rotation();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CIS_RHF") result = test_cis_rhf();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CIS_UHF") result = test_cis_uhf();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CIS_RHF_DL") result = test_cis_rhf_dl();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CIS_UHF_DL") result = test_cis_uhf_dl();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CPHF"   ) result = test_cphf   ();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL") result = test_dmatPol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL_X") result = test_dmatPolX();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_1_1") result = test_eri_1_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_3_1") result = test_eri_3_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER") result = test_unitaryOptimizer();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER_4_2") result = test_unitaryOptimizer_4_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="SCF_PERTURB") result = test_scf_perturb();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="QUAMBO") result = test_quambo();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CAMM") result = test_camm();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMTP_POT_FIELD") result = test_dmtp_pot_field();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMTP_ENERGY") result = test_dmtp_energy();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="EFP2_ENERGY") result = test_efp2_energy();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="OEP_EFP2_ENERGY") result = test_oep_efp2_energy();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="KABSCH_SUPERIMPOSITION") result = test_kabsch_superimposition();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMTP_SUPERIMPOSITION") result = test_dmtp_superimposition();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ESP_SOLVER") result = test_esp_solver();  
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CUSTOM") result = test_custom();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="POINTS_COLLECTION3D") result = test_points_collection3d();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CT_ENERGY_BENCHMARK_OL") result = test_ct_energy_benchmark_ol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CT_ENERGY_OEP_BASED_OL") result = test_ct_energy_oep_based_ol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_BENCHMARK_HS") result = test_rep_energy_benchmark_hs();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_BENCHMARK_DDS") result = test_rep_energy_benchmark_dds();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_BENCHMARK_MURRELL_ETAL") result = test_rep_energy_benchmark_murrell_etal();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_OEP_BASED_MURRELL_ETAL") result = test_rep_energy_oep_based_murrell_etal();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_BENCHMARK_OL") result = test_rep_energy_benchmark_ol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="REP_ENERGY_BENCHMARK_EFP2") result = test_rep_energy_benchmark_efp2();
  else throw psi::PSIEXCEPTION("Incorrect test name specified!");
  return result;
}
