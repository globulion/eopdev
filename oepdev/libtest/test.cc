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
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CPHF"   ) result = test_cphf   ();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL") result = test_dmatPol();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMATPOL_X") result = test_dmatPolX();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_1_1") result = test_eri_1_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_2_2") result = test_eri_2_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="ERI_3_1") result = test_eri_3_1();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER") result = test_unitaryOptimizer();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="UNITARY_OPTIMIZER_4_2") result = test_unitaryOptimizer_4_2();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="SCF_PERTURB") result = test_scf_perturb();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CAMM") result = test_camm();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="DMTP_ENERGY") result = test_dmtp_energy();
  else if (options_.get_str("OEPDEV_TEST_NAME")=="CUSTOM") result = test_custom();
  else throw psi::PSIEXCEPTION("Incorrect test name specified!");
  return result;
}
