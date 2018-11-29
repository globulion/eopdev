#include <iostream>
#include "test.h"
#include "../lib3d/dmtp.h"
#include "../libutil/wavefunction_union.h"

using namespace std;

double oepdev::test::Test::test_dmtp_energy(void) {
  // Sanity check for multimer test
  if (options_.get_str("OEPDEV_TEST_MODE") != "DIMER") 
     throw psi::PSIEXCEPTION("Monomer test mode cannot be used for this test. Set the OEPDEV_TEST_MODE to DIMER");

  // This test is for H2O dimer at HF/STO-3G (PyQuante-mod basis)
  std::shared_ptr<psi::Vector> conv_ref = std::make_shared<psi::Vector>("Convergence - reference", 5);
  conv_ref->set(0, -0.0024280728478111128);
  conv_ref->set(1, -0.0052836835555816492);
  conv_ref->set(2, -0.0074585162898145605);
  conv_ref->set(3, -0.0073177088182281823);
  conv_ref->set(4, -0.0067736731230370502);

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute CAMM's for each monomer
  std::shared_ptr<DMTPole> dmtp_1 = oepdev::DMTPole::build("CAMM", wfn_union->l_wfn(0));
  std::shared_ptr<DMTPole> dmtp_2 = oepdev::DMTPole::build("CAMM", wfn_union->l_wfn(1));
  dmtp_1->compute();
  dmtp_2->compute();

  // Compute interaction energy
  std::shared_ptr<oepdev::MultipoleConvergence> energy = dmtp_1->energy(dmtp_2, oepdev::MultipoleConvergence::ConvergenceLevel::R5);
  double conv_R1 = energy->level(oepdev::MultipoleConvergence::ConvergenceLevel::R1)->get(0,0);
  double conv_R2 = energy->level(oepdev::MultipoleConvergence::ConvergenceLevel::R2)->get(0,0);
  double conv_R3 = energy->level(oepdev::MultipoleConvergence::ConvergenceLevel::R3)->get(0,0);
  double conv_R4 = energy->level(oepdev::MultipoleConvergence::ConvergenceLevel::R4)->get(0,0);
  double conv_R5 = energy->level(oepdev::MultipoleConvergence::ConvergenceLevel::R5)->get(0,0);
  std::shared_ptr<psi::Vector> conv = std::make_shared<psi::Vector>("Convergence", 5);
  conv->set(0, conv_R1);  conv->set(1, conv_R2);  conv->set(2, conv_R3);  conv->set(3, conv_R4); conv->set(4, conv_R5);
  conv->print();
  conv_ref->print();
  conv->subtract(conv_ref);
  double result = conv->rms();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
