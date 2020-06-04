#include <iostream>
#include "test.h"
#include "../libgefp/gefp.h"
#include "../libutil/wavefunction_union.h"

using namespace std;

double oepdev::test::Test::test_efp2_energy(void) {
  /* 
     This test is between two interacting molecules whose internal coordinates are identical.
     Therefore only one EFP2 fragment factory is needed. Second parameters are obtained
     through superimposition.

     Basis sets are taken from the WavefunctionUnion object.
   */

  // Sanity check for multimer test
  if (options_.get_str("OEPDEV_TEST_MODE") != "DIMER") 
     throw psi::PSIEXCEPTION("Monomer test mode cannot be used for this test. Set the OEPDEV_TEST_MODE to DIMER");

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute EFP2 parameters for fragment 1
  std::shared_ptr<GenEffParFactory> factory = oepdev::GenEffParFactory::build("EFP2", wfn_union->l_wfn(0), options_);
  std::shared_ptr<GenEffPar> parameters = factory->compute();

  // Compute interaction energy
  //TODO

  double result = 0.0;

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
