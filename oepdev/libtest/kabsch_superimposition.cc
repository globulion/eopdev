#include <iostream>
#include "test.h"
#include "../lib3d/dmtp.h"
#include "../libutil/kabsch_superimposer.h"
#include "../libutil/wavefunction_union.h"

using namespace std;

double oepdev::test::Test::test_kabsch_superimposition(void) {
  // Sanity check for multimer test
  if (options_.get_str("OEPDEV_TEST_MODE") != "DIMER") 
     throw psi::PSIEXCEPTION("Monomer test mode cannot be used for this test. Set the OEPDEV_TEST_MODE to DIMER");

  // Reference geometry
  //TODO

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Extract monomers
  psi::SharedMolecule mol_1 = wfn_union->l_molecule(0);
  psi::SharedMolecule mol_2 = wfn_union->l_molecule(1);
  
  // Superimpose
  oepdev::KabschSuperimposer sup = oepdev::KabschSuperimposer();
  sup.compute(mol_1, mol_2);
  psi::SharedMatrix xyz = sup.get_transformed();

  psi::outfile->Printf("\n RMS of superimposition = %14.6f [A.U]\n\n", sup.rms());

  // xyz->subtract(xyz_ref);
  double result = xyz->rms();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
