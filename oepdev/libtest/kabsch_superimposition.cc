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
  psi::SharedMatrix xyz_ref = std::make_shared<psi::Matrix>("Exact Superimposed Geometry", 5, 3);
  xyz_ref->set(0, 0, -6.22131056); xyz_ref->set(0, 1, 7.72994707); xyz_ref->set(0, 2, 2.62315247);
  xyz_ref->set(1, 0, -5.59812308); xyz_ref->set(1, 1, 5.32842843); xyz_ref->set(1, 2, 1.99991792);
  xyz_ref->set(2, 0, -6.28591633); xyz_ref->set(2, 1, 7.73552723); xyz_ref->set(2, 2, 4.49230148);
  xyz_ref->set(3, 0, -4.20409110); xyz_ref->set(3, 1, 5.42309066); xyz_ref->set(3, 2, 0.60483965);
  xyz_ref->set(4, 0, -7.20037104); xyz_ref->set(4, 1, 4.48584063); xyz_ref->set(4, 2, 1.21261005);
  // Reference rotation matrix
  psi::SharedMatrix rot_ref = std::make_shared<psi::Matrix>("Exact Rotation Matrix", 3, 3);
  rot_ref->set(0, 0, 0.45104044); rot_ref->set(0, 1, 0.84150199); rot_ref->set(0, 2, 0.29738345);
  rot_ref->set(1, 0, 0.50712679); rot_ref->set(1, 1, 0.03254751); rot_ref->set(1, 2,-0.86125668);
  rot_ref->set(2, 0,-0.73442831); rot_ref->set(2, 1, 0.53927271); rot_ref->set(2, 2,-0.41206796);
  // Reference translation vector
  psi::SharedVector tra_ref = std::make_shared<psi::Vector>("Exact Translation Vector", 3);
  tra_ref->set(0, -7.34107192); tra_ref->set(1, 12.98138065); tra_ref->set(2, 4.26911244);

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Extract monomers
  psi::SharedMolecule mol_1 = wfn_union->l_molecule(0);
  psi::SharedMolecule mol_2 = wfn_union->l_molecule(1);
  
  // Superimpose
  oepdev::KabschSuperimposer sup = oepdev::KabschSuperimposer();
  sup.compute(mol_1, mol_2); // superimposition mol_1 ---> mol_2
  psi::SharedMatrix xyz = sup.get_transformed();

  psi::outfile->Printf("\n Kabsch Superimposition\n\n");

  xyz_ref->print();
  xyz->set_name("Computed geometry"); xyz->print();

  rot_ref->print();
  sup.rotation->print();

  tra_ref->print();
  sup.translation->print();
 
  psi::outfile->Printf("\n RMS of superimposition = %14.6f [A.U]\n\n", sup.rms());
 

  xyz->subtract(xyz_ref);
  double result = xyz->rms();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
