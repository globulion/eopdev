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

  // Reference interaction energy
  const double eint_ref = 0.0;

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute EFP2 parameters for fragment 1
  std::shared_ptr<GenEffParFactory> factory = oepdev::GenEffParFactory::build("EFP2", wfn_union->l_wfn(0), options_);
  std::shared_ptr<GenEffPar> parameters = factory->compute();

  // Create EFP2 Fragments
  std::shared_ptr<GenEffFrag> frag_1 = std::make_shared<oepdev::GenEffFrag>("Fragment 1");
  std::shared_ptr<GenEffFrag> frag_2 = std::make_shared<oepdev::GenEffFrag>("Fragment 2");

  frag_1->set_molecule(wfn_union->l_molecule(0));
  frag_2->set_molecule(wfn_union->l_molecule(1));

  frag_1->parameters["efp2"] = parameters->clone();
  frag_2->parameters["efp2"] = parameters->clone();
  frag_2->superimpose();

  frag_1->basissets["primary"] = wfn_union->l_primary(0);
  frag_2->basissets["primary"] = wfn_union->l_primary(1);

  // Compute interaction energy
  double eint = frag_1->energy("EFP2", frag_2);
  psi::outfile->Printf("\n EFP2 Interaction Energy= %14.6f [kcal/mol]\n\n", eint*OEPDEV_AU_KcalPerMole);

  double result = eint - eint_ref;

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
