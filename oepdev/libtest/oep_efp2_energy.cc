#include <iostream>
#include "test.h"
#include "../libgefp/gefp.h"
#include "../libutil/wavefunction_union.h"

using namespace std;

double oepdev::test::Test::test_oep_efp2_energy(void) {
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
  const double ref_eint_coul = 0.0; //TODO
  const double ref_eint_exrep= 0.0; //TODO
  const double ref_eint_ind  = 0.0; //TODO
  const double ref_eint_ct   = 0.0; //TODO
  const double ref_eint_disp = 0.0; //TODO
  const double ref_eint      = ref_eint_coul + ref_eint_exrep + ref_eint_ind + ref_eint_ct + ref_eint_disp;

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute EFP2 parameters for fragment 1
  std::shared_ptr<GenEffParFactory> factory = oepdev::GenEffParFactory::build("OEP-EFP2", wfn_union->l_wfn(0), options_,
                                                               wfn_union->l_auxiliary(0), wfn_union->l_intermediate(0));
  std::shared_ptr<GenEffPar> parameters = factory->compute();

  // Create EFP2 Fragments
  std::shared_ptr<GenEffFrag> frag_1 = std::make_shared<oepdev::GenEffFrag>("Fragment 1");
  std::shared_ptr<GenEffFrag> frag_2 = std::make_shared<oepdev::GenEffFrag>("Fragment 2");

  wfn_union->l_molecule(0)->print();
  wfn_union->l_molecule(1)->print();

  std::shared_ptr<GenEffPar> par_1 = parameters->clone();
  std::shared_ptr<GenEffPar> par_2 = parameters->clone();
  frag_1->parameters["efp2"] = par_1;
  frag_2->parameters["efp2"] = par_2;

  frag_1->set_molecule(wfn_union->l_molecule(0));
  frag_2->set_molecule(wfn_union->l_molecule(1));
  frag_1->set_basisset("primary", wfn_union->l_primary(0)); // set_basisset has to be invoked after parameters exist
  frag_2->set_basisset("primary", wfn_union->l_primary(1));
  frag_1->set_basisset("auxiliary", wfn_union->l_auxiliary(0)); 
  frag_2->set_basisset("auxiliary", wfn_union->l_auxiliary(1));

  frag_1->set_ndocc(wfn_union->l_ndocc(0));
  frag_2->set_ndocc(wfn_union->l_ndocc(1));
  frag_1->set_nbf(wfn_union->l_nbf(0));
  frag_2->set_nbf(wfn_union->l_nbf(1));

//frag_2->basissets["primary"]->print_detail();

  frag_2->superimpose();

  // Compute interaction energy
  double eint_coul = frag_1->energy("EFP2:COUL", frag_2);
  double eint_exrep= frag_1->energy("OEPb-EFP2:EXREP",frag_2);
  double eint_ind  = frag_1->energy("EFP2:IND" , frag_2);
  double eint_ct   = frag_1->energy("EFP2:CT"  , frag_2);
  double eint_disp = frag_1->energy("EFP2:DISP", frag_2);

  double eint = eint_coul + eint_ind + eint_exrep + eint_ct + eint_disp;

  psi::outfile->Printf("\n EFP2 Interaction Energy Components [a.u.] [kcal/mol]\n\n");
  psi::outfile->Printf("  COUL= %14.6f%14.6f\n", eint_coul , eint_coul *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  EXRP= %14.6f%14.6f\n", eint_exrep, eint_exrep*OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  IND = %14.6f%14.6f\n", eint_ind  , eint_ind  *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  CT  = %14.6f%14.6f\n", eint_ct   , eint_ct   *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  DISP= %14.6f%14.6f\n", eint_disp , eint_disp *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  TOT = %14.6f%14.6f\n", eint      , eint      *OEPDEV_AU_KcalPerMole);

  psi::outfile->Printf("\n EFP2 Interaction Energy Components from GAMESS-US [a.u.] [kcal/mol]\n\n");
  psi::outfile->Printf("  COUL= %14.6f%14.6f\n", ref_eint_coul , ref_eint_coul *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  EXRP= %14.6f%14.6f\n", ref_eint_exrep, ref_eint_exrep*OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  IND = %14.6f%14.6f\n", ref_eint_ind  , ref_eint_ind  *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  CT  = %14.6f%14.6f\n", ref_eint_ct   , ref_eint_ct   *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  DISP= %14.6f%14.6f\n", ref_eint_disp , ref_eint_disp *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  TOT = %14.6f%14.6f\n", ref_eint      , ref_eint      *OEPDEV_AU_KcalPerMole);


  double result = pow(eint_coul -ref_eint_coul , 2.0) +
                  pow(eint_exrep-ref_eint_exrep, 2.0) + 
                  pow(eint_ind  -ref_eint_ind  , 2.0) +
                  pow(eint_ct   -ref_eint_ct   , 2.0) +
                  pow(eint_disp -ref_eint_disp , 2.0);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
