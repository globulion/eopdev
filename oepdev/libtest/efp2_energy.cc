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

  // Reference interaction energy: GAMESS-US Ver. 30 Sept 2017(R2)
  //                           6-31++G** (options tuned to oepdev)           6-311++G(2df) (options not tuned)
  const double ref_eint_coul =-0.0171513398;                             // -0.0142777350;  
  const double ref_eint_exrep= 0.0211206100;                             //  0.0205637907; 
  const double ref_eint_ind  =-0.0023666178;                             // -0.0024764792; 
  const double ref_eint_ct   =-0.0028507410;                             // -0.0037775538; 
  const double ref_eint_disp = 0.0; //TODO
  const double ref_eint      = ref_eint_coul + ref_eint_exrep + ref_eint_ind + ref_eint_ct + ref_eint_disp;

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute EFP2 parameters for fragment 1
  std::shared_ptr<GenEffParFactory> factory = oepdev::GenEffParFactory::build("EFP2", wfn_union->l_wfn(0), options_);
  std::shared_ptr<GenEffPar> parameters = factory->compute();

  // Create EFP2 Fragments
  oepdev::SharedGenEffFrag frag_1 = oepdev::GenEffFrag::build("Fragment 1");
  oepdev::SharedGenEffFrag frag_2 = oepdev::GenEffFrag::build("Fragment 2");

  // Set the parameters first
  std::shared_ptr<GenEffPar> par_1 = parameters->clone();
  std::shared_ptr<GenEffPar> par_2 = parameters->clone();
  frag_1->parameters["efp2"] = par_1;
  frag_2->parameters["efp2"] = par_2;

  // Set necessary sizing data
  frag_1->set_ndocc(wfn_union->l_ndocc(0));
  frag_2->set_ndocc(wfn_union->l_ndocc(1));
  frag_1->set_nbf(wfn_union->l_nbf(0));
  frag_2->set_nbf(wfn_union->l_nbf(1));

  // Set the molecule and basis sets
  frag_1->set_molecule(wfn_union->l_molecule(0));
  frag_2->set_molecule(wfn_union->l_molecule(1));
  frag_1->set_basisset("primary", wfn_union->l_primary(0)); // set_basisset has to be invoked after parameters exist
  frag_2->set_basisset("primary", wfn_union->l_primary(1));

//frag_2->basissets["primary"]->print_detail();

  frag_2->superimpose();

  // Compute interaction energy
  double eint_coul = frag_1->energy_term("EFP2:COUL", frag_2);
  double eint_exrep= frag_1->energy_term("EFP2:EXREP",frag_2);
  double eint_ind  = frag_1->energy_term("EFP2:IND" , frag_2);
  double eint_ct   = frag_1->energy_term("EFP2:CT"  , frag_2);
  double eint_disp = frag_1->energy_term("EFP2:DISP", frag_2);

  double eint = eint_coul + eint_ind + eint_exrep + eint_ct + eint_disp;

  // New test: Compute from fragmented system instance
  oepdev::SharedGenEffFrag f = frag_1->clone();

  std::vector<int> ind = {0,0};
  std::vector<oepdev::SharedGenEffFrag> bsm;
  std::vector<psi::SharedMolecule> list_mol;
  std::vector<psi::SharedBasisSet> list_prim;

  bsm.push_back(f);
  list_mol.push_back(wfn_union->l_molecule(0));
  list_mol.push_back(wfn_union->l_molecule(1));
  list_prim.push_back(wfn_union->l_primary(0));
  list_prim.push_back(wfn_union->l_primary(1));

  oepdev::SharedFragmentedSystem system = oepdev::FragmentedSystem::build(bsm, ind);
  system->set_geometry(list_mol);
  system->set_primary(list_prim);

  double eint_tt = system->compute_energy("EFP2");

  double eint_coul_t = system->compute_energy_term("EFP2:COUL" , false);
  double eint_exrep_t= system->compute_energy_term("EFP2:EXREP", false);
  double eint_ind_t  = system->compute_energy_term("EFP2:IND"  , true );
  double eint_ct_t   = system->compute_energy_term("EFP2:CT"   , false);
  double eint_disp_t = system->compute_energy_term("EFP2:DISP" , false);
  double eint_t = eint_coul_t + eint_ind_t + eint_exrep_t + eint_ct_t + eint_disp_t;


  psi::outfile->Printf("\n EFP2 Interaction Energy Components [a.u.] [kcal/mol]\n\n");
  psi::outfile->Printf("  COUL= %14.6f%14.6f\n", eint_coul , eint_coul *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  EXRP= %14.6f%14.6f\n", eint_exrep, eint_exrep*OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  IND = %14.6f%14.6f\n", eint_ind  , eint_ind  *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  CT  = %14.6f%14.6f\n", eint_ct   , eint_ct   *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  DISP= %14.6f%14.6f\n", eint_disp , eint_disp *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  TOT = %14.6f%14.6f\n", eint      , eint      *OEPDEV_AU_KcalPerMole);

  psi::outfile->Printf("\n EFP2 Interaction Energy Components [a.u.] [kcal/mol] - From FragmentedSystem\n\n");
  psi::outfile->Printf("  COUL= %14.6f%14.6f\n", eint_coul_t , eint_coul_t *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  EXRP= %14.6f%14.6f\n", eint_exrep_t, eint_exrep_t*OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  IND = %14.6f%14.6f\n", eint_ind_t  , eint_ind_t  *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  CT  = %14.6f%14.6f\n", eint_ct_t   , eint_ct_t   *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  DISP= %14.6f%14.6f\n", eint_disp_t , eint_disp_t *OEPDEV_AU_KcalPerMole);
  psi::outfile->Printf("  TOT = %14.6f%14.6f\n", eint_t      , eint_t      *OEPDEV_AU_KcalPerMole);

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
         result+= pow(eint_coul_t -eint_coul , 2.0) +
                  pow(eint_exrep_t-eint_exrep, 2.0) + 
                  pow(eint_ind_t  -eint_ind  , 2.0) +
                  pow(eint_ct_t   -eint_ct   , 2.0) +
                  pow(eint_disp_t -eint_disp , 2.0);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
