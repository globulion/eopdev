#include <iostream>
#include "test.h"
#include "../libutil/wavefunction_union.h"
#include "../libutil/kabsch_superimposer.h"
#include "../libutil/basis_rotation.h"

using namespace std;

double oepdev::test::Test::test_basis_rotation(void) {
  // Sanity check for multimer test
  if (options_.get_str("OEPDEV_TEST_MODE") != "DIMER") 
     throw psi::PSIEXCEPTION("Monomer test mode cannot be used for this test. Set the OEPDEV_TEST_MODE to DIMER");

  // This test is for H2O dimer at HF/6-311++G(2df) level of theory
  // Both monomers have the same internal geometry

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Print Basis Set info
  psi::SharedBasisSet primary = wfn_union->l_primary(0);
  std::map<int,char> type = {};
  type[0] = 'S'; type[1] = 'P'; type[2] = 'D'; type[3] = 'F';
  for (int i=0; i<primary->nbf(); ++i) {
       int s = primary->ao_to_shell(i);
       int am= primary->shell(s).am();
       psi::outfile->Printf(" Basis Signature %2d   --->  %c\n", i+1, type[am]);
  }
  

  // Calculate Cartesian rotation matrix first
  oepdev::KabschSuperimposer sup = oepdev::KabschSuperimposer();
  sup.compute(wfn_union->l_molecule(0), wfn_union->l_molecule(1));
  psi::SharedMatrix r = sup.rotation;

  // Calculate AO rotation matrix 
  psi::SharedMatrix R = oepdev::ao_rotation_matrix(r, wfn_union->l_primary(0)); /* normal matrices */
  psi::SharedMatrix Ri= R->clone(); Ri->invert(); Ri->transpose_this();         /* density and molecular orbitals */

  /* S matrix */
  psi::SharedMatrix S_ref = wfn_union->l_wfn(1)->S();
  psi::SharedMatrix S     = wfn_union->l_wfn(0)->S();

  /* H-core matrix */
  psi::SharedMatrix H_ref = wfn_union->l_wfn(1)->H();
  psi::SharedMatrix H     = wfn_union->l_wfn(0)->H();

  /* Fock matrix */
  psi::SharedMatrix F_ref = wfn_union->l_wfn(1)->Fa();
  psi::SharedMatrix F     = wfn_union->l_wfn(0)->Fa();

  /* Density matrix */
  psi::SharedMatrix D_ref = wfn_union->l_wfn(1)->Da();
  psi::SharedMatrix D     = wfn_union->l_wfn(0)->Da();

  /* LCAO-MO matrix */
  psi::SharedMatrix C_ref = wfn_union->l_wfn(1)->Ca();
  psi::SharedMatrix C     = wfn_union->l_wfn(0)->Ca();

  // Rotate AO matrices
  psi::SharedMatrix S_tst = psi::Matrix::triplet(R, S, R, true, false, false);
  psi::SharedMatrix H_tst = psi::Matrix::triplet(R, H, R, true, false, false);
  psi::SharedMatrix F_tst = psi::Matrix::triplet(R, F, R, true, false, false);
  psi::SharedMatrix D_tst = psi::Matrix::triplet(Ri, D, Ri, true, false, false);
  psi::SharedMatrix C_tst = psi::Matrix::doublet(Ri, C   , true, false       );

  // Test energies
  const double e_ref = H_ref->vector_dot(D_ref) + F_ref->vector_dot(D_ref);
  const double e_tst = H_tst->vector_dot(D_tst) + F_tst->vector_dot(D_tst);
  const double e     = H    ->vector_dot(D    ) + F    ->vector_dot(D    );

  D_ref->set_name("D Ref");
//D_ref->print();
  D_tst->set_name("D Com");
//D_tst->print();

  S_tst->subtract(S_ref);
  H_tst->subtract(H_ref);
  F_tst->subtract(F_ref);
  D_tst->subtract(D_ref);
//C_tst->subtract(C_ref); -> Ca might have different phases. Therefore use different rms test (below)

  double rms_s = S_tst->rms();
  double rms_h = H_tst->rms();
  double rms_f = F_tst->rms();
  double rms_d = D_tst->rms();
  double rms_c = C_tst->sum_of_squares() - C_ref->sum_of_squares();

//C_tst->subtract(C_ref); C_tst->print();
//D_tst->set_name("D DIFF");D_tst->apply_denominator(D_ref);D_tst->scale(100.0);D_tst->print();
//F_tst->set_name("F DIFF");F_tst->apply_denominator(F_ref);D_tst->scale(100.0);F_tst->print();

  psi::outfile->Printf(" RMS of Superimposition (Kabsch) = %14.3f [a.u.]\n\n", sup.rms());
  psi::outfile->Printf(" Electronic Energy (REF) = %14.6f [a.u.]\n", e_ref);
  psi::outfile->Printf(" Electronic Energy (TST) = %14.6f [a.u.]\n", e_tst);
  psi::outfile->Printf(" Electronic Energy (   ) = %14.6f [a.u.]\n", e    );
  psi::outfile->Printf(" RMS S = %14.6f\n", rms_s);
  psi::outfile->Printf(" RMS H = %14.6f\n", rms_h);
  psi::outfile->Printf(" RMS F = %14.6f\n", rms_f);
  psi::outfile->Printf(" RMS D = %14.6f\n", rms_d);
  psi::outfile->Printf(" RMS C = %14.6f\n", rms_c);


  // Accumulate errors into result
  double result = rms_s + rms_h + rms_f + rms_d + rms_c;

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
