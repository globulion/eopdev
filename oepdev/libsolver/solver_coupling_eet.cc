#include "solver.h"
#include "psi4/libpsi4util/process.h"
#include <utility>
#include "../libutil/integrals_iter.h"
#include "../libutil/util.h"


using namespace std;
using namespace psi;
using namespace oepdev;

using SharedMTPConv = std::shared_ptr<oepdev::MultipoleConvergence>;
using SharedDMTPole = std::shared_ptr<oepdev::DMTPole>;

EETCouplingSolver::EETCouplingSolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Benchmarks
  methods_benchmark_.push_back("FUJIMOTO_TI_CIS");

  // OEP-based
  methods_oepBased_ .push_back("FUJIMOTO_TI_CIS");
}

EETCouplingSolver::~EETCouplingSolver() {}

double EETCouplingSolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "FUJIMOTO_TI_CIS")  e = compute_oep_based_fujimoto_ti_cis();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for EET coupling calculations!\n");
  }
  return e;
}

double EETCouplingSolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" ||
           method == "FUJIMOTO_TI_CIS") e = compute_benchmark_fujimoto_ti_cis();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for EET coupling calculations!\n");
  }
  return e;
}

void EETCouplingSolver::determine_electronic_state(std::shared_ptr<CISComputer> cis, int& I) {
 if (I<1) {
   int count = 1;
   const double ft = options_.get_double("OSCILLATOR_STRENGTH_THRESHOLD");
   for (int i=0; i<cis->nstates(); ++i) {
        if (cis->oscillator_strength(i) > ft) {
            if (count == -I) {
                I = i+1;
                break;
            }
            count += 1;
        }
   }
 } 
 I -= 1; // transform to C++ indexing (from 0)
}
double EETCouplingSolver::compute_benchmark_fujimoto_ti_cis() { //TODO
  double e = 0.0;
  psi::timer_on("Solver EET TI/CIS               ");

  // ---> Allocate <--- //
  int nbf   = wfn_union_->basisset()->nbf();
  int nbf_A = wfn_union_->l_nbf(0);
  int nbf_B = wfn_union_->l_nbf(1);

  SharedMatrix Sao_AB = std::make_shared<psi::Matrix>("Sao(A,B)" , nbf_A, nbf_B);
  SharedMatrix Tao_AB = std::make_shared<psi::Matrix>("Tao(A,B)" , nbf_A, nbf_B);
  SharedMatrix VaoB_AB = std::make_shared<psi::Matrix>("VaoB(A,B)" , nbf_A, nbf_B);
  SharedMatrix VaoA_BA = std::make_shared<psi::Matrix>("VaoA(B,A)" , nbf_B, nbf_A);
  SharedMatrix VaoB_AA = std::make_shared<psi::Matrix>("VaoB(A,A)" , nbf_A, nbf_A);
  SharedMatrix VaoA_BB = std::make_shared<psi::Matrix>("VaoA(B,B)" , nbf_B, nbf_B);

  psi::IntegralFactory fact_AAAA(wfn_union_->l_primary(0), wfn_union_->l_primary(0), 
                                 wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_BBBB(wfn_union_->l_primary(1), wfn_union_->l_primary(1), 
                                 wfn_union_->l_primary(1), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_ABAB(wfn_union_->l_primary(0), wfn_union_->l_primary(1), 
                                 wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_BABA(wfn_union_->l_primary(1), wfn_union_->l_primary(0), 
                                 wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_ABBB(wfn_union_->l_primary(0), wfn_union_->l_primary(1), 
                                 wfn_union_->l_primary(1), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_AAAB(wfn_union_->l_primary(0), wfn_union_->l_primary(0), 
                                 wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_AABB(wfn_union_->l_primary(0), wfn_union_->l_primary(0), 
                                 wfn_union_->l_primary(1), wfn_union_->l_primary(1));


  // Compute one-electron integrals
  std::shared_ptr<psi::OneBodyAOInt> oneInt, ovlInt(fact_ABAB.ao_overlap()), kinInt(fact_ABAB.ao_kinetic());
  std::shared_ptr<psi::PotentialInt> potInt_B_AB= std::make_shared<psi::PotentialInt>(fact_ABAB.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_A_BA= std::make_shared<psi::PotentialInt>(fact_BABA.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_B_AA= std::make_shared<psi::PotentialInt>(fact_AAAA.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_A_BB= std::make_shared<psi::PotentialInt>(fact_BBBB.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(1));


  SharedMatrix Zxyz_A = std::make_shared<psi::Matrix>(potInt_B_AB->charge_field());
  SharedMatrix Zxyz_B = std::make_shared<psi::Matrix>(potInt_A_BA->charge_field());

  potInt_B_AB->set_charge_field(Zxyz_B);
  potInt_A_BA->set_charge_field(Zxyz_A);
  potInt_B_AA->set_charge_field(Zxyz_B);
  potInt_A_BB->set_charge_field(Zxyz_A);
  oneInt = potInt_B_AB;
  oneInt->compute(VaoB_AB);
  oneInt = potInt_A_BA;
  oneInt->compute(VaoA_BA);
  oneInt = potInt_A_BB;
  oneInt->compute(VaoA_BB);
  oneInt = potInt_B_AA;
  oneInt->compute(VaoB_AA);
  ovlInt->compute(Sao_AB);
  kinInt->compute(Tao_AB);

  // Fock matrix of entire system in AB space: 1-electron contribution
  VaoA_BA->transpose_this();
  SharedMatrix Fao_AB = Tao_AB->clone();
  Fao_AB->add(VaoB_AB);
  Fao_AB->add(VaoA_BA);


  // Auxiliary data
  int I = options_.get_int("EXCITED_STATE_A");
  int J = options_.get_int("EXCITED_STATE_B");

  const int Ne = 2.0 * (wfn_union_->l_wfn(0)->nalpha() + wfn_union_->l_wfn(1)->nalpha());
  const int homo_A = wfn_union_->l_wfn(0)->nalpha() - 1;
  const int homo_B = wfn_union_->l_wfn(1)->nalpha() - 1;
  const int lumo_A = 0; //wfn_union_->l_wfn(0)->nalpha();
  const int lumo_B = 0; //wfn_union_->l_wfn(1)->nalpha();
  const double na_A = (double)wfn_union_->l_wfn(0)->nalpha();
  const double na_B = (double)wfn_union_->l_wfn(1)->nalpha();
  const double na_AB= na_A + na_B;

  // 
  SharedMatrix Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  SharedMatrix Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");
  SharedMatrix Ca_vir_A = wfn_union_->l_wfn(0)->Ca_subset("AO","VIR");
  SharedMatrix Ca_vir_B = wfn_union_->l_wfn(1)->Ca_subset("AO","VIR");
  //
  SharedVector eps_a_occ_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","OCC");
  SharedVector eps_a_occ_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","OCC");
  SharedVector eps_a_vir_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","VIR");
  SharedVector eps_a_vir_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","VIR");
  // 
  SharedMatrix Da_A = wfn_union_->l_wfn(0)->Da();
  SharedMatrix Da_B = wfn_union_->l_wfn(1)->Da();

  // Density matrices (AO)
  SharedMatrix Peg_A, Peg_B, Pe_A, Pe_B;
  // TrCAMM's
  SharedDMTPole trcamm_A, trcamm_B;


  // Obtain CIS HOMO-LUMO amplitudes and unperturbed excitation energies
  psi::outfile->Printf(" --> Running CIS calculations on the monomers <--\n");
  double t_A, t_B, E_ex_A, E_ex_B;
  const bool symm = options_.get_bool("TrCAMM_SYMMETRIZE");
  {
      std::shared_ptr<CISComputer> cis_A = CISComputer::build("RESTRICTED", wfn_union_->l_wfn(0), options_); 
      cis_A->compute();
      this->determine_electronic_state(cis_A, I);            // Excited state ID in C++ convention
      E_ex_A = cis_A->eigenvalues()->get(I);                 // Excitation energy wrn ground state
      t_A = cis_A->U_homo_lumo(I).first * sqrt(na_A/na_AB);  // CIS amplitude scaled to the dimer
      trcamm_A = cis_A->trcamm(I, symm);                     // TrCAMM moments
      Pe_A  = cis_A->Da_ao(I); Pe_A ->add(cis_A->Db_ao(I));  // Excited state bond order matrices of monomers
      Peg_A = cis_A->Ta_ao(I); Peg_A->add(cis_A->Tb_ao(I));  // Transition density matrices of monomers
      psi::outfile->Printf("     State I= %2d, f= %9.6f [a.u.] E= %9.3f [EV] t(H->L)= %9.6f\n", 
                                 I+1, cis_A->oscillator_strength(I), E_ex_A*OEPDEV_AU_EV, t_A);
  }{
      std::shared_ptr<CISComputer> cis_B = CISComputer::build("RESTRICTED", wfn_union_->l_wfn(1), options_); 
      cis_B->compute();
      this->determine_electronic_state(cis_B, J);
      E_ex_B = cis_B->eigenvalues()->get(J);
      t_B = cis_B->U_homo_lumo(J).first * sqrt(na_B/na_AB);
      trcamm_B = cis_B->trcamm(J, symm);
      Pe_B  = cis_B->Da_ao(J); Pe_B ->add(cis_B->Db_ao(J));
      Peg_B = cis_B->Ta_ao(J); Peg_B->add(cis_B->Tb_ao(J));
      psi::outfile->Printf("     State J= %2d, f= %9.6f [a.u.] E= %9.3f [EV] t(H->L)= %9.6f\n", 
                                 J+1, cis_B->oscillator_strength(J), E_ex_B*OEPDEV_AU_EV, t_B);
  }
  psi::outfile->Printf("\n");

  // Check the populations of states
  if (wfn_union_->options().get_int("PRINT") > -1) {
     psi::SharedMatrix Sao_AA = std::make_shared<psi::Matrix>("", nbf_A, nbf_A);
     psi::SharedMatrix Sao_BB = std::make_shared<psi::Matrix>("", nbf_B, nbf_B);
     std::shared_ptr<psi::OneBodyAOInt> sAA(fact_AAAA.ao_overlap()); sAA->compute(Sao_AA);
     std::shared_ptr<psi::OneBodyAOInt> sBB(fact_BBBB.ao_overlap()); sBB->compute(Sao_BB);
     const double na_g1 = psi::Matrix::doublet(Da_A, Sao_AA)->trace();
     const double na_g2 = psi::Matrix::doublet(Da_B, Sao_BB)->trace();
     const double n_e1 = psi::Matrix::doublet(Pe_A, Sao_AA)->trace();
     const double n_e2 = psi::Matrix::doublet(Pe_B, Sao_BB)->trace();
     const double n_eg1 = psi::Matrix::doublet(Peg_A, Sao_AA)->trace();
     const double n_eg2 = psi::Matrix::doublet(Peg_B, Sao_BB)->trace();
     //
     psi::outfile->Printf(" ===> Populations of Electronic states <===\n\n");
     psi::outfile->Printf("         Total electron number\n");
     psi::outfile->Printf("    g(1) %13.5f\n", na_g1*2.0);
     psi::outfile->Printf("    g(2) %13.5f\n", na_g2*2.0);
     psi::outfile->Printf("    e(1) %13.5f\n", n_e1);
     psi::outfile->Printf("    e(2) %13.5f\n", n_e2);
     psi::outfile->Printf("    eg(1)%13.5f\n", n_eg1);
     psi::outfile->Printf("    eg(2)%13.5f\n", n_eg2);
     psi::outfile->Printf("\n");
  }

  // Change the phases of the SCF HOMO and LUMO orbitals of monomers
  //if (options_.get_int("PHASE_HOMO_A") == -1) Ca_occ_A->scale_column(0, homo_A, -1.0);
  //if (options_.get_int("PHASE_HOMO_B") == -1) Ca_occ_B->scale_column(0, homo_B, -1.0);
  //if (options_.get_int("PHASE_LUMO_A") == -1) Ca_vir_A->scale_column(0, lumo_A, -1.0);
  //if (options_.get_int("PHASE_LUMO_B") == -1) Ca_vir_B->scale_column(0, lumo_B, -1.0);

  // [1.1] Compute Overlap integrals between basis functions

  // S12
  SharedMatrix PSP = psi::Matrix::triplet(Peg_A, Sao_AB, Peg_B, false, false, false);
  double S12 = psi::Matrix::doublet(PSP, Sao_AB, false, true)->trace();
  S12*= -(1.0)/(double)Ne;
  PSP.reset();

  // S13
  double shh = 0;
  for (int i=0; i<nbf_A; ++i) {
       for (int j=0; j<nbf_B; ++j) {
            shh += Sao_AB->get(i,j) * Ca_occ_A->get(i,homo_A) * Ca_occ_B->get(j,homo_B);
       }
  }
  double sll = 0;
  for (int i=0; i<nbf_A; ++i) {
       for (int j=0; j<nbf_B; ++j) {
            sll += Sao_AB->get(i,j) * Ca_vir_A->get(i,0) * Ca_vir_B->get(j,0);
       }
  }
  double S13 = -t_A * sll/(double)Ne;
  // S42
  double S42 = -t_B * sll/(double)Ne;
  // S14
  double S14 = +t_A * shh/(double)Ne;
  // S32
  double S32 = +t_B * shh/(double)Ne;
  // S34
  double S34 = -sll * shh/(double)Ne;

  if (wfn_union_->options().get_int("PRINT") > -1) {
     psi::outfile->Printf(" ===> Overlap matrix between basis states <===\n\n");
     psi::outfile->Printf("         1.       2.       3.       4.\n");
     psi::outfile->Printf("  1. %9.4f  %9.4f  %9.4f  %9.4f\n", 1.0, S12, S13, S14);
     psi::outfile->Printf("  2. %9.4f  %9.4f  %9.4f  %9.4f\n", S12, 1.0, S32, S42);
     psi::outfile->Printf("  3. %9.4f  %9.4f  %9.4f  %9.4f\n", S13, S32, 1.0, S34);
     psi::outfile->Printf("  4. %9.4f  %9.4f  %9.4f  %9.4f\n", S14, S42, S34, 1.0);
     psi::outfile->Printf("\n");
  }


  // Compute Hamiltonian diagonal elements
  double E01= E_ex_A;
  double E02= E_ex_B;
  double E03=-eps_a_occ_A->get(homo_A) + eps_a_vir_B->get(lumo_B);
  double E04= eps_a_vir_A->get(lumo_A) - eps_a_occ_B->get(homo_B);

  double E1 = E01;
  double E2 = E02;
  double E3 = E03;
  double E4 = E04;

  // Compute direct coupling constants
  double V0_Coul = 0.0; 
  double V0_Exch = 0.0;

  // Compute indirect coupling constant
  double V0_ET1 = 0.0;
  double V0_ET2 = 0.0;
  double V0_HT1 = 0.0;
  double V0_HT2 = 0.0;
  double V0_CT  = 0.0;
  double V0_CT_M= 0.0; // Mulliken-approximated V0_CT

  // ===> Accumulate (ab|cd) contributions from AO ERI's <=== //
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aaaa = oepdev::ShellCombinationsIterator::build(fact_AAAA, "ALL", 4);
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_bbbb = oepdev::ShellCombinationsIterator::build(fact_BBBB, "ALL", 4);
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aaab = oepdev::ShellCombinationsIterator::build(fact_AAAB, "ALL", 4);
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_abbb = oepdev::ShellCombinationsIterator::build(fact_ABBB, "ALL", 4);
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_abab = oepdev::ShellCombinationsIterator::build(fact_ABAB, "ALL", 4);
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aabb = oepdev::ShellCombinationsIterator::build(fact_AABB, "ALL", 4);

  std::shared_ptr<psi::TwoBodyAOInt> t_aaaa(fact_AAAA.eri()); const double* b_aaaa = t_aaaa->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_bbbb(fact_BBBB.eri()); const double* b_bbbb = t_bbbb->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_aaab(fact_AAAB.eri()); const double* b_aaab = t_aaab->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_abbb(fact_ABBB.eri()); const double* b_abbb = t_abbb->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_abab(fact_ABAB.eri()); const double* b_abab = t_abab->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_aabb(fact_AABB.eri()); const double* b_aabb = t_aabb->buffer();

  double** dA = Da_A->pointer();
  double** dB = Da_B->pointer();
  double** pA = Pe_A->pointer();
  double** pB = Pe_B->pointer();
  double** pAt= Peg_A->pointer();
  double** pBt= Peg_B->pointer();
  double** f  = Fao_AB->pointer();
  double** coA= Ca_occ_A->pointer();
  double** coB= Ca_occ_B->pointer();
  double** cvA= Ca_vir_A->pointer();
  double** cvB= Ca_vir_B->pointer();

  psi::SharedMatrix Smo_oAoB = psi::Matrix::triplet(Ca_occ_A, Sao_AB, Ca_occ_B, true, false, false);
  psi::SharedMatrix Smo_vAvB = psi::Matrix::triplet(Ca_vir_A, Sao_AB, Ca_vir_B, true, false, false);
  psi::SharedMatrix Smo_oAvB = psi::Matrix::triplet(Ca_occ_A, Sao_AB, Ca_vir_B, true, false, false);
  psi::SharedMatrix Smo_vAoB = psi::Matrix::triplet(Ca_vir_A, Sao_AB, Ca_occ_B, true, false, false);
  double s_HL_AB = Smo_oAvB->get(homo_A, lumo_B);
  double s_LH_AB = Smo_vAoB->get(lumo_A, homo_B);
  double s_HH_AB = Smo_oAoB->get(homo_A, homo_B);
  double s_LL_AB = Smo_vAvB->get(lumo_A, lumo_B);
  double Q1 = s_HL_AB * s_LH_AB / 2.000;
  double Q2 =-s_HH_AB * s_LL_AB / 4.000;
  double Q3 = Q1 + Q2;

  // ----> (AA|AB) : ET1, HT1, Fock <---- // 
  for (s_aaab->first(); s_aaab->is_done() == false; s_aaab->next()) 
  {
    s_aaab->compute_shell(t_aaab);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_aaab->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // A
         int j = ii->j(); // A
         int k = ii->k(); // A
         int l = ii->l(); // B
         double eri = b_aaab[ii->index()];

         // Fock matrix
         f[k][l] += 2.0 * dA[i][j] * eri;
         f[i][l] -=       dA[k][j] * eri;

         // V0_ET1
         V0_ET1 += 2.0 * eri * coA[j][homo_A] * cvA[i][lumo_A] * coA[k][homo_A] * cvB[l][lumo_B];
         V0_ET1 -=       eri * coA[j][homo_A] * cvA[k][lumo_A] * coA[i][homo_A] * cvB[l][lumo_B];

         // V0_HT1
         V0_HT1 += 2.0 * eri * cvA[j][lumo_A] * coA[i][homo_A] * cvA[k][lumo_A] * coB[l][homo_B];
         V0_HT1 -=       eri * cvA[j][lumo_A] * coA[k][homo_A] * cvA[i][lumo_A] * coB[l][homo_B];
    }
  }
 
  // ----> (AB|BB) : ET2, HT2, Fock <---- // 
  for (s_abbb->first(); s_abbb->is_done() == false; s_abbb->next()) 
  {
    s_abbb->compute_shell(t_abbb);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_abbb->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // A
         int j = ii->j(); // B
         int k = ii->k(); // B
         int l = ii->l(); // B
         double eri = b_abbb[ii->index()];
       //psi::outfile->Printf("(A%2d B%2d | B%2d B%2d) = %14.6f\n", i, j, k, l, eri);

         // Fock matrix
         f[i][j] += 2.0 * dB[k][l] * eri;
         f[i][l] -=       dB[k][j] * eri;

         // V0_ET2
         V0_ET2 += 2.0 * eri * cvA[i][lumo_A] * coB[j][homo_B] * cvB[l][lumo_B] * coB[k][homo_B];
         V0_ET2 -=       eri * cvA[i][lumo_A] * coB[l][homo_B] * cvB[j][lumo_B] * coB[k][homo_B];

         // V0_HT2
         V0_HT2 += 2.0 * eri * coA[i][homo_A] * cvB[j][lumo_B] * coB[l][homo_B] * cvB[k][lumo_B];
         V0_HT2 -=       eri * coA[i][homo_A] * cvB[l][lumo_B] * coB[j][homo_B] * cvB[k][lumo_B];

    }
  }

  // ----> (AB|AB) : V0_Exch, CT, E1, E2 <---- //
  for (s_abab->first(); s_abab->is_done() == false; s_abab->next()) 
  {
    s_abab->compute_shell(t_abab);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_abab->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // A
         int j = ii->j(); // B
         int k = ii->k(); // A
         int l = ii->l(); // B
         double eri = b_abab[ii->index()];

         // V0_Exch
         V0_Exch -= 0.5 * pAt[k][i] * pBt[j][l] * eri;

         // V0_CT
         V0_CT  += 2.0 * eri * coA[i][homo_A] * cvB[j][lumo_B] * cvA[k][lumo_A] * coB[l][homo_B];
         V0_CT  -=       eri * coA[i][homo_A] * cvB[l][lumo_B] * cvA[k][lumo_A] * coB[j][homo_B];

         // E1
         E1 -= eri * (pA[k][i] - 2.0 * dA[k][i]) * dB[j][l];
         // E2
         E2 -= eri * (pB[l][j] - 2.0 * dB[l][j]) * dA[i][k];

    }
  }

  // ----> (AA|BB) : V0_Coul, V0_CT_M, E1, E2, E3, E4 <---- //
  for (s_aabb->first(); s_aabb->is_done() == false; s_aabb->next()) 
  {
    s_aabb->compute_shell(t_aabb);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_aabb->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // A
         int j = ii->j(); // A
         int k = ii->k(); // B
         int l = ii->l(); // B
         double eri = b_aabb[ii->index()];

         // V0_Coul
         V0_Coul += pAt[j][i] * pBt[l][k] * eri;

         // V0_CT_M
         V0_CT_M += Q1 * eri * coA[i][homo_A] * coA[j][homo_A] * coB[k][homo_B] * coB[l][homo_B];
         V0_CT_M += Q1 * eri * cvA[i][lumo_A] * cvA[j][lumo_A] * cvB[k][lumo_B] * cvB[l][lumo_B];
         V0_CT_M += Q2 * eri * coA[i][homo_A] * coA[j][homo_A] * cvB[k][lumo_B] * cvB[l][lumo_B];
         V0_CT_M += Q2 * eri * cvA[i][lumo_A] * cvA[j][lumo_A] * coB[k][homo_B] * coB[l][homo_B];

         // E1
         E1 += 2.0 * eri * (pA[j][i] - 2.0 * dA[j][i]) * dB[l][k];
         // E2
         E2 += 2.0 * eri * (pB[l][k] - 2.0 * dB[l][k]) * dA[j][i];
         // E3
         E3 -= coA[j][homo_A] * coA[i][homo_A] * cvB[l][lumo_B] * cvB[k][lumo_B] * eri;
         // E4
         E4 -= cvA[j][lumo_A] * cvA[i][lumo_A] * coB[l][homo_B] * coB[k][homo_B] * eri;

    }
  }

  // ----> (AA|AA) : V0_CT_M <---- // 
  for (s_aaaa->first(); s_aaaa->is_done() == false; s_aaaa->next()) 
  {
    s_aaaa->compute_shell(t_aaaa);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_aaaa->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // A
         int j = ii->j(); // A
         int k = ii->k(); // A
         int l = ii->l(); // A
         double eri = b_aaaa[ii->index()];

         // V0_CT_M
         V0_CT_M += Q3 * eri * coA[i][homo_A] * coA[j][homo_A] * cvA[k][lumo_A] * cvA[l][lumo_A];
    }
  }

  // ----> (BB|BB) : V0_CT_M <---- // 
  for (s_bbbb->first(); s_bbbb->is_done() == false; s_bbbb->next()) 
  {
    s_bbbb->compute_shell(t_bbbb);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_bbbb->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) 
    {
         int i = ii->i(); // B
         int j = ii->j(); // B
         int k = ii->k(); // B
         int l = ii->l(); // B
         double eri = b_bbbb[ii->index()];

         // V0_CT_M
         V0_CT_M += Q3 * eri * coB[i][homo_B] * coB[j][homo_B] * cvB[k][lumo_B] * cvB[l][lumo_B];
    }
  }



  // Finish E1 and E2
  for (int i=0; i<nbf_A; ++i) {
  for (int j=0; j<nbf_A; ++j) {
       E1 += (pA[j][i] - 2.0 * dA[j][i])*VaoB_AA->get(i,j);
  }}
  for (int i=0; i<nbf_B; ++i) {
  for (int j=0; j<nbf_B; ++j) {
       E2 += (pB[j][i] - 2.0 * dB[j][i])*VaoA_BB->get(i,j);
  }}


  if (options_.get_bool("TI_CIS_PRINT_FOCK_MATRIX")) {
      Fao_AB->set_name("Fock matrix (A,B): unperturbed - density matrix as Hadamard sum of monomer");
      Fao_AB->print();
  }
 
  // Compute the full Fock matrix from SCF wavefunction of entire dimer
  if (options_.get_bool("TI_CIS_SCF_FOCK_MATRIX")) {
     const int nbf = wfn_union_->basisset()->nbf();
     psi::SharedWavefunction wfn_dimer = oepdev::solve_scf(wfn_union_->molecule(), wfn_union_->basisset(), 
                                wfn_union_->get_basisset("BASIS_DF_SCF"),
                                oepdev::create_superfunctional("HF", options_), options_, psi::PSIO::shared_object(), false);
     for (int i=0; i<nbf_A; ++i) {
          for (int j=0; j<nbf_B; ++j) {
               f[i][j] = wfn_dimer->Fa()->get(i, nbf_A+j);
          }
     }
     if (options_.get_int("PRINT") > -1) {                                                               
         Fao_AB->set_name("Fock matrix (A,B): perturbed - SCF");
         Fao_AB->print();
     }
  }

  // ----> Add Fock matrix contributions to all V0 <---- //
  for (int i=0; i<nbf_A; ++i) {
  for (int j=0; j<nbf_B; ++j) {
       double v = f[i][j];
       V0_ET1 += cvA[i][lumo_A] * cvB[j][lumo_B] * v;
       V0_ET2 += cvA[i][lumo_A] * cvB[j][lumo_B] * v;
       V0_HT1 -= coA[i][homo_A] * coB[j][homo_B] * v;
       V0_HT2 -= coA[i][homo_A] * coB[j][homo_B] * v;
  }}
  // ----> Multiply V0 by CIS amplitudes
  V0_ET1 *= t_A;
  V0_ET2 *= t_B;
  V0_HT1 *= t_A;
  V0_HT2 *= t_B;

  if (wfn_union_->options().get_int("PRINT") > -1) {
    psi::outfile->Printf(" ===> Basis Energies [EV] <===\n\n");
    psi::outfile->Printf("   E01 = %9.3f  E1 = %9.3f\n", E01*OEPDEV_AU_EV, E1*OEPDEV_AU_EV);
    psi::outfile->Printf("   E02 = %9.3f  E2 = %9.3f\n", E02*OEPDEV_AU_EV, E2*OEPDEV_AU_EV);
    psi::outfile->Printf("   E03 = %9.3f  E3 = %9.3f\n", E03*OEPDEV_AU_EV, E3*OEPDEV_AU_EV);
    psi::outfile->Printf("   E04 = %9.3f  E4 = %9.3f\n", E04*OEPDEV_AU_EV, E4*OEPDEV_AU_EV);
    psi::outfile->Printf("\n");
  }


  // Compute overlap-corrected indirect coupling matrix elements
  double V_ET1 = (V0_ET1 - 0.5*(S13*(E1+E2)))/(1.0 - S13*S13);
  double V_ET2 = (V0_ET2 - 0.5*(S42*(E1+E2)))/(1.0 - S42*S42);
  double V_HT1 = (V0_HT1 - 0.5*(S14*(E1+E2)))/(1.0 - S14*S14);
  double V_HT2 = (V0_HT2 - 0.5*(S32*(E1+E2)))/(1.0 - S32*S32);
  double V_CT  = (V0_CT  - 0.5*(S34*(E1+E2)))/(1.0 - S34*S34);
  double V_CT_M= (V0_CT_M- 0.5*(S34*(E1+E2)))/(1.0 - S34*S34);

  // Compute final coupling contributions
  double V_Coul = V0_Coul / (1.0 - S12*S12);
  double V_Exch = V0_Exch / (1.0 - S12*S12);
  double V_Ovrl =-(E1 + E2)*S12/(2.0 * (1.0 - S12*S12));

  //E3 = E03; E4 = E04;
  double V_TI_2 =-(V_ET1*V_HT2)/(E3-E1) -(V_HT1*V_ET2)/(E4-E1);
  double V_TI_3 = (V_ET1*V_ET2 + V_HT1*V_HT2) * V_CT / ((E3-E1)*(E4-E1));

  double V0_TI_2 =-(V0_ET1*V0_HT2)/(E3-E01) -(V0_HT1*V0_ET2)/(E4-E01);
  double V0_TI_3 = (V0_ET1*V0_ET2 + V0_HT1*V0_HT2) * V0_CT / ((E3-E01)*(E4-E01));

  double V_TI_3_M = (V_ET1*V_ET2 + V_HT1*V_HT2) * V_CT_M / ((E3-E1)*(E4-E1));
  double V0_TI_3_M = (V0_ET1*V0_ET2 + V0_HT1*V0_HT2) * V0_CT_M / ((E3-E01)*(E4-E01));

  double V_direct = V_Coul + V_Exch + V_Ovrl;
  double V_indirect = V_TI_2 + V_TI_3;

  double V0_direct = V0_Coul + V0_Exch;
  double V0_indirect = V0_TI_2 + V0_TI_3;

  double V_TI_CIS = V_direct + V_indirect;
  double V0_TI_CIS = V0_direct + V0_indirect;

  psi::timer_off("Solver EET TI/CIS               ");


  // ===> Compute TrCAMM coupling <=== //
  psi::timer_on("Solver EET TrCAMM               ");
  SharedMTPConv V0_TrCAMM = trcamm_A->energy(trcamm_B);
  double V0_TrCAMM_R1 = V0_TrCAMM->level(oepdev::MultipoleConvergence::R1)->get(0,0);
  double V0_TrCAMM_R2 = V0_TrCAMM->level(oepdev::MultipoleConvergence::R2)->get(0,0);
  double V0_TrCAMM_R3 = V0_TrCAMM->level(oepdev::MultipoleConvergence::R3)->get(0,0);
  double V0_TrCAMM_R4 = V0_TrCAMM->level(oepdev::MultipoleConvergence::R4)->get(0,0);
  double V0_TrCAMM_R5 = V0_TrCAMM->level(oepdev::MultipoleConvergence::R5)->get(0,0);
  psi::timer_off("Solver EET TrCAMM               ");


  // ---> Save <--- //
  psi::Process::environment.globals["EET V0 COUL CM-1"      ] = V0_Coul      *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TrCAMM R1 CM-1" ] = V0_TrCAMM_R1 *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TrCAMM R2 CM-1" ] = V0_TrCAMM_R2 *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TrCAMM R3 CM-1" ] = V0_TrCAMM_R3 *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TrCAMM R4 CM-1" ] = V0_TrCAMM_R4 *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TrCAMM R5 CM-1" ] = V0_TrCAMM_R5 *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 EXCH CM-1"      ] = V0_Exch      *OEPDEV_AU_CMRec;  
  psi::Process::environment.globals["EET V COUL CM-1"       ] = V_Coul       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V EXCH CM-1"       ] = V_Exch       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V OVRL CM-1"       ] = V_Ovrl       *OEPDEV_AU_CMRec;
  //                                                                                           
  psi::Process::environment.globals["EET V0 ET1 CM-1"       ] = V0_ET1       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 ET2 CM-1"       ] = V0_ET2       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 HT1 CM-1"       ] = V0_HT1       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 HT2 CM-1"       ] = V0_HT2       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 CT CM-1"        ] = V0_CT        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 CT(MULLIKEN) CM-1"] = V0_CT_M     *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V ET1 CM-1"        ] = V_ET1        *OEPDEV_AU_CMRec;  
  psi::Process::environment.globals["EET V ET2 CM-1"        ] = V_ET2        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V HT1 CM-1"        ] = V_HT1        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V HT2 CM-1"        ] = V_HT2        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V CT CM-1"         ] = V_CT         *OEPDEV_AU_CMRec;
  //psi::Process::environment.globals["EET V CT(MULLIKEN) CM-1"      ] = V_CT_M      *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TI(2) CM-1"     ] = V0_TI_2      *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TI(3) CM-1"     ] = V0_TI_3      *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TI(2) CM-1"      ] = V_TI_2       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TI(3) CM-1"      ] = V_TI_3       *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V0 Direct CM-1"    ] = V0_direct    *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 Indirect CM-1"  ] = V0_indirect  *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V Direct CM-1"     ] = V_direct     *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V Indirect CM-1"   ] = V_indirect   *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V0 TI_CIS  CM-1"   ] = V0_TI_CIS   *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TI_CIS  CM-1"    ] = V_TI_CIS    *OEPDEV_AU_CMRec;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > -1) {
     psi::outfile->Printf("  ==> SOLVER: EET coupling constant <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (TI-CIS)    <==\n\n");
     psi::outfile->Printf("     V0 Coul   = %13.2f\n", V0_Coul *OEPDEV_AU_CMRec         );
     psi::outfile->Printf("     V0 Exch   = %13.2f\n", V0_Exch *OEPDEV_AU_CMRec         );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     V Coul    = %13.2f\n", V_Coul *OEPDEV_AU_CMRec          );
     psi::outfile->Printf("     V Exch    = %13.2f\n", V_Exch *OEPDEV_AU_CMRec          );
     psi::outfile->Printf("     V Ovrl    = %13.2f\n", V_Ovrl *OEPDEV_AU_CMRec          );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     TrCAMM-R1 = %13.2f\n", V0_TrCAMM_R1 *OEPDEV_AU_CMRec    );
     psi::outfile->Printf("     TrCAMM-R2 = %13.2f\n", V0_TrCAMM_R2 *OEPDEV_AU_CMRec    );
     psi::outfile->Printf("     TrCAMM-R3 = %13.2f\n", V0_TrCAMM_R3 *OEPDEV_AU_CMRec    );
     psi::outfile->Printf("     TrCAMM-R4 = %13.2f\n", V0_TrCAMM_R4 *OEPDEV_AU_CMRec    );
     psi::outfile->Printf("     TrCAMM-R5 = %13.2f\n", V0_TrCAMM_R5 *OEPDEV_AU_CMRec    );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     V0 Direct = %13.2f\n", V0_direct*OEPDEV_AU_CMRec        );
     psi::outfile->Printf("     V  Direct = %13.2f\n", V_direct *OEPDEV_AU_CMRec        );
     psi::outfile->Printf("     ===============================\n"                      );
     psi::outfile->Printf("     V0_ET1 = %13.2f V_ET1 = %13.2f\n", V0_ET1*OEPDEV_AU_CMRec, V_ET1*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0_ET2 = %13.2f V_ET2 = %13.2f\n", V0_ET2*OEPDEV_AU_CMRec, V_ET2*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0_HT1 = %13.2f V_HT1 = %13.2f\n", V0_HT1*OEPDEV_AU_CMRec, V_HT1*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0_HT2 = %13.2f V_HT2 = %13.2f\n", V0_HT2*OEPDEV_AU_CMRec, V_HT2*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0_CT  = %13.2f V_CT  = %13.2f\n", V0_CT *OEPDEV_AU_CMRec, V_CT *OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0_CT_M= %13.2f V_CT_M= %13.2f\n", V0_CT_M*OEPDEV_AU_CMRec, V_CT_M*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     V0 TI_2= %13.2f V TI_2= %13.2f\n", V0_TI_2*OEPDEV_AU_CMRec, V_TI_2*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0 TI_3= %13.2f V TI_3= %13.2f\n", V0_TI_3*OEPDEV_AU_CMRec, V_TI_3*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     V0 TI_3_M= %13.2f V TI_3_M= %13.2f\n", V0_TI_3_M*OEPDEV_AU_CMRec, V_TI_3_M*OEPDEV_AU_CMRec);
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     V0 Indir = %13.2f\n", V0_indirect*OEPDEV_AU_CMRec       );
     psi::outfile->Printf("     V  Indir = %13.2f\n", V_indirect *OEPDEV_AU_CMRec       );
     psi::outfile->Printf("     ===============================\n"                      );
     psi::outfile->Printf("     V0 TI/CIS= %13.2f\n", V0_TI_CIS  *OEPDEV_AU_CMRec       );
     psi::outfile->Printf("     V  TI/CIS= %13.2f\n", V_TI_CIS   *OEPDEV_AU_CMRec       );
     psi::outfile->Printf("\n");
  }


  return e;
}

double EETCouplingSolver::compute_oep_based_fujimoto_ti_cis() { //TODO
  double e = 0.0;
  psi::timer_on("Solver EET TI/CIS OEP-Based     ");
  // Return
  psi::timer_off("Solver EET TI/CIS OEP-Based     ");
  return e; 
}
