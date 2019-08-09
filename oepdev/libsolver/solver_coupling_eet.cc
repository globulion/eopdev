#include "solver.h"
#include "psi4/libpsi4util/process.h"
#include <utility>
#include "../libutil/integrals_iter.h"


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

  psi::IntegralFactory fact_ABAB(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_BABA(wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_ABBB(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_AAAB(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_AABB(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(1));


  // Compute one-electron integrals
  std::shared_ptr<psi::OneBodyAOInt> oneInt, ovlInt(fact_ABAB.ao_overlap()), kinInt(fact_ABAB.ao_kinetic());
  std::shared_ptr<psi::PotentialInt> potInt_B_AB= std::make_shared<psi::PotentialInt>(fact_ABAB.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_A_BA= std::make_shared<psi::PotentialInt>(fact_BABA.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));

  SharedMatrix Zxyz_A = std::make_shared<psi::Matrix>(potInt_B_AB->charge_field());
  SharedMatrix Zxyz_B = std::make_shared<psi::Matrix>(potInt_A_BA->charge_field());

  potInt_B_AB->set_charge_field(Zxyz_B);
  potInt_A_BA->set_charge_field(Zxyz_A);
  oneInt = potInt_B_AB;
  oneInt->compute(VaoB_AB);
  oneInt = potInt_A_BA;
  oneInt->compute(VaoA_BA);
  ovlInt->compute(Sao_AB);
  kinInt->compute(Tao_AB);

  // Fock matrix of entire system in AB space: 1-electron contribution
  VaoA_BA->transpose();
  SharedMatrix Fao_AB = Tao_AB->clone();
  Fao_AB->add(VaoB_AB);
  Fao_AB->add(VaoA_BA);


  // Auxiliary data
  int I = options_.get_int("EXCITED_STATE_A");
  int J = options_.get_int("EXCITED_STATE_B");

  const int Ne = 2.0 * (wfn_union_->l_wfn(0)->nalpha() + wfn_union_->l_wfn(1)->nalpha());
  const int homo_A = wfn_union_->l_wfn(0)->nalpha() - 1;
  const int homo_B = wfn_union_->l_wfn(1)->nalpha() - 1;
  const int lumo_A = wfn_union_->l_wfn(0)->nalpha();
  const int lumo_B = wfn_union_->l_wfn(1)->nalpha();

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
  SharedMatrix Peg_A, Peg_B;
  // TrCAMM's
  SharedDMTPole trcamm_A, trcamm_B;


  // Obtain CIS HOMO-LUMO amplitudes and unperturbed excitation energies
  double t_A, t_B, E_ex_A, E_ex_B;
  {
      std::shared_ptr<CISComputer> cis_A = CISComputer::build("RESTRICTED", wfn_union_->l_wfn(0), options_); 
      cis_A->compute();
      this->determine_electronic_state(cis_A, I);
      E_ex_A = cis_A->eigenvalues()->get(I);
      t_A = cis_A->U_homo_lumo(I).first;
      trcamm_A = cis_A->trcamm(I);
      Peg_A = cis_A->Ta_ao(I);
      Peg_A->add(cis_A->Tb_ao(I));
  }{
      std::shared_ptr<CISComputer> cis_B = CISComputer::build("RESTRICTED", wfn_union_->l_wfn(1), options_); 
      cis_B->compute();
      this->determine_electronic_state(cis_B, J);
      E_ex_B = cis_B->eigenvalues()->get(J);
      t_B = cis_B->U_homo_lumo(J).first;
      trcamm_B = cis_B->trcamm(J);
      Peg_B = cis_B->Ta_ao(J);
      Peg_B->add(cis_B->Tb_ao(J));
  }


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


  // Compute Hamiltonian eigenvalues TODO
  double E01= E_ex_A;
  double E02= E_ex_B;
  double E03=-eps_a_occ_A->get(homo_A) + eps_a_vir_B->get(0     );
  double E04= eps_a_vir_A->get(0     ) - eps_a_occ_B->get(homo_B);

  double E1 = E01;
  double E2 = E02;
  double E3 = E03;
  double E4 = E04;

  // Compute direct coupling constants TODO
  double V0_Coul = 0.0; 
  double V0_Exch = 0.0;

  // Compute indirect coupling constants TODO
  double V0_ET1 = 0.0;
  double V0_ET2 = 0.0;
  double V0_HT1 = 0.0;
  double V0_HT2 = 0.0;
  double V0_CT  = 0.0;

  // ===> Accumulate (ab|cd) contributions from AO ERI's <=== //
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aaab = oepdev::ShellCombinationsIterator::build(fact_AAAB, "ALL");
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_abbb = oepdev::ShellCombinationsIterator::build(fact_ABBB, "ALL");
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_abab = oepdev::ShellCombinationsIterator::build(fact_ABAB, "ALL");
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aabb = oepdev::ShellCombinationsIterator::build(fact_AABB, "ALL");

  std::shared_ptr<psi::TwoBodyAOInt> t_aaab(fact_AAAB.eri()); const double * b_aaab = t_aaab->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_abbb(fact_ABBB.eri()); const double * b_abbb = t_abbb->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_abab(fact_ABAB.eri()); const double * b_abab = t_abab->buffer();
  std::shared_ptr<psi::TwoBodyAOInt> t_aabb(fact_AABB.eri()); const double * b_aabb = t_aabb->buffer();



  // ----> (AA|AB) : ET1, HT1, Fock <---- // 
  double** dA = Da_A->pointer();
  double** dB = Da_B->pointer();
  double** f  = Fao_AB->pointer();
  double** coA= Ca_occ_A->pointer();
  double** coB= Ca_occ_B->pointer();
  double** cvA= Ca_vir_A->pointer();
  double** cvB= Ca_vir_B->pointer();

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
         f[i][l] -=       dA[k][j];

         // V0_ET1
         V0_ET1 += 2.0 * eri * coA[i][homo_A] * cvA[j][0] * cvA[k][homo_A] * cvB[l][0];
         V0_ET1 -=       eri * coA[k][homo_A] * cvA[j][0] * cvA[i][homo_A] * cvB[l][0];

         // V0_HT1
         V0_HT1 += 2.0 * eri * cvA[j][0] * coA[i][homo_A] * cvA[k][0] * coB[l][homo_B];
         V0_HT1 -=       eri * cvA[j][0] * coA[k][homo_A] * cvA[i][0] * coB[l][homo_B];
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

         // Fock matrix
         f[i][j] += 2.0 * dA[k][l] * eri;
         f[i][l] -=       dA[k][j];

         // V0_ET2
         V0_ET1 += 2.0 * eri * cvA[i][0] * coB[j][homo_B] * cvB[l][0] * coB[k][homo_B];
         V0_ET1 -=       eri * cvA[i][0] * coB[l][homo_B] * cvB[j][0] * coB[k][homo_B];

         // V0_HT2
         V0_ET1 += 2.0 * eri * coA[i][homo_A] * cvB[j][0] * cvB[l][homo_B] * cvB[k][0];
         V0_ET1 -=       eri * coA[i][homo_A] * cvB[l][0] * cvB[j][homo_B] * cvB[k][0];

    }
  }

  // ----> (AB|AB) : CT, E1, E2 <---- //
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

         // V0_CT
         V0_CT  += 2.0 * eri * cvB[j][0] * coA[i][homo_A] * coB[l][homo_B] * cvA[k][0];
         V0_CT  -=       eri * cvB[l][0] * coA[i][homo_A] * coB[j][homo_B] * cvA[k][0];

    }
  }

  // ----> (AA|BB) : E1, E2, E3, E4 <---- //

  // ----> Add Fock matrix contributions to all V0 <---- //

  // Compute overlap-corrected indirect coupling matrix elements
  double V_ET1 = (V0_ET1 - 0.5*(S13*(E1+E2)))/(1.0 - S13*S13);
  double V_ET2 = (V0_ET2 - 0.5*(S42*(E1+E2)))/(1.0 - S42*S42);
  double V_HT1 = (V0_HT1 - 0.5*(S14*(E1+E2)))/(1.0 - S14*S14);
  double V_HT2 = (V0_HT2 - 0.5*(S32*(E1+E2)))/(1.0 - S32*S32);
  double V_CT  = (V0_CT  - 0.5*(S34*(E1+E2)))/(1.0 - S34*S34);

  // Compute final coupling contributions
  double V_Coul = V0_Coul / (1.0 - S12*S12);
  double V_Exch = V0_Exch / (1.0 - S12*S12);
  double V_Ovrl =-(E1 + E2)*S12/(2.0 * (1.0 - S12*S12));

  double V_TI_2 = -(V_ET1*V_HT2)/(E3-E1) - (V_HT1*V_ET2)/(E4-E1);
  double V_TI_3 = (V_ET1*V_CT*V_ET2 + V_HT1*V_CT*V_HT2) / ((E3-E1)*(E4-E1));

  double V0_TI_2 = -(V0_ET1*V0_HT2)/(E3-E1) - (V0_HT1*V0_ET2)/(E4-E1);
  double V0_TI_3 = (V0_ET1*V0_CT*V0_ET2 + V0_HT1*V0_CT*V0_HT2) / ((E3-E1)*(E4-E1));
 
  double V_direct = V_Coul + V_Exch + V_Ovrl;
  double V_indirect = V_TI_2 + V_TI_3;

  double V0_direct = V0_Coul + V0_Exch;
  double V0_indirect = V0_TI_2 + V0_TI_3;

  double V_TDFI_TI = V_direct + V_indirect;
  double V0_TDFI_TI = V0_direct + V0_indirect;

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
  //psi::Process::environment.globals["EET V0 CT(MULLIKEN) CM-1"     ] = V0_CT_M     *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V ET1 CM-1"        ] = V_ET1        *OEPDEV_AU_CMRec;  
  psi::Process::environment.globals["EET V ET2 CM-1"        ] = V_ET2        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V HT1 CM-1"        ] = V_HT1        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V HT2 CM-1"        ] = V_HT2        *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V CT CM-1"         ] = V_CT         *OEPDEV_AU_CMRec;
  //psi::Process::environment.globals["EET V CT(MULLIKEN) CM-1"      ] = V_CT_M      *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TI(2) CM-1"     ] = V0_TI_2       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 TI(3) CM-1"     ] = V0_TI_3       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TI(2) CM-1"      ] = V_TI_2       *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TI(3) CM-1"      ] = V_TI_3       *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V0 Direct CM-1"    ] = V0_direct     *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V0 Indirect CM-1"  ] = V0_indirect   *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V Direct CM-1"     ] = V_direct     *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V Indirect CM-1"   ] = V_indirect   *OEPDEV_AU_CMRec;
  //
  psi::Process::environment.globals["EET V0 TDFI_TI  CM-1"  ] = V0_TDFI_TI    *OEPDEV_AU_CMRec;
  psi::Process::environment.globals["EET V TDFI_TI  CM-1"   ] = V_TDFI_TI    *OEPDEV_AU_CMRec;

  return e;
}

double EETCouplingSolver::compute_oep_based_fujimoto_ti_cis() { //TODO
  double e = 0.0;
  psi::timer_on("Solver EET TI/CIS OEP-Based     ");
  // Return
  psi::timer_off("Solver EET TI/CIS OEP-Based     ");
  return e; 
}
