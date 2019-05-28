//#include "psi4/libtrans/integraltransform.h"
//#include "psi4/libdpd/dpd.h"

#include "solver.h"

using namespace std;
using namespace psi;
using namespace oepdev;

// CT Solver//
ChargeTransferEnergySolver::ChargeTransferEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Benchmarks
  methods_benchmark_.push_back("OTTO_LADIK"      );
  methods_benchmark_.push_back("EFP2"            );
  // OEP-based
  methods_oepBased_ .push_back("MURRELL_ETAL"    );
}
ChargeTransferEnergySolver::~ChargeTransferEnergySolver() {}
double ChargeTransferEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "MURRELL_ETAL"      )  e = compute_oep_based_murrell_etal();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for repulsion energy calculations!\n");
  }
  return e;
}
double ChargeTransferEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" ||
           method == "OTTO_LADIK"   ) e = compute_benchmark_otto_ladik();
  else if (method == "EFP2"         ) e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
double ChargeTransferEnergySolver::compute_benchmark_otto_ladik(){}
double ChargeTransferEnergySolver::compute_benchmark_efp2(){}
double ChargeTransferEnergySolver::compute_oep_based_murrell_etal()
{
  // ===> Compute OEP Objects <=== //
  SharedOEPotential oep_1 = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY",
                                                       wfn_union_->l_wfn(0), 
                                                       wfn_union_->l_auxiliary(0), 
                                                       wfn_union_->l_intermediate(0), 
                                                       wfn_union_->options());
  SharedOEPotential oep_2 = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY", 
                                                       wfn_union_->l_wfn(1), 
                                                       wfn_union_->l_auxiliary(1), 
                                                       wfn_union_->l_intermediate(1), 
                                                       wfn_union_->options());
  oep_1->compute();
  oep_2->compute();

  // ===> Compute Overlap Matrices <=== //
  int nbf_p1 = wfn_union_->l_nbf(0);
  int nbf_p2 = wfn_union_->l_nbf(1);
  int nbf_a1 = wfn_union_->l_auxiliary(0)->nbf();
  int nbf_a2 = wfn_union_->l_auxiliary(1)->nbf();
  int nocc_1 = wfn_union_->l_ndocc(0);
  int nocc_2 = wfn_union_->l_ndocc(1);
  int nvir_1 = wfn_union_->l_nvir(0);
  int nvir_2 = wfn_union_->l_nvir(1);

  std::shared_ptr<psi::Matrix> Sao_1p2p     = std::make_shared<psi::Matrix>("Sao 1p2p", nbf_p1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1a2p     = std::make_shared<psi::Matrix>("Sao 1a2p", nbf_a1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1p2a     = std::make_shared<psi::Matrix>("Sao 1p2a", nbf_p1, nbf_a2);

  psi::IntegralFactory fact_1p2p(wfn_union_->l_primary  (0), wfn_union_->l_primary  (1), wfn_union_->l_primary  (0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1a2p(wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1), wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1p2a(wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1), wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1));

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2p(fact_1p2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

  ovlInt_1p2p->compute(Sao_1p2p);
  ovlInt_1a2p->compute(Sao_1a2p);
  ovlInt_1p2a->compute(Sao_1p2a);

  // ---> Canonical MO's: LCAO and energies <--- //
  std::shared_ptr<psi::Matrix> Ca_occ_1 = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_2 = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_vir_1 = wfn_union_->l_wfn(0)->Ca_subset("AO","VIR");
  std::shared_ptr<psi::Matrix> Ca_vir_2 = wfn_union_->l_wfn(1)->Ca_subset("AO","VIR");
  std::shared_ptr<psi::Vector> e_occ_1  = wfn_union_->l_wfn(0)->epsilon_a_subset("AO","OCC");
  std::shared_ptr<psi::Vector> e_occ_2  = wfn_union_->l_wfn(1)->epsilon_a_subset("AO","OCC");
  std::shared_ptr<psi::Vector> e_vir_1  = wfn_union_->l_wfn(0)->epsilon_a_subset("AO","VIR");
  std::shared_ptr<psi::Vector> e_vir_2  = wfn_union_->l_wfn(1)->epsilon_a_subset("AO","VIR");

  std::shared_ptr<psi::Matrix> S1 = psi::Matrix::doublet(Ca_occ_1, Sao_1p2a, true, false); // OCC(A) x AUX(B)
  std::shared_ptr<psi::Matrix> S2 = psi::Matrix::doublet(Ca_occ_2, Sao_1a2p, true, true ); // OCC(B) x AUX(A)

  // ---> Localize occupied orbitals <--- //
  // TODO

  // ===> Compute V1 term <=== //
  std::shared_ptr<psi::Matrix> v_ab_v1 = psi::Matrix::doublet(S1, oep_2->matrix("Murrell-etal.V1"), false, false);
  std::shared_ptr<psi::Matrix> v_ba_v1 = psi::Matrix::doublet(S2, oep_1->matrix("Murrell-etal.V1"), false, false);

  // ===> Compute V2 term <=== //
  // TODO

  // ===> Compute V3 term <=== //
  // TODO

  // ===> Compute CT Energy <=== //
  double e_ab_v1 = compute_ct_component(nocc_1, nvir_2, e_occ_1, e_vir_2, v_ab_v1);
  double e_ba_v1 = compute_ct_component(nocc_2, nvir_1, e_occ_2, e_vir_1, v_ba_v1);

  double e_ab = e_ab_v1;
  double e_ba = e_ba_v1;

  double e_tot = e_ab + e_ba;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Charge-Transfer Energy Calculations    <==\n"  );
     psi::outfile->Printf("  ==>     OEP-Based (Murrell-etal               )    <==\n\n");
     psi::outfile->Printf("     E_A-->B   = %13.6f\n", e_ab                             );
     psi::outfile->Printf("     E_B-->A   = %13.6f\n", e_ba                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E_TOT     = %13.6f\n", e_tot                            );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("\n");

}

  // Return the Total CT Energy
  return e_tot; 
}

// private
double ChargeTransferEnergySolver::compute_ct_component(int nocc_X, int nvir_Y, std::shared_ptr<psi::Vector> eps_occ_X, std::shared_ptr<psi::Vector> eps_vir_Y, std::shared_ptr<psi::Matrix> V)
{
   // requires matrix elements in CMO (canonical SCF) basis
   double e_XY = 0.0;
   for (int i=0; i<nocc_X; ++i) { /* X-->Y term */
       for (int n=0; n<nvir_Y; ++n) {
            double vin = V->get(i, n); /* V: OCC_X x VIR_Y */
            e_XY += (vin * vin) / (eps_occ_X->get(i) - eps_vir_Y->get(n));
       }
  }
  e_XY *= 2.0;
  return e_XY;
}
