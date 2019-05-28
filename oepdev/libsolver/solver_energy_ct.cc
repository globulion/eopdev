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
           method == "MURRELL_ETAL" ) e = compute_benchmark_otto_ladik(); // TODO: change to murrell_etal or revise naming!
  else if (method == "EFP2"         ) e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
// Copied from MC_OTTO_LADIK for testing
double ChargeTransferEnergySolver::compute_benchmark_otto_ladik(){

  //psi::timer_on("SOLVER: Charge-transfer Energy Calculations (Otto-Ladik)");
  psi::timer_on("Solver E(CT) Otto-Ladik       ");
  

  int nbf_1 = wfn_union_->l_nbf(0);
  int nbf_2 = wfn_union_->l_nbf(1);

  // ===> One electron part <=== //
  // Term 1//
  std::shared_ptr<psi::Matrix> VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> VaoB11    = std::make_shared<psi::Matrix>("VaoB(1,1)" , nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> VaoA22    = std::make_shared<psi::Matrix>("VaoA(2,2)" , nbf_2, nbf_2);
  std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)"  , nbf_1, nbf_2);

  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_11(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_22(wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));

  std::shared_ptr<psi::OneBodyAOInt> oneInt, ovlInt(fact_12.ao_overlap());
  std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_12.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_21.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_11= std::make_shared<psi::PotentialInt>(fact_11.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_22= std::make_shared<psi::PotentialInt>(fact_22.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(1));


  std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
  std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

  potInt_1->set_charge_field(Zxyz_2);
  potInt_11->set_charge_field(Zxyz_2);
  potInt_22->set_charge_field(Zxyz_1);
  oneInt = potInt_1;
  oneInt->compute(VaoB12);
  oneInt = potInt_11;
  oneInt->compute(VaoB11);
  oneInt = potInt_22;
  oneInt->compute(VaoA22);
  ovlInt->compute(Sao12);

  std::shared_ptr<psi::Matrix> Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_vir_B = wfn_union_->l_wfn(1)->Ca_subset("AO","VIR");

  std::shared_ptr<psi::Matrix> Smo12  = psi::Matrix::triplet(Ca_occ_A, Sao12 , Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> Smo1Y  = psi::Matrix::triplet(Ca_occ_A, Sao12 , Ca_vir_B, true, false, false);  
  std::shared_ptr<psi::Matrix> VmoB1Y = psi::Matrix::triplet(Ca_occ_A, VaoB12, Ca_vir_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoAY2 = psi::Matrix::triplet(Ca_vir_B, VaoA22, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB11 = psi::Matrix::triplet(Ca_occ_A, VaoB11, Ca_occ_A, true, false, false);

  //Potential integral II//
  double** S1 = Smo1Y->pointer();
  std::shared_ptr<psi::Matrix> PI2 = psi::Matrix::doublet(VmoB11, Smo1Y, false, false);

  //Potential integral III//
  double** S2 = Smo12->pointer();
  std::shared_ptr<psi::Matrix> PI3 = psi::Matrix::doublet(Smo12, VmoAY2, false, true);


  //Sum of one-electron integrals written in VmoB1Y
  std::shared_ptr<psi::Matrix> VmoB1Y_group_1 = VmoB1Y->clone();                     // debug for OEP
  std::shared_ptr<psi::Matrix> VmoB1Y_group_2 = std::make_shared<psi::Matrix>(PI2);  // debug for OEP 
  std::shared_ptr<psi::Matrix> VmoB1Y_group_3 = std::make_shared<psi::Matrix>(PI3);  // debug for OEP
  VmoB1Y->subtract(PI2);
  VmoB1Y->subtract(PI3);



  //Term 2//
  std::shared_ptr<psi::Matrix> VaoA21      = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_1);
  std::shared_ptr<psi::Matrix> VaoA22_2    = std::make_shared<psi::Matrix>("VaoA_2(2,2)" , nbf_2, nbf_2);
  std::shared_ptr<psi::Matrix> VaoB11_2    = std::make_shared<psi::Matrix>("VaoB_2(1,1)" , nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> Sao21       = std::make_shared<psi::Matrix>("Sao(2,1)"  , nbf_2, nbf_1);

  std::shared_ptr<psi::OneBodyAOInt> oneInt_2, ovlInt_2(fact_21.ao_overlap());
  std::shared_ptr<psi::PotentialInt> potInt_11_2 = std::make_shared<psi::PotentialInt>(fact_11.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_22_2 = std::make_shared<psi::PotentialInt>(fact_22.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(1));

  potInt_2->set_charge_field(Zxyz_1);
  potInt_22_2->set_charge_field(Zxyz_1);
  potInt_11_2->set_charge_field(Zxyz_2);
  oneInt_2 = potInt_2;
  oneInt_2->compute(VaoA21);
  oneInt_2 = potInt_22_2;
  oneInt_2->compute(VaoA22_2);
  oneInt_2 = potInt_11_2;
  oneInt_2->compute(VaoB11_2);
  ovlInt_2->compute(Sao21);
  

  std::shared_ptr<psi::Matrix> Ca_vir_A = wfn_union_->l_wfn(0)->Ca_subset("AO","VIR");

  std::shared_ptr<psi::Matrix> Smo21  = psi::Matrix::triplet(Ca_occ_B, Sao21 , Ca_occ_A, true, false, false);
  std::shared_ptr<psi::Matrix> Smo2X  = psi::Matrix::triplet(Ca_occ_B, Sao21 , Ca_vir_A, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA2X = psi::Matrix::triplet(Ca_occ_B, VaoA21, Ca_vir_A, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA22 = psi::Matrix::triplet(Ca_occ_B, VaoA22_2, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoBX1 = psi::Matrix::triplet(Ca_vir_A, VaoB11_2, Ca_occ_A, true, false, false);

  //Potential integral II//
  double** S1_2 = Smo2X->pointer();
  std::shared_ptr<psi::Matrix> PI4 = psi::Matrix::doublet(VmoA22, Smo2X, false, false);

  //Potential integral III//
  double** S2_2 = Smo21->pointer();
  std::shared_ptr<psi::Matrix> PI5 = psi::Matrix::doublet(Smo21, VmoBX1, false, true);

  //Sum of one-electron integrals written to VmoA2X
  std::shared_ptr<psi::Matrix> VmoA2X_group_1 = VmoA2X->clone();                      // debug for OEP
  std::shared_ptr<psi::Matrix> VmoA2X_group_2 = std::make_shared<psi::Matrix>(PI4);   // debug for OEP 
  std::shared_ptr<psi::Matrix> VmoA2X_group_3 = std::make_shared<psi::Matrix>(PI5);   // debug for OEP
  VmoA2X->subtract(PI4);
  VmoA2X->subtract(PI5);
  

  

  // ===> Two electron part <=== //
  std::shared_ptr<psi::IntegralTransform> integrals_ = wfn_union_->integrals();
  dpd_set_default(integrals_->get_dpd_id());
  dpdbuf4 buf_Y122, buf_1122, buf_Y211, buf_Y212, buf_X211, buf_X122, buf_X121;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  //Term 1//
  global_dpd_->buf4_init(&buf_Y122, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[Y,1]"  ), integrals_->DPD_ID("[2,2]"  ),
                         integrals_->DPD_ID("[Y,1]"  ), integrals_->DPD_ID("[2>=2]+"), 0, "MO Ints (Y1|22)");
  global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0,
                          integrals_->DPD_ID("[1,1]"  ), integrals_->DPD_ID("[2,2]"  ),
                          integrals_->DPD_ID("[1>=1]+"), integrals_->DPD_ID("[2>=2]+"), 0, "MO Ints (11|22)");
  global_dpd_->buf4_init(&buf_Y211, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[Y,2]"  ), integrals_->DPD_ID("[1,1]"  ),
                         integrals_->DPD_ID("[Y,2]"  ), integrals_->DPD_ID("[1>=1]+"), 0, "MO Ints (Y2|11)");
  global_dpd_->buf4_init(&buf_Y212, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[Y,2]"  ), integrals_->DPD_ID("[1,2]"  ),
                         integrals_->DPD_ID("[Y,2]"  ), integrals_->DPD_ID("[1,2]"  ), 0, "MO Ints (Y2|12)");

  //Term 2//
  global_dpd_->buf4_init(&buf_X211, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,2]"  ), integrals_->DPD_ID("[1,1]"  ),
                         integrals_->DPD_ID("[X,2]"  ), integrals_->DPD_ID("[1>=1]+"), 0, "MO Ints (X2|11)");
  global_dpd_->buf4_init(&buf_X122, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,1]"  ), integrals_->DPD_ID("[2,2]"  ),
                         integrals_->DPD_ID("[X,1]"  ), integrals_->DPD_ID("[2>=2]+"), 0, "MO Ints (X1|22)");
  global_dpd_->buf4_init(&buf_X121, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,1]"  ), integrals_->DPD_ID("[2,1]"  ),
                         integrals_->DPD_ID("[X,1]"  ), integrals_->DPD_ID("[2,1]"  ), 0, "MO Ints (X1|21)");
  


  //Term 1//
  //ERI I //
  double int_1 = 0.0;
  double ERI1[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_Y122, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_Y122, h);
                 for (int pq = 0; pq < buf_Y122.params->rowtot[h]; ++pq) {
                      int p = buf_Y122.params->roworb[h][pq][0];
                      int q = buf_Y122.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_Y122.params->coltot[h]; ++rs) {
                           int r = buf_Y122.params->colorb[h][rs][0];
                           int s = buf_Y122.params->colorb[h][rs][1];
                           if ((p==n) && (q==i) && (r==s)) int_1 += buf_Y122.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_Y122, h);
            }
            ERI1[i][n] = 2.0 * int_1;
            int_1 = 0.0;
       }
  }

  //ERI II //
  double int_2 = 0.0;
  double ERI2[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
                 for (int pq = 0; pq < buf_1122.params->rowtot[h]; ++pq) {
                      int p = buf_1122.params->roworb[h][pq][0];
                      int q = buf_1122.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_1122.params->coltot[h]; ++rs) {
                           int r = buf_1122.params->colorb[h][rs][0];
                           int s = buf_1122.params->colorb[h][rs][1];
                           if ((p==i) && (r==s)) int_2 += S1[q][n] * buf_1122.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
            }
            ERI2[i][n] = 2.0 * int_2; 
            int_2 = 0.0;
       }
  }

  //ERI III //
  double int_3 = 0.0;
  double ERI3[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_Y211, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_Y211, h);
                 for (int pq = 0; pq < buf_Y211.params->rowtot[h]; ++pq) {
                      int p = buf_Y211.params->roworb[h][pq][0];
                      int q = buf_Y211.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_Y211.params->coltot[h]; ++rs) {
                           int r = buf_Y211.params->colorb[h][rs][0];
                           int s = buf_Y211.params->colorb[h][rs][1];
                           if ((p==n) && (r==s) && (i!=r)) int_3 += S2[i][q] * buf_Y211.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_Y211, h);
            }
            ERI3[i][n] = 2.0 * int_3; 
            int_3 = 0.0;
       }
  }

  //ERI IV //
  double int_4 = 0.0;
  double ERI4[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_Y212, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_Y212, h);
                 for (int pq = 0; pq < buf_Y212.params->rowtot[h]; ++pq) {
                      int p = buf_Y212.params->roworb[h][pq][0];
                      int q = buf_Y212.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_Y212.params->coltot[h]; ++rs) {
                           int r = buf_Y212.params->colorb[h][rs][0];
                           int s = buf_Y212.params->colorb[h][rs][1];
                           if ((p==n) && (r==i) && (q==s)) int_4 += buf_Y212.matrix[h][pq][rs];
                      }
		 }
		 global_dpd_->buf4_mat_irrep_close(&buf_Y212, h);
            } 
            ERI4[i][n] = int_4; 
            int_4 = 0.0;
       }
  }

  //ERI V //
  double int_5 = 0.0;
  double ERI5[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_Y211, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_Y211, h);
                 for (int pq = 0; pq < buf_Y211.params->rowtot[h]; ++pq) {
                      int p = buf_Y211.params->roworb[h][pq][0];
                      int q = buf_Y211.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_Y211.params->coltot[h]; ++rs) {
                           int r = buf_Y211.params->colorb[h][rs][0];
                           int s = buf_Y211.params->colorb[h][rs][1];
                           if ((p==n) && (r==i) && (i!=s)) int_5 += S2[s][q] * buf_Y211.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_Y211, h);
            }
            ERI5[i][n] = int_5;
            int_5 = 0.0;
       }
  }


  //Term 2//
  //ERI I //
  double int_6 = 0.0;
  double ERI6[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_X211, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_X211, h);
                 for (int pq = 0; pq < buf_X211.params->rowtot[h]; ++pq) {
                      int p = buf_X211.params->roworb[h][pq][0];
                      int q = buf_X211.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_X211.params->coltot[h]; ++rs) {
                           int r = buf_X211.params->colorb[h][rs][0];
                           int s = buf_X211.params->colorb[h][rs][1];
                           if ((p==n) && (q==i) && (r==s)) int_6 += buf_X211.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_X211, h);
            }
            ERI6[i][n] = 2.0 * int_6;
            int_6 = 0.0;
       }
  }

  //ERI II //
  double int_7 = 0.0;
  double ERI7[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
                 for (int pq = 0; pq < buf_1122.params->rowtot[h]; ++pq) {
                      int p = buf_1122.params->roworb[h][pq][0];
                      int q = buf_1122.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_1122.params->coltot[h]; ++rs) {
                           int r = buf_1122.params->colorb[h][rs][0];
                           int s = buf_1122.params->colorb[h][rs][1];
                           if ((p==q) && (r==i)) int_7 += S1_2[s][n] * buf_1122.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
            }
            ERI7[i][n] = 2.0 * int_7;
            int_7 = 0.0;
       }
  }

  //ERI III //
  double int_8 = 0.0;
  double ERI8[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_X122, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_X122, h);
                 for (int pq = 0; pq < buf_X122.params->rowtot[h]; ++pq) {
                      int p = buf_X122.params->roworb[h][pq][0];
                      int q = buf_X122.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_X122.params->coltot[h]; ++rs) {
                           int r = buf_X122.params->colorb[h][rs][0];
                           int s = buf_X122.params->colorb[h][rs][1];
                           if ((p==n) && (r==s) && (i!=r)) int_8 += S2_2[i][q] * buf_X122.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_X122, h);
            }
            ERI8[i][n] = 2.0 * int_8;
            int_8 = 0.0;
       }
  }

  //ERI IV //
  double int_9 = 0.0;
  double ERI9[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_X121, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_X121, h);
                 for (int pq = 0; pq < buf_X121.params->rowtot[h]; ++pq) {
                      int p = buf_X121.params->roworb[h][pq][0];
                      int q = buf_X121.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_X121.params->coltot[h]; ++rs) {
                           int r = buf_X121.params->colorb[h][rs][0];
                           int s = buf_X121.params->colorb[h][rs][1];
                           if ((p==n) && (r==i) && (q==s)) int_9 += buf_X121.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_X121, h);
            }
            ERI9[i][n] = int_9;
            int_9 = 0.0;
       }
  }

  //ERI V //
  double int_10 = 0.0;
  double ERI10[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            for (int h = 0; h < wfn_union_->nirrep(); ++h) {
                 global_dpd_->buf4_mat_irrep_init(&buf_X122, h);
                 global_dpd_->buf4_mat_irrep_rd(&buf_X122, h);
                 for (int pq = 0; pq < buf_X122.params->rowtot[h]; ++pq) {
                      int p = buf_X122.params->roworb[h][pq][0];
                      int q = buf_X122.params->roworb[h][pq][1];
                      for (int rs = 0; rs < buf_X122.params->coltot[h]; ++rs) {
                           int r = buf_X122.params->colorb[h][rs][0];
                           int s = buf_X122.params->colorb[h][rs][1];
                           if ((p==n) && (r==i) && (i!=s)) int_10 += S2_2[s][q] * buf_X122.matrix[h][pq][rs];
                      }
                 }
                 global_dpd_->buf4_mat_irrep_close(&buf_X122, h);
            }
            ERI10[i][n] = int_10;
            int_10 = 0.0;
       }
  }
  

  global_dpd_->buf4_close(&buf_Y122);
  global_dpd_->buf4_close(&buf_1122);
  global_dpd_->buf4_close(&buf_Y211);
  global_dpd_->buf4_close(&buf_Y212);

  global_dpd_->buf4_close(&buf_X211);
  global_dpd_->buf4_close(&buf_X122);
  global_dpd_->buf4_close(&buf_X121);


  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);


  //Term 1: Sum of ERIs//
  double ERI_1[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];
  double ERI_1_group_1[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];   // debug for OEP
  double ERI_1_group_2[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];   // debug for OEP
  double ERI_1_group_3[wfn_union_->l_ndocc(0)][wfn_union_->l_nvir(1)];   // debug for OEP
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
            ERI_1[i][n] = ERI1[i][n] - ERI2[i][n] - ERI3[i][n] - ERI4[i][n] + ERI5[i][n];  // All Groups
            ERI_1_group_1[i][n] = ERI1[i][n] - ERI4[i][n]; // Group I: debug for OEP
       }
  }


  //Term 2: Sum of ERIs//
  double ERI_2[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];
  double ERI_2_group_1[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];   // debug for OEP
  double ERI_2_group_2[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];   // debug for OEP
  double ERI_2_group_3[wfn_union_->l_ndocc(1)][wfn_union_->l_nvir(0)];   // debug for OEP
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            ERI_2[i][n] = ERI6[i][n] - ERI7[i][n] - ERI8[i][n] - ERI9[i][n] + ERI10[i][n];  // All Groups
            ERI_2_group_1[i][n] = ERI6[i][n] - ERI9[i][n]; // Group I: debug for OEP
       }
  }


  //Term 1//
  std::shared_ptr<psi::Vector> Eps_occ_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> Eps_vir_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","VIR");

  //Term 2//
  std::shared_ptr<psi::Vector> Eps_occ_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> Eps_vir_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","VIR");


  //Term 1: CT Energy//
  double E_ct_1 = 0.0;
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
	    double v = ERI_1_group_1[i][n] + VmoB1Y_group_1->get(i,n);
            E_ct_1 += (v*v)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
       }
  }
  E_ct_1 = 2.0 * E_ct_1; 


  //Term 2: CT Energy//
  double E_ct_2 = 0.0;
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            double v = ERI_2_group_1[i][n] + VmoA2X_group_1->get(i,n);
            E_ct_2 += (v*v)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
       }
  }
  E_ct_2 = 2.0 * E_ct_2;



  //Whole CT Energy E_CT(A+B-) + E_CT(A-B+)//
  double E_ct = E_ct_1 + E_ct_2;
  


  //psi::timer_off("SOLVER: Charge-transfer Energy Calculations (Otto-Ladik)");
  psi::timer_off("Solver E(CT) Otto-Ladik       ");


  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Charge-transfer energy calculations <==\n");
     psi::outfile->Printf("  ==>         Benchmark (Otto-Ladik)           <==\n\n");
     psi::outfile->Printf("     E_CT (A+B-)   = %13.10f au\n", E_ct_1        );
     psi::outfile->Printf("     E_CT (A-B+)   = %13.10f au\n", E_ct_2        );
     psi::outfile->Printf("     E_CT          = %13.10f au\n", E_ct        ); 
     psi::outfile->Printf("\n");
  }

  return E_ct;


}
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
  oep_1->compute("Murrell-etal.V1");
  oep_2->compute("Murrell-etal.V1");

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
  std::shared_ptr<psi::Vector> e_occ_1  = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> e_occ_2  = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> e_vir_1  = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","VIR");
  std::shared_ptr<psi::Vector> e_vir_2  = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","VIR");

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
  double e_ab_v1 = compute_ct_component(e_occ_1, e_vir_2, v_ab_v1);
  double e_ba_v1 = compute_ct_component(e_occ_2, e_vir_1, v_ba_v1);

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
double ChargeTransferEnergySolver::compute_ct_component(
 std::shared_ptr<psi::Vector> eps_occ_X, 
 std::shared_ptr<psi::Vector> eps_vir_Y, 
 std::shared_ptr<psi::Matrix> V)
{
   // requires matrix elements in CMO (canonical SCF) basis
   const int nocc_X = eps_occ_X->dim();
   const int nvir_Y = eps_vir_Y->dim();
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
