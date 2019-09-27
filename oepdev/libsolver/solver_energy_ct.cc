#include "solver.h"
#include "psi4/libpsi4util/process.h"
#include <ctime>

using namespace std;
using namespace psi;
using namespace oepdev;

// CT Solver//
ChargeTransferEnergySolver::ChargeTransferEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Sanity check
  if (options_.get_bool("OEPDEV_LOCALIZE")) 
      throw psi::PSIEXCEPTION("Error. OEPDEV_LOCALIZE must be set to False because CMO's are necessary");
  // Benchmarks
  methods_benchmark_.push_back("OTTO_LADIK");
  methods_benchmark_.push_back("EFP2"      );
  // OEP-based
  methods_oepBased_ .push_back("OTTO_LADIK");
}
ChargeTransferEnergySolver::~ChargeTransferEnergySolver() {}
double ChargeTransferEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "OTTO_LADIK"      )  e = compute_oep_based_murrell_etal();
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
           method == "OTTO_LADIK" ) e = compute_benchmark_murrell_etal();
  else if (method == "EFP2"         ) e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
double ChargeTransferEnergySolver::compute_benchmark_murrell_etal(){

  //psi::timer_on("SOLVER: Charge-transfer Energy Calculations (Otto-Ladik)");
  //psi::timer_on("Solver E(CT) Otto-Ladik       ");
  

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
  std::shared_ptr<psi::Matrix> VmoB1Y_group_1 = VmoB1Y->clone();  // debug for OEP
  std::shared_ptr<psi::Matrix> VmoB1Y_group_2 = PI2   ->clone();  // debug for OEP 
  std::shared_ptr<psi::Matrix> VmoB1Y_group_3 = PI3   ->clone();  // debug for OEP
  VmoB1Y_group_2->scale(-1.0); VmoB1Y_group_3->scale(-1.0);
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
  std::shared_ptr<psi::Matrix> VmoA2X_group_1 = VmoA2X->clone();   // debug for OEP
  std::shared_ptr<psi::Matrix> VmoA2X_group_2 = PI4   ->clone();   // debug for OEP 
  std::shared_ptr<psi::Matrix> VmoA2X_group_3 = PI5   ->clone();   // debug for OEP
  VmoA2X_group_2->scale(-1.0); VmoA2X_group_3->scale(-1.0);
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
                         integrals_->DPD_ID("[Y,I]"  ), integrals_->DPD_ID("[J,J]"  ),
                         integrals_->DPD_ID("[Y,I]"  ), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (YI|JJ)");
  global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0,
                          integrals_->DPD_ID("[I,I]"  ), integrals_->DPD_ID("[J,J]"  ),
                          integrals_->DPD_ID("[I>=I]+"), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (II|JJ)");
  global_dpd_->buf4_init(&buf_Y211, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I,I]"  ),
                         integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I>=I]+"), 0, "MO Ints (YJ|II)");
  global_dpd_->buf4_init(&buf_Y212, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I,J]"  ),
                         integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I,J]"  ), 0, "MO Ints (YJ|IJ)");

  //Term 2//
  global_dpd_->buf4_init(&buf_X211, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,J]"  ), integrals_->DPD_ID("[I,I]"  ),
                         integrals_->DPD_ID("[X,J]"  ), integrals_->DPD_ID("[I>=I]+"), 0, "MO Ints (XJ|II)");
  global_dpd_->buf4_init(&buf_X122, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J,J]"  ),
                         integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (XI|JJ)");
  global_dpd_->buf4_init(&buf_X121, PSIF_LIBTRANS_DPD, 0,
                         integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J,I]"  ),
                         integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J,I]"  ), 0, "MO Ints (XI|JI)");
  


  //Term 1//
  //ERI I // Group: I
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

  //ERI II // Group II
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

  //ERI III // Group III
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

  //ERI IV // Group I
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

  //ERI V // Group III
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
  //ERI I // Group I
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

  //ERI II // Group II
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

  //ERI III // Group III
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

  //ERI IV // Group I
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

  //ERI V // Group III
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
            ERI_1_group_2[i][n] =-ERI2[i][n]             ; // Group II: debug for OEP
            ERI_1_group_3[i][n] =-ERI3[i][n] + ERI5[i][n]; // Group III: debug for OEP
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
            ERI_2_group_2[i][n] =-ERI7[i][n]             ; // Group II: debug for OEP
            ERI_2_group_3[i][n] =-ERI8[i][n] + ERI10[i][n];// Group III: debug for OEP
       }
  }


  //Term 1//
  std::shared_ptr<psi::Vector> Eps_occ_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> Eps_vir_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","VIR");

  //Term 2//
  std::shared_ptr<psi::Vector> Eps_occ_B = wfn_union_->l_wfn(1)->epsilon_a_subset("MO","OCC");
  std::shared_ptr<psi::Vector> Eps_vir_A = wfn_union_->l_wfn(0)->epsilon_a_subset("MO","VIR");


  //Term 1: CT Energy A ---> B //
  double E_ct_1 = 0.0;
  double E_ct_1_12 = 0.0, E_ct_1_13 = 0.0, E_ct_1_23 = 0.0;
  double E_ct_1_1 = 0.0, E_ct_1_2 = 0.0, E_ct_1_3 = 0.0;
  for (int i=0; i<wfn_union_->l_ndocc(0); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(1); ++n) {
	    double v = ERI_1[i][n] + VmoB1Y->get(i,n);
            E_ct_1 += (v*v)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            //
	    double v1= ERI_1_group_1[i][n] + VmoB1Y_group_1->get(i,n);
	    double v2= ERI_1_group_2[i][n] + VmoB1Y_group_2->get(i,n);
	    double v3= ERI_1_group_3[i][n] + VmoB1Y_group_3->get(i,n);
            double v123 = v1 + v2 + v3;
            double v12  = v1 + v2;
            double v13  = v1 + v3;
            double v23  = v2 + v3;
            E_ct_1_1 += (v1*v1)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            E_ct_1_2 += (v2*v2)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            E_ct_1_3 += (v3*v3)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            E_ct_1_12 += (v12*v12)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            E_ct_1_13 += (v13*v13)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
            E_ct_1_23 += (v23*v23)/(Eps_occ_A->get(i) - Eps_vir_B->get(n));
       }
  }
  E_ct_1 *= 2.0; 
  E_ct_1_1 *= 2.0;
  E_ct_1_2 *= 2.0;
  E_ct_1_3 *= 2.0;
  E_ct_1_12 *= 2.0;
  E_ct_1_13 *= 2.0;
  E_ct_1_23 *= 2.0;


  //Term 2: CT Energy B ---> A //
  double E_ct_2 = 0.0;
  double E_ct_2_12 = 0.0, E_ct_2_13 = 0.0, E_ct_2_23 = 0.0;
  double E_ct_2_1 = 0.0, E_ct_2_2 = 0.0, E_ct_2_3 = 0.0;
  for (int i=0; i<wfn_union_->l_ndocc(1); ++i) {
       for (int n=0; n<wfn_union_->l_nvir(0); ++n) {
            double v = ERI_2[i][n] + VmoA2X->get(i,n);
            E_ct_2 += (v*v)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            //
            double v1= ERI_2_group_1[i][n] + VmoA2X_group_1->get(i,n);
            double v2= ERI_2_group_2[i][n] + VmoA2X_group_2->get(i,n);
            double v3= ERI_2_group_3[i][n] + VmoA2X_group_3->get(i,n);
            double v123 = v1 + v2 + v3;
            double v12  = v1 + v2;
            double v13  = v1 + v3;
            double v23  = v2 + v3;
            E_ct_2_1 += (v1*v1)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            E_ct_2_2 += (v2*v2)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            E_ct_2_3 += (v3*v3)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            E_ct_2_12 += (v12*v12)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            E_ct_2_13 += (v13*v13)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
            E_ct_2_23 += (v23*v23)/(Eps_occ_B->get(i) - Eps_vir_A->get(n));
       }
  }
  E_ct_2 *= 2.0;
  E_ct_2_1 *= 2.0;
  E_ct_2_2 *= 2.0;
  E_ct_2_3 *= 2.0;
  E_ct_2_12 *= 2.0;
  E_ct_2_13 *= 2.0;
  E_ct_2_23 *= 2.0;


  //Whole CT Energy E_CT(A+B-) + E_CT(A-B+)//
  double E_ct = E_ct_1 + E_ct_2;
  psi::Process::environment.globals["EINT CT OTTO-LADIK KCAL"] = E_ct*OEPDEV_AU_KcalPerMole;

  //psi::timer_off("SOLVER: Charge-transfer Energy Calculations (Otto-Ladik)");
  //psi::timer_off("Solver E(CT) Otto-Ladik       ");


  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > -1) {
     psi::outfile->Printf("  ==> SOLVER: Charge-transfer energy calculations <==\n");
     psi::outfile->Printf("  ==>         Benchmark (Otto-Ladik)           <==\n\n");
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("     Group I\n"                                   );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_1      );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_1      );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("     Group II\n"                                  );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_2      );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_2      );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("     Group III\n"                                 );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_3      );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_3      );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     ===============================\n"           );
     psi::outfile->Printf("     Group I+II\n"                                );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_12     );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_12     );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("     Group I+III\n"                               );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_13     );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_13     );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("     Group II+III\n"                              );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1_23     );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2_23     );
   //psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     ===============================\n"           );
     psi::outfile->Printf("     Total\n"                                     );
     psi::outfile->Printf("     E_CT (A-->B)   = %13.6f au\n", E_ct_1        );
     psi::outfile->Printf("     E_CT (B-->A)   = %13.6f au\n", E_ct_2        );
     psi::outfile->Printf("     E_CT           = %13.6f au\n", E_ct          ); 
     psi::outfile->Printf("     -------------------------------\n"           );
     psi::outfile->Printf("\n");
  }

  return E_ct;


}
double ChargeTransferEnergySolver::compute_benchmark_efp2()
  {
  double e_ct  = 0.0;
  

  // ---> Timer-on <--- //
  psi::timer_on("Solver E(CT) EFP2 MO-Expanded");
  
  int nbf_1     = wfn_union_->l_nbf(0);
  int nbf_2	= wfn_union_->l_nbf(1);  
  int ndocc_1	= wfn_union_->l_ndocc(0);
  int ndocc_2   = wfn_union_->l_ndocc(1);
  int nvir_1	= wfn_union_->l_nvir(0);
  int nvir_2	= wfn_union_->l_nvir(1);

  // ===> ONE-ELECTRON PART <=== //

  // V matrices //
  std::shared_ptr<psi::Matrix> VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> VaoB11    = std::make_shared<psi::Matrix>("VaoB(1,1)" , nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> VaoA21    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_1);
  std::shared_ptr<psi::Matrix> VaoA22    = std::make_shared<psi::Matrix>("VaoA(2,2)" , nbf_2, nbf_2);

  // S matrices //
  std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)"  , nbf_1, nbf_2);

  // T matrices //
  std::shared_ptr<psi::Matrix> Tao12     = std::make_shared<psi::Matrix>("Tao(1,2)", nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> Tao11     = std::make_shared<psi::Matrix>("Tao(1,1)", nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> Tao22     = std::make_shared<psi::Matrix>("Tao(2,2)", nbf_2, nbf_2);

  // F matrices //
  std::shared_ptr<psi::Matrix> Fao11     = wfn_union_->l_wfn(0)->Fa();
  std::shared_ptr<psi::Matrix> Fao22     = wfn_union_->l_wfn(1)->Fa();

  // Ca matrices //
  std::shared_ptr<psi::Matrix> Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_vir_A = wfn_union_->l_wfn(0)->Ca_subset("AO","VIR");
  std::shared_ptr<psi::Matrix> Ca_vir_B = wfn_union_->l_wfn(1)->Ca_subset("AO","VIR");

  // IntegralFactory //
  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_11(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_22(wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));

  // ---> Type of effective potential matrix elements <--- //
  enum t_potential_type {ERI_Potential_Type, DMTP_Potential_Type};
  t_potential_type potential_type;
  if (options_.get_str("EFP2_CT_POTENTIAL_INTS") == "DMTP") 
       potential_type = DMTP_Potential_Type;
  else potential_type =  ERI_Potential_Type;

  clock_t t_time = -clock(); // Clock BEGIN

  if (potential_type == ERI_Potential_Type) {
     // PotentialInt - nuclear part //
     std::shared_ptr<psi::PotentialInt> potInt_12 = std::make_shared<psi::PotentialInt>(fact_12.spherical_transform(), 
                                                                                       wfn_union_->l_primary(0),
                                                                                       wfn_union_->l_primary(1));
     std::shared_ptr<psi::PotentialInt> potInt_21 = std::make_shared<psi::PotentialInt>(fact_21.spherical_transform(),
                                                                                       wfn_union_->l_primary(1),
                                                                                       wfn_union_->l_primary(0));
     std::shared_ptr<psi::PotentialInt> potInt_11 = std::make_shared<psi::PotentialInt>(fact_11.spherical_transform(),
                                                                                       wfn_union_->l_primary(0),
                                                                                       wfn_union_->l_primary(0));
     std::shared_ptr<psi::PotentialInt> potInt_22 = std::make_shared<psi::PotentialInt>(fact_22.spherical_transform(),
                                                                                       wfn_union_->l_primary(1),
                                                                                       wfn_union_->l_primary(1));
                                                                                                                       
                                                                                                                       
     // Set charge field //
     std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_11->charge_field());
     std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_22->charge_field());
                                                                                                                       
     potInt_12->set_charge_field(Zxyz_2);
     potInt_11->set_charge_field(Zxyz_2);
     potInt_21->set_charge_field(Zxyz_1);
     potInt_22->set_charge_field(Zxyz_1);
                                                                                                                       
     // Potential integrals (nuclear contribution) //
     std::shared_ptr<psi::OneBodyAOInt> oneInt;
     oneInt = potInt_12;
     oneInt->compute(VaoB12);
     oneInt = potInt_11;
     oneInt->compute(VaoB11);
     oneInt = potInt_21;
     oneInt->compute(VaoA21);
     oneInt = potInt_22;
     oneInt->compute(VaoA22);
  }
  t_time += clock(); // Clock END

  // Overlap integrals //
  std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_12.ao_overlap());
  ovlInt->compute(Sao12);

  // Kinetic energy integrals //
  std::shared_ptr<psi::OneBodyAOInt> kinInt12(fact_12.ao_kinetic());
  kinInt12->compute(Tao12); 

  t_time -= clock(); // Clock BEGIN
  std::shared_ptr<psi::OneBodyAOInt> kinInt11(fact_11.ao_kinetic());
  kinInt11->compute(Tao11);

  std::shared_ptr<psi::OneBodyAOInt> kinInt22(fact_22.ao_kinetic());
  kinInt22->compute(Tao22);
  t_time += clock(); // Clock END

  if (potential_type == DMTP_Potential_Type) {
      // Compute CAMM on monomers 1 and 2
      std::shared_ptr<oepdev::DMTPole> camm_1 = oepdev::DMTPole::build("CAMM", wfn_union_->l_wfn(0));
      std::shared_ptr<oepdev::DMTPole> camm_2 = oepdev::DMTPole::build("CAMM", wfn_union_->l_wfn(1));
      camm_1->compute();
      camm_2->compute();
      //
      const double p1 = 1.0 / 3.0;
      const double p2 = 2.0 / 3.0;
      const double p3 = 1.0 /15.0;
      const double p4 = 3.0 /15.0;
      const double p5 = 6.0 /15.0;
      const double prefacs[20] = {
     /* 0    X    Y    Z    XX  YY  ZZ  XY  XZ  YZ */
        1.0, 1.0, 1.0, 1.0, p1, p1, p1, p2, p2, p2, 
     /*    XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ */
           p3, p3, p3, p4, p4, p4, p4, p4, p4, p5};
      //
      t_time -= clock(); // Clock BEGIN

      /* __12  ->   V from molecule 2 
         __11  ->   V from molecule 2
         __21  ->   V from molecule 1
         __22  ->   V from molecule 1 */

      // 
      size_t n_multipole_1 = wfn_union_->l_molecule(0)->natom();
      size_t n_multipole_2 = wfn_union_->l_molecule(1)->natom();
      assert (n_multipole_1 == camm_1->n_sites());
      assert (n_multipole_2 == camm_2->n_sites());
      auto xyz_1 = this->extract_xyz(wfn_union_->l_molecule(0));
      auto xyz_2 = this->extract_xyz(wfn_union_->l_molecule(1));
      auto mult_1= this->extract_dmtp(camm_1);
      auto mult_2= this->extract_dmtp(camm_2);
      //
      std::shared_ptr<psi::OneBodyAOInt> efp_ints_12(fact_12.ao_efp_multipole_potential());
      std::shared_ptr<psi::OneBodyAOInt> efp_ints_11(fact_11.ao_efp_multipole_potential());
      std::shared_ptr<psi::OneBodyAOInt> efp_ints_21(fact_21.ao_efp_multipole_potential());
      std::shared_ptr<psi::OneBodyAOInt> efp_ints_22(fact_22.ao_efp_multipole_potential());
      //
      std::vector<psi::SharedMatrix> mats_12, mats_11, mats_21, mats_22;
      for (int i = 0; i < 20; ++i) {
           mats_12.push_back(std::make_shared<psi::Matrix>("", nbf_1, nbf_2));
           mats_11.push_back(std::make_shared<psi::Matrix>("", nbf_1, nbf_1));
           mats_21.push_back(std::make_shared<psi::Matrix>("", nbf_2, nbf_1));
           mats_22.push_back(std::make_shared<psi::Matrix>("", nbf_2, nbf_2));
      }
      //
      double *xyz_1_p = xyz_1->pointer();
      double *xyz_2_p = xyz_2->pointer();
      double *mult_1_p = mult_1->pointer();
      double *mult_2_p = mult_2->pointer();

      // Molecule B CAMM
      for (size_t n2 = 0; n2 < n_multipole_2; n2++) {

           for (int i = 0; i < 20; ++i) {
                mats_12[i]->zero();
                mats_11[i]->zero();
           }
           psi::Vector3 coords_2(xyz_2_p[n2 * 3], xyz_2_p[n2 * 3 + 1], xyz_2_p[n2 * 3 + 2]);
           efp_ints_12->set_origin(coords_2);
           efp_ints_11->set_origin(coords_2);
           efp_ints_12->compute(mats_12);
           efp_ints_11->compute(mats_11);

           for (int i = 0; i < 20; ++i) {                            
                mats_12[i]->scale(-prefacs[i] * mult_2_p[20 * n2 + i]);
                mats_11[i]->scale(-prefacs[i] * mult_2_p[20 * n2 + i]);
                VaoB12->add(mats_12[i]);
                VaoB11->add(mats_11[i]);
           }
      }

      // Molecule A CAMM
      for (size_t n1 = 0; n1 < n_multipole_1; n1++) {

           for (int i = 0; i < 20; ++i) {
                mats_21[i]->zero();
                mats_22[i]->zero();
           }
           psi::Vector3 coords_1(xyz_1_p[n1 * 3], xyz_1_p[n1 * 3 + 1], xyz_1_p[n1 * 3 + 2]);
           efp_ints_21->set_origin(coords_1);
           efp_ints_22->set_origin(coords_1);
           efp_ints_21->compute(mats_21);
           efp_ints_22->compute(mats_22);

           for (int i = 0; i < 20; ++i) {                            
                mats_21[i]->scale(-prefacs[i] * mult_1_p[20 * n1 + i]);
                mats_22[i]->scale(-prefacs[i] * mult_1_p[20 * n1 + i]);
                VaoA21->add(mats_21[i]);
                VaoA22->add(mats_22[i]);
           }
      }
      //
      t_time += clock(); // Clock END
  }


  // ---> Transform one electron contributions to MO basis <--- //
   

  // Transform S matrices //
  t_time -= clock(); // Clock BEGIN
  std::shared_ptr<psi::Matrix> Smoij    = psi::Matrix::triplet(Ca_occ_A, Sao12,  Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> Smoxn    = psi::Matrix::triplet(Ca_occ_A, Sao12,  Ca_vir_B, true, false, false); // x \in OCC_A
  std::shared_ptr<psi::Matrix> Smoyn    = psi::Matrix::triplet(Ca_vir_A, Sao12,  Ca_vir_B, true, false, false); // y \in VIR_A
  
  std::shared_ptr<psi::Matrix> Smoxm    = psi::Matrix::triplet(Ca_occ_B, Sao12,  Ca_vir_A, true, true , false); // x \in OCC_B
  std::shared_ptr<psi::Matrix> Smoym    = psi::Matrix::triplet(Ca_vir_B, Sao12,  Ca_vir_A, true, true , false); // y \in VIR_B
                                                                                                                                 

  // Transform T matrices //
  std::shared_ptr<psi::Matrix> Tmonn    = psi::Matrix::triplet(Ca_vir_B, Tao22,  Ca_vir_B, true, false, false);
  std::shared_ptr<psi::Matrix> Tmonj    = psi::Matrix::triplet(Ca_vir_B, Tao22,  Ca_occ_B, true, false, false);                            
  std::shared_ptr<psi::Matrix> Tmoxj    = psi::Matrix::triplet(Ca_occ_A, Tao12,  Ca_occ_B, true, false, false);  // x \in OCC_A (T_mj)
  std::shared_ptr<psi::Matrix> Tmoyj    = psi::Matrix::triplet(Ca_vir_A, Tao12,  Ca_occ_B, true, false, false);  // y \in VIR_A (T_mj)          
  std::shared_ptr<psi::Matrix> Tmomm    = psi::Matrix::triplet(Ca_vir_A, Tao11,  Ca_vir_A, true, false, false);
  std::shared_ptr<psi::Matrix> Tmomi    = psi::Matrix::triplet(Ca_vir_A, Tao11,  Ca_occ_A, true, false, false);  
  std::shared_ptr<psi::Matrix> Tmoxi    = psi::Matrix::triplet(Ca_occ_B, Tao12,  Ca_occ_A, true, true , false);  // x \in OCC_B (T_ni)
  std::shared_ptr<psi::Matrix> Tmoyi    = psi::Matrix::triplet(Ca_vir_B, Tao12,  Ca_occ_A, true, true , false);  // y \in VIR_B (T_ni)

                                                                                                                                 
  // Transform F matrices //
  std::shared_ptr<psi::Matrix> Fmoii    = psi::Matrix::triplet(Ca_occ_A, Fao11,  Ca_occ_A, true, false, false);
  std::shared_ptr<psi::Matrix> Fmojj    = psi::Matrix::triplet(Ca_occ_B, Fao22,  Ca_occ_B, true, false, false);


   // Transform V matrices //                                                                                                   
  std::shared_ptr<psi::Matrix> VmoBin   = psi::Matrix::triplet(Ca_occ_A, VaoB12, Ca_vir_B, true, false, false); 
  std::shared_ptr<psi::Matrix> VmoBix   = psi::Matrix::triplet(Ca_occ_A, VaoB11, Ca_occ_A, true, false, false); // x \in OCC_A (V_im)
  std::shared_ptr<psi::Matrix> VmoBiy   = psi::Matrix::triplet(Ca_occ_A, VaoB11, Ca_vir_A, true, false, false); // y \in VIR_A (V_im)

  std::shared_ptr<psi::Matrix> VmoAjm   = psi::Matrix::triplet(Ca_occ_B, VaoA21, Ca_vir_A, true, false, false); 
  std::shared_ptr<psi::Matrix> VmoAjx   = psi::Matrix::triplet(Ca_occ_B, VaoA22, Ca_occ_B, true, false, false); // x \in OCC_B (V_jn)
  std::shared_ptr<psi::Matrix> VmoAjy   = psi::Matrix::triplet(Ca_occ_B, VaoA22, Ca_vir_B, true, false, false); // y \in VIR_B (V_jn)


  // ==> TWO-ELECTRON PART (electronic contibutions to potentials) <== //
  // --> add electron contributions to VmoA** and VmoB** (6 matrix elements) <-- //

  if (potential_type == ERI_Potential_Type) { /* Switch off 2-electron part for debug purposes */

      std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals();                                                                                                      
      dpd_set_default(integrals->get_dpd_id());
      dpdbuf4 buf_Y122, buf_X122, buf_1122; // components of electronic contribution to V^B
      dpdbuf4 buf_X211, buf_Y211; 		// components of electronic contribution to V^A
      //dpdbuf4 buf_2211; // BBB-T
      std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
                                                                                                                                                                                        
      // -- Open DPD library to calculate ERI <--
      psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
                                                                                                                                                                                        
      // --> Generate integral buffers <-- 
      // Y122 = (V_2 O_1 | O_2 O_2) = (ni|jj)
      global_dpd_->buf4_init(&buf_Y122, PSIF_LIBTRANS_DPD, 0,
                             integrals->DPD_ID("[Y,I]"  ), integrals->DPD_ID("[J,J]"  ),
                             integrals->DPD_ID("[Y,I]"  ), integrals->DPD_ID("[J>=J]+"), 0, "MO Ints (YI|JJ)");
      // X122 = (V_1 O_1 | O_2 O_2) = (mi|jj)
      global_dpd_->buf4_init(&buf_X122, PSIF_LIBTRANS_DPD, 0,
            	         integrals->DPD_ID("[X,I]"  ), integrals->DPD_ID("[J,J]"  ),
            		 integrals->DPD_ID("[X,I]"  ), integrals->DPD_ID("[J>=J]+"), 0, "MO Ints (XI|JJ)");
      // 1122 = (O_1 O_1 | O_2 O_2) = (11|22)
      global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0,
            	         integrals->DPD_ID("[I,I]"  ), integrals->DPD_ID("[J,J]"  ),
            		 integrals->DPD_ID("[I>=I]+"), integrals->DPD_ID("[J>=J]+"), 0, "MO Ints (II|JJ)");
      // Y211 = (V_2 O_2 | O_1 O_1) = (nj|kk)
      global_dpd_->buf4_init(&buf_Y211, PSIF_LIBTRANS_DPD, 0,
                             integrals->DPD_ID("[Y,J]"  ), integrals->DPD_ID("[I,I]"  ),
                             integrals->DPD_ID("[Y,J]"  ), integrals->DPD_ID("[I>=I]+"), 0, "MO Ints (YJ|II)");
      // X211 = (V_1 O_2 | O_1 O_1) = (mj|kk)
      global_dpd_->buf4_init(&buf_X211, PSIF_LIBTRANS_DPD, 0,
            	         integrals->DPD_ID("[X,J]"  ), integrals->DPD_ID("[I,I]"  ),
            		 integrals->DPD_ID("[X,J]"  ), integrals->DPD_ID("[I>=I]+"), 0, "MO Ints (XJ|II)");

      // Potentials (nuc. + elec.) included in expression of E_CT(A->B)
      double** vmoBin = VmoBin->pointer();
      double** vmoBix = VmoBix->pointer();
      double** vmoBiy = VmoBiy->pointer();

      // Potentials (nuc. + elec.)included in expression of E_CT(B->A)
      double** vmoAjm = VmoAjm->pointer();
      double** vmoAjx = VmoAjx->pointer();
      double** vmoAjy = VmoAjy->pointer();

      psi::outfile->Printf("ndocc_1: %2d, ndocc_2: %2d, nvir_1: %2d, nvir_2:  %2d \n", ndocc_1, ndocc_2, nvir_1, nvir_2);

      // vmoBin (elec.) = (in|jj) = (1Y|22) = (Y1|22)
      for (int h = 0; h < wfn_union_->nirrep(); ++h) 
           {
           global_dpd_->buf4_mat_irrep_init(&buf_Y122, h);
           global_dpd_->buf4_mat_irrep_rd(&buf_Y122, h);
           for (int pq = 0; pq < buf_Y122.params->rowtot[h]; ++pq) 
                {
                // pq = max p x max q 
                int p = buf_Y122.params->roworb[h][pq][0]; // indx (0;nvir_2-1) Y
                int q = buf_Y122.params->roworb[h][pq][1]; // indx (0;nocc_1-1) 1
//    	    psi::outfile->Printf("pq: %d, p: %2d, q: %2d \n", pq, p, q);
                for (int rs = 0; rs < buf_Y122.params->coltot[h]; ++rs) 
                	 {
            	 int r = buf_Y122.params->colorb[h][rs][0]; // indx (0;nocc_2-1) 2
            	 int s = buf_Y122.params->colorb[h][rs][1]; // indx (0;nocc_2-1) 2
                  //   psi::outfile->Printf("rs: %d, r: %2d, s: %2d \n", rs r, s); 
                     if (r == s) 
                     	{
            		vmoBin[q][p] += 2.0*buf_Y122.matrix[h][pq][rs]; 
                     	//psi::outfile->Printf("To V^B_in[%d][%d] added: 2x(%2d %2d | %2d %2d) = %16.10f\n", q, p, p, q, r, s, buf_Y122.matrix[h][pq][rs]);
            	 	}
            	 }
                }
           global_dpd_->buf4_mat_irrep_close(&buf_Y122, h);
           }
      global_dpd_->buf4_close(&buf_Y122);
                                                                                                                                                                                        
      // vmoBix (elec.) = 2(ix|jj) = 2(11|22)  BB
      // vmoAjx (elec.) = 2(jx|kk) = 2(22|11)  BB
      for (int h = 0; h < wfn_union_->nirrep(); ++h) 
           {
           global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
           global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
           for (int pq = 0; pq < buf_1122.params->rowtot[h]; ++pq) 
                {
                int p = buf_1122.params->roworb[h][pq][0]; // ndocc_1
                int q = buf_1122.params->roworb[h][pq][1]; // ndocc_1
                                                                                                                                                                                        
                //if (p == q) {
                //    for (int rs = 0; rs < buf_1122.params->coltot[h]; ++rs) {
                //         int r = buf_1122.params->colorb[h][rs][0]; // ndocc_2
                //         int s = buf_1122.params->colorb[h][rs][1]; // ndocc_2
                //         vmoAjx[r][s] += 2.0 * buf_1122.matrix[h][pq][rs];
                //    }
                //}
//    	    psi::outfile->Printf("pq: %d, p: %2d, q: %2d \n", pq, p, q);
                for (int rs = 0; rs < buf_1122.params->coltot[h]; ++rs) 
            	 {
            	 int r = buf_1122.params->colorb[h][rs][0]; // ndocc_2
            	 int s = buf_1122.params->colorb[h][rs][1]; // ndocc_2
            	 // psi::outfile->Printf("rs: %d, r: %2d, s: %2d \n", rs, r, s); 
//    		 if (p==q)
            	 //if ((p == q) && (r == s))  BB
                     if (r == s)
            		  {
            		  vmoBix[p][q] += 2.0 * buf_1122.matrix[h][pq][rs];
            		  //psi::outfile->Printf("To V^B_ix[%d][%d] added: 2x(%2d %2d | %2d %2d) = %16.10f\n", p, q, p, q, r, s, buf_1122.matrix[h][pq][rs]);
            		  } 
            	 //if ((p == q) && (r == s)) BB
      	         if (p == q)
            	          {
            	          vmoAjx[r][s] += 2.0 * buf_1122.matrix[h][pq][rs];
            	          //psi::outfile->Printf("To V^A_jx[%d][%d] added: 2x(%2d %2d | %2d %2d) = %16.10f\n", r, s, p, q, r, s, buf_1122.matrix[h][pq][rs]);			   
            	          }
                     }
                }
           global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
           }
      global_dpd_->buf4_close(&buf_1122);

      // vmoBiy (elec.) = 2(im|jj) = 2(1X|22) = 2(X1|22) = 2(mi|jj)
      for (int h = 0; h < wfn_union_->nirrep(); ++h) 
           {
           global_dpd_->buf4_mat_irrep_init(&buf_X122, h);
           global_dpd_->buf4_mat_irrep_rd(&buf_X122, h);
           for (int pq = 0; pq < buf_X122.params->rowtot[h]; ++pq) 
                {
                int p = buf_X122.params->roworb[h][pq][0]; // nvir_1 X
                int q = buf_X122.params->roworb[h][pq][1]; // nocc_1 1
//    	    psi::outfile->Printf("pq: %d, p: %2d, q: %2d \n", pq, p, q);
                for (int rs = 0; rs < buf_X122.params->coltot[h]; ++rs) 
                     {
            	 int r = buf_X122.params->colorb[h][rs][0]; // nocc_2 2
            	 int s = buf_X122.params->colorb[h][rs][1]; // nocc_2 2
//    		 psi::outfile->Printf("rs: %d, r: %2d, s: %2d \n", rs, r, s); 
            	 if (r == s) 
            		{
            		vmoBiy[q][p] += 2.0 * buf_X122.matrix[h][pq][rs];
            		//psi::outfile->Printf("To V^B_iy[%d][%d] added: (%2d %2d | %2d %2d) = %16.10f\n", q, p, p, q, r, s, buf_X122.matrix[h][pq][rs]);
            		}
                     }
                }
           global_dpd_->buf4_mat_irrep_close(&buf_X122, h);
           }
      global_dpd_->buf4_close(&buf_X122);

      // vmoAjm (elec.) = 2(jm|kk) = 2(2X|11) = 2(X2|11) = 2(mj|kk)
      for (int h = 0; h < wfn_union_->nirrep(); ++h) 
           {
           global_dpd_->buf4_mat_irrep_init(&buf_X211, h);
           global_dpd_->buf4_mat_irrep_rd(&buf_X211, h);
           for (int pq = 0; pq < buf_X211.params->rowtot[h]; ++pq) 
                {
                int p = buf_X211.params->roworb[h][pq][0]; // nvir_1 X
                int q = buf_X211.params->roworb[h][pq][1]; // nocc_2 2
                for (int rs = 0; rs < buf_X211.params->coltot[h]; ++rs) 
                     {
            	 int r = buf_X211.params->colorb[h][rs][0]; // nocc_1 1
            	 int s = buf_X211.params->colorb[h][rs][1]; // nocc_1 1
            	 if (r == s) 
            	 	{
            	 	vmoAjm[q][p] += 2.0 * buf_X211.matrix[h][pq][rs];
            	 	//psi::outfile->Printf("To V^A_jm[%d][%d] added: (%2d %2d | %2d %2d) = %16.10f\n", q, p, p, q, r, s, buf_X211.matrix[h][pq][rs]);
            	 	}
                	 }
                }
           global_dpd_->buf4_mat_irrep_close(&buf_X211, h);
           }
      global_dpd_->buf4_close(&buf_X211);

      // vmoAjy (elec.) = 2(jn|kk) =  2(2Y|11) = 2(Y2|11) = 2(nj|kk)
      for (int h = 0; h < wfn_union_->nirrep(); ++h) 
           {
           global_dpd_->buf4_mat_irrep_init(&buf_Y211, h);
           global_dpd_->buf4_mat_irrep_rd(&buf_Y211, h);
           for (int pq = 0; pq < buf_Y211.params->rowtot[h]; ++pq) 
                {
                int p = buf_Y211.params->roworb[h][pq][0]; // nvir_2 Y
                int q = buf_Y211.params->roworb[h][pq][1]; // nocc_2 2
                for (int rs = 0; rs < buf_Y211.params->coltot[h]; ++rs) 
                     {
            	 int r = buf_Y211.params->colorb[h][rs][0]; // nocc_1
            	 int s = buf_Y211.params->colorb[h][rs][1]; // nocc_1
            	 if (r == s)
            		{ 
            		vmoAjy[q][p] += 2.0 * buf_Y211.matrix[h][pq][rs]; 
            		//psi::outfile->Printf("To V^A_jy[%d][%d] added: 2x(%2d %2d | %2d %2d) = %16.10f\n", q, p, p, q, r, s, buf_Y211.matrix[h][pq][rs]);			   
            		}
                     }
                }
           global_dpd_->buf4_mat_irrep_close(&buf_Y211, h);
           }
      global_dpd_->buf4_close(&buf_Y211);

      // Close DPD library
      psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  } /* End of 2-electron part */


 
  // ==> Vin_I and Wjm_I <== // 

  // Vin_I = VBin - VBix*Sxn - VBiy*Syn // 
  std::shared_ptr<psi::Matrix> Vin_I = VmoBin->clone();
  Vin_I->gemm(false, false, -1.0, VmoBix, Smoxn, 1.0);
  Vin_I->gemm(false, false, -1.0, VmoBiy, Smoyn, 1.0);

  // Wjm_I = VAjm - Vjx*Sxm - Vjy*Sym 
  std::shared_ptr<psi::Matrix> Wjm_I = VmoAjm->clone();
  Wjm_I->gemm(false, false, -1.0, VmoAjx, Smoxm, 1.0);
  Wjm_I->gemm(false, false, -1.0, VmoAjy, Smoym, 1.0);

  // ==> Vin_II and Wjm_II <== //

  // Vin_II = Vin_I + Sij*(Tnj - Sxn*Txj - Syn*Tyj) //
  Tmonj->gemm(true, false, -1.0, Smoxn, Tmoxj, 1.0);
  Tmonj->gemm(true, false, -1.0, Smoyn, Tmoyj, 1.0);
  std::shared_ptr<psi::Matrix> Vin_II = psi::Matrix::doublet(Smoij, Tmonj, false, true);
  Vin_II->add(Vin_I);

  // Wjm_II = Wjm_I + Sij*(Tmi - Sxm*Txi - Sym*Tyi)
  Tmomi->gemm(true, false, -1.0, Smoxm, Tmoxi, 1.0);
  Tmomi->gemm(true, false, -1.0, Smoym, Tmoyi, 1.0);
  std::shared_ptr<psi::Matrix> Wjm_II = psi::Matrix::doublet(Smoij, Tmomi, true, true);
  Wjm_II->add(Wjm_I);


  // ==> Final calculation of CT energy <== //

  double e_ct_AB = 0.0; // A ---> B
  for (int i=0; i<ndocc_1; ++i)
       {
       for (int n=0; n<nvir_2; ++n)
            {
            // Vin_I/(F_ii - T_nn)
            double value = Vin_I->get(i, n) / (Fmoii->get(i, i) - Tmonn->get(n, n));

            // (1 - S2_mn) -> normalizing factor
            double s2 = 0.0;
            for (int x=0; x<ndocc_1; ++x)
                 {
                 double s = Smoxn->get(x, n);
                 s2 += s*s;
                 }
            for (int y=0; y<nvir_1; ++y)
                 {
                 double s = Smoyn->get(y, n);
                 s2 += s*s;
                 }
            value /= (1.0 - s2);

            // Vin_I * Vin_II/(F_ii-T_nn)/(1-S2mn)
            value *= Vin_II->get(i, n);

            e_ct_AB += value;
            }
        }

  double e_ct_BA = 0.0; // B ---> A
  for (int j=0; j<ndocc_2; ++j)
       {
       for (int m=0; m<nvir_1; ++m)
            {
            // Wjm_I/(F_jj - T_mm)
            double value = Wjm_I->get(j, m) / (Fmojj->get(j, j) - Tmomm->get(m, m));

            // 1 - S2mn -> normalizing factor
            double s2 = 0.0;
            for (int x=0; x<ndocc_2; ++x)
                 {
                 double s = Smoxm->get(x, m);
                 s2 += s*s;
                 }
            for (int y=0; y<nvir_2; ++y)
                 {
                 double s = Smoym->get(y, m);
                 s2 += s*s;
                 }
            value /= (1.0 - s2);

            // Wjm_I * Wjm_II/(F_jj-T_mm)/(1-S2mn)
            value *= Wjm_II->get(j, m);

            e_ct_BA += value;
            }
        }

  // --> Compute total CT energy <-- //
  e_ct_AB *= 2.0;
  e_ct_BA *= 2.0;
  e_ct = e_ct_AB + e_ct_BA;

  psi::Process::environment.globals["EINT CT EFP2 KCAL"] = e_ct * OEPDEV_AU_KcalPerMole;

  // ---> Timer-off <--- //
  psi::timer_off("Solver E(CT) EFP2 MO-Expanded");
  t_time += clock(); // Clock END
  cout << " o TIME EFP2: " << ((double)t_time/CLOCKS_PER_SEC) << endl;
    
  // ---> Print results <--- //
  if (wfn_union_->options().get_int("PRINT") > -1) 
     {
     psi::outfile->Printf("  ==> SOLVER: Charge transfer energy calculations <==\n");
     psi::outfile->Printf("  ==> EFP2 Model (MO-based) <==\n\n");
     psi::outfile->Printf("     E_CT_A->B  = %13.6f au\n", e_ct_AB);
     psi::outfile->Printf("     E_CT_B->A  = %13.6f au\n", e_ct_BA);
     psi::outfile->Printf("     E_CT       = %13.6f au\n", e_ct   );
     psi::outfile->Printf("\n");
     }
  return e_ct;
  }

double ChargeTransferEnergySolver::compute_oep_based_murrell_etal()
{
  // ===> Set-up OEP Objects <=== //
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
  oep_1->compute("Otto-Ladik.V1.GDF");
  oep_2->compute("Otto-Ladik.V1.GDF");
  oep_1->localize();
  oep_2->localize();
  oep_1->compute("Otto-Ladik.V3.CAMM-nj");
  oep_2->compute("Otto-Ladik.V3.CAMM-nj");


  // ===> Molecules <=== //
  std::shared_ptr<psi::Molecule> mol_1 = wfn_union_->l_wfn(0)->molecule();
  std::shared_ptr<psi::Molecule> mol_2 = wfn_union_->l_wfn(1)->molecule();

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
  ovlInt_1p2p->compute(Sao_1p2p);

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

  clock_t t_time = -clock(); // Clock BEGIN
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



  // ---> Localized occupied orbitals <--- //
  std::shared_ptr<psi::Matrix> La_occ_1 = oep_1->localizer()->L();
  std::shared_ptr<psi::Matrix> La_occ_2 = oep_2->localizer()->L();
  std::shared_ptr<psi::Matrix> S12= psi::Matrix::triplet(La_occ_1, Sao_1p2p, La_occ_2, true, false, false); // LOCC(A) x LOCC(B)
  std::shared_ptr<psi::Matrix> S1Y= psi::Matrix::triplet(La_occ_1, Sao_1p2p, Ca_vir_2, true, false, false); // LOCC(A) x  VIR(B)
  std::shared_ptr<psi::Matrix> S2X= psi::Matrix::triplet(La_occ_2, Sao_1p2p, Ca_vir_1, true, true , false); // LOCC(B) x  VIR(A)
  
  // ---> Get LMO centroids <--- //
  std::vector<std::shared_ptr<psi::Vector>> rmo_1 = oep_1->lmoc();
  std::vector<std::shared_ptr<psi::Vector>> rmo_2 = oep_2->lmoc();

  // ---> Compute auxiliary tensors u and w <--- //
  std::shared_ptr<psi::Vector> u_1 = this->compute_u_vector(rmo_1, rmo_2, mol_2);
  std::shared_ptr<psi::Vector> u_2 = this->compute_u_vector(rmo_2, rmo_1, mol_1);
  std::shared_ptr<psi::Matrix> w_1 = this->compute_w_matrix(mol_1, mol_2, rmo_1);
  std::shared_ptr<psi::Matrix> w_2 = this->compute_w_matrix(mol_2, mol_1, rmo_2);

  // ---> Get distributed effective charges <--- //
  std::vector<std::shared_ptr<psi::Matrix>> q_1 = oep_1->oep("Otto-Ladik.V3.CAMM-nj").dmtp->charges();
  std::vector<std::shared_ptr<psi::Matrix>> q_2 = oep_2->oep("Otto-Ladik.V3.CAMM-nj").dmtp->charges();

  // ===> Compute V1 term <=== //
  std::shared_ptr<psi::Matrix> v_ab_v1 = psi::Matrix::doublet(S1, oep_2->matrix("Otto-Ladik.V1.GDF"), false, false);
  std::shared_ptr<psi::Matrix> v_ba_v1 = psi::Matrix::doublet(S2, oep_1->matrix("Otto-Ladik.V1.GDF"), false, false);
  //const double sc = 0.25;
  //v_ab_v1->scale(sc);
  //v_ba_v1->scale(sc);

  // ===> Compute V2 term <=== //
  std::shared_ptr<psi::Matrix> v_ab_v2 = std::make_shared<psi::Matrix>("", nocc_1, nvir_2);
  std::shared_ptr<psi::Matrix> v_ba_v2 = std::make_shared<psi::Matrix>("", nocc_2, nvir_1);
  for (int i=0; i<nocc_1; ++i) {
       for (int n=0; n<nvir_2; ++n) {
            double v = S1Y->get(i,n) * u_1->get(i);
            v_ab_v2->set(i, n, v); 
       }
  }
  for (int i=0; i<nocc_2; ++i) {
       for (int n=0; n<nvir_1; ++n) {
            double v = S2X->get(i,n) * u_2->get(i);
            v_ba_v2->set(i, n, v); 
       }
  }
  v_ab_v2->gemm(false, false, 1.0, oep_1->localizer()->U(), v_ab_v2->clone(), 0.0);
  v_ba_v2->gemm(false, false, 1.0, oep_2->localizer()->U(), v_ba_v2->clone(), 0.0);

  // ===> Compute V3 term <=== //
  std::shared_ptr<psi::Matrix> v_ab_v3 = std::make_shared<psi::Matrix>("", nocc_1, nvir_2);
  std::shared_ptr<psi::Matrix> v_ba_v3 = std::make_shared<psi::Matrix>("", nocc_2, nvir_1);
  for (int i=0; i<nocc_1; ++i) {
       for (int n=0; n<nvir_2; ++n) {
            double v = 0;
            for (int j=0; j<nocc_2; ++j) {              
                 for (int y=0; y<mol_2->natom(); ++y) {
                      v += S12->get(i,j) * w_1->get(i, y) * q_2[nvir_2*j+n]->get(y, 0);
                 }
            }
            v_ab_v3->set(i, n, v);
       }
  }
  //
  for (int i=0; i<nocc_2; ++i) {
       for (int n=0; n<nvir_1; ++n) {
            double v = 0;
            for (int j=0; j<nocc_1; ++j) {              
                 for (int y=0; y<mol_1->natom(); ++y) {
                      v += S12->get(j,i) * w_2->get(i, y) * q_1[nvir_1*j+n]->get(y, 0);
                 }
            }
            v_ba_v3->set(i, n, v);
       }
  }
  v_ab_v3->gemm(false, false, 1.0, oep_1->localizer()->U(), v_ab_v3->clone(), 0.0);
  v_ba_v3->gemm(false, false, 1.0, oep_2->localizer()->U(), v_ba_v3->clone(), 0.0);


  // ---> Add coupling constant contributions <--- //
  std::shared_ptr<psi::Matrix> v_ab_v12 = v_ab_v1->clone(); v_ab_v12->add(v_ab_v2);
  t_time += clock(); // Clock END
  std::shared_ptr<psi::Matrix> v_ba_v12 = v_ba_v1->clone(); v_ba_v12->add(v_ba_v2);
  std::shared_ptr<psi::Matrix> v_ab_v13 = v_ab_v1->clone(); v_ab_v13->add(v_ab_v3);
  std::shared_ptr<psi::Matrix> v_ba_v13 = v_ba_v1->clone(); v_ba_v13->add(v_ba_v3);
  std::shared_ptr<psi::Matrix> v_ab_v23 = v_ab_v2->clone(); v_ab_v23->add(v_ab_v3);
  std::shared_ptr<psi::Matrix> v_ba_v23 = v_ba_v2->clone(); v_ba_v23->add(v_ba_v3);
  std::shared_ptr<psi::Matrix> v_ab_v123= v_ab_v12->clone(); v_ab_v123->add(v_ab_v3);
  std::shared_ptr<psi::Matrix> v_ba_v123= v_ba_v12->clone(); v_ba_v123->add(v_ba_v3);

  // ===> Compute CT Energy <=== //
  t_time -= clock(); // Clock BEGIN
  double e_ab_v1  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v1);
  t_time += clock(); // Clock END
  double e_ba_v1  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v1);
  double e_ab_v2  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v2);
  double e_ba_v2  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v2);
  double e_ab_v3  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v3);
  double e_ba_v3  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v3);
  double e_ab_v12 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v12);
  double e_ba_v12 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v12);
  double e_ab_v13 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v13);
  double e_ba_v13 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v13);
  double e_ab_v23 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v23);
  double e_ba_v23 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v23);
  double e_ab_v123= this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v123);
  double e_ba_v123= this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v123);

  double e_ab = e_ab_v123;
  double e_ba = e_ba_v123;

  double e_tot = e_ab + e_ba;
  psi::Process::environment.globals["EINT CT OEP-OTTO-LADIK KCAL"] = e_tot*OEPDEV_AU_KcalPerMole;

  cout << " o TIME OEP : " << ((double)t_time/CLOCKS_PER_SEC) << endl;
  //cout << " Time OEP : " << t_end - t_begin - (t2 - t1) << endl;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > -1) {
     psi::outfile->Printf("  ==> SOLVER: Charge-Transfer Energy Calculations    <==\n"  );
     psi::outfile->Printf("  ==>     OEP-Based (Otto-Ladik                 )    <==\n\n");
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     Group I\n"                                              );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v1                        );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v1                        );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     Group II\n"                                             );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v2                        );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v2                        );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     Group III\n"                                            );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v3                        );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v3                        );
     psi::outfile->Printf("     ===============================\n"                      );
     psi::outfile->Printf("     Group I+II\n"                                           );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v12                       );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v12                       );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     Group I+III\n"                                          );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v13                       );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v13                       );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     Group II+III\n"                                         );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab_v23                       );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba_v23                       );
     psi::outfile->Printf("     ===============================\n"                      );
     psi::outfile->Printf("     Total\n"                                                );
     psi::outfile->Printf("     E (A-->B)   = %13.6f\n", e_ab                           );
     psi::outfile->Printf("     E (B-->A)   = %13.6f\n", e_ba                           );
   //psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E_TOT       = %13.6f\n", e_tot                          );
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
std::shared_ptr<psi::Vector> ChargeTransferEnergySolver::compute_u_vector(
  std::vector<std::shared_ptr<psi::Vector>> rmo_1, 
  std::vector<std::shared_ptr<psi::Vector>> rmo_2, 
  std::shared_ptr<psi::Molecule> mol_2
)
{
  const int nocc_1 = rmo_1[0]->dim();
  const int nocc_2 = rmo_2[0]->dim();
  const int nat_2  = mol_2->natom();

  std::shared_ptr<psi::Vector> u = std::make_shared<psi::Vector>("", nocc_1);

  for (int i=0; i<nocc_1; ++i) {
       double v = 0.0;
       //
       for (int y=0; y<nat_2; ++y) {
            double ryi = sqrt(pow(mol_2->x(y) - rmo_1[0]->get(i), 2.0) +
                              pow(mol_2->y(y) - rmo_1[1]->get(i), 2.0) +
                              pow(mol_2->z(y) - rmo_1[2]->get(i), 2.0) );
            v += (double)mol_2->Z(y) / ryi;
       }
       //
       for (int j=0; j<nocc_2; ++j) {
            double rji = sqrt(pow(rmo_2[0]->get(j) - rmo_1[0]->get(i), 2.0) +
                              pow(rmo_2[1]->get(j) - rmo_1[1]->get(i), 2.0) +
                              pow(rmo_2[2]->get(j) - rmo_1[2]->get(i), 2.0) );
            v -= 2.0 / rji;
       }
       u->set(i, v);
  }
  return u;
}
std::shared_ptr<psi::Matrix> ChargeTransferEnergySolver::compute_w_matrix(
  std::shared_ptr<psi::Molecule> mol_1,
  std::shared_ptr<psi::Molecule> mol_2,
  std::vector<std::shared_ptr<psi::Vector>> rmo_1
)
{
  const int nat_1 = mol_1->natom();
  const int nat_2 = mol_2->natom();
  const int nocc_1= rmo_1[0]->dim();

  std::shared_ptr<psi::Matrix> w = std::make_shared<psi::Matrix>("", nocc_1, nat_2);

  for (int y=0; y<nat_2; ++y) {
       double vy = 0.0;
       //
       for (int x=0; x<nat_1; ++x) { 
            double rxy = sqrt(pow(mol_1->x(x) - mol_2->x(y), 2.0) +
                              pow(mol_1->y(x) - mol_2->y(y), 2.0) +
                              pow(mol_1->z(x) - mol_2->z(y), 2.0) );
            vy += (double)mol_1->Z(x) / rxy;
       }
       //
       for (int k=0; k<nocc_1; ++k) {
            double rky = sqrt(pow(rmo_1[0]->get(k) - mol_2->x(y), 2.0) +
                              pow(rmo_1[1]->get(k) - mol_2->y(y), 2.0) +
                              pow(rmo_1[2]->get(k) - mol_2->z(y), 2.0) );
            vy-= 2.0 / rky;
       }
       //
       for (int i=0; i<nocc_1; ++i) {
            //
            double riy = sqrt(pow(rmo_1[0]->get(i) - mol_2->x(y), 2.0) +
                              pow(rmo_1[1]->get(i) - mol_2->y(y), 2.0) +
                              pow(rmo_1[2]->get(i) - mol_2->z(y), 2.0) );
            w->set(i, y, vy + 2.0/riy );
       }
  }
  w->scale(-1.0);
  return w;
}
//
std::shared_ptr<psi::Vector> ChargeTransferEnergySolver::extract_xyz(std::shared_ptr<psi::Molecule> mol)
{
  auto xyz = std::make_shared<psi::Vector>(3 * mol->natom());
  double* xyz_p = xyz->pointer();
  for (int i = 0; i < mol->natom(); ++i) {
       *xyz_p++ = mol->x(i);
       *xyz_p++ = mol->y(i);
       *xyz_p++ = mol->z(i); 
  }
  return xyz;
}
//
std::shared_ptr<psi::Vector> ChargeTransferEnergySolver::extract_dmtp(std::shared_ptr<oepdev::DMTPole> camm)
{
  auto mult= std::make_shared<psi::Vector>((1 + 3 + 6 + 10) * camm->n_sites());
  psi::SharedMatrix m_0 = camm->charges(0);
  psi::SharedMatrix m_1 = camm->dipoles(0);
  psi::SharedMatrix m_2 = camm->quadrupoles(0);
  psi::SharedMatrix m_3 = camm->octupoles(0);
  //m_1->zero();
  //m_2->zero(); 
  if (options_.get_bool("EFP2_CT_NO_OCTUPOLES")) m_3->zero();

  double* mult_p = mult->pointer();
  double** p_0 = m_0->pointer();
  double** p_1 = m_1->pointer();
  double** p_2 = m_2->pointer();
  double** p_3 = m_3->pointer();

  for (int i = 0; i < camm->n_sites(); ++i) {
       *mult_p++ = p_0[i][0];   // 0

       *mult_p++ = p_1[i][0];   // X
       *mult_p++ = p_1[i][1];   // Y
       *mult_p++ = p_1[i][2];   // Z

        double t = 0.5 * (p_2[i][0] + p_2[i][3] + p_2[i][5]);
       *mult_p++ = p_2[i][0] * 1.5 - t;   // XX
       *mult_p++ = p_2[i][3] * 1.5 - t;   // YY
       *mult_p++ = p_2[i][5] * 1.5 - t;   // ZZ
       *mult_p++ = p_2[i][1] * 1.5;   // XY
       *mult_p++ = p_2[i][2] * 1.5;   // XZ
       *mult_p++ = p_2[i][4] * 1.5;   // YZ

        double tx = 0.5 * (p_3[i][0] + p_3[i][3] + p_3[i][5]);
        double ty = 0.5 * (p_3[i][6] + p_3[i][1] + p_3[i][8]);
        double tz = 0.5 * (p_3[i][9] + p_3[i][2] + p_3[i][7]);
       *mult_p++ = p_3[i][0] * 2.5 - 3.0 * tx;   // XXX
       *mult_p++ = p_3[i][6] * 2.5 - 3.0 * ty;   // YYY
       *mult_p++ = p_3[i][9] * 2.5 - 3.0 * tz;   // ZZZ
       *mult_p++ = p_3[i][1] * 2.5 -       ty;   // XXY
       *mult_p++ = p_3[i][2] * 2.5 -       tz;   // XXZ
       *mult_p++ = p_3[i][3] * 2.5 -       tx;   // XYY
       *mult_p++ = p_3[i][7] * 2.5 -       tz;   // YYZ
       *mult_p++ = p_3[i][5] * 2.5 -       tx;   // XZZ
       *mult_p++ = p_3[i][8] * 2.5 -       ty;   // YZZ
       *mult_p++ = p_3[i][4] * 2.5           ;   // XYZ
  }
  return mult;
}
