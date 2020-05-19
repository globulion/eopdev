#include "cis.h"
#include <iostream>

namespace oepdev{

U_CISComputer_Explicit::U_CISComputer_Explicit(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 U_CISComputer(wfn, opt)
{
}

U_CISComputer_Explicit::~U_CISComputer_Explicit() {}

void U_CISComputer_Explicit::set_beta_(void) {
// Fb_oo_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","OCC"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","OCC"), true, false, false);
// Fb_vv_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","VIR"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","VIR"), true, false, false);
 eps_b_o_ = ref_wfn_->epsilon_b_subset("MO", "ACTIVE_OCC");
 eps_b_v_ = ref_wfn_->epsilon_b_subset("MO", "ACTIVE_VIR");
}

void U_CISComputer_Explicit::build_hamiltonian_(void) {

 std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
 psi::dpd_set_default(inttrans_->get_dpd_id());
 psi::dpdbuf4 buf_OOVV, buf_OVOV, buf_OVov, buf_oovv, buf_ovov;
 psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
 //psio->tocprint(PSIF_LIBTRANS_DPD);

 psi::global_dpd_->buf4_init(&buf_OOVV, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[O,O]"  ), inttrans_->DPD_ID("[V,V]"  ),
                        inttrans_->DPD_ID("[O>=O]+"), inttrans_->DPD_ID("[V>=V]+"), 0, "MO Ints (OO|VV)");

 psi::global_dpd_->buf4_init(&buf_OVOV, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[O,V]"  ),
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[O,V]"  ), 0, "MO Ints (OV|OV)");

 psi::global_dpd_->buf4_init(&buf_OVov, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[o,v]"  ),
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[o,v]"  ), 0, "MO Ints (OV|ov)");

 psi::global_dpd_->buf4_init(&buf_oovv, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[o,o]"  ), inttrans_->DPD_ID("[v,v]"  ),
                        inttrans_->DPD_ID("[o>=o]+"), inttrans_->DPD_ID("[v>=v]+"), 0, "MO Ints (oo|vv)");

 psi::global_dpd_->buf4_init(&buf_ovov, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[o,v]"  ), inttrans_->DPD_ID("[o,v]"  ),
                        inttrans_->DPD_ID("[o,v]"  ), inttrans_->DPD_ID("[o,v]"  ), 0, "MO Ints (ov|ov)");


 double** H = this->H_->pointer();
 //double** Fa_oo = this->Fa_oo_->pointer();
 //double** Fa_vv = this->Fa_vv_->pointer();
 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
 double* eps_b_o = this->eps_b_o_->pointer();
 double* eps_b_v = this->eps_b_v_->pointer();
 const int off = this->naocc_ * this->navir_;


 // (OV|OV) integral contributions and Fock matrix contributions
 for (int h = 0; h < ref_wfn_->nirrep(); ++h) {
      psi::global_dpd_->buf4_mat_irrep_init(&buf_OVOV, h);
      psi::global_dpd_->buf4_mat_irrep_rd(&buf_OVOV, h);
      for (int ia = 0; ia < buf_OVOV.params->rowtot[h]; ++ia) {
           int i = buf_OVOV.params->roworb[h][ia][0];
           int a = buf_OVOV.params->roworb[h][ia][1];
           for (int jb = 0; jb < buf_OVOV.params->coltot[h]; ++jb) {
                int j = buf_OVOV.params->colorb[h][jb][0];
                int b = buf_OVOV.params->colorb[h][jb][1];
                int ia_= this->navir_ * i + a;
                int jb_= this->navir_ * j + b;
                double ia_jb = buf_OVOV.matrix[h][ia][jb];
                double v = ia_jb;
                if ((i==j) && (a==b)) v+= eps_a_v[a] - eps_a_o[i];
                H[ia_][jb_    ] += v    ; // block AA 
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_OVOV, h);
 }

 // (ov|ov) integral contributions and Fock matrix contributions
 for (int h = 0; h < ref_wfn_->nirrep(); ++h) {
      psi::global_dpd_->buf4_mat_irrep_init(&buf_ovov, h);
      psi::global_dpd_->buf4_mat_irrep_rd(&buf_ovov, h);
      for (int ia = 0; ia < buf_ovov.params->rowtot[h]; ++ia) {
           int i = buf_ovov.params->roworb[h][ia][0];
           int a = buf_ovov.params->roworb[h][ia][1];
           for (int jb = 0; jb < buf_ovov.params->coltot[h]; ++jb) {
                int j = buf_ovov.params->colorb[h][jb][0];
                int b = buf_ovov.params->colorb[h][jb][1];
                int ia_= off + this->nbvir_ * i + a;
                int jb_= off + this->nbvir_ * j + b;
                double ia_jb = buf_ovov.matrix[h][ia][jb];
                double v = ia_jb;
                if ((i==j) && (a==b)) v+= eps_b_v[a] - eps_b_o[i];
                H[ia_][jb_    ] += v    ; // block BB
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_ovov, h);
 }

 // (OO|VV) integral contributions
 for (int h = 0; h < ref_wfn_->nirrep(); ++h) {
      psi::global_dpd_->buf4_mat_irrep_init(&buf_OOVV, h);
      psi::global_dpd_->buf4_mat_irrep_rd(&buf_OOVV, h);
      for (int ij = 0; ij < buf_OOVV.params->rowtot[h]; ++ij) {
           int i = buf_OOVV.params->roworb[h][ij][0];
           int j = buf_OOVV.params->roworb[h][ij][1];
           for (int ab = 0; ab < buf_OOVV.params->coltot[h]; ++ab) {
                int a = buf_OOVV.params->colorb[h][ab][0];
                int b = buf_OOVV.params->colorb[h][ab][1];
                int ia = this->navir_ * i + a;
                int jb = this->navir_ * j + b;
                H[ia][jb] -= buf_OOVV.matrix[h][ij][ab]; // block AA
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_OOVV, h);
 }

 // (oo|vv) integral contributions
 for (int h = 0; h < ref_wfn_->nirrep(); ++h) {
      psi::global_dpd_->buf4_mat_irrep_init(&buf_oovv, h);
      psi::global_dpd_->buf4_mat_irrep_rd(&buf_oovv, h);
      for (int ij = 0; ij < buf_oovv.params->rowtot[h]; ++ij) {
           int i = buf_oovv.params->roworb[h][ij][0];
           int j = buf_oovv.params->roworb[h][ij][1];
           for (int ab = 0; ab < buf_oovv.params->coltot[h]; ++ab) {
                int a = buf_oovv.params->colorb[h][ab][0];
                int b = buf_oovv.params->colorb[h][ab][1];
                int ia = off + this->nbvir_ * i + a;
                int jb = off + this->nbvir_ * j + b;
                H[ia][jb] -= buf_oovv.matrix[h][ij][ab]; // block BB
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_oovv, h);
 }

 // (OV|ov) integral contributions
 for (int h = 0; h < ref_wfn_->nirrep(); ++h) {
      psi::global_dpd_->buf4_mat_irrep_init(&buf_OVov, h);
      psi::global_dpd_->buf4_mat_irrep_rd(&buf_OVov, h);
      for (int ia = 0; ia < buf_OVov.params->rowtot[h]; ++ia) {
           int i = buf_OVov.params->roworb[h][ia][0];
           int a = buf_OVov.params->roworb[h][ia][1];
           for (int jb = 0; jb < buf_OVov.params->coltot[h]; ++jb) {
                int j = buf_OVov.params->colorb[h][jb][0];
                int b = buf_OVov.params->colorb[h][jb][1];
                int ia_ = this->navir_ * i + a;
                int jb_ = this->nbvir_ * j + b;
                double ia_jb = buf_OVov.matrix[h][ia][jb];
                H[ia_][jb_+off] += ia_jb; // block AB
                H[jb_+off][ia_] += ia_jb; // block BA
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_OVov, h);
 }

 //
 psi::global_dpd_->buf4_close(&buf_OOVV);
 psi::global_dpd_->buf4_close(&buf_OVOV);
 psi::global_dpd_->buf4_close(&buf_oovv);
 psi::global_dpd_->buf4_close(&buf_ovov);
 psi::global_dpd_->buf4_close(&buf_OVov);
 psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
}


} // EndNameSpace oepdev
