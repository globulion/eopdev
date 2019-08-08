#include "cis.h"
#include <iostream>

namespace oepdev{

R_CISComputer::R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt, psi::IntegralTransform::TransformationType::Restricted)
{
}

R_CISComputer::~R_CISComputer() {}

void R_CISComputer::set_beta_(void) {
 //Fb_oo_ = Fa_oo_;
 //Fb_vv_ = Fa_vv_;
 eps_b_o_ = eps_a_o_;
 eps_b_v_ = eps_a_v_;
 // They are not used anyway
}

void R_CISComputer::build_hamiltonian_(void) {

 std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
 psi::dpd_set_default(inttrans_->get_dpd_id());
 psi::dpdbuf4 buf_OOVV, buf_OVOV;
 psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
 psio->tocprint(PSIF_LIBTRANS_DPD);

 psi::global_dpd_->buf4_init(&buf_OOVV, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[O,O]"  ), inttrans_->DPD_ID("[V,V]"  ),
                        inttrans_->DPD_ID("[O>=O]+"), inttrans_->DPD_ID("[V>=V]+"), 0, "MO Ints (OO|VV)");

 psi::global_dpd_->buf4_init(&buf_OVOV, PSIF_LIBTRANS_DPD, 0,
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[O,V]"  ),
                        inttrans_->DPD_ID("[O,V]"  ), inttrans_->DPD_ID("[O,V]"  ), 0, "MO Ints (OV|OV)");


 double** H = this->H_->pointer();
 //double** Fa_oo = this->Fa_oo_->pointer();
 //double** Fa_vv = this->Fa_vv_->pointer();
 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
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
                //if (i==j) v+= Fa_vv[a][b];
                //if (a==b) v-= Fa_oo[i][j];
                H[ia_][jb_    ] += v    ; // block AA 
                H[ia_][jb_+off] += ia_jb; // block AB 
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_OVOV, h);
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
                                                         // block AB -> no contribution
           }
      }
      psi::global_dpd_->buf4_mat_irrep_close(&buf_OOVV, h);
 }

 //
 psi::global_dpd_->buf4_close(&buf_OOVV);
 psi::global_dpd_->buf4_close(&buf_OVOV);
 psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

 // Copy the AA and AB blocks into BB and BA blocks, respectively
 for (int i=0; i<nbocc_; ++i) {
 for (int a=0; a<nbvir_; ++a) {
      int ia = nbvir_*i + a;

      for (int j=0; j<nbocc_; ++j) {
      for (int b=0; b<nbvir_; ++b) {
           int jb = nbvir_*j + b;
           // block BB --> AA
           H[ia+off][jb+off] = H[ia][jb];
           // block BA --> AB.T
           H[ia+off][jb] = H[ia][jb+off];
      }
      }
 }
 }

}



} // EndNameSpace oepdev
