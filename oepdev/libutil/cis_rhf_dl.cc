#include "cis.h"
#include "integrals_iter.h"
#include <iostream>

namespace oepdev{

R_CISComputer_DL::R_CISComputer_DL(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 R_CISComputer(wfn, opt),
 Ca_occ__(wfn->Ca_subset("AO","ACTIVE_OCC")),
 Ca_vir__(wfn->Ca_subset("AO","ACTIVE_VIR"))
{
}

void R_CISComputer_DL::set_nstates_() {
  // Set the number of states as the number of roots in Davidson-Liu method
  this->nstates_ = options_.get_int("DAVIDSON_LIU_NROOTS");
}

R_CISComputer_DL::~R_CISComputer_DL() {}

void R_CISComputer_DL::transform_integrals_(void) 
{
 // Do nothing here since it is direct CIS
}

void R_CISComputer_DL::allocate_hamiltonian_(void) 
{
  // Set the number of starting guess vectors
  int L_start = this->nstates_;
  const int i = options_.get_int("DAVIDSON_LIU_SPACE_START");
  if (i>0) {
     if (i < L_start) throw psi::PSIEXCEPTION("The initial space size in Davidson-Liu calculation has to be at least equal to the number of roots!");
     L_start = i;
  }

  // Initialize Davidson-Liu solver
  this->davidson_liu_initialize(this->ndets_, L_start, this->nstates_);

  // Assign the pointers to the Davidson-Liu eigenpairs (do double memory allocation)
  this->E_ = this->E_davidson_liu_; this->E_->set_name("CIS Eigenvalues");
  this->U_ = this->U_davidson_liu_; this->U_->set_name("CIS Eigenvectors");
}

void R_CISComputer_DL::build_hamiltonian_(void) 
{
 // Do nothing here since it is Davidson-Liu CIS
}

void R_CISComputer_DL::diagonalize_hamiltonian_(void) 
{
 this->run_davidson_liu();
}


void R_CISComputer_DL::davidson_liu_compute_diagonal_hamiltonian(void) {

 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
 double** cao = Ca_occ__->pointer();
 double* h = this->H_diag_davidson_liu_->pointer();
 const int off = this->naocc_ * this->navir_;
 const int nbf = this->ref_wfn_->basisset()->nbf();

 // Fock matrix and ERI contribution directly in MO basis
 std::vector<psi::SharedMatrix>& C_left = jk_->C_left();
 std::vector<psi::SharedMatrix>& C_right= jk_->C_right();
 const std::vector<psi::SharedMatrix>& J = jk_->J();
 const std::vector<psi::SharedMatrix>& K = jk_->K();
 C_left.clear(); C_right.clear();

 psi::SharedMatrix identity = std::make_shared<psi::Matrix>("", nbf, nbf);
 identity->identity();

 for (int i=0; i<this->naocc_; ++i) {
      psi::SharedMatrix Wia = std::make_shared<psi::Matrix>("", nbf, nbf); 
      psi::SharedVector Cai = Ca_occ__->get_column(0, i);
      double*  cai = Cai->pointer();
      double** wai = Wia->pointer();
      for (int ii=0; ii<nbf; ++ii) {
           double cai_ii = cai[ii];
           wai[ii][ii] = cai_ii * cai_ii;
           for (int jj=0; jj<ii; ++jj) {
                double v = cai_ii * cai[jj];
                wai[ii][jj] = v; 
                wai[jj][ii] = v;
           }
      }

      psi::SharedMatrix I1 = identity->clone();
      C_left.push_back(Wia);
      C_right.push_back(I1);
 }

 this->jk_->compute();

 // Add to Hamiltonian
 for (int i=0; i<this->naocc_; ++i) {
      psi::SharedMatrix Ji_aa = psi::Matrix::triplet(Ca_vir__, J[i], Ca_vir__, true, false, false);
      psi::SharedMatrix Ki_aa = psi::Matrix::triplet(Ca_vir__, K[i], Ca_vir__, true, false, false);
      for (int a=0; a<this->navir_; ++a) {
           int ia = navir_*i + a;
           double iaia = Ki_aa->get(a,a);
           double iiaa = Ji_aa->get(a,a);
           double v  = eps_a_v[a] - eps_a_o[i] + iaia - iiaa;
           h[ia    ] = v; // A block
           h[ia+off] = v; // B block
      }
 }


}

void R_CISComputer_DL::davidson_liu_compute_sigma(void) {
 // Loop over all uncomputed sigma vectors
 // to_compute = this->L_davidson_liu_ - this->davidson_liu_n_sigma_computed_;

 // Sizing
 const int nbf = this->ref_wfn_->basisset()->nbf();
 const int off = this->naocc_ * this->navir_;

 // Prepare LCAO-MO matrices and Fock matrix in MO basis
 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();

 // Access the JK memory
 std::vector<psi::SharedMatrix>& C_left = jk_->C_left();
 std::vector<psi::SharedMatrix>& C_right= jk_->C_right();
 const std::vector<psi::SharedMatrix>& J = jk_->J();
 const std::vector<psi::SharedMatrix>& K = jk_->K();

 // Clear JK buffers
 C_left.clear(); C_right.clear();

 // Compute generalized density matrices from guess vectors
 std::vector<psi::SharedMatrix> Ta_set, Tb_set;
 for (int k=this->davidson_liu_n_sigma_computed_; k<this->L_davidson_liu_; ++k) {
      psi::SharedMatrix Ba = std::make_shared<psi::Matrix>("", this->naocc_, this->navir_); 
      psi::SharedMatrix Bb = std::make_shared<psi::Matrix>("", this->nbocc_, this->nbvir_);
                                                                                            
      for (int i=0; i<naocc_; ++i) {            
      for (int a=0; a<navir_; ++a) {
           int ia = navir_*i + a;
                                                                                            
           Ba->set(i, a, this->guess_vectors_davidson_liu_->V(k)->get(ia    ));
           Bb->set(i, a, this->guess_vectors_davidson_liu_->V(k)->get(ia+off));
      }}
      psi::SharedMatrix Ta = psi::Matrix::triplet(Ca_occ__, Ba, Ca_vir__, false, false, true);
      psi::SharedMatrix Tb = psi::Matrix::triplet(Ca_occ__, Bb, Ca_vir__, false, false, true);

      Ta_set.push_back(Ta);
      Tb_set.push_back(Tb);
 }

 // Compute generalized J and K matrices for A and B blocks of guess vectors
 psi::SharedMatrix identity = std::make_shared<psi::Matrix>("", nbf, nbf);
 identity->identity();

 for (int k=this->davidson_liu_n_sigma_computed_; k<this->L_davidson_liu_; ++k) {
      psi::SharedMatrix I1 = identity->clone();
      psi::SharedMatrix I2 = identity->clone();

      C_left.push_back(Ta_set[k-this->davidson_liu_n_sigma_computed_]);
      C_right.push_back(I1);
      C_left.push_back(Tb_set[k-this->davidson_liu_n_sigma_computed_]);
      C_right.push_back(I2);
 }
 this->jk_->compute();

 // Transform J and K matrices from AO to MO basis and compute sigma vectors
 for (int k=this->davidson_liu_n_sigma_computed_; k<this->L_davidson_liu_; ++k) {

     const int kk = k-this->davidson_liu_n_sigma_computed_;
     psi::SharedMatrix Ja_OV = psi::Matrix::triplet(Ca_occ__, J[2*kk+0], Ca_vir__, true, false, false);
     psi::SharedMatrix Jb_OV = psi::Matrix::triplet(Ca_occ__, J[2*kk+1], Ca_vir__, true, false, false);
     psi::SharedMatrix Ka_OV = psi::Matrix::triplet(Ca_occ__, K[2*kk+0], Ca_vir__, true, false, false);
     psi::SharedMatrix Kb_OV = psi::Matrix::triplet(Ca_occ__, K[2*kk+1], Ca_vir__, true, false, false);
     Ja_OV->add(Jb_OV);

     psi::SharedMatrix S_OV = Ja_OV->clone();
     psi::SharedMatrix S_ov = Ja_OV->clone();
     S_OV->subtract(Ka_OV);
     S_ov->subtract(Kb_OV);

     psi::SharedVector Sigma = std::make_shared<psi::Vector>("", this->N_davidson_liu_);
     double* sigma = Sigma->pointer();
     double* bvec  = this->guess_vectors_davidson_liu_->V(k)->pointer();

     for (int i=0; i<naocc_; ++i) {            
     for (int a=0; a<navir_; ++a) {
          int ia = navir_*i + a;
                                               
          double v  = eps_a_v[a] - eps_a_o[i];
          sigma[ia    ] = v * bvec[ia    ];
          sigma[ia+off] = v * bvec[ia+off];

          sigma[ia    ]+= S_OV->get(i,a);
          sigma[ia+off]+= S_ov->get(i,a);
     }
     }

     this->sigma_vectors_davidson_liu_.push_back(Sigma);

 }

}


} // EndNameSpace oepdev
