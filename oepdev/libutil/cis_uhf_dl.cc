#include "cis.h"
#include "integrals_iter.h"
#include <iostream>

namespace oepdev{

U_CISComputer_DL::U_CISComputer_DL(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 U_CISComputer(wfn, opt),
 Ca_occ__(wfn->Ca_subset("AO","OCC")),
 Ca_vir__(wfn->Ca_subset("AO","VIR")),
 Cb_occ__(wfn->Cb_subset("AO","OCC")),
 Cb_vir__(wfn->Cb_subset("AO","VIR"))
{
}

void U_CISComputer_DL::set_nstates_() {
  // Set the number of states as the number of roots in Davidson-Liu method
  this->nstates_ = options_.get_int("DAVIDSON_LIU_NROOTS");
}

U_CISComputer_DL::~U_CISComputer_DL() {}

void U_CISComputer_DL::transform_integrals_(void) 
{
 // Do nothing here since it is direct CIS
}

void U_CISComputer_DL::allocate_hamiltonian_(void) 
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

void U_CISComputer_DL::build_hamiltonian_(void) 
{
 // Do nothing here since it is Davidson-Liu CIS
}

void U_CISComputer_DL::diagonalize_hamiltonian_(void) 
{
 this->run_davidson_liu();
}


void U_CISComputer_DL::davidson_liu_compute_diagonal_hamiltonian(void) {

 /* 
    Loop over ERI's explicitly. This is more efficient in terms of memory because
    too many generalized AO density matrices needs to be generated for JK object.
 */

 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
 double* eps_b_o = this->eps_b_o_->pointer();
 double* eps_b_v = this->eps_b_v_->pointer();
 double* h = this->H_diag_davidson_liu_->pointer();
 const int off = this->naocc_ * this->navir_;

 // Fock matrix contribution directly in MO basis
 for (int i=0; i<naocc_; ++i) {
 for (int a=0; a<navir_; ++a) {
      int ia = navir_*i + a;

      double v  = eps_a_v[a] - eps_a_o[i];
      h[ia    ] = v;  // A block
 }
 }
 for (int i=0; i<nbocc_; ++i) {
 for (int a=0; a<nbvir_; ++a) {
      int ia = nbvir_*i + a;

      double v  = eps_b_v[a] - eps_b_o[i];
      h[ia+off] = v;  // B block
 }
 }

 // ERI contribution in AO basis on the fly
 psi::IntegralFactory fact(ref_wfn_->basisset());
 std::shared_ptr<psi::TwoBodyAOInt> tei(fact.eri());                                           
 const double * buffer = tei->buffer();

 std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact, "ALL");

 double** cao = Ca_occ__->pointer();
 double** cav = Ca_vir__->pointer();
 double** cbo = Cb_occ__->pointer();
 double** cbv = Cb_vir__->pointer();

 const double eri_cutoff = options_.get_double("CIS_SCHWARTZ_CUTOFF");

 for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
 {
      shellIter->compute_shell(tei);
      std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
      for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
      {
           int I = intsIter->i();  // \alpha  
           int J = intsIter->j();  // \beta   
           int K = intsIter->k();  // \gamma  
           int L = intsIter->l();  // \delta  

           double eri = buffer[intsIter->index()];

           if (std::abs(eri) > eri_cutoff) {

               // A block
               for (int i=0; i<naocc_; ++i) {                
                    double cIi = cao[I][i];
                    double cJi = cao[J][i];
                    double cKi = cao[K][i];
               for (int a=0; a<navir_; ++a) {
                    int ia = navir_*i + a;

                    double cJa = cav[J][a];
                    double cKa = cav[K][a];
                    double cLa = cav[L][a];

                    double v = eri * cIi * cLa * (cJa * cKi - cJi * cKa);
                    h[ia    ] += v;
               }
               }

               // B block
               for (int i=0; i<nbocc_; ++i) {                
                    double cIi = cbo[I][i];
                    double cJi = cbo[J][i];
                    double cKi = cbo[K][i];
               for (int a=0; a<nbvir_; ++a) {
                    int ia = nbvir_*i + a;

                    double cJa = cbv[J][a];
                    double cKa = cbv[K][a];
                    double cLa = cbv[L][a];

                    double v = eri * cIi * cLa * (cJa * cKi - cJi * cKa);
                    h[ia+off] += v;
               }
               }

           }
      }
 }

}

void U_CISComputer_DL::davidson_liu_compute_sigma(void) {
 // Loop over all uncomputed sigma vectors
 // to_compute = this->L_davidson_liu_ - this->davidson_liu_n_sigma_computed_;

 // Sizing
 const int nbf = this->ref_wfn_->basisset()->nbf();
 const int off = this->naocc_ * this->navir_;

 // Prepare LCAO-MO matrices and Fock matrix in MO basis
 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
 double* eps_b_o = this->eps_b_o_->pointer();
 double* eps_b_v = this->eps_b_v_->pointer();

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

      // A block                                                                                            
      for (int i=0; i<naocc_; ++i) {            
      for (int a=0; a<navir_; ++a) {
           int ia = navir_*i + a;
                                                                                            
           Ba->set(i, a, this->guess_vectors_davidson_liu_->V(k)->get(ia    ));
      }}
      // B block
      for (int i=0; i<nbocc_; ++i) {            
      for (int a=0; a<nbvir_; ++a) {
           int ia = nbvir_*i + a;
                                                                                            
           Bb->set(i, a, this->guess_vectors_davidson_liu_->V(k)->get(ia+off));
      }}

      psi::SharedMatrix Ta = psi::Matrix::triplet(Ca_occ__, Ba, Ca_vir__, false, false, true);
      psi::SharedMatrix Tb = psi::Matrix::triplet(Cb_occ__, Bb, Cb_vir__, false, false, true);

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
      C_left.push_back(Tb_set[k-this->davidson_liu_n_sigma_computed_]);
      C_right.push_back(I1);
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

     psi::SharedMatrix Ja_ov = psi::Matrix::triplet(Cb_occ__, J[2*kk+0], Cb_vir__, true, false, false);
     psi::SharedMatrix Jb_ov = psi::Matrix::triplet(Cb_occ__, J[2*kk+1], Cb_vir__, true, false, false);
     psi::SharedMatrix Ka_ov = psi::Matrix::triplet(Cb_occ__, K[2*kk+0], Cb_vir__, true, false, false);
     psi::SharedMatrix Kb_ov = psi::Matrix::triplet(Cb_occ__, K[2*kk+1], Cb_vir__, true, false, false);

     Ja_OV->add(Jb_OV);
     Ja_ov->add(Jb_ov);

     psi::SharedMatrix S_OV = Ja_OV->clone();
     psi::SharedMatrix S_ov = Ja_ov->clone();
     S_OV->subtract(Ka_OV);
     S_ov->subtract(Kb_ov);

     psi::SharedVector Sigma = std::make_shared<psi::Vector>("", this->N_davidson_liu_);
     double* sigma = Sigma->pointer();
     double* bvec  = this->guess_vectors_davidson_liu_->V(k)->pointer();

     // A block
     for (int i=0; i<naocc_; ++i) {            
     for (int a=0; a<navir_; ++a) {
          int ia = navir_*i + a;
                                               
          double v  = eps_a_v[a] - eps_a_o[i];

          sigma[ia    ] = v * bvec[ia    ] + S_OV->get(i,a); 
     }
     }

     // B block
     for (int i=0; i<nbocc_; ++i) {            
     for (int a=0; a<nbvir_; ++a) {
          int ia = nbvir_*i + a;
                                               
          double v  = eps_b_v[a] - eps_b_o[i];

          sigma[ia+off] = v * bvec[ia+off] + S_ov->get(i,a);
     }
     }

     this->sigma_vectors_davidson_liu_.push_back(Sigma);

 }

}


} // EndNameSpace oepdev
