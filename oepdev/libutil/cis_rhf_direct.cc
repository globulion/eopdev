#include "cis.h"
#include "integrals_iter.h"
#include <iostream>

namespace oepdev{

R_CISComputer_Direct::R_CISComputer_Direct(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 R_CISComputer(wfn, opt)
{
}

R_CISComputer_Direct::~R_CISComputer_Direct() {}

void R_CISComputer_Direct::transform_integrals_(void) 
{
 // Do nothing here since it is direct CIS
}

void R_CISComputer_Direct::build_hamiltonian_(void) {

 this->H_->zero();
 double** H = this->H_->pointer();
 //double** Fa_oo = this->Fa_oo_->pointer();
 //double** Fa_vv = this->Fa_vv_->pointer();
 double* eps_a_o = this->eps_a_o_->pointer();
 double* eps_a_v = this->eps_a_v_->pointer();
 const int off = this->naocc_ * this->navir_;

 // Fock matrix contribution directly in MO basis
 for (int i=0; i<naocc_; ++i) {
 for (int a=0; a<navir_; ++a) {
      int ia = navir_*i + a;

      H[ia][ia] = eps_a_v[a] - eps_a_o[i];
 }
 }

 // ERI contribution in AO basis on the fly
 psi::IntegralFactory fact(ref_wfn_->basisset());
 std::shared_ptr<psi::TwoBodyAOInt> tei(fact.eri());                                           
 const double * buffer = tei->buffer();

 std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact, "ALL");
 
 double** ca_occ = ref_wfn_->Ca_subset("AO","OCC")->pointer();
 double** ca_vir = ref_wfn_->Ca_subset("AO","VIR")->pointer();

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

               for (int i=0; i<naocc_; ++i) {                
               for (int a=0; a<navir_; ++a) {
                    int ia = navir_*i + a;
                                                            
	            double cIi = ca_occ[I][i];
	            double cJa = ca_vir[J][a];
                    double cKa = ca_vir[K][a];
                                                             
                    for (int j=0; j<naocc_; ++j) {
                    for (int b=0; b<navir_; ++b) {
      	                 int jb = navir_*j + b;
                         
	                 double cJj = ca_occ[J][j];
	                 double cKj = ca_occ[K][j];
	                 double cLb = ca_vir[L][b];
                                                             
      	                 double c_1 = cIi * cJa * cKj * cLb;
      	                 double c_2 = cIi * cJj * cKa * cLb;
                                                             
                         H[ia][jb    ] += (c_1 - c_2) * eri;
      	                 H[ia][jb+off] +=  c_1        * eri;
      	            }
      	            }
               }
               }

           }

      }
 }
 
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
