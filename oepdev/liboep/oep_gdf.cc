#include "oep_gdf.h"
#include <iostream>

#include "psi4/libmints/vector.h"
#include "psi4/libmints/mintshelper.h"

using namespace oepdev;

// GeneralizedDensityFit - Abstract Base
GeneralizedDensityFit::GeneralizedDensityFit()
 : G_(nullptr), 
   H_(nullptr),
   bs_a_(nullptr),
   bs_i_(nullptr),
   ints_aa_(nullptr),
   ints_ai_(nullptr),
   ints_ii_(nullptr),
   V_(nullptr),
   n_a_(0),
   n_i_(0),
   n_o_(0)
{
}
GeneralizedDensityFit::~GeneralizedDensityFit() {}
std::shared_ptr<psi::Matrix> GeneralizedDensityFit::compute(void) {}
void GeneralizedDensityFit::invert_matrix(std::shared_ptr<psi::Matrix>& M) {
  // Invert the M matrix
  std::shared_ptr<psi::Matrix> m = std::make_shared<psi::Matrix>(M);
  M->invert();
  // Perform the Identity Test
  std::shared_ptr<psi::Matrix> I = psi::Matrix::doublet(m, M);
  double s = 0.0;
  for (int a=0; a<m->ncol(); ++a) {
       for (int b=0; b<m->ncol(); ++b) {
            s += sqrt( I->get(a,b) * I->get(a,b) );
       }
  }
  s /= (double)(m->ncol());
  // std::cout << " Hessian Matrix Identity Test = " << s << std::endl; // for debugging
  if (std::abs(s-1.0)>0.000001) std::cout << " ----> Warning!! Hessian inverse has numerical error! <----\n";
}
// SingleGeneralizedDensityFit
SingleGeneralizedDensityFit::SingleGeneralizedDensityFit(
  std::shared_ptr<psi::BasisSet> bs_auxiliary,
  std::shared_ptr<psi::Matrix> v_vector)
 : GeneralizedDensityFit()
{
  n_a_ = v_vector->nrow();
  n_o_ = v_vector->ncol();
  bs_a_= bs_auxiliary;
  if (n_a_ != bs_a_->nbf()) throw psi::PSIEXCEPTION("The size of auxiliary basis set must be equal to the number of rows of V matrix!");
  V_   = std::make_shared<psi::Matrix>(v_vector);
}
SingleGeneralizedDensityFit::~SingleGeneralizedDensityFit() {}
std::shared_ptr<psi::Matrix> SingleGeneralizedDensityFit::compute(void)
{
  // Compute overlap integrals between auxiliary basis functions
  ints_aa_ = std::make_shared<oepdev::IntegralFactory>(bs_a_, bs_a_);
  std::shared_ptr<psi::Matrix> S = std::make_shared<psi::Matrix>("S", n_a_, n_a_);
  std::shared_ptr<psi::OneBodyAOInt> p(ints_aa_->ao_overlap()); p->compute(S); // ints_aa_->ao_overlap()->compute(S);--> mem leak?
  this->invert_matrix(S);
  // Perform GDF
  G_ = psi::Matrix::doublet(S, V_, false, false); 
  // Return
  return G_;
}
// DoubleGeneralizedDensityFit
DoubleGeneralizedDensityFit::DoubleGeneralizedDensityFit(
  std::shared_ptr<psi::BasisSet> bs_auxiliary,
  std::shared_ptr<psi::BasisSet> bs_intermediate,
  std::shared_ptr<psi::Matrix> v_vector)
 : GeneralizedDensityFit()
{
  n_i_ = v_vector->nrow();
  n_o_ = v_vector->ncol();
  n_a_ = bs_auxiliary->nbf();
  bs_a_= bs_auxiliary;
  bs_i_= bs_intermediate;
  if (n_i_ != bs_i_->nbf()) throw psi::PSIEXCEPTION("The size of intermediate basis set must be equal to the number of rows of V matrix!");
  V_   = std::make_shared<psi::Matrix>(v_vector);
}
DoubleGeneralizedDensityFit::~DoubleGeneralizedDensityFit() {}
std::shared_ptr<psi::Matrix> DoubleGeneralizedDensityFit::compute(void)
{
  // Compute overlap integrals between intermediate basis functions
  ints_ii_ = std::make_shared<oepdev::IntegralFactory>(bs_i_);
  std::shared_ptr<psi::Matrix> S = std::make_shared<psi::Matrix>("S", n_i_, n_i_);
  std::shared_ptr<psi::OneBodyAOInt> p(ints_ii_->ao_overlap()) ; p->compute(S); // ints_ii_->ao_overlap()->compute(S);--> mem leak?
  this->invert_matrix(S);

  // Perform 1st GDF
  H_ = psi::Matrix::doublet(S, V_, false, false); 

  // Compute \f$ \left( \xi \vert\vert \epsilon \right) \f$ integrals
  ints_ai_ = std::make_shared<oepdev::IntegralFactory>(bs_a_, bs_i_);
  std::shared_ptr<psi::Matrix> R = std::make_shared<psi::Matrix>("R Matrix", n_a_, n_i_);
  std::shared_ptr<oepdev::TwoBodyAOInt> pai(ints_ai_->eri_1_1()); pai->compute(R, 0, 1); // ints_ai_->eri_1_1()->compute(R, 0, 1);

  // Compute \f$ \left( \xi \vert\vert \xi' \right) \f$ integrals
  ints_aa_ = std::make_shared<oepdev::IntegralFactory>(bs_a_);
  std::shared_ptr<psi::Matrix> A = std::make_shared<psi::Matrix>("A", n_a_, n_a_);
  std::shared_ptr<oepdev::TwoBodyAOInt> paa(ints_aa_->eri_1_1()); paa->compute(A); // ints_aa_->eri_1_1()->compute(A);
  this->invert_matrix(A);

  // Perform 2nd GDF
  G_ = psi::Matrix::triplet(A, R, H_, false, false, false); 

  // Return
  return G_;
}
// OverlapGeneralizedDensityFit
OverlapGeneralizedDensityFit::OverlapGeneralizedDensityFit(
  std::shared_ptr<psi::BasisSet> bs_auxiliary,
  std::shared_ptr<psi::BasisSet> bs_intermediate,
  std::shared_ptr<psi::Matrix> v_vector)
 : GeneralizedDensityFit()
{
  n_i_ = v_vector->nrow();
  n_o_ = v_vector->ncol();
  n_a_ = bs_auxiliary->nbf();
  bs_a_= bs_auxiliary;
  bs_i_= bs_intermediate;
  if (n_i_ != bs_i_->nbf()) throw psi::PSIEXCEPTION("The size of intermediate basis set must be equal to the number of rows of V matrix!");
  V_   = std::make_shared<psi::Matrix>(v_vector);
}
OverlapGeneralizedDensityFit::~OverlapGeneralizedDensityFit() {}
std::shared_ptr<psi::Matrix> OverlapGeneralizedDensityFit::compute(void)
{
  // Compute overlap integrals between intermediate basis functions
  ints_ii_ = std::make_shared<oepdev::IntegralFactory>(bs_i_);
  std::shared_ptr<psi::Matrix> S_ii = std::make_shared<psi::Matrix>("S", n_i_, n_i_);
  std::shared_ptr<psi::OneBodyAOInt> p(ints_ii_->ao_overlap()) ; p->compute(S_ii); 
  psi::SharedMatrix S_ii_copy = S_ii->clone();
  this->invert_matrix(S_ii);

  // Perform 1st GDF
  H_ = psi::Matrix::doublet(S_ii, V_, false, false); 

  // Find matrix T_i
  psi::SharedMatrix x = S_ii_copy->clone(); x->power(-0.5000000);
  psi::SharedMatrix y = S_ii_copy->clone(); y->power( 0.5000000);
  psi::SharedMatrix h = psi::Matrix::doublet(y, H_, false, false);
  psi::SharedMatrix c = psi::Matrix::doublet(h, h, false, true);
  const int n = c->nrow();
  psi::SharedVector l = std::make_shared<psi::Vector>("", n);
  psi::SharedMatrix u = std::make_shared<psi::Matrix>("", n, n);
  c->diagonalize(u, l, psi::diagonalize_order::descending);
  const double eps = 1.0e-7;
  int I = 0;
  for (int i=0; i<n; ++i) {
       if (l->get(i) < eps) break;
       I += 1;
  }
  psi::SharedMatrix T_i_ = std::make_shared<psi::Matrix>("", n, I);
  for (int i=0; i<I; ++i) {
       psi::SharedVector col = u->get_column(0, i);
       T_i_->set_column(0, i, col);
  }
  psi::SharedMatrix T_i = psi::Matrix::doublet(x, T_i_, false, false);
  T_i_.reset(); l.reset(); u.reset();

  // Find matrix T_m
  std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(bs_i_);
  psi::SharedMatrix s_im = mints->ao_overlap(bs_i_, bs_a_);                 
  psi::SharedMatrix s_mm = mints->ao_overlap(bs_a_, bs_a_); s_mm->invert();
  psi::SharedMatrix A = psi::Matrix::doublet(T_i, s_im, true, false);
  psi::SharedMatrix B = psi::Matrix::doublet(A, s_mm, false, false);
  psi::SharedMatrix S_BBm = psi::Matrix::doublet(B, A, false, true);
  S_BBm->power(-0.50000000000000);
  psi::SharedMatrix Y = psi::Matrix::triplet(s_mm, s_im, T_i, false, true, false);
  psi::SharedMatrix T_m = psi::Matrix::doublet(Y, S_BBm, false, false);
  S_BBm->identity();

  // Find approximate G
  psi::SharedMatrix C = psi::Matrix::triplet(T_m, S_BBm, T_m, false, false, true);
  G_ = psi::Matrix::triplet(C, s_im, H_, false, true, false);

  // Return
  return G_;
}

// Static Factory Methods
std::shared_ptr<GeneralizedDensityFit> 
GeneralizedDensityFit::build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                             std::shared_ptr<psi::Matrix> v_vector)
{
  return std::make_shared<SingleGeneralizedDensityFit>(bs_auxiliary, v_vector);
}
std::shared_ptr<GeneralizedDensityFit> 
GeneralizedDensityFit::build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                             std::shared_ptr<psi::BasisSet> bs_intermediate,
                             std::shared_ptr<psi::Matrix> v_vector)
{
  return std::make_shared<DoubleGeneralizedDensityFit>(bs_auxiliary, bs_intermediate, v_vector);
}
std::shared_ptr<GeneralizedDensityFit>
GeneralizedDensityFit::build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                             std::shared_ptr<psi::BasisSet> bs_intermediate,
                             std::shared_ptr<psi::Matrix> v_vector,
                             int dummy)
{
  return std::make_shared<OverlapGeneralizedDensityFit>(bs_auxiliary, bs_intermediate, v_vector);
}
