#include "oep_gdf.h"

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
  // cout << " Hessian Matrix Identity Test = " << s << endl; // for debugging
  if (std::abs(s-1.0)>0.000001) cout << " ----> Warning!! Hessian inverse has numerical error! <----\n";
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
  ints_aa_->ao_overlap()->compute(S);
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
  ints_ii_->ao_overlap()->compute(S);
  this->invert_matrix(S);

  // Perform 1st GDF
  H_ = psi::Matrix::doublet(S, V_, false, false); 

  // Compute \f$ \left( \xi \vert\vert \epsilon \right) \f$ integrals
  ints_ai_ = std::make_shared<oepdev::IntegralFactory>(bs_a_, bs_i_);
  std::shared_ptr<psi::Matrix> R = std::make_shared<psi::Matrix>("R Matrix", n_a_, n_i_);
  ints_ai_->eri_1_1()->compute(R, 0, 1);

  // Compute \f$ \left( \xi \vert\vert \xi' \right) \f$ integrals
  ints_aa_ = std::make_shared<oepdev::IntegralFactory>(bs_a_);
  std::shared_ptr<psi::Matrix> A = std::make_shared<psi::Matrix>("A", n_a_, n_a_);
  ints_aa_->eri_1_1()->compute(A);
  this->invert_matrix(A);

  // Perform 2nd GDF
  G_ = psi::Matrix::triplet(A, R, H_, false, false, false); 

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

