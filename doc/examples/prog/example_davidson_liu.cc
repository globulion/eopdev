// Define a class that inherits from DavidsonLiu
class Diagonalize : public DavidsonLiu {

 public:
  Diagonalize(psi::SharedMatrix H, psi::Options& opt, int M);
  virtual void ~Diagonalize() {};

 protected: 
  // Desired number of roots to be found
  int M_;
  // Matrix to be diagonalized (explicitly stored)
  psi::SharedMatrix matrix_;

  // Implementation of pure methods must be declared
  void davidson_liu_compute_diagonal_hamiltonian();
  void davidson_liu_compute_sigma();
};

Diagonalize::Diagonalize(psi::SharedMatrix H, psi::Options& opt, int M) :
  DavidsonLiu(opt) , matrix_(nullptr), M_(M)
{
  matrix_ = std::make_shared<psi::Matrix>(H);
  int N = H->ncol();
  int L = M_;
  
  // Must be run in order to allocate memory
  this->davidson_liu_initialize(N, L, M_);
}

// Implementation of pure methods
void Diagonalize::davidson_liu_compute_diagonal_hamiltonian() {
 for (int i=0; i<this->matrix_->ncol(); ++i) {
      double v = matrix_->get(i, i);
      this->H_diag_davidson_liu_->set(i, v);
 }
}

void Diagonalize::davidson_liu_compute_sigma() {
 for (int k=this->davidson_liu_n_sigma_computed_; k<this->L_davidson_liu_; ++k) {
      psi::SharedVector Sigma = std::make_shared<psi::Vector>("", this->N_davidson_liu_);
      Sigma->gemv(false, 1.0, *this->matrix_, *Sigma, 0.0);
      this->sigma_vectors_davidson_liu_.push_back(Sigma);
 }
}

// Testing function
void example_davidson_liu(psi::SharedMatrix H, int M, psi::Options& opt){

  // Construct the solver object
  Diagonalize solver(H, opt, M);

  // Find *M* lowest eigenpairs of a given matrix **H**
  solver.run_davidson_liu();
};
