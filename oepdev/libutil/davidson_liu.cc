#include "davidson_liu.h"
#include <random>
#include <functional>

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
//#include "psi4/libmints/writer.h"
//#include "psi4/libmints/writer_file_prefix.h"


oepdev::DavidsonLiu::DavidsonLiu(psi::Options& opt) : 
  davidson_liu_initialized_(false), 
  davidson_liu_finalized_(false), 
  davidson_liu_n_sigma_computed_(0),
  options_(opt),
  gs_(std::make_shared<oepdev::GramSchmidt>()),
  H_diag_(nullptr),
  E_old_(nullptr),
  E_(nullptr),
  U_(nullptr)
{
}

oepdev::DavidsonLiu::~DavidsonLiu() {}

// Algorithm draw implementation is given here
void oepdev::DavidsonLiu::run_davidson_liu() {
  if (!davidson_liu_initialized_) throw psi::PSIEXCEPTION("Davidson-Liu must be initialized first!");

  double conv = 1.0e+10;
  const double eps = this->options_.get_double("DAVIDSON_LIU_CONVER");
  const int maxit = this->options_.get_int("DAVIDSON_LIU_MAXITER");
  int iter = 1;

  psi::outfile->Printf("\n ===> Starting Davidson-Liu Iterations <===\n\n");

  psi::outfile->Printf(" @Davidson-Liu: Initializing guess vectors.\n");
  this->davidson_liu_initialize_guess_vectors();

  psi::outfile->Printf(" @Davidson-Liu: Computing diagonal Hamiltonian.\n");
  this->davidson_liu_compute_diagonal_hamiltonian();

  psi::outfile->Printf("\n @Davidson-Liu: Starting iteration process.\n");
  while (conv > eps) {

     this->davidson_liu_compute_sigma();
     this->davidson_liu_add_guess_vectors();
     conv = this->davidson_liu_compute_convergence();

     psi::outfile->Printf(" @Davidson-Liu Iter=%4d Conv=%18.8f\n", iter, conv);

     iter++;
     if (iter > maxit) {
         psi::outfile->Printf(" @Davidson-Liu: Maximum iterations %d exceeded!\n", maxit);
         break;}

     if (conv < eps  ) {
         psi::outfile->Printf(" @Davidson-Liu: Maximum iterations converged successfully!\n");
         break;}

  }

  this->davidson_liu_finalize();
  psi::outfile->Printf("\n @Davidson-Liu: Done.\n");
}

// Helper interface
void oepdev::DavidsonLiu::davidson_liu_initialize(int N, int L, int M) {
 this->N_ = N;
 this->L_ = L;
 this->M_ = M;
 this->H_diag_ = std::make_shared<psi::Vector>("", N);
 this->E_old_->copy(*this->E_);
 this->davidson_liu_initialized_ = true;
}


void oepdev::DavidsonLiu::davidson_liu_initialize_guess_vectors()
{
  if (this->options_.get_str("DAVIDSON_LIU_GUESS") == "RANDOM") 
       {this->davidson_liu_initialize_guess_vectors_by_random();}
  else {this->davidson_liu_initialize_guess_vectors_by_custom();}
}

void oepdev::DavidsonLiu::davidson_liu_compute_diagonal_hamiltonian() 
{
 // abstract method: nothing to do here
}

void oepdev::DavidsonLiu::davidson_liu_initialize_guess_vectors_by_random()
{
  // Setup the random number generator
  std::default_random_engine randomNumberGenerator;
  std::uniform_real_distribution<double> randomDistribution(-1.0, 1.0);
  auto draw = std::bind( randomDistribution, randomNumberGenerator );

  // Draw random guess (non-orthonormal) vectors
  std::vector<psi::SharedVector> guess_vectors;
  for (int i=0; i<this->L_; ++i) {
       psi::SharedVector vec = std::make_shared<psi::Vector>("", this->N_);
       double* v = vec->pointer();
       for (int k=0; k<this->N_; ++k) v[k] = draw();
       guess_vectors.push_back(vec);
  }

  // Orthonormalize guess vectors
  this->gs_->reset(guess_vectors);
  this->gs_->orthonormalize();
}

void oepdev::DavidsonLiu::davidson_liu_initialize_guess_vectors_by_custom()
{
  throw psi::PSIEXCEPTION("No custom initialization of the guess vectors for Davidson-Liu algorithm implemented!");
}

void oepdev::DavidsonLiu::davidson_liu_compute_sigma() 
{
 // abstract method: nothing to do here
}

void oepdev::DavidsonLiu::davidson_liu_add_guess_vectors()
{
  // Actualize the number of sigma vectors already computed
  this->davidson_liu_n_sigma_computed_ = this->gs_->L();

  // Read the data from options
  const double threshold_large = options_.get_double("DAVIDSON_LIU_THRESH_LARGE");
  const double threshold_small = options_.get_double("DAVIDSON_LIU_THRESH_SMALL");
  const int L_max = options_.get_int("DAVIDSON_LIU_SPACE_MAX");
  this->E_old_->copy(*this->E_);

  // Form and diagonalize the sub-space Hamiltonian
  psi::SharedMatrix G = std::make_shared<psi::Matrix>("", L_, L_);   // subspace Hamiltonian
  psi::SharedVector E = std::make_shared<psi::Vector>("", L_);       // contains E_
  psi::SharedMatrix U = std::make_shared<psi::Matrix>("", N_, L_);   // contains A matrix
  double** g = G->pointer();
  double** u = U->pointer();
  double** w =U_->pointer();
  double* e = E_->pointer();
  double* h = H_diag_->pointer();

  for (int i=0; i<L_; ++i) {
       g[i][i] = gs_->V(i)->vector_dot(sigma_vectors_[i]);
       for (int j=0; j<i; ++j) {
            double v = gs_->V(i)->vector_dot(sigma_vectors_[j]);
            g[i][j] = v;
            g[j][i] = v;
       }
  }

  G->diagonalize(U, E);

  // Save current eigenpairs
  for (int k=0; k<M_; ++k) {
       e[k] = E->get(k);
       for (int n=0; n<N_; ++n) {
            double v = 0.0;
            for (int l=0; l<L_; ++l) {
                 v += u[l][k] * gs_->V(l)->get(n);
            }
            w[n][k] = v;
       }
  }

  // Enlarge the guess space
  for (int k=0; k<M_; ++k) {

       // Compute correction vector
       psi::SharedVector d = std::make_shared<psi::Vector>("", N_);

       for (int n=0; n<N_; ++n) {
            double v = -e[k] * w[n][k];
            for (int l=0; l<L_; ++l) {
                 v += sigma_vectors_[l]->get(n) * u[l][k];
            }
            d->set(n, v/(e[k] - h[n]) );
       }

       // Add vectors to the guess space
       this->gs_->orthogonalize_vector(d, false);
       double dn = d->norm();

       double threshold = threshold_small;
       if (k==0) threshold = threshold_large;
       
       if (std::abs(dn) > threshold) {
           d->scale(1.0/dn); 
           this->gs_->append(d);
           this->L_ += 1;

           // Check if the size of guess space exceeds the imposed limits
           if (gs_->L() > L_max) {
               std::vector<psi::SharedVector> new_vec;
               for (int q=0; q<M_; ++q) new_vec.push_back(this->U_->get_column(0, q));
               gs_->reset(new_vec);
               this->L_ = gs_->L();
               this->davidson_liu_n_sigma_computed_ = 0;
               break;
           }
       }
  }
}

double oepdev::DavidsonLiu::davidson_liu_compute_convergence()
{
  this->E_old_->subtract(this->E_);
  return E_old_->rms();
}

void oepdev::DavidsonLiu::davidson_liu_finalize()
{
  this->sigma_vectors_.clear();
  this->gs_->reset();
  this->davidson_liu_finalized_ = true;
}
