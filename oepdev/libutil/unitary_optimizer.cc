#include <cmath>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include "psi4/psi4-dec.h"
#include "unitary_optimizer.h"

namespace oepdev {

UnitaryOptimizer::UnitaryOptimizer(int n, double conv, int maxiter, bool verbose) :
 n_(n),
 R_(nullptr),
 P_(nullptr),
 R0_(nullptr),
 P0_(nullptr),
 X_(nullptr),
 W_(nullptr),
 Xold_(nullptr),
 Xnew_(nullptr),
 conv_(conv),
 maxiter_(maxiter),
 verbose_(verbose),
 niter_(0),
 conv_current_(1.0e+18),
 success_(false)
{
}
UnitaryOptimizer::UnitaryOptimizer(double* R, double* P, int n, double conv, int maxiter, bool verbose) 
 : UnitaryOptimizer::UnitaryOptimizer(n, conv, maxiter, verbose)
{
  this->common_init_();

  for (int i = 0; i < n_   ; ++i) {
      P0_[i] = P[i];
      P_ [i] = P[i];
  }
  for (int i = 0; i < n_*n_; ++i) {
      R0_[i] = R[i];
      R_ [i] = R[i];
  }

  this->Zinit_ = this->eval_Z_(X_, R0_, P0_);
  this->Zold_  = this->Zinit_;
}
UnitaryOptimizer::UnitaryOptimizer(std::shared_ptr<psi::Matrix> R, std::shared_ptr<psi::Vector> P, double conv, int maxiter, bool verbose)
 : UnitaryOptimizer::UnitaryOptimizer(R->ncol(), conv, maxiter, verbose)
{
  this->common_init_();

  for (int i = 0; i < n_; ++i) {
       P0_[i] = P->get(i);
       P_ [i] = P->get(i);
       for (int j = 0; j < n_; ++j) {
            R0_[IDX(i,j,n_)] = R->get(i,j);
            R_ [IDX(i,j,n_)] = R->get(i,j);
       }
  }

  this->Zinit_ = this->eval_Z_(X_, R0_, P0_);
  this->Zold_  = this->Zinit_;
}
void UnitaryOptimizer::common_init_()
{
  // Allocate and initialize
  std::memset(S_, 0, 4*sizeof(double));
  P_  = new double[n_];
  R_  = new double[n_*n_];
  P0_ = new double[n_];
  R0_ = new double[n_*n_];
  X_  = new double[n_*n_];
  W_  = new double[n_*n_];
  Xold_ = new double[n_*n_];
  Xnew_ = new double[n_*n_];

  for (int i = 0; i < n_; ++i) {
       X_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
            if (i!=j) X_[IDX(i,j,n_)] = 0.0;
       }
  }
}
UnitaryOptimizer::~UnitaryOptimizer() 
{
  delete[] R_;
  delete[] P_;
  delete[] R0_;
  delete[] P0_;
  delete[] X_;
  delete[] W_;
  delete[] Xold_;
  delete[] Xnew_;
}
bool UnitaryOptimizer::maximize() {this->run_("max"); return success_;}
bool UnitaryOptimizer::minimize() {this->run_("min"); return success_;}
void UnitaryOptimizer::run_(const std::string& opt) 
{
  this->refresh_();
  this->optimize_(opt);
}
void UnitaryOptimizer::optimize_(const std::string& opt) 
{
  if (verbose_) 
     psi::outfile->Printf(" Start  : Z[1] = %15.6f\n", Zold_);
  while (conv_current_ > conv_) {
     this->form_next_X_(opt);
     this->update_RP_();
     this->update_Z_();
     this->update_X_();
     this->update_conv_();
     this->update_iter_();
     if (verbose_) 
         psi::outfile->Printf(" Iter %2d: Z[X] = %15.6f  Conv= %15.6f\n", niter_, Znew_, conv_current_);
     if (niter_ > maxiter_) {
         psi::outfile->Printf(" Optimization unsuccesfull! Maximum iteration number %d exceeded!\n", maxiter_);
         success_ = false;
         break;
     }
  }
  success_ = ((niter_ <= maxiter_) ? true : false );
  if (verbose_ && success_) {
      psi::outfile->Printf(" Optimization succesfull!\n");
      psi::outfile->Printf(" Optimized Z[X] value: %15.6f\n", this->Z());
  }
  
}
void UnitaryOptimizer::refresh_()
{
   for (int i = 0; i < n_*n_; ++i) R_[i] = R0_[i];
   for (int i = 0; i < n_   ; ++i) {
        P_[i] = P0_[i];
        X_[IDX(i,i,n_)] = 1.0; 
        for (int j = 0; j < n_; ++j) {
             if (i!=j) X_[IDX(i,j,n_)] = 0.0;
        }
   }
   std::memset(S_, 0, 4*sizeof(double));
   // Compute initial Z value
   this->Zinit_ = this->eval_Z_(X_, R0_, P0_);
   this->Zold_  = this->Zinit_;
   this->niter_ = 0;
   this->conv_current_ = 1.0e+18;
   this->success_ = false;
}
void UnitaryOptimizer::update_conv_()
{
  conv_current_ = std::abs(Znew_-Zold_);
}
void UnitaryOptimizer::update_iter_()
{
  for (int i = 0; i < n_*n_; ++i) Xold_[i] = Xnew_[i];
  Zold_ = Znew_;
  niter_ += 1;
}
void UnitaryOptimizer::update_Z_()
{
   double Z = 0.0;
   for (int i = 0; i < n_; ++i) {
        Z -= P_[i];
        for (int j = 0; j < n_; ++j) {
             Z += R_[IDX(i,j,n_)];
        }
   }
   Znew_ = Z;
}
void UnitaryOptimizer::update_RP_()
{
   // Update P vector
   for (int i = 0; i < n_; ++i) {
        double v = 0.0;
        for (int j = 0; j < n_; ++j) {
             v += Xnew_[IDX(i,j,n_)] * P_[j];
        }
        W_[i] = v;
   }
   for (int i = 0; i < n_; ++i) P_[i] = W_[i];

   // Uptade R matrix
   for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < n_; ++j) {
             double v = 0.0;
             for (int k = 0; k < n_; ++k) {
                  for (int l = 0; l < n_; ++l) {
                       v += Xnew_[IDX(i,k,n_)] * Xnew_[IDX(j,l,n_)] * R_[IDX(k,l,n_)];
                  }
             }
             W_[IDX(i,j,n_)] = v;
        }
   }
   for (int i = 0; i < n_*n_; ++i) R_[i] = W_[i];
}
void UnitaryOptimizer::update_X_()
{
   for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < n_; ++j) {
             double v = 0.0;
             for (int k = 0; k < n_; ++k) {
                  v += Xnew_[IDX(i,k,n_)] * X_[IDX(k,j,n_)];
             }
             W_[IDX(i,j,n_)] = v;
        }
   }
   for (int i = 0; i < n_*n_; ++i) X_[i] = W_[i];
}
double UnitaryOptimizer::eval_Z_(double* X, double* R, double* P)
{
   double Z = 0.0;
   for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < n_; ++j) {
             Z -= X[IDX(i,j,n_)] * P[j];
             for (int k = 0; k < n_; ++k) {
                  for (int l = 0; l < n_; ++l) {
                       Z += X[IDX(i,j,n_)] * X[IDX(k,l,n_)] * R[IDX(j,l,n_)];
                  }
             }
        }
   }
   return Z;
}
double UnitaryOptimizer::eval_Z_()
{
   return this->eval_Z_(X_, R0_, P0_);
}
double UnitaryOptimizer::eval_dZ_(double g, double* R, double* P, int i, int j)
{
  double c = cos(g); 
  double s = sin(g);
  double s2= sin(2.0*g);
  double c2= cos(2.0*g);
  double dZ = (1.0 - c) * (P[i] + P[j]) 
                    + s * (P[i] - P[j])
                    + s2* (R[IDX(j,j,n_)] - R[IDX(i,i,n_)])
             - 2.0 * s*s* (R[IDX(j,i,n_)] + R[IDX(i,j,n_)]); 
  
  for (int ii = 0; ii < n_; ++ii) {
       if ((ii != i) && (ii != j)) {
           dZ -= (1.0 - c) * ( R_[IDX(j,ii,n_)] + R_[IDX(i,ii,n_)] + R_[IDX(ii,j,n_)] + R_[IDX(ii,i,n_)] ); 
           dZ +=        s  * ( R_[IDX(j,ii,n_)] - R_[IDX(i,ii,n_)] + R_[IDX(ii,j,n_)] - R_[IDX(ii,i,n_)] );
       }
  }
  return dZ;
}
double UnitaryOptimizer::eval_Z_trial_(int ii, int jj, double gamma)
{
  double Z = 0.0; 
  for (int i = 0; i < n_; ++i) {
       W_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
           if (i!=j) W_[IDX(i,j,n_)] = 0.0;
       }
  }
  double c = cos(gamma);
  double s = sin(gamma);
  W_[IDX(ii,ii,n_)] = c;
  W_[IDX(jj,jj,n_)] = c;
  W_[IDX(ii,jj,n_)] = s;
  W_[IDX(jj,ii,n_)] =-s;

  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            Z -= W_[IDX(i,j,n_)] * P_[j];
            for (int k = 0; k < n_; ++k) {
                 for (int l = 0; l < n_; ++l) {
                      Z += W_[IDX(i,j,n_)] * W_[IDX(k,l,n_)] * R_[IDX(j,l,n_)];
                 }
            }
       }
  }
  return Z;
}
void UnitaryOptimizer::form_X0_()
{
   for (int i = 0; i < n_; ++i) {
        Xold_[IDX(i,i,n_)] = 1.0;
        for (int j = 0; j < n_; ++j) {
             if (i != j) Xold_[IDX(i,j,n_)] = 0.0;
        }
   }
}
void UnitaryOptimizer::form_X_(int i, int j, double gamma)
{
  for (int ii = 0; ii < n_; ++ii) {
       Xnew_[IDX(ii,ii,n_)] = 1.0;
       for (int jj = 0; jj < n_; ++jj) {
            if (ii != jj) Xnew_[IDX(ii,jj,n_)] = 0.0;
       }
  }
  double c = cos(gamma);
  double s = sin(gamma);
  Xnew_[IDX(i,i,n_)] = c;
  Xnew_[IDX(j,j,n_)] = c;
  Xnew_[IDX(i,j,n_)] = s;
  Xnew_[IDX(j,i,n_)] =-s;
}
void UnitaryOptimizer::form_next_X_(const std::string& opt)
{
  bool (UnitaryOptimizer::*optfunc)(double, double) = (opt == "min" ? &UnitaryOptimizer::lt_ : &UnitaryOptimizer::gt_);

  int I = 0;
  int J = 1;
  double Gamma = 0.0;
  double dZold = 1.0e+18; if (opt == "max") dZold = -1.0e+18;
  for (int j = 0; j < n_; ++j) {
       for (int i = 0; i < j; ++i) {
            ABCD abcd = this->get_ABCD_(i, j);
            this->find_roots_boyd_(abcd);
            double gamma = this->find_gamma_(abcd, i, j, opt);
            double dZ = this->eval_dZ_(gamma, R_, P_, i, j);
            if ((this->*optfunc)(dZ, dZold)) {
                Gamma = gamma;
                I = i; J = j;
                dZold = dZ;
            }
       }
  }
  this->form_X_(I, J, Gamma);
}
ABCD UnitaryOptimizer::get_ABCD_(int i, int j)
{
  double A = P_[i] + P_[j];
  double B = P_[i] - P_[j];
  double C =-2.0 * (R_[IDX(i,j,n_)] + R_[IDX(j,i,n_)]);
  double D = 2.0 * (R_[IDX(j,j,n_)] - R_[IDX(i,i,n_)]);
  for (int ii = 0; ii < n_; ++ii) {
       if ((ii != i) && (ii != j)) {
           A -= R_[IDX(j,ii,n_)] + R_[IDX(i,ii,n_)] + R_[IDX(ii,j,n_)] + R_[IDX(ii,i,n_)];
           B += R_[IDX(j,ii,n_)] - R_[IDX(i,ii,n_)] + R_[IDX(ii,j,n_)] - R_[IDX(ii,i,n_)];
       }
  }
  ABCD abcd = {A, B, C, D};
  return abcd;
}
void UnitaryOptimizer::find_roots_boyd_(const ABCD& abcd)
{
    // Allocate
    Eigen::Matrix4cd B;
    Eigen::ComplexEigenSolver<Eigen::Matrix4cd> eigen;

    // Build up B matrix
    std::complex<double> a1= -(abcd.D + abcd.C * 1_i);
    std::complex<double> a2= -(abcd.B + abcd.A * 1_i);
    std::complex<double> a3=   0.0    + 0.0    * 1_i ;
    std::complex<double> a4= -(abcd.B - abcd.A * 1_i);
    std::complex<double> d =   abcd.D - abcd.C * 1_i ;

    B(0,0) = 0.0; B(0,1) = 1.0; B(0,2) = 0.0; B(0,3) = 0.0;
    B(1,0) = 0.0; B(1,1) = 0.0; B(1,2) = 1.0; B(1,3) = 0.0;
    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 0.0; B(2,3) = 1.0;
    B(3,0) =  a1/ d;
    B(3,1) =  a2/ d;
    B(3,2) =  a3   ;
    B(3,3) =  a4/ d;

    // Compute all the roots
    eigen.compute(B);
    Eigen::Vector4cd roots = -1_i * eigen.eigenvalues().array().log();

    // Save in work place
    for (int i = 0; i < 4; ++i) {
         //if (std::abs(roots.imag()(i)) > 1.0e-4) std::cout << " Warning! Imaginary root is non-zero = " << roots.imag()(i) << std::endl;
         double r = roots.real()(i);
         if (r < 0.0) r += 2.0 * M_PI;
         S_[i] = r;
    }
}
double UnitaryOptimizer::find_root_halley_(double x0, const ABCD& abcd)
{
  double xold = x0;
  double xnew;
  double conv = 1.0e+8;
  int niter = 0;
  while (conv > 1.0e-8) {
     double f0 = func_0_(xold, abcd);
     double f1 = func_1_(xold, abcd);
     double f2 = func_2_(xold, abcd);
     xnew = xold - 2.0 * f1 / (2.0*f1*f1 - f0*f2);
     conv = std::abs(xold - xnew);
     xold = xnew;
     niter+= 1;
     if (niter > 3000) break;
  }
  return xnew;
}
double UnitaryOptimizer::find_gamma_(const ABCD& abcd, int i, int j, const std::string& opt)
{
   double gamma = 0.0;
   if (opt == "min") {
       double Zold = 1.0e+18;
       for (int ii = 0; ii < 4; ++ii) {
            if (this->func_1_(S_[ii], abcd) > 0.0) {
                double Z = this->eval_Z_trial_(i, j, S_[ii]);
                if (Z < Zold) {
                    gamma = S_[ii];
                    Zold = Z;
                }
            }
       }
   } else { 
       double Zold =-1.0e+18;
       for (int ii = 0; ii < 4; ++ii) {
            if (this->func_1_(S_[ii], abcd) < 0.0) {
                double Z = this->eval_Z_trial_(i, j, S_[ii]);
                if (Z > Zold) {
                    gamma = S_[ii];
                    Zold = Z;
                }
            }
       }
   }
   return gamma;
}
bool UnitaryOptimizer::lt_(double a, double b)
{
   if (a < b) return true;
   return false;
}
bool UnitaryOptimizer::gt_(double a, double b)
{
   if (a < b) return false;
   return true;
}
double UnitaryOptimizer::func_0_(double g, const ABCD& abcd)
{
  double s = sin(g);
  double c = cos(g);
  double s2= sin(2.0*g);
  double c2= cos(2.0*g);
  return abcd.A * s + abcd.B * c + abcd.C * s2 + abcd.D * c2;
}
double UnitaryOptimizer::func_1_(double g, const ABCD& abcd)
{
  double s = sin(g);
  double c = cos(g);
  double s2= sin(2.0*g);
  double c2= cos(2.0*g);
  return abcd.A * c - abcd.B * s + 2.0 * abcd.C * c2 - 2.0 * abcd.D * s2;
}
double UnitaryOptimizer::func_2_(double g, const ABCD& abcd)
{
  double s = sin(g);
  double c = cos(g);
  double s2= sin(2.0*g);
  double c2= cos(2.0*g);
  return -(abcd.A * s + abcd.B * c + 4.0 * abcd.C * s2 + 4.0 * abcd.D * c2);
}
std::shared_ptr<psi::Matrix> UnitaryOptimizer::psi_X_()
{
  std::shared_ptr<psi::Matrix> X = std::make_shared<psi::Matrix>("MO-MO Transformation matrix X", n_, n_);
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            X->set(i, j, X_[IDX(i,j,n_)]);
       }
  }
  return X;
}
//-------------------------------------------------------------------------------------------------------------
UnitaryOptimizer_4_2::UnitaryOptimizer_4_2(int n, double conv, int maxiter, bool verbose) : // same
 n_(n),
 R_(nullptr),
 P_(nullptr),
 R0_(nullptr),
 P0_(nullptr),
 X_(nullptr),
 W_(nullptr),
 Xold_(nullptr),
 Xnew_(nullptr),
 conv_(conv),
 maxiter_(maxiter),
 verbose_(verbose),
 niter_(0),
 conv_current_(1.0e+18),
 success_(false),
 n2_(n*n),
 n3_(n*n*n),
 n4_(n*n*n*n),
 n5_(n*n*n*n*n)
{
}
UnitaryOptimizer_4_2::UnitaryOptimizer_4_2(double* R, double* P, int n, double conv, int maxiter, bool verbose) 
 : UnitaryOptimizer_4_2::UnitaryOptimizer_4_2(n, conv, maxiter, verbose) // override
{
  this->common_init_();

  for (int i = 0; i < n_*n_*n_; ++i) {
      P0_[i] = P[i];
      P_ [i] = P[i];
  }
  for (int i = 0; i < n_*n_*n_*n_*n_*n_; ++i) {
      R0_[i] = R[i];
      R_ [i] = R[i];
  }

  this->Zinit_ = this->eval_Z_(X_, R0_, P0_);
  this->Zold_  = this->Zinit_;
}
void UnitaryOptimizer_4_2::common_init_() // override
{
  // Allocate and initialize
  std::memset(S_, 0, 8*sizeof(double));
  size_t n2 = n_*n_;
  size_t n3 = n_*n_*n_;
  size_t n6 = n3*n3;
  try {
     R_  = new double[n6];
     R0_ = new double[n6];
     W_  = new double[n6];
  } catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating 6-rank tensors \n%s\n", e.what());
        exit(EXIT_FAILURE);
  }
  P_  = new double[n3];
  P0_ = new double[n3];
  X_  = new double[n2];
  Xold_ = new double[n2];
  Xnew_ = new double[n2];

  for (int i = 0; i < n_; ++i) {
       X_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
            if (i!=j) X_[IDX(i,j,n_)] = 0.0;
       }
  }
}
UnitaryOptimizer_4_2::~UnitaryOptimizer_4_2() // same
{
  delete[] R_;
  delete[] P_;
  delete[] R0_;
  delete[] P0_;
  delete[] X_;
  delete[] W_;
  delete[] Xold_;
  delete[] Xnew_;
}
bool UnitaryOptimizer_4_2::maximize() {this->run_("max"); return success_;} // same
bool UnitaryOptimizer_4_2::minimize() {this->run_("min"); return success_;} // same
void UnitaryOptimizer_4_2::run_(const std::string& opt) // same
{
  this->refresh_();
  this->optimize_(opt);
}
void UnitaryOptimizer_4_2::optimize_(const std::string& opt) // same
{
  if (verbose_) 
     psi::outfile->Printf(" Start  : Z[1] = %15.6f\n", Zold_);
  while (conv_current_ > conv_) {
     this->form_next_X_(opt);
     this->update_RP_();
     this->update_Z_();
     this->update_X_();
     this->update_conv_();
     this->update_iter_();
     if (verbose_) 
         psi::outfile->Printf(" Iter %2d: Z[X] = %15.6f  Conv= %15.6f\n", niter_, Znew_, conv_current_);
     if (niter_ > maxiter_) {
         psi::outfile->Printf(" Optimization unsuccesfull! Maximum iteration number %d exceeded!\n", maxiter_);
         success_ = false;
         break;
     }
  }
  success_ = ((niter_ <= maxiter_) ? true : false );
  if (verbose_ && success_) {
      psi::outfile->Printf(" Optimization succesfull!\n");
      psi::outfile->Printf(" Optimized Z[X] value: %15.6f\n", this->Z());
  }
  
}
void UnitaryOptimizer_4_2::refresh_() // override (can be same if introducing nr_, np_ and ns_ dimensions)
{
   for (int i = 0; i < n5_*n_; ++i) R_[i] = R0_[i];
   for (int i = 0; i < n3_   ; ++i) P_[i] = P0_[i];
   for (int i = 0; i < n_   ; ++i) {
        X_[IDX(i,i,n_)] = 1.0; 
        for (int j = 0; j < n_; ++j) {
             if (i!=j) X_[IDX(i,j,n_)] = 0.0;
        }
   }
   std::memset(S_, 0, 8*sizeof(double));
   // Compute initial Z value
   this->Zinit_ = this->eval_Z_(X_, R0_, P0_);
   this->Zold_  = this->Zinit_;
   this->niter_ = 0;
   this->conv_current_ = 1.0e+18;
   this->success_ = false;
}
void UnitaryOptimizer_4_2::update_conv_() // same
{
  conv_current_ = std::abs(Znew_-Zold_);
}
void UnitaryOptimizer_4_2::update_iter_() // same
{
  for (int i = 0; i < n_*n_; ++i) Xold_[i] = Xnew_[i];
  Zold_ = Znew_;
  niter_ += 1;
}
void UnitaryOptimizer_4_2::update_Z_() // override
{
   double Z = 0.0;
   for (int i = 0; i < n_; ++i) {
        Z += P_[IDX3(i,i,i)];
        for (int j = 0; j < n_; ++j) {
             Z += R_[IDX6(i,j,i,j,i,j)];
        }
   }
   Znew_ = Z;
}
void UnitaryOptimizer_4_2::update_RP_() // override
{
  // P tensor
  for (int i = 0; i < n_; ++i) {
       for (int J = 0; J < n_; ++J) {
            for (int K = 0; K < n_; ++K) {
                 double v = 0.0;
                 for (int j = 0; j < n_; ++j) {
                      for (int k = 0; k < n_; ++k) {
                           v += Xnew_[IDX(j,J,n_)] * Xnew_[IDX(k,K,n_)] * P_[IDX3(i,j,k)];
                      }
                 }
                 W_[IDX3(i,J,K)] = v;
            }
       }
  }
  for (int i = 0; i < n3_; ++i) P_[i] = W_[i];

  // R tensor
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            for (int K = 0; K < n_; ++K) {
            for (int L = 0; L < n_; ++L) {
            for (int M = 0; M < n_; ++M) {
            for (int N = 0; N < n_; ++N) {
                 double v = 0.0;
                 for (int k = 0; k < n_; ++k) {
                      for (int l = 0; l < n_; ++l) {
                           for (int m = 0; m < n_; ++m) {
                                for (int n = 0; n < n_; ++n) {
                                     v += Xnew_[IDX(k,K,n_)] * Xnew_[IDX(l,L,n_)] * Xnew_[IDX(m,M,n_)] * Xnew_[IDX(n,N,n_)] * R_[IDX6(i,j,k,l,m,n)];
                                }
                           }
                      }
                 }
                 W_[IDX6(i,j,K,L,M,N)] = v;
            }}}}
       }
  }
  for (int i = 0; i < n5_*n_; ++i) R_[i] = W_[i];
}
void UnitaryOptimizer_4_2::update_X_() // override (can be same if adding argument 'left'==bool true/false)
{
   for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < n_; ++j) {
             double v = 0.0;
             for (int k = 0; k < n_; ++k) {
                  v += X_[IDX(i,k,n_)] * Xnew_[IDX(k,j,n_)]; // here - left or right multiplication (depends on the problem)
             }
             W_[IDX(i,j,n_)] = v;
        }
   }
   for (int i = 0; i < n_*n_; ++i) X_[i] = W_[i];
}
double UnitaryOptimizer_4_2::eval_Z_(double* X, double* R, double* P) // override
{
  double Z = 0.0;
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            for (int k = 0; k < n_; ++k) {
                Z += X[IDX(j,i,n_)] * X[IDX(k,i,n_)] * P[IDX3(i,j,k)];
                for (int l = 0; l < n_; ++l) {
                     for (int m = 0; m < n_; ++m) {
                          for (int n = 0; n < n_; ++n) {
                               Z += X[IDX(k,i,n_)] * X[IDX(l,j,n_)] * X[IDX(m,i,n_)] * X[IDX(n,j,n_)] * R[IDX6(i,j,k,l,m,n)];
                          }
                     }
                }
            }
       }
  }
  return Z;
}
double UnitaryOptimizer_4_2::eval_Z_() // same
{
   return this->eval_Z_(X_, R0_, P0_);
}
double UnitaryOptimizer_4_2::eval_dZ_(double gamma, double* R, double* P, int I, int J) // override
{
  // New Z
  for (int i = 0; i < n_; ++i) {
       W_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
           if (i!=j) W_[IDX(i,j,n_)] = 0.0;
       }
  }
  double c = cos(gamma);
  double s = sin(gamma);
  W_[IDX(I,I,n_)] = c;
  W_[IDX(J,J,n_)] = c;
  W_[IDX(I,J,n_)] = s;
  W_[IDX(J,I,n_)] =-s;

  double Znew = this->eval_Z_(W_, R, P);

  // Old Z
  for (int i = 0; i < n_; ++i) {
       W_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
           if (i!=j) W_[IDX(i,j,n_)] = 0.0;
       }
  }

  double Zold = this->eval_Z_(W_, R, P);

  double dZ = Znew - Zold;
  return dZ;
}
double UnitaryOptimizer_4_2::eval_Z_trial_(int I, int J, double gamma) // override
{
  double Z = 0.0; 
  for (int i = 0; i < n_; ++i) {
       W_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
           if (i!=j) W_[IDX(i,j,n_)] = 0.0;
       }
  }
  double c = cos(gamma);
  double s = sin(gamma);
  W_[IDX(I,I,n_)] = c;
  W_[IDX(J,J,n_)] = c;
  W_[IDX(I,J,n_)] = s;
  W_[IDX(J,I,n_)] =-s;

  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            for (int k = 0; k < n_; ++k) {
                Z += W_[IDX(j,i,n_)] * W_[IDX(k,i,n_)] * P_[IDX3(i,j,k)];
                for (int l = 0; l < n_; ++l) {
                     for (int m = 0; m < n_; ++m) {
                          for (int n = 0; n < n_; ++n) {
                               Z += W_[IDX(k,i,n_)] * W_[IDX(l,j,n_)] * W_[IDX(m,i,n_)] * W_[IDX(n,j,n_)] * R_[IDX6(i,j,k,l,m,n)];
                          }
                     }
                }
            }
       }
  }
  return Z;
}
void UnitaryOptimizer_4_2::form_X0_() // same
{
   for (int i = 0; i < n_; ++i) {
        Xold_[IDX(i,i,n_)] = 1.0;
        for (int j = 0; j < n_; ++j) {
             if (i != j) Xold_[IDX(i,j,n_)] = 0.0;
        }
   }
}
void UnitaryOptimizer_4_2::form_X_(int I, int J, double gamma) // same
{
  for (int i = 0; i < n_; ++i) {
       Xnew_[IDX(i,i,n_)] = 1.0;
       for (int j = 0; j < n_; ++j) {
            if (i != j) Xnew_[IDX(i,j,n_)] = 0.0;
       }
  }
  double c = cos(gamma);
  double s = sin(gamma);
  Xnew_[IDX(I,I,n_)] = c;
  Xnew_[IDX(J,J,n_)] = c;
  Xnew_[IDX(I,J,n_)] = s;
  Xnew_[IDX(J,I,n_)] =-s;
}
void UnitaryOptimizer_4_2::form_next_X_(const std::string& opt) // override (can be same if generalized ABCD to FourierN template)
{
  bool (UnitaryOptimizer_4_2::*optfunc)(double, double) = (opt == "min" ? &UnitaryOptimizer_4_2::lt_ : &UnitaryOptimizer_4_2::gt_);

  int I = 0;
  int J = 1;
  double Gamma = 0.0;
  double dZold = 1.0e+18; if (opt == "max") dZold = -1.0e+18;
  for (int j = 0; j < n_; ++j) {
       for (int i = 0; i < j; ++i) {
            Fourier9 abcd = this->get_fourier_(i, j);
            this->find_roots_boyd_(abcd);
            double gamma = this->find_gamma_(abcd, i, j, opt);
            double dZ = this->eval_dZ_(gamma, R_, P_, i, j);
            if ((this->*optfunc)(dZ, dZold)) {
                Gamma = gamma;
                I = i; J = j;
                dZold = dZ;
            }
       }
  }
  this->form_X_(I, J, Gamma);
}
Fourier9 UnitaryOptimizer_4_2::get_fourier_(int I, int J) // override
{
  // Initialize
  double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0,
                   b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;
  // P-contribution
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            for (int k = 0; k < n_; ++k) {
                 double p = P_[IDX3(i,j,k)];
                 double A =   (KroneckerDelta_(I,j)*KroneckerDelta_(J,i) - KroneckerDelta_(I,i)*KroneckerDelta_(J,j)) * (1.0-KroneckerDelta_(i,j));
                 double B =   -KroneckerDelta_(i,j)*(KroneckerDelta_(I,j)+KroneckerDelta_(J,i));
                 double C =    KroneckerDelta_(i,k)*(1.0-KroneckerDelta_(I,k))*(1.0-KroneckerDelta_(J,i));
                 double D =    KroneckerDelta_(i,k)*(KroneckerDelta_(I,k)+KroneckerDelta_(J,i));
                 double E =   -KroneckerDelta_(I,i)*KroneckerDelta_(J,k)*(1.0-KroneckerDelta_(i,k));
                 double F =    KroneckerDelta_(I,k)*KroneckerDelta_(J,i)*(1.0-KroneckerDelta_(i,k));
                 double G =    (1.0-KroneckerDelta_(i,k))*(KroneckerDelta_(I,k)*KroneckerDelta_(J,i) - KroneckerDelta_(I,i)*KroneckerDelta_(J,k));
                 double H =   -KroneckerDelta_(i,k)*(KroneckerDelta_(I,k)+KroneckerDelta_(J,i));
                 double II=    KroneckerDelta_(i,j)*(1.0-KroneckerDelta_(I,j))*(1.0-KroneckerDelta_(J,i));
                 double JJ=    KroneckerDelta_(i,j)*(KroneckerDelta_(I,j)+KroneckerDelta_(J,i));
                 double K =   -KroneckerDelta_(I,i)*KroneckerDelta_(J,j)*(1.0-KroneckerDelta_(i,j));
                 double L =    KroneckerDelta_(I,j)*KroneckerDelta_(J,i)*(1.0-KroneckerDelta_(i,j));

                 a0 += p * (A*D+G*JJ+B*(E+F)+H*(K+L)) / 2.0;
                 a1 += p * (A*C+G*II);
                 a2 += p * (A*D+G*JJ-B*(E+F)-H*(K+L)) / 2.0;
                 b1 += p * (B*C+H*II);
                 b2 += p * (H*JJ+B*D+G*(K+L)+A*(E+F)) / 2.0;
            }
       }
  }

  // R-contribution
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            for (int k = 0; k < n_; ++k) {
                 for (int l = 0; l < n_; ++l) {
                      for (int m = 0; m < n_; ++m) {
                           for (int n = 0; n < n_; ++n) {
                                double r = R_[IDX6(i,j,k,l,m,n)];

                                // First R batch
                                double A = this->a_(i,k,I,J);
                                double B =-this->b_(i,k,I,J);
                                double C = this->a_(i,m,I,J);
                                double D = this->c_(i,m,I,J);
                                double E = this->b_(i,m,I,J);
                                double F = this->a_(j,l,I,J);
                                double G = this->c_(j,l,I,J);
                                double H = this->b_(j,l,I,J);
                                double II= this->a_(j,n,I,J);
                                double JJ= this->c_(j,n,I,J);
                                double K = this->b_(j,n,I,J);
                                a0 += r * (A*C*F*K+A*C*H*II+4.0*A*D*G*K+4.0*A*D*H*JJ+A*E*F*II+4.0*A*E*G*JJ+3.0*A*E*H*K
                                          +3.0*B*C*F*II+4.0*B*C*G*JJ+B*C*H*K+4.0*B*D*F*JJ+4.0*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0;
                                a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4.0*A*D*G*JJ+3.0*A*D*H*K+3.0*A*E*G*K+3.0*A*E*H*JJ
                                          +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0;
                                a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0;
                                a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ
                                           -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0;
                                a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0;
                                b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3.0*B*C*F*JJ+3.0*B*C*G*II
                                           +3.0*B*D*F*II+4.0*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b2 += r * (A*C*F*II+2.0*A*C*G*JJ+A*C*H*K+2.0*A*D*F*JJ+2.0*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K
                                          +B*C*H*II+2.0*B*D*G*K+2.0*B*D*H*JJ+B*E*F*II+2.0*B*E*G*JJ+B*E*H*K) / 4.0;
                                b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II
                                          -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0;
                                // Second R batch
                                A = this->a_(i,m,I,J);
                                B =-this->b_(i,m,I,J);
                                C = this->a_(i,k,I,J);
                                D = this->c_(i,k,I,J);
                                E = this->b_(i,k,I,J);
                                F = this->a_(j,l,I,J);
                                G = this->c_(j,l,I,J);
                                H = this->b_(j,l,I,J);
                                II= this->a_(j,n,I,J);
                                JJ= this->c_(j,n,I,J);
                                K = this->b_(j,n,I,J);
                                a0 += r * (A*C*F*K+A*C*H*II+4.0*A*D*G*K+4.0*A*D*H*JJ+A*E*F*II+4.0*A*E*G*JJ+3.0*A*E*H*K
                                          +3.0*B*C*F*II+4.0*B*C*G*JJ+B*C*H*K+4.0*B*D*F*JJ+4.0*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0;
                                a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4.0*A*D*G*JJ+3.0*A*D*H*K+3.0*A*E*G*K+3.0*A*E*H*JJ
                                          +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0;
                                a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0;
                                a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ
                                           -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0;
                                a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0;
                                b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3.0*B*C*F*JJ+3.0*B*C*G*II
                                           +3.0*B*D*F*II+4.0*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b2 += r * (A*C*F*II+2.0*A*C*G*JJ+A*C*H*K+2.0*A*D*F*JJ+2.0*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K
                                          +B*C*H*II+2.0*B*D*G*K+2.0*B*D*H*JJ+B*E*F*II+2.0*B*E*G*JJ+B*E*H*K) / 4.0;
                                b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II
                                          -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0;
                                // Third R batch
                                A = this->a_(j,l,I,J);
                                B =-this->b_(j,l,I,J);
                                C = this->a_(i,k,I,J);
                                D = this->c_(i,k,I,J);
                                E = this->b_(i,k,I,J);
                                F = this->a_(i,m,I,J);
                                G = this->c_(i,m,I,J);
                                H = this->b_(i,m,I,J);
                                II= this->a_(j,n,I,J);
                                JJ= this->c_(j,n,I,J);
                                K = this->b_(j,n,I,J);
                                a0 += r * (A*C*F*K+A*C*H*II+4.0*A*D*G*K+4.0*A*D*H*JJ+A*E*F*II+4.0*A*E*G*JJ+3.0*A*E*H*K
                                          +3.0*B*C*F*II+4.0*B*C*G*JJ+B*C*H*K+4.0*B*D*F*JJ+4.0*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0;
                                a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4.0*A*D*G*JJ+3.0*A*D*H*K+3.0*A*E*G*K+3.0*A*E*H*JJ
                                          +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0;
                                a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0;
                                a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ
                                           -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0;
                                a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0;
                                b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3.0*B*C*F*JJ+3.0*B*C*G*II
                                           +3.0*B*D*F*II+4.0*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b2 += r * (A*C*F*II+2.0*A*C*G*JJ+A*C*H*K+2.0*A*D*F*JJ+2.0*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K
                                          +B*C*H*II+2.0*B*D*G*K+2.0*B*D*H*JJ+B*E*F*II+2.0*B*E*G*JJ+B*E*H*K) / 4.0;
                                b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II
                                          -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0;
                                // Fourth R batch
                                A = this->a_(j,n,I,J);
                                B =-this->b_(j,n,I,J);
                                C = this->a_(i,k,I,J);
                                D = this->c_(i,k,I,J);
                                E = this->b_(i,k,I,J);
                                F = this->a_(i,m,I,J);
                                G = this->c_(i,m,I,J);
                                H = this->b_(i,m,I,J);
                                II= this->a_(j,l,I,J);
                                JJ= this->c_(j,l,I,J);
                                K = this->b_(j,l,I,J);
                                a0 += r * (A*C*F*K+A*C*H*II+4.0*A*D*G*K+4.0*A*D*H*JJ+A*E*F*II+4.0*A*E*G*JJ+3.0*A*E*H*K
                                          +3.0*B*C*F*II+4.0*B*C*G*JJ+B*C*H*K+4.0*B*D*F*JJ+4.0*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0;
                                a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4.0*A*D*G*JJ+3.0*A*D*H*K+3.0*A*E*G*K+3.0*A*E*H*JJ
                                          +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0;
                                a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0;
                                a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ
                                           -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0;
                                a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0;
                                b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3.0*B*C*F*JJ+3.0*B*C*G*II
                                           +3.0*B*D*F*II+4.0*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b2 += r * (A*C*F*II+2.0*A*C*G*JJ+A*C*H*K+2.0*A*D*F*JJ+2.0*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K
                                          +B*C*H*II+2.0*B*D*G*K+2.0*B*D*H*JJ+B*E*F*II+2.0*B*E*G*JJ+B*E*H*K) / 4.0;
                                b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II
                                          -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0;
                                b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0;
                           }
                      }
                 }
            }
       }
  }

  // Save
  Fourier9 fourier = {a0, a1, a2, a3, a4, b1, b2, b3, b4};

  // Return
  return fourier;
}
void UnitaryOptimizer_4_2::find_roots_boyd_(const Fourier9& abcd) // override (can be same if generalized to FourierN)
{
    // Allocate
    
    Eigen::Matrix<std::complex<double>,8,8> B;
    Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>,8,8>> eigen;

    // Build up B matrix
    std::complex<double> c1= -(abcd.a4 + abcd.b4 * 1_i);
    std::complex<double> c2= -(abcd.a3 + abcd.b3 * 1_i);
    std::complex<double> c3= -(abcd.a2 + abcd.b2 * 1_i);
    std::complex<double> c4= -(abcd.a1 + abcd.b1 * 1_i);
    std::complex<double> c5= -(abcd.a0 * 2.0          );
    std::complex<double> c6= -(abcd.a1 - abcd.b1 * 1_i);
    std::complex<double> c7= -(abcd.a2 - abcd.b2 * 1_i);
    std::complex<double> c8= -(abcd.a3 - abcd.b3 * 1_i);

    std::complex<double> d =   abcd.a4 - abcd.b4 * 1_i ;

    B(0,0) = 0.0; B(0,1) = 1.0; B(0,2) = 0.0; B(0,3) = 0.0; B(0,4) = 0.0; B(0,5) = 0.0; B(0,6) = 0.0; B(0,7) = 0.0;
    B(1,0) = 0.0; B(1,1) = 0.0; B(1,2) = 1.0; B(1,3) = 0.0; B(1,4) = 0.0; B(1,5) = 0.0; B(1,6) = 0.0; B(1,7) = 0.0;
    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 0.0; B(2,3) = 1.0; B(2,4) = 0.0; B(2,5) = 0.0; B(2,6) = 0.0; B(2,7) = 0.0;
    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0; B(3,4) = 1.0; B(3,5) = 0.0; B(3,6) = 0.0; B(3,7) = 0.0;
    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0; B(4,4) = 0.0; B(4,5) = 1.0; B(4,6) = 0.0; B(4,7) = 0.0;
    B(5,0) = 0.0; B(5,1) = 0.0; B(5,2) = 0.0; B(5,3) = 0.0; B(5,4) = 0.0; B(5,5) = 0.0; B(5,6) = 1.0; B(5,7) = 0.0;
    B(6,0) = 0.0; B(6,1) = 0.0; B(6,2) = 0.0; B(6,3) = 0.0; B(6,4) = 0.0; B(6,5) = 0.0; B(6,6) = 0.0; B(6,7) = 1.0;
    B(7,0) =  c1/ d;
    B(7,1) =  c2/ d;
    B(7,2) =  c3/ d;
    B(7,3) =  c4/ d;
    B(7,4) =  c5/ d;
    B(7,5) =  c6/ d;
    B(7,6) =  c7/ d;
    B(7,7) =  c8/ d;
   
    // Compute all the roots
    eigen.compute(B);
    Eigen::Matrix<std::complex<double>,8,1> roots = -1_i * eigen.eigenvalues().array().log();

    // Save in work place
    for (int i = 0; i < 8; ++i) {
         //if (std::abs(roots.imag()(i)) > 1.0e-4) std::cout << " Warning! Imaginary root is non-zero = " << roots.imag()(i) << std::endl;
         double r = roots.real()(i);
         if (r < 0.0) r += M_PI;
         S_[i] = r;
    }
}
double UnitaryOptimizer_4_2::find_gamma_(const Fourier9& abcd, int I, int J, const std::string& opt) // override (can be same if adding bool for checking the Hessian or not). In this case, Fourier9 is not necessary because Hessian is not computed
{
   double gamma = 0.0;
   if (opt == "min") {
       double Zold = 1.0e+18;
       for (int i = 0; i < 8; ++i) {
            //if (this->func_1_(S_[i], abcd) > 0.0) {
                double Z = this->eval_Z_trial_(I, J, S_[i]);
                if (Z < Zold) {
                    gamma = S_[i];
                    Zold = Z;
                }
            //}
       }
   } else { 
       double Zold =-1.0e+18;
       for (int i = 0; i < 8; ++i) {
            //if (this->func_1_(S_[i], abcd) < 0.0) {
                double Z = this->eval_Z_trial_(I, J, S_[i]);
                if (Z > Zold) {
                    gamma = S_[i];
                    Zold = Z;
                }
            //}
       }
   }
   return gamma;
}
bool UnitaryOptimizer_4_2::lt_(double a, double b) // same
{
   if (a < b) return true;
   return false;
}
bool UnitaryOptimizer_4_2::gt_(double a, double b) // same
{
   if (a < b) return false;
   return true;
}
std::shared_ptr<psi::Matrix> UnitaryOptimizer_4_2::psi_X_() // same
{
  std::shared_ptr<psi::Matrix> X = std::make_shared<psi::Matrix>("MO-MO Transformation matrix X", n_, n_);
  for (int i = 0; i < n_; ++i) {
       for (int j = 0; j < n_; ++j) {
            X->set(i, j, X_[IDX(i,j,n_)]);
       }
  }
  return X;
}

} // EndNameSpace oepdev
