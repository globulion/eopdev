#include <cmath>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include "unitary_optimizer.h"
#include "psi4/psi4-dec.h"

using namespace oepdev;

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
void UnitaryOptimizer::update_conv_()
{
  conv_current_ = std::abs(Znew_-Zold_);
}
void UnitaryOptimizer::update_iter_()
{
  for (int i = 0; i < n_; ++i) Xold_[i] = Xnew_[i];
  Zold_ = Znew_;
  niter_ += 1;
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
void UnitaryOptimizer::form_X0_()
{
   for (int i = 0; i < n_; ++i) {
        Xold_[IDX(i,i,n_)] = 1.0;
        for (int j = 0; j < n_; ++j) {
             if (i != j) Xold_[IDX(i,j,n_)] = 0.0;
        }
   }
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

