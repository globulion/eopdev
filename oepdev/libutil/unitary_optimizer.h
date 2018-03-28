#ifndef _oepdev_libutil_unitary_optimizer_h
#define _oepdev_libutil_unitary_optimizer_h
/** @file unitary_optimizer.h */
#include <string>

namespace oepdev {

using namespace std;

#define IDX(i,j,n) ((n)*(i)+(j))
#define SDIM 100

/** \addtogroup OEPDEV_UTILITIES
 * @{
 */

struct ABCD 
{
  double A;
  double B;
  double C;
  double D;
};

class UnitaryOptimizer
{
   public:
     /// Create from R and P matrices and optimization options
     UnitaryOptimizer(double* R, double* P, int n, double conv = 1.0e-6, int maxiter = 100, bool verbose = true);

     /// Clear memory
    ~UnitaryOptimizer();

     /// Run the minimization
     bool maximize();

     /// Run the maximization
     bool minimize();

     /// Get the unitary matrix (solution)
     double* X() const {return X_;}

     /// Get the actual value of Z function
     double  Z() {return this->eval_Z_();}

     /// Get the status of the optimization
     bool success() const {return success_;}

   protected:
     /// Dimension of the problem
     const int n_;
     /// Convergence
     const double conv_;
     /// Maximum number of iterations
     const int maxiter_;
     /// Verbose mode
     const bool verbose_;

     /// R matrix
     double*  R_;
     /// P vector
     double*  P_;
     /// Reference R matrix
     double*  R0_;
     /// Reference P vector
     double*  P0_;
     /// X Matrix (accumulated solution)
     double*  X_;
     /// Work place
     double*  W_;
     /// Temporary X matrix
     double*  Xold_;
     /// Temporary X matrix
     double*  Xnew_;
     /// Current number of iterations
     int niter_;
     /// Current solutions
     double S_[4];
     /// Initial Z value
     double Zinit_;
     /// Old Z value
     double Zold_;
     /// New Z value
     double Znew_;
     /// Current convergence
     double conv_current_;
     /// Status of optimization
     bool success_;

     /// Prepare the optimizer
     void common_init_();
     /// Run the optimization (intermediate interface)
     void run_(const std::string& opt);
     /// Run the optimization (inner interface)
     void optimize_(const std::string& opt);
     /// Restore the initial state of the optimizer
     void refresh_();

     /// Update the convergence
     void update_conv_();
     /// Update the iterates
     void update_iter_();
     /// Update Z value
     void update_Z_();
     /// Uptade R and P matrices
     void update_RP_();
     /// Update the solution matrix X
     void update_X_();

     /// Evaluate the objective Z function
     double eval_Z_(double* X, double* R, double* P);
     double eval_Z_();
     /// Evaluate the change in Z
     double eval_dZ_(double g, double* R, double* P, int i, int j);
     /// Evaluate the trial Z value
     double eval_Z_trial_(int i, int j, double gamma);

     /// Create identity matrix
     void form_X0_();
     /// Form unitary matrix X (store in buffer Xnew_)
     void form_X_(int i, int j, double gamma);
     /// Form the next unitary matrix X
     void form_next_X_(const std::string& opt);

     /// Retrieve ABCD parameters for root search
     ABCD get_ABCD_(int i, int j);
     /// Solve for all roots of equation A*sin(g) + B*cos(g) + C*sin(2*g) + D*cos(2*g) = 0 -> implements Boyd's method
     void find_roots_boyd_(const ABCD& abcd);
     /// Solve for root of equation A*sin(g) + B*cos(g) + C*sin(2*g) + D*cos(2*g) = 0 -> implements Halley's method
     double find_root_halley_(double x0, const ABCD& abcd);
     /// Compute gamma from roots of base equations
     double find_gamma_(const ABCD& abcd, int i, int j, const std::string& opt);

     /// less-than function
     bool lt_(double a, double b);
     /// greater-than function
     bool gt_(double a, double b);

     /// Function f(gamma) = d(dZ)/dgamma
     inline double func_0_(double g, const ABCD& abcd);
     /// Gradient of f(gamma)
     inline double func_1_(double g, const ABCD& abcd);
     /// Hessian of f(gamma)
     inline double func_2_(double g, const ABCD& abcd);

};

/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_unitary_optimizer_h

