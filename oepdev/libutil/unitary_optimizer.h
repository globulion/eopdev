#ifndef _oepdev_libutil_unitary_optimizer_h
#define _oepdev_libutil_unitary_optimizer_h
/** @file unitary_optimizer.h */
#include <string>
#include <complex>
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"

namespace oepdev {

using namespace std;

#define IDX(i,j,n) ((n)*(i)+(j))

constexpr std::complex<double> operator""_i(unsigned long long d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}
constexpr std::complex<double> operator""_i(long double d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}

/** \addtogroup OEPDEV_UTILITIES
 * @{
 */

/**\brief Simple structure to hold the Fourier series expansion coefficients.
 *
 */
struct ABCD 
{
  double A;
  double B;
  double C;
  double D;
};

/**\brief Find the optimim unitary matrix of quadratic matrix equation.
 *
 * The objective function of the orthogonal matrix \f$ {\bf X} \f$ 
 * \f[
 *   Z({\bf X}) \equiv \sum_{ijkl} X_{ij} X_{kl} R_{jl} - \sum_{ij} X_{ij} P_j
 * \f]
 * is optimized by using the Jacobi iteration algorithm.
 * In the above equation, \f$ {\bf R} \f$ is a square, general real matrix 
 * of size \f$ N\times N\f$ whereas \f$ {\bf P} \f$ is a real vector of length
 * \f$ N \f$.
 *
 * ## Algorithm.
 * Optimization of \f$ {\bf X} \f$ is factorized into a sequence of 
 * 2-dimensional rotations with one real parameter \f$ \gamma \f$:
 * \f[
 *   {\bf X}^{\rm New} = {\bf U}(\gamma) \cdot {\bf X}^{\rm Old}
 * \f]
 * where
 * \f[
 * {\bf U}(\gamma) \equiv
 * \begin{pmatrix}
 * \ddots  &  &  & & \\ 
 *  & \cos(\gamma) & \cdots  & \sin(\gamma)& \\ 
 *  & \vdots  & \ddots & \vdots & \\ 
 *  & -\sin(\gamma) & \cdots  & \cos(\gamma) & \\
 *  &  &   &  &  & \ddots 
 * \end{pmatrix}
 * \f]
 * is the Jacobi transformation matrix constructed for the \f$ I\f$th and \f$ J\f$th
 * element from the entire \f$ N\f$-dimensional set.
 * For the sake of algirithmic simplicity,
 * every iteration after \f$ {\bf U}(\gamma) \f$ has been formed, \f$ {\bf X}^{\rm Old} \f$ is for a while assumed to be an 
 * identity matrix and
 * the \f$ {\bf R} \f$ matrix and \f$ {\bf P} \f$ vector
 * are transformed according to the following formulae
 * \f{align*}{
 *   {\bf R} &\rightarrow {\bf U} {\bf R} {\bf U}^T \\
 *   {\bf P} &\rightarrow {\bf U} {\bf P}
 * \f}
 * The full transformation matrix is accumulated in the memory buffer until convergence.
 *
 * In each iteration, the optimum angle \f$ \gamma \f$ is found as follows: First, the roots of
 * the finite Fourier series 
 * \f[
 *  A \sin(\gamma) + B \cos(\gamma) + C \sin(2\gamma) + D \cos(2\gamma) = 0
 * \f]
 * are found. In the above equations, the expansion coefficients are given as
 * \f{align*}{
 *   A &= P_I + P_J - \sum_{k \neq I,J} \left( R_{Ik} + R_{Jk} + R_{kI} + R_{kJ} \right) \\
 *   B &= P_I - P_J - \sum_{k \neq I,J} \left( R_{Ik} - R_{Jk} + R_{kI} - R_{kJ} \right) \\
 *   C &= -2(R_{IJ} + R_{JI}) \\
 *   D &= -2(R_{II} - R_{JJ})
 * \f}
 * and \f$ I,J \f$ are the chosen indices in the Jacobi iteration subspace.
 * The roots are evaluated by applying the Boyd's method[1], in which they are given as 
 * \f[
 *   \gamma_n = \Re\left[-i\ln(\lambda_n)\right]
 * \f]
 * where \f$ \lambda_n \f$ is an eivenvalue of the following 4 by 4 complex matrix:
 * \f[
 *  \begin{pmatrix}
 * 0 & 1 & 0 & 0 \\ 
 * 0 & 0 & 1 & 0 \\ 
 * 0 & 0 & 0 & 1 \\ 
 * -\frac{D+iC}{D-iC} &
 * -\frac{B+iC}{A-iC} &
 * 0 & 
 * -\frac{B-iC}{A-iC}
 * \end{pmatrix}
 * \f]
 * Once the four roots of the Fourier series equation are found, one solution out of four
 * is chosen which satisfies the global optimum condition, i.e., the largest increase/decrease
 * in the objective function given by
 * \f[
 *  \delta Z = A(1-\cos(\gamma)) + B\sin(\gamma) + C\sin^2(\gamma) + \frac{D}{2}\sin(2\gamma)
 * \f]
 * The discrimination between the minimae/maximae is performed based on the evaluation of the Hessian
 * of \f$ Z \f$ with respect to \f$ \gamma \f$,
 * \f[
 *  \frac{\partial^2 Z}{\partial \gamma^2} = A \cos(\gamma) - B \sin(\gamma) + 2C \cos(2\gamma) - 2D \sin(2\gamma)
 * \f]
 * All the \f$ N(N-1)/2 \f$ unique pairs of molecular orbitals are checked
 * and the optimal set of \f$ \gamma, I, J\f$ is chosen to construct \f$ {\bf X}^{\rm New} \f$.
 *  
 * ## References:
 * [1] Boyd, J.P.; J. Eng. Math. (2006) 56, pp. 203-219
 */
class UnitaryOptimizer
{
   public:
     /**\brief Create from R and P matrices and optimization options
      * @param R - \f$ {\bf R} \f$ matrix
      * @param P - \f$ {\bf P} \f$ vector
      * @param n - dimensionality of the problem (\f$ N \f$)
      * @param conv - convergence in the \f$ Z \f$ function
      * @param maxiter - maximum number of iterations
      * @param verbose - whether print information of iteration process or not
      * Sets up the optimizer.
      */
     UnitaryOptimizer(double* R, double* P, int n, double conv = 1.0e-6, int maxiter = 100, bool verbose = true);

     /**\brief Create from R and P matrices and optimization options
      * @param R - \f$ {\bf R} \f$ matrix
      * @param P - \f$ {\bf P} \f$ vector
      * @param conv - convergence in the \f$ Z \f$ function
      * @param maxiter - maximum number of iterations
      * @param verbose - whether print information of iteration process or not
      * Sets up the optimizer.
      */
     UnitaryOptimizer(std::shared_ptr<psi::Matrix> R, std::shared_ptr<psi::Vector> P, double conv = 1.0e-6, int maxiter = 100, bool verbose = true);

     /// Clear memory
    ~UnitaryOptimizer();

     /// Run the minimization
     bool maximize();

     /// Run the maximization
     bool minimize();

     /// Get the unitary matrix (solution)
     std::shared_ptr<psi::Matrix> X() {return this->psi_X_();}

     /// Get the unitary matrix (pointer to solution)
     double* get_X() const {return this->X_;}

     /// Get the actual value of Z function
     double  Z() {return this->eval_Z_();}

     /// Get the status of the optimization
     bool success() const {return success_;}

   protected:
     /**\brief Initialize the basic memory
      * @param n - dimensionality of the problem (\f$ N \f$)
      * @param conv - convergence in the \f$ Z \f$ function
      * @param maxiter - maximum number of iterations
      * @param verbose - whether print information of iteration process or not
      * Sets up the optimizer.
      */
     UnitaryOptimizer(int n, double conv, int maxiter, bool verbose);

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
     /// Hessian of f(gamma) - used only for Halley method (not implemented since Boyd method is more suitable here)
     inline double func_2_(double g, const ABCD& abcd);

     /// Form the Psi4 matrix with the transformation matrix
     std::shared_ptr<psi::Matrix> psi_X_();

};

/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_unitary_optimizer_h

