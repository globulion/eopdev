#ifndef _oepdev_libutil_davidson_liu_h
#define _oepdev_libutil_davidson_liu_h
/** @file davidson_liu.h */

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "gram_schmidt.h"

namespace oepdev{

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */


/**
 *  \brief Davidson-Liu diagonalization method.
 *
 * Find the lowest *M* eigenvalues and associated eigenvectors 
 * of the real, symmetric (square) matrix **H**.
 *
 * Associated options:
 *  * `DAVIDSON_LIU_NROOTS`       - number of roots of interest. Default: 1.
 *  * `DAVIDSON_LIU_CONVER`       - convergence of the iterative procedure as RMS of old and current eigenvalues. Default: 1.0E-10.
 *  * `DAVIDSON_LIU_MAXITER`      - maximum number of iterations. Default: 500.
 *  * `DAVIDSON_LIU_GUESS`        - Type of guess vectors. Default: `RANDOM`, which is constructing ranrom vectors.
 *  * `DAVIDSON_LIU_THRESH_LARGE` - Small correction vector threshold (see description below). Default: 1.0E-03. 
 *  * `DAVIDSON_LIU_THRESH_SMALL` - Small correction vector threshold (see description below). Default: 1.0E-06. 
 *  * `DAVIDSON_LIU_SPACE_MAX`    - Maximum number of guess vectors. Default: 200.     
 *  * `DAVIDSON_LIU_SPACE_START`  - Starting amount of guess vectors. Must be larger or equal to number of roots. 
 *                                  Default: -1, which means that number of roots is taken.
 *  * `DAVIDSON_LIU_STOP_WHEN_UNCONVERGED` - Raise error when iterations do not converge. Default: True.
 *
 *
 * # Usage in C++ programming
 *
 * This class is an abstract base. In order to use the Davidson-Liu method 
 * fully implemented here, one must define a child class inheriting from oepdev::DavidsonLiu
 * and implementing two of the pure methods:
 *  * `davidson_liu_compute_diagonal_hamiltonian` - method specifying
 *     the calculation of the \f$ \boldsymbol{\sigma} \f$ vectors, which are stored
 *     in the `std::vector<psi::SharedVector> sigma_vectors_davidson_liu_`;
 *  * `davidson_liu_compute_diagonal_hamiltonian` - method specifying
 *     the calculation of the diagonal elements of the Hamiltonian,
 *     stored in the `psi::SharedVector H_diag_davidson_liu_`.
 *
 * \see Examples for demo use.
 *
 * # Implementation
 * 
 * The implementation follows Figure 5, Section 3.2.1 in Ref.[1].
 * Dimensionality:
 *  * `N` - number of rows/collumns of matrix to diagonalize
 *  * `L` - current number of guess vectors
 *  * `M` - number of roots of interest
 *
 * Sigma vectors are defined to be 
 * \f[
 *   {\bf S} = {\bf H} {\bf B}
 * \f]
 * where **B** are the guess vectors stored as a matrix of size (N, L) in core memory.
 * Subspace Hamiltonain is then given by
 * \f[
 *   {\bf G} = {\bf B}^{\rm T} {\bf S}
 * \f]
 * and is diagonalized using standard diagonalization technique,
 * \f[
 *   {\bf G} = {\bf U} {\bf z} {\bf U}^{\rm T}
 * \f]
 * where **z** are the eigenvalues. First *M* lowest eigenvalues and associated eigenvectors
 * are saved in **E** and **A**, respectively (with the latter having size of (L, M)).
 * The current eigenvector matrix **C** containing roots
 * is given by
 * \f[
 *   {\bf C} = {\bf B} {\bf A}
 * \f]
 * Once this step is completed, the correction vectors are computed for each eigenvalue according to
 * \f[
 *   \delta_{Ik} = \frac{1}{E_k - H_{II}} 
 *    \left[ -E_k C_{Ik} + \sum_l^L \sigma_{Il} A_{lk} \right]
 * \f]
 * and they are orthonormalized against all the collumns of **B** by using the Gram-Schmidt procedure.
 * If the norm of such orthonormalized correction vector is larger than threshold value,
 * it is appended to **B** as new guess vector.
 *
 * \note
 *   Note that the current implementation uses the original Davidson's preconditioner, 
 *   which might have problems with breaking spin symmetry of the solution.
 * 
 * ## Treatment of correction vector threshold.
 * In the current implementation, two threshold values are defined:
 *  * larger threshold, controlled by `DAVIDSON_LIU_THRESH_LARGE` Psi4 option, 
 *    is used for the first lowest eigenvalue.
 *  * smaller threshold, controlled by `DAVIDSON_LIU_THRESH_SMALL` Psi4 option, 
 *    is used for the next eigenvalues if *M* > 1.
 * 
 * 
 * 
 * 
 * 
 * # References
 *  [1] C. David Sherrillt and Henry F. Schaefer III, *Adv. Quant. Chem.* **1999** (34), pp. 94720-1460.
 */
class DavidsonLiu {

  // --> public interface <-- //
  public:

   /// Constructor
   DavidsonLiu(psi::Options& opt);

   /// Destructor
   virtual ~DavidsonLiu();

   /// Run the Davidson-Liu solver
   virtual void run_davidson_liu();

   /// Get the eigenvalues
   psi::SharedVector eigenvalues_davidson_liu() const {return E_davidson_liu_;}
   psi::SharedVector E_davidson_liu() const {return E_davidson_liu_;}

   /// Get the eigenvectors
   psi::SharedMatrix eigenvectors_davidson_liu() const {return U_davidson_liu_;}
   psi::SharedMatrix U_davidson_liu() const {return U_davidson_liu_;}


  protected:

    /// Dimensionality of Hamiltonian                               
    int N_davidson_liu_;
                                                                    
    /// Number of guess vectors
    int L_davidson_liu_;
                                                                    
    /// Number of roots of interest
    int M_davidson_liu_;

    /// Psi4 options
    psi::Options& options_;
                                                                    
    /// Eigenvalues
    psi::SharedVector E_davidson_liu_;
                                                                    
    /// Eigenvectors
    psi::SharedMatrix U_davidson_liu_;

    /// Diagonal elements of the matrix to diagonalize
    psi::SharedVector H_diag_davidson_liu_;
                                                                    
    /// Old estimation of eigenvalues
    psi::SharedVector E_old_davidson_liu_;
                                                                    
    /// Is Davidson-Liu computer ready for calculations?
    bool davidson_liu_initialized_;
                                                                    
    /// Is Davidson-Liu computer finished with calculations?
    bool davidson_liu_finalized_;

    // Nmber of sigma vectors already computed
    int davidson_liu_n_sigma_computed_;
                                                                    
    /// Sigma vectors stored
    std::vector<psi::SharedVector> sigma_vectors_davidson_liu_;
                                                                    
    /// Object storing guess vectors
    std::shared_ptr<oepdev::GramSchmidt> guess_vectors_davidson_liu_;
                                                                    

    /// Helper interface
    virtual void davidson_liu_initialize(int N, int L, int M);
    virtual void davidson_liu_initialize_guess_vectors();
    virtual void davidson_liu_initialize_guess_vectors_by_random();
    virtual void davidson_liu_initialize_guess_vectors_by_custom();
    virtual void davidson_liu_compute_diagonal_hamiltonian() = 0;
    virtual void davidson_liu_compute_sigma() = 0;
    virtual void davidson_liu_add_guess_vectors();
    virtual double davidson_liu_compute_convergence();
    virtual void davidson_liu_finalize(bool);


};



/** @}*/

} // EndNameSpace oepdev

/** \example example_davidson_liu.cc
 *
 * This example is a trivial demo to use `oepdev::DavidsonLiu`
 * in order to diagonalize a real, symmetric matrix **H**,
 * stored in a `psi::SharedMatrix H`. This can help you to construct
 * more complicated classes that need to solve eigenpairs for very large, sparse
 * matrices such as CI Hamiltonians.
 *
 * \note This example might need compile properly (it's only a draft). Debug if necessary.
 */



#endif // _oepdev_libutil_davidson_liu_h
