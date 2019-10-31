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
 * # Usage
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
 * # Implementation
 * 
 * Th
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
    virtual void davidson_liu_finalize();


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
