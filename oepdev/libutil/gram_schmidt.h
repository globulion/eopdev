#ifndef _oepdev_libutil_gram_schmidt_h
#define _oepdev_libutil_gram_schmidt_h
/** @file gram_schmidt.h */

#include "psi4/libmints/vector.h"


namespace oepdev{

using SharedVector             = std::shared_ptr<psi::Vector>;

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/**
 *  \brief Gram-Schmidt orthogonalization method.
 *
 * Orthonormalize a set of *L* vectors, i.e.,
 * \f[
 *   \left\{ {\bf v}_k \right\} \rightarrow \left\{ {\bf u}_k \right\} \text{ for } k=1,2,\ldots,L
 * \f]
 *
 * # Implementation
 * 
 * The orthogonalized vectors are generated according to
 * \f[
 *    {\bf u}_k = \left[ 1 - \sum_{i=1}^{k-1} \hat{P}_{{\bf u}_i}  \right] {\bf v}_k
 * \f]
 * 
 * where the projection operator is given by
 * \f[
 *   \hat{P}_{{\bf u}} = \frac{1}{u^2} {\bf u} [\square \cdot {\bf u}]
 * \f]
 */
class GramSchmidt {

  // --> public interface <-- //
  public:


   /** \brief Construct the blank Gram-Schmidt Orthonormalizer.
    *
    */
   GramSchmidt();

   /** \brief Construct the Gram-Schmidt Orthonormalizer.
    *
    * @param vectors - list of vectors to be orthogonalized.
    *
    */
   GramSchmidt(std::vector<psi::SharedVector> vectors);

   /// Destructor
   virtual ~GramSchmidt();

   /// Retrieve all the vectors
   virtual std::vector<psi::SharedVector> V(void) const {return V_;}

   /// Retrieve the number of vectors
   virtual int L(void) const {return L_;}

   /// Retrieve the *i*th vector
   virtual psi::SharedVector V(int i) const {return V_.at(i);}

   /// Normalize all the vectors
   void normalize(void);

   /// Orthonormalize all the vectors
   void orthonormalize(void);

   /// Orthogonalize all the vectors
   void orthogonalize(void);

   /// Orthogonalize vector with respect to the vector set. Modifies **d**.
   void orthogonalize_vector(psi::SharedVector& d, bool normalize = false) const;

   /** Compute the projection vector.
    *
    * @param u - projected direction
    * @param v - projected vector
    * @returns a new vector \f$ {\bf v}' \f$ such that
    *
    *  \f[
    *  {\bf v}' = \hat{P}_{{\bf u}} {\bf v} 
    *  \f]
    */
   psi::SharedVector projection(psi::SharedVector u, psi::SharedVector v) const;

   /// Append new vector to the list
   void append(psi::SharedVector d);

   /// Reset by providing new vectors
   void reset(std::vector<psi::SharedVector> V);

   /// Reset to empty state
   void reset(void);


  protected:

  /// Vectors stored
  std::vector<psi::SharedVector> V_;

  /// Number of vectors
  int L_;

};



/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_gram_schmidt_h
