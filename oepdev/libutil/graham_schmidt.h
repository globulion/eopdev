#ifndef _oepdev_libutil_graham_schmidt_h
#define _oepdev_libutil_graham_schmidt_h
/** @file graham_schmidt.h */

#include "psi4/libmints/vector.h"


namespace oepdev{

using SharedVector             = std::shared_ptr<psi::Vector>;

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/**
 *  \brief Graham-Schmidt orthogonalization method.
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
class GrahamSchmidt {

  // --> public interface <-- //
  public:

   /** \brief Construct the Graham-Schmidt Orthonormalizer.
    *
    * @param vectors - list of vectors to be orthogonalized.
    *
    */
   GrahamSchmidt(std::vector<psi::SharedVector> vectors);
   /// Destructor
   virtual ~GrahamSchmidt();


   /// Normalize all the vectors
   void normalize(void);

   /// Orthonormalize all the vectors
   void orthonormalize(void);

   /// Orthogonalize all the vectors
   void orthogonalize(void);

   /// Orthogonalize vector with respect to the vector set
   psi::SharedVector orthogonalize_vector(psi::SharedVector d, bool normalize = false);

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
   psi::SharedVector projection(psi::SharedVector u, psi::SharedVector v);

   /// Append new vector to the list
   void append(psi::SharedVector d);

   /// Reset by providing other vectors
   void reset(std::vector<psi::SharedVector> V);


  protected:

  /// Vectors stored
  std::vector<psi::SharedVector> V_;

  /// Number of vectors
  int L_;

};



/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_graham_schmidt_h
