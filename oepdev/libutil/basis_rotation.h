#ifndef _oepdev_libutil_basis_rotation_h
#define _oepdev_libutil_basis_rotation_h
/** @file basis_rotation.h */

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"

namespace psi{
  using SharedBasisSet = std::shared_ptr<BasisSet>;
}

namespace oepdev{


/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \section rot Theory
 *
 *   The objective is to find the formulae for rotation matrices
 *   of the AO spaces as functions of the Cartesian 3 x 3 rotation matrices.
 *   It is obvious that *p*-type functions transform as a usual Cartesian vectors.
 *   However, higher angular momentum functions transform in a more complex way.
 *
 *   ## Problem
 *   
 *   Define a vectorized AO space *M* of rank *r*>1 that is constructed from unique
 *   tensor components of fully symmetric *r*-th rank AO tensor populated in standard order,
 *   \f[
 *     M_{\{ab\ldots k\}} = M_{ab\ldots k} \quad {\rm for} \quad a\leq b\leq \ldots \leq k
 *   \f]
 *   Given a general rotation of Cartesian tensors
 *   \f[
 *     M_{ab\ldots k} = \sum_{a'b'\ldots k'} M_{a'b'\ldots k'} r_{a'a} r_{b'b} \cdots r_{k'k}
 *   \f]
 *   find closed expressions for the rotation matrix in reduced composite AO space
 *   obeying
 *   \f[
 *     M_{[ab\ldots k]} = \sum_{\{a'b'\ldots k'\}} M_{\{a'b'\ldots k'\}} R_{\{a'b'\ldots k'\},[ab\ldots k]}
 *   \f]
 *   In the derivations below the following identity of first-order partitioning will be of use:
 *   \f[
 *     \sum_{ab} M_{ab} \hat{s}_{ab} = \sum_{\{ab\}} M_{\{ab\}} 
 *          \left( \hat{s}_{ab} + \Delta_{ab} \hat{s}_{ba} \right)
 *   \f]
 *   where
 *   \f[
 *      \Delta_{ab} \equiv 1 - \delta_{ab}
 *   \f]
 *   and the operator *s* of rank *r* acts as follows
 *   \f[
 *    s^{ab\ldots k}_{a'b'\ldots k'} \equiv \hat{s}_{a'b'\ldots k'}
 *    \underbrace{ {\bf r} \otimes {\bf r} \otimes \cdots \otimes {\bf r} }_{r}
 *         = r_{a'a} r_{b'b} \cdots r_{k'k}
 *   \f]
 *
 *   ## Rotation of 6D functions
 *
 *   The rotation of the 6D vector full tensor AO space of rank 2 and dimensions (3,3) 
 *   is given by
 *   \f[
 *      M_{ab} = \sum_{a'b'} M_{a'b'} r_{a'b} r_{b'b}
 *   \f]
 *   Applying the identity of first-order partitioning directly leads to
 *   the formula for a reduced 6D tensor rotation of rank 1 and dimension (6),
 *   \f[
 *      M_{[ab]} = \sum_{\{a'b'\}} M_{\{a'b'\}} R_{\{a'b'\},[ab]} 
 *   \f]
 *   where the 6 x 6 rotation matrix is given by
 *   \f[
 *      R_{\{a'b'\},[ab]} = r_{a'a} r_{b'b} + \Delta_{a'b'} r_{b'a} r_{a'b} 
 *   \f]
 *   
 *
 *   ## Rotation of 10F functions
 *
 *   The rotation of the 10D vector full tensor AO space of rank 3 and dimensions (3,3,3) 
 *   is given by
 *   \f[
 *      M_{abc} = \sum_{a'b'c'} M_{a'b'c'} r_{a'b} r_{b'b} r_{c'c}
 *   \f]
 *   First of all, notice that one can perform the following partitioning
 *   \f[
 *     \sum_a \sum_{b\neq a} \sum_{c\neq b \neq a} M_{abc} \hat{s}_{abc} = 
 *     \sum_{\{abc\}} M_{\{abc\}} \left( 
 *       \hat{s}_{abc} + \hat{s}_{acb} + \hat{s}_{bac} + \hat{s}_{bca} + \hat{s}_{cab} + \hat{s}_{cba}
 *     \right)
 *   \f]
 *   Then, perform a partitioning of the triple sum,
 *   \f{multline*}{
 *    \sum_{abc} M_{abc} \hat{s}_{abc} = 
 *                 \sum_a \sum_{b\neq a} \sum_{c\neq b \neq a} M_{abc} \hat{s}_{abc} \\
 *               + \sum_a \sum_{b\geq a} M_{abb} \hat{s}_{abb} + \sum_a \sum_{b < a} M_{abb} \hat{s}_{abb} \\
 *               + \sum_a \sum_{b > a  } M_{aba} \hat{s}_{aba} + \sum_a \sum_{b < a} M_{aba} \hat{s}_{aba} \\
 *               + \sum_a \sum_{b > a  } M_{bba} \hat{s}_{bba} + \sum_a \sum_{b < a} M_{bba} \hat{s}_{bba} \\
 *   \f}
 *   Using the first-order partitioning theorem and interchanging the dummy indices
 *   one finds that
 *   \f[
 *      M_{[abc]} = \sum_{\{a'b'c'\}} M_{\{a'b'c'\}} R_{\{a'b'c'\},[abc]} 
 *   \f]
 *   where the 10 x 10 rotation matrix is given by
 *   \f{multline*}{
 *      R_{\{a'b'c'\},[abc]} = \delta_{b'c'} 
 *            \left( s^{abc}_{a'b'b'} + \Delta_{a'b'} \left\{ s^{abc}_{b'a'b'} + s^{abc}_{b'b'a'} \right\} \right) \\
 *                           + \delta_{a'b'} \Delta_{b'c'}
 *            \left( s^{abc}_{c'a'a'} + s^{abc}_{a'c'a'} + s^{abc}_{a'a'c'}  \right)  \\
 *          +\Delta_{a'b'}\Delta_{b'c'} 
 *          \left( s^{abc}_{'a'b'c} + s^{abc}_{a'c'b'} + s^{abc}_{b'a'c'} + 
                   s^{abc}_{b'c'a'} + s^{abc}_{c'a'b'} + s^{abc}_{c'b'a'} \right)
 *   \f}
 *   where
 *   \f[
 *    s^{abc}_{a'b'c'} \equiv \hat{s}_{a'b'c'} {\bf r} \otimes {\bf r} \otimes {\bf r} = r_{a'a} r_{b'b} r_{c'c}
 *   \f]
 *
 *
 */

/** \name Rotation of AO Space */

//@{


/** \brief Compute the 6 x 6 rotation matrix of the 6D orbitals. 
 * 
 *  @param r - Cartesian 3 x 3 rotation matrix
 *  @return 6 x 6 rotation matrix of the 6D orbitals
 */
psi::SharedMatrix r6(psi::SharedMatrix r);

/** \brief Compute the 10 x 10 rotation matrix of the 10F orbitals. 
 * 
 *  @param r - Cartesian 3 x 3 rotation matrix
 *  @return 10 x 10 rotation matrix of the 10F orbitals
 */
psi::SharedMatrix r6(psi::SharedMatrix r);


void populate(double** R, double** r, std::vector<int> idx_am, const int& nam);



/** \brief Compute the full rotation matrix of AO orbital space. 
 * 
 *  @param r - Cartesian 3 x 3 rotation matrix
 *  @param b - Basis set 
 */
psi::SharedMatrix ao_rotation_matrix(psi::SharedMatrix r, psi::SharedBasisSet b);

//@}

/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_basis_rotation_h
