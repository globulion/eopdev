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



/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_basis_rotation_h
