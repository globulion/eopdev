#ifndef _oepdev_libutil_kabsch_superimposer_h
#define _oepdev_libutil_kabsch_superimposer_h
/** @file kabsch_superimposer.h */

#include <string>

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"


namespace oepdev{

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief Compute the Cartesian rotation matrix between two structures. 
 *
 *  The superimposition is defined as:
 *  \f[
 *    {\bf X}' = {\bf t} + {\bf X} \cdot {\bf r} \approx {\bf X}_0
 *  \f]
 *  where \f$ X_{iu} \f$ is the *u*-th Cartesian component of the *i*-th atom's position,
 *  \f$ {\bf t} \f$ is the superimposition translation vector, \f$ {\bf r} \f$ 
 *  is the superimposition rotation matrix, and prime denotes transformed coordinates.
 *
 *  The superimposition uses the Kabsch algorithm.
 *  
 *  # The Kabsch Algorithm.
 *  
 *  Rotation matrix is calculated from
 *  \f[
 *   {\bf r} = {\bf U} \cdot {\bf V}^{\rm T}
 *  \f]
 *  where 
 *  \f[
 *   {\bf A} = {\bf U} \cdot {\bf S} \cdot {\bf V}^{\rm T}
 *  \f]
 *  is the singular value decomposition of the covariance matrix
 *  \f[
 *    {\bf A} = \left[ {\bf X} - \left< {\bf X} \right> \right]^{\rm T} \cdot \left[ {\bf X}_0 - \left< {\bf X}_0 \right> \right]
 *  \f]
 *  The average of position is given by
 *  \f[
 *    \left< {\bf X} \right>_{u} = \frac{1}{N} \sum_i X_{iu}
 *  \f]
 *  where *N* is the number of atoms.
 *  If determinant of rotation matrix is negative (indicating inversion), 
 *  rotation matrix is recomputed by inverting the sign of the third column of \f$ {\bf V} \f$.
 *  
 *  The translation vector is then calculated by
 *  \f[
 *   {\bf t} = \left< {\bf X}_0 \right> - \left< {\bf X} \right> \cdot {\bf r}
 *  \f]
 */
class KabschSuperimposer
{
  public:
   /// Constructor
   KabschSuperimposer() : rotation(std::make_shared<psi::Matrix>("",3,3)), 
                          translation(std::make_shared<psi::Vector>("",3)),
                          initial_xyz(nullptr), final_xyz(nullptr) {rotation->identity();translation->zero();};
   /// Destructor
  ~KabschSuperimposer() {};
   /** \brief Run the Kabsch algorithm.
    *  
    *  @param initial_xyz - position vectors \f$ {\bf X} \f$
    *  @param final_xyz   - position vectors \f$ {\bf X}_0 \f$
    */
   void compute(psi::SharedMatrix initial_xyz, psi::SharedMatrix final_xyz);

   /** \brief Run the Kabsch algorithm.
    *  
    *  @param initial_mol - molecule with atomic positions at \f$ {\bf X} \f$
    *  @param final_mol   - molecule with atomic positions at \f$ {\bf X}_0 \f$
    */
   void compute(psi::SharedMolecule initial_mol, psi::SharedMolecule final_mol) {
        psi::SharedMatrix initial_xyz = std::make_shared<psi::Matrix>(initial_mol->geometry());
        psi::SharedMatrix final_xyz = std::make_shared<psi::Matrix>(final_mol->geometry());
        compute(initial_xyz, final_xyz);}
   /// Rotation matrix \f$ {\bf r} \f$
   psi::SharedMatrix rotation;
   /// Translation vector \f$ {\bf t} \f$
   psi::SharedVector translation;
   /// Initial xyz \f$ {\bf X} \f$
   psi::SharedMatrix initial_xyz;
   /// Final xyz \f$ {\bf X}_0 \f$
   psi::SharedMatrix final_xyz;
   /// Return transformed coordinates \f$ {\bf X}' \f$
   psi::SharedMatrix get_transformed(void);
   /// Compute RMS or superimposition
   double rms(void);
   /// Clear all previous calculations
   void clear(void);
};




/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_kabsch_superimposer_h
