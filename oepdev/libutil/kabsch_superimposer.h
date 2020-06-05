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
 *  Uses the Kabsch algorithm.
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
   /// Run the superimposition
   void compute(psi::SharedMatrix initial_xyz, psi::SharedMatrix final_xyz);
   /// Rotation matrix
   psi::SharedMatrix rotation;
   /// Translation vector
   psi::SharedVector translation;
   /// Initial xyz
   psi::SharedMatrix initial_xyz;
   /// Final xyz
   psi::SharedMatrix final_xyz;
   /// Return transformed coordinates
   psi::SharedMatrix get_transformed(void);
   /// Compute RMS or superimposition
   double rms(void);
   /// Clear all previous calculations
   void clear(void);
};




/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_kabsch_superimposer_h
