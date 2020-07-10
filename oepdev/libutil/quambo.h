#ifndef _oepdev_libutil_quambo_h
#define _oepdev_libutil_quambo_h
/** @file quambo.h */

#include <string>

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"


namespace oepdev{

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief The Quasiatomic Minimal Basis Set Molecular Orbitals (QUAMBO)
 *
 *  
 *  # Calculation Algorithm.
 *  TODO
 */
class QUAMBO
{
  public:
   /// Constructor
   QUAMBO();
   /// Destructor
  ~QUAMBO();

   void compute(void);

  protected:

};




/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_quambo_h
