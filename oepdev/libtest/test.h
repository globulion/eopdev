#ifndef _oepdev_libtest_test_h_
#define _oepdev_libtest_test_h_
/** @file test.h */

#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/integral.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"

#include "../libpsi/integral.h"
#include "../libutil/integrals_iter.h"

namespace oepdev{
namespace test{

using namespace std;
/** \addtogroup OEPDEV_TEST
 * @{
 */

/** 
 *  \class Test
 *  \brief Manages test routines.
 */
class Test 
{
  protected: 
   /// Wavefunction object
   std::shared_ptr<psi::Wavefunction> wfn_;
   /// Psi4 Options
   psi::Options options_;

   /// Test the oepdev::ERI_2_2 class against psi::ERI
   double test_eri_2_2(void);

  public:
   Test(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options);
  ~Test();
   /// Pefrorm the test
   double run(void);
};

/** @}*/
} // EndNameSpace test
} // EndNameSpace oepdev

#endif
