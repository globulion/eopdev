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
/** \addtogroup OEPDEV_TESTS
 * @{
 */

/** 
 *  \class Test
 *  \brief Manages test routines.
 */
class Test 
{
  public:
   /// Construct the tester
   Test(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options);

   /// Destructor
  ~Test();

   /// Pefrorm the test
   double run(void);

  protected: 
   /// Wavefunction object
   std::shared_ptr<psi::Wavefunction> wfn_;

   /// Psi4 Options
   psi::Options& options_;

   // ---> Tests <--- //

   /// Test the basic functionalities of OEPDev
   double test_basic(void);

   /// Test the CPHF method
   double test_cphf(void);

   /// Test the density matrix susceptibility (X = 1)
   double test_dmatPol(void);

   /// Test the density matrix susceptibility
   double test_dmatPolX(void);

   /// Test the oepdev::ERI_1_1 class against psi::ERI
   double test_eri_1_1(void);

   /// Test the oepdev::ERI_2_2 class against psi::ERI
   double test_eri_2_2(void);

   /// Test the oepdev::ERI_3_1 class against psi::ERI
   double test_eri_3_1(void);

   /// Test the oepdev::UnitaryOptimizer class
   double test_unitaryOptimizer(void);

   /// Test the oepdev::UnitaryOptimizer_4_2 class
   double test_unitaryOptimizer_4_2(void);

   /// Test the oepdev::RHFPerturbed class
   double test_scf_perturb(void);

   /// Test the oepdev::CAMM class
   double test_camm(void);

   /// Test the oepdev::DMTP class for energy calculations
   double test_dmtp_energy(void);

   /// Test the custom code
   double test_custom(void);

};

/** @}*/
} // EndNameSpace test
} // EndNameSpace oepdev

#endif
