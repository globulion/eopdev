#ifndef _oepdev_libtest_test_h_
#define _oepdev_libtest_test_h_
/** @file test.h */

#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"

#include "../libpsi/integral.h"
#include "../libutil/integrals_iter.h"

#include "../../include/oepdev_files.h"

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

   /// Test the AO basis set rotation from oepdev::ao_rotation_matrix
   double test_basis_rotation(void);

   /// Test the CIS(RHF) method
   double test_cis_rhf(void);

   /// Test the CIS(UHF) method
   double test_cis_uhf(void);

   /// Test the CIS(RHF) method with Davidson-Liu algorithm
   double test_cis_rhf_dl(void);

   /// Test the CIS(UHF) method with Davidson-Liu algorithm
   double test_cis_uhf_dl(void);

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

   /// Test the oepdev::UnitaryOptimizer_2 class
   double test_unitaryOptimizer_2(void);

   /// Test the oepdev::UnitaryOptimizer_4_2 class
   double test_unitaryOptimizer_4_2(void);

   /// Test the oepdev::RHFPerturbed class
   double test_scf_perturb(void);

   /// Test the oepdev::QUAMBO class
   double test_quambo(void);

   /// Test the oepdev::CAMM class
   double test_camm(void);

   /// Test the oepdev::MultipoleConvergence class: potential and field calculations
   double test_dmtp_pot_field(void);

   /// Test the oepdev::DMTP class for energy calculations
   double test_dmtp_energy(void);

   /// Test the oepdev::EFP2_GenEffPar and oepdev::EFP2_Computer classes
   double test_efp2_energy(void);

   /// Test the oepdev::EFP2_GenEffPar and oepdev::EFP2_Computer classes
   double test_oep_efp2_energy(void);

   /// Test the oepdev::KabschSuperimposer
   double test_kabsch_superimposition(void);

   /// Test the oepdev::DMTP class for superimposition
   double test_dmtp_superimposition(void);

   /// Test the oepdev::ESPSolver
   double test_esp_solver(void);

   /// Test the cube file generation (oepdev::Field3D electrostatic potential and oepdev::Points3DIterator for cube collection)
   double test_points_collection3d(void);
   
   ///Test the Charge-transfer Energy Solver (benchmark method Otto-Ladik)
   double test_ct_energy_benchmark_ol(void);

   ///Test the Charge-transfer Energy Solver (oep-based method Otto-Ladik)
   double test_ct_energy_oep_based_ol(void);

   ///Test the Repulsion Energy Solver: (benchmark method Hayes-Stone)
   double test_rep_energy_benchmark_hs(void);

   ///Test the Repulsion Energy Solver: (benchmark method Density-Based - DDS/HF)
   double test_rep_energy_benchmark_dds(void);

   ///Test the Repulsion Energy Solver: (benchmark method Murrell-etal)
   double test_rep_energy_benchmark_murrell_etal(void);

   ///Test the Repulsion Energy Solver: (OEP-based method Murrell-etal)
   double test_rep_energy_oep_based_murrell_etal(void);

   ///Test the Repulsion Energy Solver: (benchmark method Otto-Ladik)
   double test_rep_energy_benchmark_ol(void);

   ///Test the Repulsion Energy Solver: (benchmark method EFP2)
   double test_rep_energy_benchmark_efp2(void);

   /// Test the custom code (to be deprecated)
   double test_custom(void);

};

/** @}*/
} // EndNameSpace test
} // EndNameSpace oepdev

#endif
