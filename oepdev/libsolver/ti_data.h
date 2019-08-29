#ifndef _oepdev_libutil_ti_data_h
#define _oepdev_libutil_ti_data_h
/** @file ti_data.h */

#include<cstdio>
#include<string>
#include<map>

#include "../lib3d/dmtp.h"

namespace oepdev{

using namespace std;

using SharedDMTPConvergence = std::shared_ptr<oepdev::MultipoleConvergence>;

/** \addtogroup OEPDEV_SOLVERS
 * @{
 */

/**\brief Solver of properties of molecular aggregates. Abstract base.
 *
 * Uses only a wavefunction union object to initialize.
 *
 */
class TIData 
{
 public:

   /// Constructor
   TIData();
   /// Destroctor
   virtual ~TIData();


   // <--- Computers ---> //

   /**\brief Compute overlap corrected matrix elements
    * 
    */
   virtual double coupling_direct(void);
   virtual double coupling_direct_coul(void);
   virtual double coupling_direct_exch(void);
   virtual double coupling_indirect(void);
   virtual double coupling_indirect_ti2(void);
   virtual double coupling_indirect_ti3(void);
   virtual double coupling_total(void);

   virtual double overlap_corrected(const std::string& type);
   virtual double overlap_corrected_direct(void);
   virtual double overlap_corrected_direct(double v);
   virtual double overlap_corrected_indirect(double v, double s);

 public:

   /// Overlap matrix elements between basis functions
   double s12, s13, s14, s23, s34, s24;

   /// Diagonal Hamiltonian matrix elements
   double e1, e2, e3, e4;

   /// Environmental corrections to the E1 and E2
   double de1, de2;
   bool diagonal_correction;
   bool mulliken_approximation;
   bool overlap_correction;
   bool trcamm_approximation;
   oepdev::MultipoleConvergence::ConvergenceLevel trcamm_convergence;

   /// Dictionary of all zeroth-order off-diagonal matrix elements
   std::map<std::string, double> v0;

   /// V0_Coul multipole convergence
   oepdev::SharedDMTPConvergence v0_trcamm;
};




/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_ti_data_h
