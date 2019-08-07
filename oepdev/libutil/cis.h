#ifndef _oepdev_libutil_cis_h
#define _oepdev_libutil_cis_h
/** @file cis.h */

#include <string>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
//#include "psi4/libmints/integral.h"
#include "psi4/libtrans/integraltransform.h"
//#include "psi4/libdpd/dpd.h"


namespace oepdev{

using SharedMolecule           = std::shared_ptr<psi::Molecule>;
using SharedMatrix             = std::shared_ptr<psi::Matrix>;
using SharedVector             = std::shared_ptr<psi::Vector>;
using SharedIntegralTransform  = std::shared_ptr<psi::IntegralTransform>;

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief CISComputer
 *
 */
class CISComputer {

  // --> public interface <-- //
  public:

   /** \brief Build CIS Computer.
    *
    * @param type          - Type of computer
    * @param wfn           - Psi4 wavefunction
    * @param opt           - Psi4 options
    * @param reference     - Reference Slater determinant (`RHF`, `UHF` available).
    *
    * Available computer types:
    *  - `RESTRICTED` or `RCIS` - RHF wavefunction is used as reference state
    *  - `UNRESTRICTED` or `UCIS` - UHF wavefunction is used as reference state
    *
    * \note Useful options:
    *   - `CIS_NSTATES` - Number of lowest-energy excited states to include. Default: `-1` (means all states are saved)
    */
   static std::shared_ptr<CISComputer> build(const std::string& type, 
                                             std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt,
                                             const std::string& reference = "");

   /// Destructor
   virtual ~CISComputer();

   /// Solve the CIS problem
   virtual void compute(void);

   static const std::vector<std::string> reference_types;
 
  protected:
   /// Reference wavefunction
   std::shared_ptr<psi::Wavefunction> ref_wfn_;
   /// Psi4 Options
   psi::Options options_;

   /// Number of MO's
   const int nmo_;
   /// Number of alpha occupied MO's
   const int naocc_;
   /// Number of beta occupied MO's
   const int nbocc_;
   /// Number of alpha virtual MO's
   const int navir_;
   /// Number of beta virtual MO's
   const int nbvir_;
   /// Number of excited determinants
   int ndets_;

   /// CIS Excited State Hamiltonian in Slater determinantal basis
   SharedMatrix H_;
   /// CIS Coefficients \f$ U_{uI} \f$ for each excited state *I* and basis Slater determinant *u*
   SharedMatrix U_;
   /// Electronic excitation energies \f$ E_{I} \f$ wrt ground state
   SharedVector E_;

   /// Fock matrices: OO, oo, VV and vv blocks
   SharedMatrix Fa_oo_, Fb_oo_, Fa_vv_, Fb_vv_;

  // --> protected interface <-- //
  protected:

   CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

   virtual void prepare_for_cis_(void);
   virtual void build_hamiltonian_(void);
   virtual void diagonalize_hamiltonian_(void);

   virtual void set_beta_(void) = 0;

  // --> private interface <-- //
  private:
   void common_init(void);
};

class R_CISComputer: public CISComputer {
  public:
   R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~R_CISComputer(); 
  protected:
   void set_beta_(void);
};

class U_CISComputer: public CISComputer {
  public:
   U_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~U_CISComputer(); 
  protected:
   void set_beta_(void);
};


/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_cis_h
