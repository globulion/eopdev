#ifndef _oepdev_liboep_gefp_h_ 
#define _oepdev_liboep_gefp_h_ 
/** @file gefp.h */

#include <vector>
#include <string>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "../liboep/oep.h"
#include "../libutil/cphf.h"

namespace oepdev{

using namespace std;

/** \addtogroup OEPDEV_GEFP
 * @{
 */

/** \brief Generalized Effective Fragment. Container Class.
 *
 * Describes the GEFP fragment that is in principle designed to work
 * at correlated levels of theory.
 */
class GenEffFrag
{
  public:
   /// Initialize
   GenEffFrag();
   /// Destruct
  ~GenEffFrag();

   /// Rotate
   void rotate(std::shared_ptr<psi::Matrix> R);

   /// Translate
   void translate(std::shared_ptr<psi::Vector> T);

   /// Superimpose
   void superimpose(std::shared_ptr<psi::Matrix> targetXYZ, std::vector<int> supList);

  protected:
   // ===> Generalized Fragment Parameters <=== //

   /// Density Matrix Susceptibility Tensor
   std::shared_ptr<GenEffFrag> densityMatrixSusceptibilityGEF_;

   /// Electrostatic Energy Effective One-Electron Potential
   std::shared_ptr<GenEffFrag> electrostaticEnergyGEF_;

   /// Exchange-Repulsion Effective One-Electron Potential
   std::shared_ptr<GenEffFrag> repulsionEnergyGEF_;

   /// Charge-Transfer Effective One-Electron Potential
   std::shared_ptr<GenEffFrag> chargeTransferEnergyGEF_;

   /// EET Coupling Effective One-Electron Potential
   std::shared_ptr<GenEffFrag> EETCouplingConstantGEF_;
};

/** \brief Generalized Effective Fragment Factory. Abstract Base.
 *
 * Describes the GEFP fragment that is in principle designed to work
 * at correlated levels of theory.
 */
class GenEffFragFactory
{
  public: 
   /// Construct from wavefunction and Psi4 options
   GenEffFragFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

   /// Destruct
   virtual ~GenEffFragFactory();

   /// Compute the fragment parameters
   virtual void compute(void) = 0;

   /// Grab wavefunction
   virtual std::shared_ptr<psi::Wavefunction> wfn(void) const {return wfn_;}

   /// Grab options
   virtual psi::Options& options(void) const {return options_;}

  protected:
   /// Wavefunction
   std::shared_ptr<psi::Wavefunction> wfn_;

   /// Psi4 Options
   psi::Options& options_;


  
  private:
};

/** \brief Polarization GEFP Factory.
 *
 * Implements creation of the density matrix susceptibility tensors.
 */
class PolarGEFactory : public GenEffFragFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   PolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   PolarGEFactory(std::shared_ptr<CPHF> cphf); 

   /// Destruct
   virtual ~PolarGEFactory();

   /// Compute the density matrix susceptibility tensors
   void compute(void);

   /// Grab the density matrix susceptibility tensor for the *i*-th LMO
   std::vector<std::shared_ptr<psi::Matrix>> susceptibility(int i) const {return densityMatrixSusceptibility_[i];}

   /// Grab the density matrix susceptibility tensor *x*-th component for the *i*-th LMO
   std::shared_ptr<psi::Matrix> susceptibility(int i, int x) const {return densityMatrixSusceptibility_[i][x];}

   /// Grab the CPHF solver
   std::shared_ptr<CPHF> cphf_solver(void) const {return cphfSolver_;}

  protected:
   /// The CPHF object
   std::shared_ptr<CPHF> cphfSolver_;

   /// The density matrix susceptibility tensors
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixSusceptibility_;
};

/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_gefp_h_ 
