#ifndef _oepdev_liboep_gefp_h_ 
#define _oepdev_liboep_gefp_h_ 
/** @file gefp.h */

#include <vector>
#include <string>
#include <random>
#include <cmath>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector3.h"
#include "../liboep/oep.h"
#include "../libutil/cphf.h"

namespace oepdev{

using namespace std;

/** \addtogroup OEPDEV_GEFP
 * @{
 */

/** \brief Generalized Effective Fragment Parameters. Container Class.
 *
 */
class GenEffPar
{
  public:
   /// Create with name of this parameter type
   GenEffPar(std::string name) : name_(name) {};
   /// Destruct
  ~GenEffPar() {};

   /// Set The Density Matrix Susceptibility Tensors  
   void set_susceptibility(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixSusceptibility_=susc;};

   /// Grab the density matrix susceptibility tensor for the *i*-th LMO
   std::vector<std::shared_ptr<psi::Matrix>> susceptibility(int i) const {return densityMatrixSusceptibility_[i];}

   /// Grab the density matrix susceptibility tensor *x*-th component for the *i*-th LMO
   std::shared_ptr<psi::Matrix> susceptibility(int i, int x) const {return densityMatrixSusceptibility_[i][x];}


  protected:
   /// The Name of Parameter Type
   std::string name_;

   /// The Density Matrix Susceptibility Tensors
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixSusceptibility_;
};

/** \brief Generalized Effective Fragment. Container Class.
 *
 * Describes the GEFP fragment that is in principle designed to work
 * at correlated levels of theory.
 */
class GenEffFrag
{
  protected: 
   /// Name of GEFP
   std::string name_;

  public:
   /// Initialize with default name of GEFP (Default)
   GenEffFrag();
   /// Initialize with custom name of GEFP
   GenEffFrag(std::string name);
   /// Destruct
  ~GenEffFrag();

   /// Dictionary of All GEF Parameters
   std::map<std::string, std::shared_ptr<GenEffPar>> parameters;


   // ---> Mutators <--- //

   /// Rotate
   void rotate(std::shared_ptr<psi::Matrix> R);

   /// Translate
   void translate(std::shared_ptr<psi::Vector> T);

   /// Superimpose
   void superimpose(std::shared_ptr<psi::Matrix> targetXYZ, std::vector<int> supList);

   /// Set the Density Matrix Susceptibility Tensors
   void set_gefp_polarization(const std::shared_ptr<GenEffPar>& par) {densityMatrixSusceptibilityGEF_=par;}

   void set_dmat_susceptibility(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) 
        { if (densityMatrixSusceptibilityGEF_) 
              densityMatrixSusceptibilityGEF_->set_susceptibility(susc); }

   //void set_lmo_centroids();

   /// Grab the density matrix susceptibility tensor for the *i*-th LMO
   std::vector<std::shared_ptr<psi::Matrix>> susceptibility(int i) const {return densityMatrixSusceptibilityGEF_->susceptibility(i);}

   /// Grab the density matrix susceptibility tensor *x*-th component for the *i*-th LMO
   std::shared_ptr<psi::Matrix> susceptibility(int i, int x) const {return densityMatrixSusceptibilityGEF_->susceptibility(i,x);}


  protected:
   // ===> Generalized Fragment Parameters <=== //

   /// Density Matrix Susceptibility Tensor
   std::shared_ptr<GenEffPar> densityMatrixSusceptibilityGEF_;

   /// Electrostatic Energy Effective One-Electron Potential
   std::shared_ptr<GenEffPar> electrostaticEnergyGEF_;

   /// Exchange-Repulsion Effective One-Electron Potential
   std::shared_ptr<GenEffPar> repulsionEnergyGEF_;

   /// Charge-Transfer Effective One-Electron Potential
   std::shared_ptr<GenEffPar> chargeTransferEnergyGEF_;

   /// EET Coupling Effective One-Electron Potential
   std::shared_ptr<GenEffPar> EETCouplingConstantGEF_;
};


/** \brief Generalized Effective Fragment Factory. Abstract Base.
 *
 * Describes the GEFP fragment that is in principle designed to work
 * at correlated levels of theory.
 */
class GenEffParFactory //: public std::enable_shared_from_this<GenEffParFactory>
{
  public: 
   /// Construct from wavefunction and Psi4 options
   GenEffParFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

   /// Destruct
   virtual ~GenEffParFactory();

   /// Compute the fragment parameters
   virtual std::shared_ptr<GenEffPar> compute(void) = 0;

   /// Grab wavefunction
   virtual std::shared_ptr<psi::Wavefunction> wfn(void) const {return wfn_;}

   /// Grab options
   virtual psi::Options& options(void) const {return options_;}

  protected:
   /// Wavefunction
   std::shared_ptr<psi::Wavefunction> wfn_;

   /// Psi4 Options
   psi::Options& options_;

   /// Random number generators
   std::default_random_engine randomNumberGenerator_;
   std::uniform_real_distribution<double> randomDistribution_;

   /// Draw random number
   virtual double random_double() {return randomDistribution_(randomNumberGenerator_);};

   /// Draw random point in 3D space, excluding the vdW region
   virtual std::shared_ptr<psi::Vector> draw_random_point();

   /// Is the point inside a vdW region?
   virtual bool is_in_vdWsphere(double x, double y, double z) const;

   /// Matrix with vdW sphere information
   std::shared_ptr<psi::Matrix> excludeSpheres_;

   /// Map with vdW radii
   std::map<std::string, double> vdwRadius_;

   /// Centre-of-mass coordinates
   double cx_, cy_, cz_;

   /// Radius of padding sphere around the molecule
   double radius_;
};

/** \brief Polarization GEFP Factory.
 *
 * Implements creation of the density matrix susceptibility tensors for which \f$ {\bf X} = {\bf 1}\f$.
 * Guarantees the idempotency of the density matrix up to first-order in LCAO-MO variation.
 */
class PolarGEFactory : public GenEffParFactory, public std::enable_shared_from_this<PolarGEFactory>
{
  public:
   /// Construct from CPHF object and Psi4 options
   PolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   PolarGEFactory(std::shared_ptr<CPHF> cphf); 

   /// Destruct
   virtual ~PolarGEFactory();

   /// Compute the density matrix susceptibility tensors
   virtual std::shared_ptr<GenEffPar> compute(void);

  protected:
   /// The CPHF object
   std::shared_ptr<CPHF> cphfSolver_;

   /// Randomly draw electric field value
   std::shared_ptr<psi::Vector> draw_field();

   /// Solve SCF equations to find perturbed one-particle density matrix due to uniform electric field
   std::shared_ptr<psi::Matrix> perturbed_dmat(const std::shared_ptr<psi::Vector>& field);

   /// Solve SCF equations to find perturbed one-particle density matrix due to point charge
   std::shared_ptr<psi::Matrix> perturbed_dmat(const std::shared_ptr<psi::Vector>& pos, const double& charge);

};

/** \brief Polarization GEFP Factory with Least-Squares Scaling of MO Space.
 *
 * Implements creation of the density matrix susceptibility tensors for which \f$ {\bf X} \neq {\bf 1}\f$.
 * Guarantees the idempotency of the density matrix up to first-order in LCAO-MO variation.
 */
class MOScaledPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   MOScaledPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   MOScaledPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~MOScaledPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};

/** \brief Polarization GEFP Factory with Least-Squares Scaling of Cartesian Degrees of freedom.
 *
 *  The resulting density matrix does not guarantee idempotency.
 */
class FieldScaledPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   FieldScaledPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   FieldScaledPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~FieldScaledPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};


/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_gefp_h_ 
