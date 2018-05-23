#ifndef _oepdev_liboep_gefp_h_ 
#define _oepdev_liboep_gefp_h_ 
/** @file gefp.h */

#include <vector>
#include <string>
#include <random>
#include <cmath>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/vector3.h"
#include "../liboep/oep.h"
#include "../libutil/util.h"
#include "../libutil/cphf.h"
#include "../libutil/scf_perturb.h"

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
   GenEffPar(std::string name) : name_(name), hasDensityMatrixDipolePolarizability_(false), 
                                              hasDensityMatrixDipoleDipoleHyperpolarizability_(false),
                                              hasDensityMatrixQuadrupolePolarizability_(false) {};
   /// Destruct
  ~GenEffPar() {};


   // ---> Mutators <--- //


   /** \brief Set the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field \f$ {\bf F} \f$
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient \f$ \nabla \otimes {\bf F} \f$
    *  @param susc              - the susceptibility tensor
    *
    *  The following susceptibilities are supported (fieldRank, fieldGradientRank):
    *   - (1, 0) - dipole polarizability, interacts with \f$ {\bf F} \f$
    *   - (2, 0) - dipole-dipole hyperpolarizability, interacts with \f$ {\bf F} \otimes {\bf F} \f$
    *   - (0, 1) - quadrupole polarizability, interacts with \f$ \nabla \otimes {\bf F} \f$
    */
   void set_susceptibility(int fieldRank, int fieldGradientRank, const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc)
   {
      std::string notsupported = "Susceptibilities for this rank are not supported yet.";
      if        (fieldRank == 0) { // Not dependent on electric field
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");}
           else if (fieldGradientRank == 1){set_quadrupole_polarizability(susc);}
      } else if (fieldRank == 1) { // Linear wrt electric field
           if (fieldGradientRank == 0)     {set_dipole_polarizability(susc);} 
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else if (fieldRank == 2) {  // Quadratic wrt electric field
           if (fieldGradientRank == 0)     {set_dipole_dipole_hyperpolarizability(susc);}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else {
           throw psi::PSIEXCEPTION(notsupported);
      }
   }

   /// Set The Density Matrix Dipole Polarizability
   void set_dipole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixDipolePolarizability_=susc;hasDensityMatrixDipolePolarizability_=true;}

   /// Set The Density Matrix Dipole-Dipole Hyperpolarizability
   void set_dipole_dipole_hyperpolarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixDipoleDipoleHyperpolarizability_=susc;hasDensityMatrixDipoleDipoleHyperpolarizability_=true;}

   /// Set The Density Matrix Quadrupole Polarizability
   void set_quadrupole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixQuadrupolePolarizability_=susc;hasDensityMatrixQuadrupolePolarizability_=true;}

   /// Set the distributed centres' positions
   void set_centres(const std::vector<std::shared_ptr<psi::Vector>>& centres) {distributedCentres_=centres;}

   // ---> Allocators <--- //


   /** \brief Allocate the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field \f$ {\bf F} \f$
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient \f$ \nabla \otimes {\bf F} \f$
    *  @param nsites            - number of distributed sites
    *  @param nbf               - number of basis functions in the basis set
    *
    *  The following susceptibilities are supported (fieldRank, fieldGradientRank):
    *   - (1, 0) - dipole polarizability, interacts with \f$ {\bf F} \f$
    *   - (2, 0) - dipole-dipole hyperpolarizability, interacts with \f$ {\bf F} \otimes {\bf F} \f$
    *   - (0, 1) - quadrupole polarizability, interacts with \f$ \nabla \otimes {\bf F} \f$
    */
   void allocate(int fieldRank, int fieldGradientRank, int nsites, int nbf)
   {
      std::string notsupported = "Susceptibilities for this rank are not supported yet.";
      if        (fieldRank == 0) { // Not dependent on electric field
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");}
           else if (fieldGradientRank == 1){allocate_quadrupole_polarizability(nsites, nbf);}
      } else if (fieldRank == 1) { // Linear wrt electric field
           if (fieldGradientRank == 0)     {allocate_dipole_polarizability(nsites, nbf);} 
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else if (fieldRank == 2) {  // Quadratic wrt electric field
           if (fieldGradientRank == 0)     {allocate_dipole_dipole_hyperpolarizability(nsites, nbf);}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else {
           throw psi::PSIEXCEPTION(notsupported);
      }
   }

   /// Allocate The Density Matrix Dipole Polarizability
   void allocate_dipole_polarizability(int nsites, int nbf);

   /// Allocate The Density Matrix Dipole-Dipole Hyperpolarizability
   void allocate_dipole_dipole_hyperpolarizability(int nsites, int nbf);

   /// Allocate The Density Matrix Quadrupole Polarizability
   void allocate_quadrupole_polarizability(int nsites, int nbf);


   // ---> Descriptors <--- //
   bool hasDensityMatrixDipolePolarizability() const {return hasDensityMatrixDipolePolarizability_;}
   bool hasDensityMatrixDipoleDipoleHyperpolarizability () const {return hasDensityMatrixDipoleDipoleHyperpolarizability_;}
   bool hasDensityMatrixQuadrupolePolarizability() const {return hasDensityMatrixQuadrupolePolarizability_;}

   // ---> Accessors <--- //


   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    *  @param i                 - id of the distributed site
    *  @param x                 - id of the composite Cartesian component
    *
    *  The following susceptibilities are supported (fieldRank, fieldGradientRank):
    *   - (1, 0) - dipole polarizability, interacts with \f$ {\bf F} \f$
    *   - (2, 0) - dipole-dipole hyperpolarizability, interacts with \f$ {\bf F} \otimes {\bf F} \f$
    *   - (0, 1) - quadrupole polarizability, interacts with \f$ \nabla \otimes {\bf F} \f$
    *   
    *  The distributed sites are assumed to be atomic sites or molecular orbital centroids (depending on the polarization
    *  factory used).
    *  For the electric field, the composite Cartesian index is just an ordinary Cartesian index.
    *  For the electric field gradient and electric field squared, the composite Cartesian index is
    *  given as
    *  \f[
    *    I(x, y) = 3x + y
    *  \f]
    *  where the values of 0, 1 and 2 correspond to *x*, *y* and *z* Cartesian components, respectively.
    *  Therefore, in the latter case, there is 9 distinct composite Cartesian components.
    */
   std::shared_ptr<psi::Matrix> susceptibility(int fieldRank, int fieldGradientRank, int i, int x) const
   {
      std::string notsupported = "Susceptibilities for this rank are not supported yet.";
      if        (fieldRank == 0) { // Not dependent on electric field
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");}
           else if (fieldGradientRank == 1){return quadrupole_polarizability(i, x);}
      } else if (fieldRank == 1) { // Linear wrt electric field
           if (fieldGradientRank == 0)     {return dipole_polarizability(i, x);}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else if (fieldRank == 2) {  // Quadratic wrt electric field
           if (fieldGradientRank == 0)     {return dipole_dipole_hyperpolarizability(i, x);}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else {
           throw psi::PSIEXCEPTION(notsupported);
      }
   }

   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    *  @param i                 - id of the distributed site
    *
    *  The following susceptibilities are supported (fieldRank, fieldGradientRank):
    *   - (1, 0) - dipole polarizability, interacts with \f$ {\bf F} \f$
    *   - (2, 0) - dipole-dipole hyperpolarizability, interacts with \f$ {\bf F} \otimes {\bf F} \f$
    *   - (0, 1) - quadrupole polarizability, interacts with \f$ \nabla \otimes {\bf F} \f$
    *   
    *  The distributed sites are assumed to be atomic sites or molecular orbital centroids (depending on the polarization
    *  factory used).
    */
   std::vector<std::shared_ptr<psi::Matrix>> susceptibility(int fieldRank, int fieldGradientRank, int i) const
   {
      std::string notsupported = "Susceptibilities for this rank are not supported yet.";
      if        (fieldRank == 0) { // Not dependent on electric field
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");}
           else if (fieldGradientRank == 1){return quadrupole_polarizability(i);}
      } else if (fieldRank == 1) { // Linear wrt electric field
           if (fieldGradientRank == 0)     {return dipole_polarizability(i);}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else if (fieldRank == 2) {  // Quadratic wrt electric field
           if (fieldGradientRank == 0)     {return dipole_dipole_hyperpolarizability(i);} 
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else {
           throw psi::PSIEXCEPTION(notsupported);
      }
   }

   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    *
    *  The following susceptibilities are supported (fieldRank, fieldGradientRank):
    *   - (1, 0) - dipole polarizability, interacts with \f$ {\bf F} \f$
    *   - (2, 0) - dipole-dipole hyperpolarizability, interacts with \f$ {\bf F} \otimes {\bf F} \f$
    *   - (0, 1) - quadrupole polarizability, interacts with \f$ \nabla \otimes {\bf F} \f$
    */
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susceptibility(int fieldRank, int fieldGradientRank) const
   {
      std::string notsupported = "Susceptibilities for this rank are not supported yet.";
      if        (fieldRank == 0) { // Not dependent on electric field
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");}
           else if (fieldGradientRank == 1){return quadrupole_polarizability();}
      } else if (fieldRank == 1) { // Linear wrt electric field
           if (fieldGradientRank == 0)     {return dipole_polarizability();}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else if (fieldRank == 2) {  // Quadratic wrt electric field
           if (fieldGradientRank == 0)     {return dipole_dipole_hyperpolarizability();}
           else                            {throw psi::PSIEXCEPTION(notsupported);}
      } else {
           throw psi::PSIEXCEPTION(notsupported);
      }
   }

   /// Grab the density matrix dipole polarizability tensor
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> dipole_polarizability() const {return densityMatrixDipolePolarizability_;}
   /// Grab the density matrix dipole polarizability tensor's *x*-th component
   std::vector<std::shared_ptr<psi::Matrix>> dipole_polarizability(int i) const {return densityMatrixDipolePolarizability_[i];}
   /// Grab the density matrix dipole polarizability tensor's *x*-th component of the *i*-th distributed site
   std::shared_ptr<psi::Matrix> dipole_polarizability(int i, int x) const {return densityMatrixDipolePolarizability_[i][x];}

   /// Grab the density matrix dipole-dipole hyperpolarizability tensor
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> dipole_dipole_hyperpolarizability() const {return densityMatrixDipoleDipoleHyperpolarizability_;}
   /// Grab the density matrix dipole-dipole hyperpolarizability tensor's *x*-th component
   std::vector<std::shared_ptr<psi::Matrix>> dipole_dipole_hyperpolarizability(int i) const {return densityMatrixDipoleDipoleHyperpolarizability_[i];}
   /// Grab the density matrix dipole-dipole hyperpolarizability tensor's *x*-th component of the *i*-th distributed site
   std::shared_ptr<psi::Matrix> dipole_dipole_hyperpolarizability(int i, int x) const {return densityMatrixDipoleDipoleHyperpolarizability_[i][x];}

   /// Grab the density matrix quadrupole polarizability tensor
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> quadrupole_polarizability() const {return densityMatrixQuadrupolePolarizability_;}
   /// Grab the density matrix quadrupole polarizability tensor's *x*-th component
   std::vector<std::shared_ptr<psi::Matrix>> quadrupole_polarizability(int i) const {return densityMatrixQuadrupolePolarizability_[i];}
   /// Grab the density matrix quadrupole polarizability tensor's *x*-th component of the *i*-th distributed site
   std::shared_ptr<psi::Matrix> quadrupole_polarizability(int i, int x) const {return densityMatrixQuadrupolePolarizability_[i][x];}

   /// Grab the centres' positions
   std::vector<std::shared_ptr<psi::Vector>> centres() const {return distributedCentres_;}
   /// Grab the position of the *i*-th distributed site
   std::shared_ptr<psi::Vector> centre(int i) const {return distributedCentres_[i];}

   // ---> Computers <--- //


   /** \brief Compute the density matrix due to the uniform electric field perturbation.
   *
   *  @param field - the uniform electric field vector (A.U.)
   */
   std::shared_ptr<psi::Matrix> compute_density_matrix(std::shared_ptr<psi::Vector> field);
   /** \brief Compute the density matrix due to the uniform electric field perturbation.
   *
   *  @param fx - *x*-th Cartesian component of the uniform electric field vector (A.U.)
   *  @param fy - *y*-th Cartesian component of the uniform electric field vector (A.U.)
   *  @param fz - *z*-th Cartesian component of the uniform electric field vector (A.U.)
   */
   std::shared_ptr<psi::Matrix> compute_density_matrix(double fx, double fy, double fz);
   /** \brief Compute the density matrix due to the non-uniform electric field perturbation.
   *
   *  @param fields - the list of non-uniform electric field vector (A.U.) evaluated at the distributed DMatPol sites
   */
   std::shared_ptr<psi::Matrix> compute_density_matrix(std::vector<std::shared_ptr<psi::Vector>> fields);
   /** \brief Compute the density matrix due to the non-uniform electric field perturbation.
   *
   *  @param fields - the list of electric field vectors (A.U.) evaluated at the distributed DMatPol sites
   *  @param grads  - the list of electric field gradient matrices (A.U.) evaluated at the distributed DMatPol sites
   */
   std::shared_ptr<psi::Matrix> compute_density_matrix(std::vector<std::shared_ptr<psi::Vector>> fields,
                                                       std::vector<std::shared_ptr<psi::Matrix>> grads);



  protected:
   /// The Name of Parameter Type
   std::string name_;

   /// The Density Matrix Dipole Polarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipolePolarizability_;

   /// The Density Matrix Dipole-Dipole Hyperpolarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipoleDipoleHyperpolarizability_;

   /// The Density Matrix Quadrupole Polarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixQuadrupolePolarizability_;

   /// The Positions of the Distributed Centres
   std::vector<std::shared_ptr<psi::Vector>> distributedCentres_;

   /// 
   bool hasDensityMatrixDipolePolarizability_;
   bool hasDensityMatrixDipoleDipoleHyperpolarizability_;
   bool hasDensityMatrixQuadrupolePolarizability_;
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

   /// Set the Density Matrix Susceptibility Tensor Object
   void set_gefp_polarization(const std::shared_ptr<GenEffPar>& par) {densityMatrixSusceptibilityGEF_=par;}

   /// Set the Density Matrix Dipole Polarizability
   void set_dmat_dipole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc)
        { if (densityMatrixSusceptibilityGEF_) 
              densityMatrixSusceptibilityGEF_->set_dipole_polarizability(susc); }

   /// Set the Density Matrix Dipole-Dipole Hyperpolarizability
   void set_dmat_dipole_dipole_hyperpolarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc)
        { if (densityMatrixSusceptibilityGEF_)
              densityMatrixSusceptibilityGEF_->set_dipole_dipole_hyperpolarizability(susc); }

   /// Set the Density Matrix Quadrupole Polarizability
   void set_dmat_quadrupole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc)
        { if (densityMatrixSusceptibilityGEF_)
              densityMatrixSusceptibilityGEF_->set_quadrupole_polarizability(susc); }


   //void set_lmo_centroids();


   // --> Accessors <-- //


   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    *  @param i                 - id of the distributed site
    *  @param x                 - id of the composite Cartesian component
    */
   std::shared_ptr<psi::Matrix> susceptibility(int fieldRank, int fieldGradientRank, int i, int x) const
   {
       return densityMatrixSusceptibilityGEF_->susceptibility(fieldRank, fieldGradientRank, i, x);
   }

   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    *  @param i                 - id of the distributed site
    */
   std::vector<std::shared_ptr<psi::Matrix>> susceptibility(int fieldRank, int fieldGradientRank, int i) const
   {
       return densityMatrixSusceptibilityGEF_->susceptibility(fieldRank, fieldGradientRank, i);
   }

   /** \brief Grab the Density Matrix Susceptibility
    *
    *  @param fieldRank         - power dependency with respect to the electric field
    *  @param fieldGradientRank - power dependency with respect to the electric field gradient
    */
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susceptibility(int fieldRank, int fieldGradientRank) const
   {
       return densityMatrixSusceptibilityGEF_->susceptibility(fieldRank, fieldGradientRank);
   }



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
class GenEffParFactory
{
  public: 
   /** \brief Build Density Matrix Susceptibility Generalized Factory.
    *
    * @param type          - Type of factory
    * @param wfn           - Psi4 wavefunction
    * @param opt           - Psi4 options
    *
    * Available factory types:
    *  - `POLARIZATION` - creates the polarization generalized effective fragment parameters' factory
    * Factory subtype is specified in Psi4 options (input file).
    *
    * \note Useful options:
    *  - `POLARIZATION` factory type:
    *   - `DMATPOL_TRAINING_MODE` - training mode. Default: `EFIELD`
    *   - `DMATPOL_NSAMPLES`      - number of random samples (field or test charges sets). Default: `30`
    *   - `DMATPOL_FIELD_SCALE`   - electric field scale factor (relevant if training mode is `EFIELD`). Default: `0.01` [au]
    *   - `DMATPOL_NTEST_CHARGE`  - number of test charges per sample (relevant if training mode is `CHARGES`). Default: `1`
    *   - `DMATPOL_TEST_CHARGE`   - test charge value (relevant if training mode is `CHARGES`). Default: `0.001` [au]
    *   - `DMATPOL_FIELD_RANK`    - electric field rank. Default: `1`
    *   - `DMATPOL_GRADIENT_RANK` - electric field gradient rank. Default: `0`
    *   - `DMATPOL_TEST_FIELD_X`  - test electric field in X direction. Default: `0.000` [au]
    *   - `DMATPOL_TEST_FIELD_Y`  - test electric field in Y direction. Default: `0.000` [au]
    *   - `DMATPOL_TEST_FIELD_Z`  - test electric field in Z direction. Default: `0.008` [au]
    *   - `DMATPOL_OUT_STATS`     - output file name for statistical evaluation results. Default: `dmatpol.stats.dat`
    */
   static std::shared_ptr<GenEffParFactory> build(const std::string& type, 
                                                  std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

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

   /// Grab the CPHF object
   std::shared_ptr<CPHF> cphf_solver() const {return cphfSolver_;}

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

   /// Number of basis functions
   const int nbf_;

   /// The CPHF object
   std::shared_ptr<CPHF> cphfSolver_;
};

/** \brief Polarization GEFP Factory. Abstract Base.
 * 
 *  Basic interface for the polarization density matrix susceptibility parameters.
 */
class PolarGEFactory : public GenEffParFactory
{
  public:
   /// Construct from Psi4 options
   PolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

   /// Destruct
   virtual ~PolarGEFactory();

   /// Compute the density matrix susceptibility tensors
   virtual std::shared_ptr<GenEffPar> compute(void) = 0;

  protected:
   /// The CPHF object
   ///std::shared_ptr<CPHF> cphfSolver_;

   /// Randomly draw electric field value
   std::shared_ptr<psi::Vector> draw_field();

   /// Randomly draw charge value
   double draw_charge();

   /// Solve SCF equations to find perturbed state due to uniform electric field
   std::shared_ptr<oepdev::RHFPerturbed> perturbed_state(const std::shared_ptr<psi::Vector>& field);

   /// Solve SCF equations to find perturbed state due to point charge
   std::shared_ptr<oepdev::RHFPerturbed> perturbed_state(const std::shared_ptr<psi::Vector>& pos, const double& charge);

   /// Solve SCF equations to find perturbed state due to set of point charges
   std::shared_ptr<oepdev::RHFPerturbed> perturbed_state(const std::shared_ptr<psi::Matrix>& charges);

   /// Evaluate electric field at point (x,y,z) due to point charges 
   std::shared_ptr<psi::Vector> field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                     const double& x, const double& y, const double& z);
   std::shared_ptr<psi::Vector> field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                     const std::shared_ptr<psi::Vector>& pos);

   /// Evaluate electric field gradient at point (x,y,z) due to point charges 
   std::shared_ptr<psi::Matrix> field_gradient_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                              const double& x, const double& y, const double& z);
   std::shared_ptr<psi::Matrix> field_gradient_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                              const std::shared_ptr<psi::Vector>& pos);
};

/** \brief Polarization GEFP Factory from First Principles. Hartree-Fock Approximation.
 *
 * Implements creation of the density matrix susceptibility tensors for which \f$ {\bf X} = {\bf 1}\f$.
 * Guarantees the idempotency of the density matrix up to first-order in LCAO-MO variation.
 * The density matrix susceptibility tensor is represented by:
 * \f[
 *   \delta D_{\alpha\beta} = \sum_i 
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 * \f]
 * where \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 * defined for the distributed LMO site at \f$ {\bf r}_i \f$.
 */
class AbInitioPolarGEFactory : public PolarGEFactory
{
  public:
   AbInitioPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~AbInitioPolarGEFactory();
   virtual std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements a general class of methods for the density matrix susceptibility tensors represented by:
 * \f[
 *   \delta D_{\alpha\beta} = \sum_i 
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 *        +  {\bf B}_{i;\alpha\beta}^{(20)} : {\bf F}({\bf r}_i) \otimes {\bf F}({\bf r}_i) 
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla_i \otimes {\bf F}({\bf r}_i) 
 *        +  \ldots
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(20)} \f$ is the density matrix dipole-dipole hyperpolarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(01)} \f$ is the density matrix quadrupole polarizability
 *
 * all defined for the generalized distributed site at \f$ {\bf r}_i \f$.
 * 
 * Available models:
 *
 *  1. Training against uniform electric fields
 *      - oepdev::LinearUniformEFieldPolarGEFactory - linear with respect to electric field
 *      - oepdev::QuadraticUniformEFieldPolarGEFactory - quadratic with respect to electric field
 *
 *  2. Training against non-uniform electric fields
 *      - oepdev::LinearNonUniformEFieldPolarGEFactory - linear with respect to electric field, distributed site model
 *      - oepdev::QuadraticNonUniformEFieldPolarGEFactory - quadratic with respect to electric field, distributed site model
 *      - oepdev::LinearGradientNonUniformEFieldPolarGEFactory - linear with respect to electric field 
 *        and linear with respect to electric field gradient, distributed site model
 *      - oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory - linear with respect to electric field 
 *        and linear with respect to electric field gradient, distributed site model
 *
 * For the non-linear field training, a set of point charges in each training sample is assumed.
 * Distributed models use atomic centers as expansion points.
 */
class GeneralizedPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from Psi4 wavefunction and options
   GeneralizedPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   /// Destruct
   virtual ~GeneralizedPolarGEFactory();
   /// Pefrorm Least-Squares Fit
   virtual std::shared_ptr<GenEffPar> compute(void);

   /// Dipole Polarizability (interacting with \f$ {\bf F} \f$)
   bool has_dipole_polarizability            () const {return hasDipolePolarizability_           ;}
   /// Dipole-Dipole Hyperpolarizability (interacting with \f$ {\bf F}^2 \f$)
   bool has_dipole_dipole_hyperpolarizability() const {return hasDipoleDipoleHyperpolarizability_;}
   /// Quadrupole Polarizability (interacting with \f$ \nabla \otimes {\bf F} \f$)
   bool has_quadrupole_polarizability        () const {return hasQuadrupolePolarizability_       ;}
   /// Ab Initio Dipole Polarizability (interacting with \f$ {\bf F} \f$)
   bool has_ab_initio_dipole_polarizability  () const {return hasAbInitioDipolePolarizability_   ;}

   /// Grab initial summaric Z value
   double Zinit() const {return Zinit_;}  // Not implemented
   /// Grab final summaric Z value
   double Z() const {return Z_;}          // Not implemented


  protected:

   // --> Dimensions <-- //

   /// Number of parameter blocks
   int nBlocks_;
   /// Number of distributed sites
   int nSites_;
   /// Dimensionality of entire parameter space
   int nParameters_;
   /// Dimensionality of parameter space per block
   std::vector<int> nParametersBlock_;
   /// Number of statistical samples
   const int nSamples_;


   // --> Parameters <-- //

   /// Gradient
   std::shared_ptr<psi::Matrix> Gradient_;
   /// Hessian
   std::shared_ptr<psi::Matrix> Hessian_;
   /// Parameters
   std::shared_ptr<psi::Matrix> Parameters_;
   /// Density Matrix Susceptibility Tensors Object
   std::shared_ptr<oepdev::GenEffPar> PolarizationSusceptibilities_;
   /// Density Matrix Susceptibility Tensors Object for Ab Initio Model
   std::shared_ptr<oepdev::GenEffPar> abInitioPolarizationSusceptibilities_;
   /// Allocate memory
   void allocate(void);
   /// Invert Hessian (do also the identity test)
   void invert_hessian(void);


   // --> Qualifiers <-- //

   /// Has Dipole Polarizability?
   bool hasDipolePolarizability_;
   /// Has Dipole-Dipole Hyperpolarizability?
   bool hasDipoleDipoleHyperpolarizability_;
   /// Has Quadrupole Polarizability?
   bool hasQuadrupolePolarizability_;
   /// Has Ab Initio Dipole Polarizability?
   bool hasAbInitioDipolePolarizability_;


   // --> Sets of statistical data <-- //

   /// A structure to handle statistical data
   struct StatisticalSet {
      /// Interaction energy set
      std::vector<double> InducedInteractionEnergySet;
      /// Density matrix set
      std::vector<std::shared_ptr<psi::Matrix>> DensityMatrixSet;
      /// Induced dipole moment set
      std::vector<std::shared_ptr<psi::Vector>> InducedDipoleSet;
      /// Induced quadrupole moment set
      std::vector<std::shared_ptr<psi::Matrix>> InducedQuadrupoleSet;
      /// Sum of J and K matrix set
      std::vector<std::shared_ptr<psi::Matrix>> JKMatrixSet;
   };

   /// Reference statistical data
   StatisticalSet referenceStatisticalSet_;
   /// Model statistical data
   StatisticalSet modelStatisticalSet_;
   /// Ab Initio Model statistical data
   StatisticalSet abInitioModelStatisticalSet_;
   
   /// Potential matrix set
   std::vector<std::shared_ptr<psi::Matrix>> VMatrixSet_;
   /// Electric field set
   std::vector<std::vector<std::shared_ptr<Vector>>> electricFieldSet_;
   /// Electric field gradient set
   std::vector<std::vector<std::shared_ptr<Matrix>>> electricFieldGradientSet_;
   /// Electric field sum set
   std::vector<std::vector<double>> electricFieldSumSet_;
   /// Electric field gradient sum set
   std::vector<std::vector<std::shared_ptr<psi::Vector>>> electricFieldGradientSumSet_;
   /// Electric field set for Ab Initio Model (LMO-distributed)
   std::vector<std::vector<std::shared_ptr<Vector>>> abInitioModelElectricFieldSet_;


   /// Compute electric field sum set
   void compute_electric_field_sums(void);
   /// Compute electric field gradient sum set
   void compute_electric_field_gradient_sums(void);
   /// Level shifters for Hessian blocks
   const double mField_;
   /// Run the statistical evaluation of results
   void compute_statistics(void);


   // --> Statistical descriptors <-- //

   /// Initial summaric Z value
   double Zinit_;
   /// Final summaric Z value
   double Z_;


   // --> Computers <-- //

   /// Computer of generalized JK objects
   std::shared_ptr<psi::JK> jk_;

   /// Set the distributed centres
   void set_distributed_centres(void);
   
   /// Compute the parameters
   void compute_parameters(void);

   /// Perform least-squares fit
   void fit(void);

   /// Compute ab initio parameters
   void compute_ab_initio(void);

   /// Save susceptibility tensors associated with the *i*-th and *j*-th basis set function
   void save(int i, int j);

   /// Compute samples of density matrices and select electric field distributions
   virtual void compute_samples(void) = 0;

   /// Compute Gradient vector associated with the *i*-th and *j*-th basis set function
   virtual void compute_gradient(int i, int j) = 0;

   /// Compute Hessian matrix (independent on the parameters)
   virtual void compute_hessian(void) = 0;

};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements a class of density matrix susceptibility models for parameterization
 * in uniform electric fields.
 */
class UniformEFieldPolarGEFactory : public GeneralizedPolarGEFactory
{
  public:
   UniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~UniformEFieldPolarGEFactory();
   void compute_samples(void);
   virtual void compute_gradient(int i, int j) = 0;
   virtual void compute_hessian(void) = 0;
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements a class of density matrix susceptibility models for parameterization
 * in non-uniform electric fields.
 */
class NonUniformEFieldPolarGEFactory : public GeneralizedPolarGEFactory
{
  public:
   NonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~NonUniformEFieldPolarGEFactory();
   void compute_samples(void);
   virtual void compute_gradient(int i, int j) = 0;
   virtual void compute_hessian(void) = 0;
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx
 *           {\bf B}_{\alpha\beta}^{(10)} \cdot {\bf F}
 * \f]
 * where:
 *  - \f$ {\bf B}_{\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 */
class LinearUniformEFieldPolarGEFactory : public UniformEFieldPolarGEFactory
{
  public:
   LinearUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~LinearUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx 
 *           {\bf B}_{\alpha\beta}^{(10)} \cdot {\bf F}
 *        +  {\bf B}_{\alpha\beta}^{(20)} : {\bf F} \otimes {\bf F}
 * \f]
 * where:
 *  - \f$ {\bf B}_{\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{\alpha\beta}^{(20)} \f$ is the density matrix dipole-dipole hyperpolarizability
 */
class QuadraticUniformEFieldPolarGEFactory : public UniformEFieldPolarGEFactory
{
  public:
   QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~QuadraticUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *    defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class LinearNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~LinearNonUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 *        +  {\bf B}_{i;\alpha\beta}^{(20)} : {\bf F}({\bf r}_i) \otimes {\bf F}({\bf r}_i) 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(20)} \f$ is the density matrix dipole-dipole hyperpolarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class QuadraticNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~QuadraticNonUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla_i \otimes {\bf F}({\bf r}_i) 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(01)} \f$ is the density matrix quadrupole polarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class LinearGradientNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~LinearGradientNonUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(10)} \cdot {\bf F}({\bf r}_i)
 *        +  {\bf B}_{i;\alpha\beta}^{(20)} : {\bf F} \otimes {\bf F}({\bf r}_i) 
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla_i \otimes {\bf F}({\bf r}_i) 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(20)} \f$ is the density matrix dipole-dipole hyperpolarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(01)} \f$ is the density matrix quadrupole polarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class QuadraticGradientNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~QuadraticGradientNonUniformEFieldPolarGEFactory();
   void compute_gradient(int i, int j);
   void compute_hessian(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Scaling of MO Space.
 *
 * Implements creation of the density matrix susceptibility tensors for which \f$ {\bf X} \neq {\bf 1}\f$.
 * Guarantees the idempotency of the density matrix up to first-order in LCAO-MO variation.
 *
 * \note
 * This method does not give better results than the X=1 method and is extremely time and memory consuming.
 * Therefore, it is placed here only for future reference about solving unitary optimization problem in case
 * it occurs.
 */
class UnitaryTransformedMOPolarGEFactory : public AbInitioPolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   UnitaryTransformedMOPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);

   /// Destruct
   virtual ~UnitaryTransformedMOPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};


/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_gefp_h_ 
