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
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Unphysical susceptibility!");}
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
   void set_dipole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixDipolePolarizability_=susc;}

   /// Set The Density Matrix Dipole-Dipole Hyperpolarizability
   void set_dipole_dipole_hyperpolarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixDipoleDipoleHyperpolarizability_=susc;}

   /// Set The Density Matrix Quadrupole Polarizability
   void set_quadrupole_polarizability(const std::vector<std::vector<std::shared_ptr<psi::Matrix>>>& susc) {densityMatrixQuadrupolePolarizability_=susc;}


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
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Unphysical susceptibility!");}
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
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Unphysical susceptibility!");}
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
           if      (fieldGradientRank == 0){throw psi::PSIEXCEPTION("Unphysical susceptibility!");}
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



  protected:
   /// The Name of Parameter Type
   std::string name_;

   /// The Density Matrix Dipole Polarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipolePolarizability_;

   /// The Density Matrix Dipole-Dipole Hyperpolarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipoleDipoleHyperpolarizability_;

   /// The Density Matrix Quadrupole Polarizability
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixQuadrupolePolarizability_;
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
class PolarGEFactory : public GenEffParFactory //, public std::enable_shared_from_this<PolarGEFactory>
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

   /// Solve SCF equations to find perturbed one-particle density matrix due set of point charges
   std::shared_ptr<psi::Matrix> perturbed_dmat(const std::shared_ptr<psi::Matrix>& charges);

   /// Evaluate electric field at point (x,y,z) due to point charges 
   std::shared_ptr<psi::Vector> field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                     const double& x, const double& y, const double& z);
   std::shared_ptr<psi::Vector> field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges, 
                                                     const std::shared_ptr<psi::Vector>& pos);

   /// Draw samples of density matrices from constant electric fields
   virtual void draw_samples(std::vector<std::shared_ptr<psi::Matrix>>& electricFieldSet,
                             std::vector<std::shared_ptr<psi::Matrix>>& densityMatrixSet);

   /// Draw samples of density matrices from sets of point charges
   virtual void draw_samples(std::vector<std::shared_ptr<psi::Matrix>>& electricFieldSet,
                             std::vector<std::shared_ptr<psi::Matrix>>& electricFieldGradientSet,
                             std::vector<std::shared_ptr<psi::Matrix>>& densityMatrixSet);

};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements a general class of methods for the density matrix susceptibility tensors represented by:
 * \f[
 *   \delta D_{\alpha\beta} = \sum_i 
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(00)} \cdot {\bf F}
 *        +  {\bf B}_{i;\alpha\beta}^{(10)} : {\bf F} \otimes {\bf F} 
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla \otimes {\bf F} 
 *        +  \ldots
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole-dipole hyperpolarizability
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

   /// Construct from CPHF object and Psi4 options
   GeneralizedPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   /// Construct from CPHF object only (options will be read from CPHF object)
   GeneralizedPolarGEFactory(std::shared_ptr<CPHF> cphf);
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

  protected:

   /// Number of parameter blocks
   int nBlocks_;
   /// Number of distributed sites
   int nSites_;
   /// Dimensionality of entire parameter space
   int nParameters_;
   /// Dimensionality of parameter space per block
   std::vector<int> nParametersBlock_;

   /// Gradient
   std::shared_ptr<psi::Matrix> Gradient_;
   /// Hessian
   std::shared_ptr<psi::Matrix> Hessian_;
   /// Parameters
   std::shared_ptr<psi::Matrix> Parameters_;
   /// Density Matrix Susceptibility Tensors Object
   std::shared_ptr<oepdev::GenEffPar> PolarizationSusceptibilities_;

   /// Dipole Polarizability
   bool hasDipolePolarizability_;
   /// Dipole-Dipole Hyperpolarizability 
   bool hasDipoleDipoleHyperpolarizability_;
   /// Quadrupole Polarizability
   bool hasQuadrupolePolarizability_;
   

   /// Compute the parameters
   void compute_parameters(void);

   /// Perform least-squares fit
   void fit(void);

   /// Save susceptibility tensors associated with the *i*-th and *j*-th basis set function
   void save(int i, int j);

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
   UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~UniformEFieldPolarGEFactory();
   //virtual std::shared_ptr<GenEffPar> compute(void) = 0;
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements a class of density matrix susceptibility models for parameterization
 * in non-uniform electric fields.
 */
class NonUniformEFieldPolarGEFactory : public GeneralizedPolarGEFactory
{
  public:
   NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~NonUniformEFieldPolarGEFactory();
   //virtual std::shared_ptr<GenEffPar> compute(void) = 0;
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx
 *           {\bf B}_{\alpha\beta}^{(00)} \cdot {\bf F}
 * \f]
 * where:
 *  - \f$ {\bf B}_{\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 */
class LinearUniformEFieldPolarGEFactory : public UniformEFieldPolarGEFactory
{
  public:
   LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~LinearUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx 
 *           {\bf B}_{\alpha\beta}^{(00)} \cdot {\bf F}
 *        +  {\bf B}_{\alpha\beta}^{(10)} : {\bf F} \otimes {\bf F} 
 * \f]
 * where:
 *  - \f$ {\bf B}_{\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{\alpha\beta}^{(10)} \f$ is the density matrix dipole-dipole hyperpolarizability
 */
class QuadraticUniformEFieldPolarGEFactory : public UniformEFieldPolarGEFactory
{
  public:
   QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   QuadraticUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~QuadraticUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *           {\bf B}_{i;\alpha\beta}^{(00)} \cdot {\bf F}
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *    defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class LinearNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~LinearNonUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(00)} \cdot {\bf F}
 *        +  {\bf B}_{i;\alpha\beta}^{(10)} : {\bf F} \otimes {\bf F} 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole-dipole hyperpolarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class QuadraticNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~QuadraticNonUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(00)} \cdot {\bf F}
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla \otimes {\bf F} 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(01)} \f$ is the density matrix quadrupole polarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class LinearGradientNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~LinearGradientNonUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Parameterization.
 *
 * Implements the density matrix susceptibility model of the form
 * \f[
 *   \delta D_{\alpha\beta} \approx \sum_i
 *              \left\{ 
 *           {\bf B}_{i;\alpha\beta}^{(00)} \cdot {\bf F}
 *        +  {\bf B}_{i;\alpha\beta}^{(10)} : {\bf F} \otimes {\bf F} 
 *        +  {\bf B}_{i;\alpha\beta}^{(01)} : \nabla \otimes {\bf F} 
 *           \right\}                
 * \f]
 * where:
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(00)} \f$ is the density matrix dipole polarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(10)} \f$ is the density matrix dipole-dipole hyperpolarizability
 *  - \f$ {\bf B}_{i;\alpha\beta}^{(01)} \f$ is the density matrix quadrupole polarizability
 * all defined for the distributed site at \f$ {\bf r}_i \f$.
 */
class QuadraticGradientNonUniformEFieldPolarGEFactory : public NonUniformEFieldPolarGEFactory
{
  public:
   QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);
   QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf);
   virtual ~QuadraticGradientNonUniformEFieldPolarGEFactory();
   //std::shared_ptr<GenEffPar> compute(void);
};

/** \brief Polarization GEFP Factory with Least-Squares Scaling of MO Space.
 *
 * Implements creation of the density matrix susceptibility tensors for which \f$ {\bf X} \neq {\bf 1}\f$.
 * Guarantees the idempotency of the density matrix up to first-order in LCAO-MO variation.
 */
class UnitaryTransformedMOPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   UnitaryTransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   UnitaryTransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~UnitaryTransformedMOPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};

/** \brief Polarization GEFP Factory with Least-Squares Scaling of Cartesian Degrees of freedom.
 *
 *  The resulting density matrix does not guarantee idempotency.
 */
class ScaledXYZPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   ScaledXYZPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   ScaledXYZPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~ScaledXYZPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};

/** \brief Polarization GEFP Factory with Least-Squares Transformation of Cartesian Degrees of freedom.
 *
 *  The resulting density matrix does not guarantee idempotency.
 */
class TransformedXYZPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   TransformedXYZPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   TransformedXYZPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~TransformedXYZPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};


/** \brief Polarization GEFP Factory with Least-Squares Scaling of AO degrees of freedom.
 *
 *  The resulting density matrix does not guarantee idempotency.
 */
class ScaledAOPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   ScaledAOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   ScaledAOPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~ScaledAOPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

};

/** \brief Polarization GEFP Factory with Least-Squares Transformation of MO degrees of freedom.
 *
 *  The resulting density matrix does not guarantee idempotency.
 */
class TransformedMOPolarGEFactory : public PolarGEFactory
{
  public:
   /// Construct from CPHF object and Psi4 options
   TransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt);

   /// Construct from CPHF object only (options will be read from CPHF object)
   TransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf);

   /// Destruct
   virtual ~TransformedMOPolarGEFactory();

   /// Pefrorm Least-Squares Fit
   std::shared_ptr<GenEffPar> compute(void);

  private:
   /// Gradient for the Newton-Raphson optimization
   void gradient(std::shared_ptr<psi::Matrix> A, std::shared_ptr<psi::Matrix> C, std::shared_ptr<psi::Matrix> S);

};



/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_gefp_h_ 
