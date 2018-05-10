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

   /// Grab the density matrix susceptibility tensors
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susceptibility() const {return densityMatrixSusceptibility_;}

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


   /// Grab the density matrix susceptibility tensors
   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susceptibility() const {return densityMatrixSusceptibilityGEF_->susceptibility();}

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