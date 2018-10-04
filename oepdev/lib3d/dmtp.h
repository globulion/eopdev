#ifndef _oepdev_libutil_dmtp_h
#define _oepdev_libutil_dmtp_h
/** @file dmtp.h */

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
namespace psi{
 using SharedBasisSet = std::shared_ptr<BasisSet>;
}

namespace oepdev{

/** \addtogroup OEPDEV_3DFIELDS
 * @{
 */

using namespace std;
using namespace psi;

class DMTPole;

/** \brief Multipole Convergence. 
 *
 * Handles the convergence of the distributed multipole expansions up to hexadecapole.
 */
class MultipoleConvergence
{
  public:
    enum ConvergenceLevel {R1, R2, R3, R4, R5};
    enum Property {Energy, Potential};
    MultipoleConvergence(std::shared_ptr<DMTPole> dmtp1, std::shared_ptr<DMTPole> dmtp2, 
                         ConvergenceLevel max_clevel = R4);
    virtual ~MultipoleConvergence();
    void compute(Property property = Energy);
    std::vector<double> level(ConvergenceLevel clevel = R4);
  protected:
    ConvergenceLevel max_clevel_;
    std::shared_ptr<DMTPole> dmtp_1_;
    std::shared_ptr<DMTPole> dmtp_2_;
    std::map<std::string, std::shared_ptr<psi::Vector>> convergenceList_;
    void compute_energy();
    void compute_potential();
};

/** \brief Distributed Multipole Analysis Container and Computer. Abstract Base.
 *
 * Handles the distributed multipole expansions up to hexadecapole.
 */
class DMTPole : public std::enable_shared_from_this<DMTPole>
{
  friend class MultipoleConvergence;

  public:

    // <--- Constructors and Destructor ---> //


    /** \brief Build an empty DMTP object from the wavefunction.
     *
     *  @param type - DMTP method. Available: `CAMM`.
     *  @param wfn - wavefunction
     *  @param n   - number of DMTP sets
     *  @return DMTP distribution
     */
    static std::shared_ptr<DMTPole> build(const std::string& type,
                                          std::shared_ptr<psi::Wavefunction> wfn, 
                                          int n = 1);

    /// Destructor
    virtual ~DMTPole();
  
 
    // <--- Accessors ---> //

    /// Has distributed charges?
    virtual bool has_charges()       const {return hasCharges_;      }
    /// Has distributed dipoles?
    virtual bool has_dipoles()       const {return hasDipoles_;      }
    /// Has distributed quadrupoles?
    virtual bool has_quadrupoles()   const {return hasQuadrupoles_;  }
    /// Has distributed octupoles?
    virtual bool has_octupoles()     const {return hasOctupoles_;    }
    /// Has distributed hexadecapoles?
    virtual bool has_hexadecapoles() const {return hasHexadecapoles_;}

    /// Get the positions of distribution centres
    virtual psi::SharedMatrix centres() const {return centres_;}
    /// Get the positions of distribution origins
    virtual psi::SharedMatrix origins() const {return origins_;}

    /// Get the distributed charges
    virtual std::vector<psi::SharedMatrix> charges() const {return charges_;}
    /// Get the distributed dipoles
    virtual std::vector<psi::SharedMatrix> dipoles() const {return dipoles_;}
    /// Get the distributed quadrupoles
    virtual std::vector<psi::SharedMatrix> quadrupoles() const {return quadrupoles_;}
    /// Get the distributed octupoles
    virtual std::vector<psi::SharedMatrix> octupoles() const {return octupoles_;}
    /// Get the distributed hexadecapoles
    virtual std::vector<psi::SharedMatrix> hexadecapoles() const {return hexadecapoles_;}
    /// Get the distributed charges for the \f$ i \f$th distribution
    virtual psi::SharedMatrix charges(int i) const {return charges_.at(i);}
    /// Get the distributed dipoles for the \f$ i \f$th distribution
    virtual psi::SharedMatrix dipoles(int i) const {return dipoles_.at(i);}
    /// Get the distributed quadrupoles for the \f$ i \f$th distribution
    virtual psi::SharedMatrix quadrupoles(int i) const {return quadrupoles_.at(i);}
    /// Get the distributed octupoles for the \f$ i \f$th distribution
    virtual psi::SharedMatrix octupoles(int i) const {return octupoles_.at(i);}
    /// Get the distributed hexadecapoles for the \f$ i \f$th distribution
    virtual psi::SharedMatrix hexadecapoles(int i) const {return hexadecapoles_.at(i);}

    /// Get the number of distributed sites
    virtual int n_sites() const {return nSites_;}
    /// Get the number of distributions
    virtual int n_dmtp() const {return nDMTPs_;}


    // <--- Mutators ---> //

    /// Set the distributed charges
    void set_charges(std::vector<psi::SharedMatrix> M) {charges_ = M;}
    /// Set the distributed dipoles
    void set_dipoles(std::vector<psi::SharedMatrix> M) {dipoles_ = M;}
    /// Set the distributed quadrupoles
    void set_quadrupoles(std::vector<psi::SharedMatrix> M) {quadrupoles_ = M;}
    /// Set the distributed octupoles
    void set_octupoles(std::vector<psi::SharedMatrix> M) {octupoles_ = M;}
    /// Set the distributed hexadecapoles
    void set_hexadecapoles(std::vector<psi::SharedMatrix> M) {hexadecapoles_ = M;}

    /// Set the distributed charges for the \f$ i \f$th distribution
    void set_charges(psi::SharedMatrix M, int i) {charges_[i] = std::make_shared<psi::Matrix>(M);}
    /// Set the distributed dipoles for the \f$ i \f$th distribution
    void set_dipoles(psi::SharedMatrix M, int i) {dipoles_[i] = std::make_shared<psi::Matrix>(M);}
    /// Set the distributed quadrupoles for the \f$ i \f$th distribution
    void set_quadrupoles(psi::SharedMatrix M, int i) {quadrupoles_[i] = std::make_shared<psi::Matrix>(M);}
    /// Set the distributed octupoles for the \f$ i \f$th distribution
    void set_octupoles(psi::SharedMatrix M, int i) {octupoles_[i] = std::make_shared<psi::Matrix>(M);}
    /// Set the distributed hexadecapoles for the \f$ i \f$th distribution
    void set_hexadecapoles(psi::SharedMatrix M, int i) {hexadecapoles_[i] = std::make_shared<psi::Matrix>(M);}

    /// Change origins of the distributed multipole moments of ith set
    virtual void recenter(psi::SharedMatrix new_origins, int i);
    /// Change origins of the distributed multipole moments of all sets
    virtual void recenter(psi::SharedMatrix new_origins);


    // <--- Computers ---> //

    /// Compute DMTP's from the one-particle density matrix
    virtual void compute(psi::SharedMatrix D, bool transition, int i) = 0;
    /// Compute DMTP's from the set of the one-particle density matrices
    void compute(std::vector<psi::SharedMatrix> D, std::vector<bool> transition);
    /// Compute DMTP's from the *sum* of the ground-state alpha and beta one-particle density matrices (transition=false, i=0)
    void compute(void);

    /** \brief Evaluate the generalized interaction energy.
     *
     *  @param other - interacting DMTP distribution. 
     *  @param type  - convergence level (see below).
     *  @return The generalized interaction energy convergence (A.U. units)
     *
     *  The following convergence levels are available:
     *    - `R-1`: includes qq terms.
     *    - `R-2`: includes dq terms and above.
     *    - `R-3`: includes qQ, dd terms and above.
     *    - `R-4`: includes qO, dQ terms and above. 
     *    - `R-5`: includes qH, dO, QQ terms and above.
     */
    std::shared_ptr<MultipoleConvergence> energy(std::shared_ptr<DMTPole> other, const std::string& type = "R-5");

    /** \brief Evaluate the generalized potential.
     *
     *  @param other - interacting DMTP distribution. 
     *  @param type  - convergence level (see below).
     *  @return The generalized potential convergence (A.U. units)
     *
     *  The following convergence levels are available:
     *    - `R-1`: includes qq terms.
     *    - `R-2`: includes dq terms and above.
     *    - `R-3`: includes qQ, dd terms and above.
     *    - `R-4`: includes qO, dQ terms and above. 
     *    - `R-5`: includes qH, dO, QQ terms and above.
     */
     std::shared_ptr<MultipoleConvergence> potential(std::shared_ptr<DMTPole> other, const std::string& type = "R-5");


  protected:

    /** \brief Construct an empty DMTP object from the wavefunction.
     *
     *  @param wfn - wavefunction
     *  @param n   - number of DMTP sets
     *  Do not use this constructor. Use the DMTPole::build method.
     */
    DMTPole(std::shared_ptr<psi::Wavefunction> wfn, int n);


    /// Compute multipole integrals
    void compute_integrals();
    /// Compute maximum order of the integrals
    void compute_order();

    /// Name of the distribution method
    std::string name_;
    /// Molecule associated with this DMTP
    psi::SharedMolecule mol_;
    /// Wavefunction associated with this DMTP
    psi::SharedWavefunction wfn_;
    /// Basis set (primary)
    psi::SharedBasisSet primary_;

    /// Number of DMTP's
    int nDMTPs_;
    /// Number of DMTP sites
    int nSites_;

    /// Maximum order of the multipole
    int order_;

    /// Multipole integrals
    std::vector<psi::SharedMatrix> mpInts_;

    /// Has distributed charges?
    bool hasCharges_;
    /// Has distributed dipoles?
    bool hasDipoles_;
    /// Has distributed quadrupoles?
    bool hasQuadrupoles_;
    /// Has distributed octupoles?
    bool hasOctupoles_;
    /// Has distributed hexadecapoles?
    bool hasHexadecapoles_;

    /// DMTP centres
    psi::SharedMatrix centres_;
    /// DMTP origins
    psi::SharedMatrix origins_;

    /// DMTP charges
    std::vector<psi::SharedMatrix> charges_;
    /// DMTP dipoles
    std::vector<psi::SharedMatrix> dipoles_;
    /// DMTP quadrupoles
    std::vector<psi::SharedMatrix> quadrupoles_;
    /// DMTP octupoles
    std::vector<psi::SharedMatrix> octupoles_;
    /// DMTP hexadecapoles
    std::vector<psi::SharedMatrix> hexadecapoles_;

    /// Initialize and allocate memory
    virtual void allocate();
};

/** \brief Cumulative Atomic Multipole Moments.
 *
 *  Cumulative atomic multipole representation of the molecular charge distribution. Method of Sokalski and Poirier.
 *  Ref.: W. A. Sokalski and R. A. Poirier, *Chem. Phys. Lett.*, 98(1) **1983**
 */
class CAMM : public DMTPole 
{
 public:
   CAMM(psi::SharedWavefunction wfn, int n);
   virtual ~CAMM();
   virtual void compute(psi::SharedMatrix D, bool transition, int n);
 private:
   /// Set the distribution sites to atoms
   void set_sites(void);
};

/** @}*/

} // EndNameSpace oepdev
#endif //_oepdev_libutil_dmtp_h
