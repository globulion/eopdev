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
 * Takes shared pointers to existing DMTPole objects and computes the generalized
 * property:
 *  - energy
 *  - potential
 *  from the DMTP sets. The results are stored in vector of length equal to the number
 *  of DMTP's in a set decribed by DMTPole objects given. 
 *
 * \note
 *  The number of DMTP's in each object has to be the same.
 */
class MultipoleConvergence
{
  public:
    /** 
     * Convergence level of the multipole expansion:
     *
     * @param R1 - qq term
     * @param R2 - qd and sum of the above
     * @param R3 - qQ, dd and sum of the above
     * @param R4 - qO, dQ and sum of the above
     * @param R5 - qH, dO, QQ and sum of the above
     *
     */
    enum ConvergenceLevel {R1, R2, R3, R4, R5};

    /** 
     * Property to be evaluated from interacting DMTP's:
     *
     * @param Energy    - generalized energy
     * @param Potential - generalized potential
     *
     */
    enum Property {Energy, Potential};

    /** \brief Construct from two shared DMTPole objects.
     *
     * @param dmtp1      - first DMTPole object
     * @param dmtp2      - second DMTPole object
     * @param max_clevel - maximul allowed convergence level
     */
    MultipoleConvergence(std::shared_ptr<DMTPole> dmtp1, std::shared_ptr<DMTPole> dmtp2, 
                         ConvergenceLevel max_clevel = R5);

    /// Destructor
    virtual ~MultipoleConvergence();

    /** Compute the generalized property
     * 
     *  @param property - generalized Property
     */
    void compute(Property property = Energy);

    /** Grab the generalized property at specified level of convergence
     * 
     * @param clevel - ConvergenceLevel
     * @return vector of results (each element corresponds to each DMTP pair in a set)
     */
    std::shared_ptr<psi::Vector> level(ConvergenceLevel clevel = R5);


  protected:

    /// Maximum allowed convergence level
    ConvergenceLevel max_clevel_;
    /// First DMTP set
    std::shared_ptr<DMTPole> dmtp_1_;
    /// Second DMTP set
    std::shared_ptr<DMTPole> dmtp_2_;
    /// Dictionary of available convergence level results
    std::map<std::string, std::shared_ptr<psi::Vector>> convergenceList_;

    /// Compute the generalized energy
    void compute_energy();
    /// Void compute the generalized potential
    void compute_potential();
};

/** \brief Distributed Multipole Analysis Container and Computer. Abstract Base.
 *
 * Handles the distributed multipole expansions up to hexadecapoles.
 * Distributed centres as well as DMTP origins 
 * are allowed to be located in arbitrary points in space.
 * The object describes a set of \f$ N \f$ DMTP's, that can be generated
 * by providing one-particle density matrices in AO basis. 
 * Nuclear contributions can be switched on or off separately 
 * for each DMTP within a set. The following operations on the DMTP sets are
 * available through the API:
 *  - translation
 *  - rotation
 *  - superimposition
 *  - recentering the origins
 *  - computing the generalized property from another DMTP set
 */
class DMTPole : public std::enable_shared_from_this<DMTPole>
{
  /**
   * Convergence of multipole moment series.
   * @relates MultipoleConvergence
   */
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

    /** 
     *  Change origins of the distributed multipole moments of all sets
     *
     *  @param new_origins - matrix with coordinates of the new origins \f$ \{ {\bf r}_{\rm new} \}\f$.
     *  \note The number of origins has to be equal to the number of distributed centres.
     *
     *  Recentering of the multipoles affects the distributed dipoles and higher moments.
     *  The moments are given as
     *  \f{align*}{
     *   q_{\rm new}  &= q_{\rm old} \\
     *   {\boldsymbol{\upmu}}_{\rm new}   &= {\boldsymbol{\upmu}}_{\rm old} - q_{\rm old} {\boldsymbol \Delta}^{(1)} \\
     *   {\boldsymbol{\Theta}}_{\rm new}  &= {\boldsymbol{\Theta}}_{\rm old} + q_{\rm old} {\boldsymbol \Delta}^{(2)}
     *       - \sum_{\mathscr{P}_2} \mathscr{P}_2 \left[
     *          \left(
     *          q_{\rm old} {\bf r}_{\rm old} +
     *          {\boldsymbol{\upmu}}_{\rm old} 
     *          \right)
     *          \otimes {\boldsymbol \Delta}^{(1)} 
     *          \right] \\
     *   {\boldsymbol{\Omega}}_{\rm new}  &= {\boldsymbol{\Omega}}_{\rm old} - q_{\rm old} {\boldsymbol \Delta}^{(3)}
     *       + \sum_{\mathscr{P}_3} \mathscr{P}_3 \left[ 
     *            \left(
     *            q_{\rm old} {\bf r}_{\rm old} +
     *            {\boldsymbol{\upmu}}_{\rm old}
     *            \right)
     *            \otimes {\boldsymbol \Delta}^{(2)}
     *              \right]
     *       - \sum_{\mathscr{P}_6} \mathscr{P}_6 \left[ 
     *            \left(
     *            q_{\rm old} {\bf r}_{\rm old}^2  +
     *            {\boldsymbol{\upmu}}_{\rm old} \otimes {\bf r}_{\rm old} +
     *            {\boldsymbol{\Theta}}_{\rm old} 
     *            \right)
     *            \otimes {\boldsymbol \Delta}^{(1)}
     *              \right] \\
     *    {\boldsymbol{\Xi}}_{\rm new}  &= {\boldsymbol{\Xi}}_{\rm old} + q_{\rm old} {\boldsymbol \Delta}^{(4)}
     *        - \sum_{\mathscr{P}_3} \mathscr{P}_3 \left[ 
     *            \left(
     *            q_{\rm old} {\bf r}_{\rm old} +
     *            {\boldsymbol{\upmu}}_{\rm old} 
     *            \right)
     *            \otimes {\boldsymbol \Delta}^{(3)}
     *              \right]
     *        + \sum_{\mathscr{P}_3} \mathscr{P}_3 \left[ 
     *            \left(
     *            q_{\rm old} {\bf r}_{\rm old}^2 +
     *            {\boldsymbol{\upmu}}_{\rm old} \otimes {\bf r}_{\rm old} +
     *            {\boldsymbol{\Theta}}_{\rm old} 
     *            \right)
     *            \otimes {\boldsymbol \Delta}^{(2)}
     *              \right]
     *        - \sum_{\mathscr{P}_3} \mathscr{P}_3 \left[ 
     *            \left(
     *            q_{\rm old} {\bf r}_{\rm old}^3 +
     *            {\boldsymbol{\upmu}}_{\rm old} \otimes {\bf r}_{\rm old}^2 +
     *            {\boldsymbol{\Theta}}_{\rm old} \otimes {\bf r}_{\rm old} +
     *            {\boldsymbol{\Omega}}_{\rm old}
     *            \right)
     *            \otimes {\boldsymbol \Delta}^{(1)}
     *              \right]
     *  \f}
     *  where 
     *  \f{align*}{
     *   {\boldsymbol \Delta}^{(1)} &\equiv {\bf r}_{\rm new}   - {\bf r}_{\rm old}   \\
     *   {\boldsymbol \Delta}^{(2)} &\equiv {\bf r}_{\rm new}^2 - {\bf r}_{\rm old}^2 \\
     *   {\boldsymbol \Delta}^{(3)} &\equiv {\bf r}_{\rm new}^3 - {\bf r}_{\rm old}^3 \\
     *   {\boldsymbol \Delta}^{(4)} &\equiv {\bf r}_{\rm new}^4 - {\bf r}_{\rm old}^4
     *  \f}
     *  In the above equations, the distributed centre label was omitted (redundant) as 
     *  each distributed site of multipoles is independent of the others.
     *  TODO - Finish for octupoles and hexadecapoles! -> define the permutation operators!
     */
    virtual void recenter(psi::SharedMatrix new_origins);

    /// Translate the DMTP sets
    void translate(psi::SharedVector transl);
    /// Rotate the DMTP sets
    void rotate(psi::SharedMatrix rotmat);
    /// Superimpose the DMTP sets
    void superimpose(psi::SharedMatrix ref_xyz, std::vector<int> suplist);


    // <--- Computers ---> //

    /// Compute DMTP's from the one-particle density matrix
    virtual void compute(psi::SharedMatrix D, bool transition, int i) = 0;
    /// Compute DMTP's from the set of the one-particle density matrices
    void compute(std::vector<psi::SharedMatrix> D, std::vector<bool> transition);
    /// Compute DMTP's from the *sum* of the ground-state alpha and beta one-particle density matrices (transition=false, i=0)
    void compute(void);

    /** \brief Evaluate the generalized interaction energy.
     *
     *  @param other       - interacting DMTP distribution. 
     *  @param max_clevel  - maximum convergence level (see below).
     *  @return The generalized interaction energy convergence (A.U. units)
     *
     *  The following convergence levels are available:
     *    - `MultipoleConvergence::R1`: includes qq terms.
     *    - `MultipoleConvergence::R2`: includes dq terms and above.
     *    - `MultipoleConvergence::R3`: includes qQ, dd terms and above.
     *    - `MultipoleConvergence::R4`: includes qO, dQ terms and above. 
     *    - `MultipoleConvergence::R5`: includes qH, dO, QQ terms and above.
     */
    std::shared_ptr<MultipoleConvergence> energy(std::shared_ptr<DMTPole> other, 
                                                 MultipoleConvergence::ConvergenceLevel max_clevel = MultipoleConvergence::R5);

    /** \brief Evaluate the generalized potential.
     *
     *  @param other       - interacting DMTP distribution. 
     *  @param max_clevel  - maximum convergence level (see below).
     *  @return The generalized potential convergence (A.U. units)
     *
     *  The following convergence levels are available:
     *    - `MultipoleConvergence::R1`: includes qq terms.
     *    - `MultipoleConvergence::R2`: includes dq terms and above.
     *    - `MultipoleConvergence::R3`: includes qQ, dd terms and above.
     *    - `MultipoleConvergence::R4`: includes qO, dQ terms and above. 
     *    - `MultipoleConvergence::R5`: includes qH, dO, QQ terms and above.
     */
     std::shared_ptr<MultipoleConvergence> potential(std::shared_ptr<DMTPole> other,
                                                 MultipoleConvergence::ConvergenceLevel max_clevel = MultipoleConvergence::R5);


    // <--- Printers ---> //

    /// Print the header
    virtual void print_header() const = 0;

    /// Print the contents
    void print() const;

  protected:

    /** \brief Construct an empty DMTP object from the wavefunction.
     *
     *  @param wfn - wavefunction
     *  @param n   - number of DMTP sets
     *
     *  Do not use this constructor. Use the DMTPole::build method.
     */
    DMTPole(std::shared_ptr<psi::Wavefunction> wfn, int n);


    /// Compute multipole integrals
    void compute_integrals();
    /// Compute maximum order of the integrals
    void compute_order();

    /// Change origins of the distributed multipole moments of ith set
    virtual void recenter(psi::SharedMatrix new_origins, int i);

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

/** 
 *  \brief Cumulative Atomic Multipole Moments.
 *
 *  Cumulative atomic multipole representation of the molecular charge distribution. Method of Sokalski and Poirier.
 *  Ref.: W. A. Sokalski and R. A. Poirier, *Chem. Phys. Lett.*, 98(1) **1983**
 *
 *  # Methodology.
 *  The distributed multipole moments are computed in the following way: 
 *   - first the atomic additive multipole moments (AAMM's) with origins set to the global
 *     coordinate system origin are computed. AO basis set partitioning is used
 *     to dostribute the AAMM's onto the atomic centres.
 *   - subsequently, the AAMM's origins are moved to the corresponding atomic site. 
 *
 *   The computation of the AAMM's is performed according to the following prescription:
 *   \f[
 *   M^{(A)}_{uw\ldots z} ({\bf 0})
 *   = \sum_{\alpha\in A} \sum_{\beta\in{\rm all AO's}} 
 *     D_{\alpha\beta}^{\rm OED} 
 *     \left< \alpha \vert \mathscr{M}_{uw\ldots z} ({\bf 0}) \vert \beta \right>
 *   \f]
 *   where \f$ M^{(A)}_{uw\ldots z} \f$ denotes the \f$ (uw\ldots z) \f$-th component
 *   of the multipole centered at atomic site \f$ A \f$, the symbol \f$ \mathscr{M}({\bf 0}) \f$
 *   is the associated quantum mechanical operator and \f$ D_{\alpha\beta}^{\rm OED} \f$ is the (generalized) 
 *   one-particle density matrx element in AO basis (Greek indices).
 *
 *   Recentering of the multipole moments is described in the documentation of oepdev::DMTPole::recenter.
 *
 */
class CAMM : public DMTPole 
{
 public:
   /// Construct CAMM DMTPole object
   CAMM(psi::SharedWavefunction wfn, int n);
   virtual ~CAMM();
   virtual void compute(psi::SharedMatrix D, bool transition, int n);
   virtual void print_header(void) const;
 private:
   /// Set the distribution sites to atoms
   void set_sites(void);
};

/** @}*/

} // EndNameSpace oepdev
#endif //_oepdev_libutil_dmtp_h
