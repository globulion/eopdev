#ifndef _oepdev_libutil_dmtp_h
#define _oepdev_libutil_dmtp_h
/** @file dmtp.h */

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
namespace psi{
 using SharedBasisSet = std::shared_ptr<BasisSet>;
}

namespace oepdev{

/** \addtogroup OEPDEV_3DFIELDS
 * @{
 */

using namespace std;
using namespace psi;

/** \brief Distributed Multipole Analysis Container and Computer. Abstract Base.
 *
 */
class DMTPole
{
  public:

    // <--- Constructors and Destructor ---> //

    /** \brief Construct an empty DMTP object from the wavefunction.
     *
     */
    DMTPole(std::shared_ptr<psi::Wavefunction> wfn, int n);

    /** \brief Build an empty DMTP object from the wavefunction.
     *
     */
    static std::shared_ptr<DMTPole> build(const std::string& type,
                                          std::shared_ptr<psi::Wavefunction> wfn, 
                                          int n = 1);

    /// Destructor
    virtual ~DMTPole();
  
 
    // <--- Accessors ---> //
    virtual bool has_charges()       const {return hasCharges_;      }
    virtual bool has_dipoles()       const {return hasDipoles_;      }
    virtual bool has_quadrupoles()   const {return hasQuadrupoles_;  }
    virtual bool has_octupoles()     const {return hasOctupoles_;    }
    virtual bool has_hexadecapoles() const {return hasHexadecapoles_;}

    /// Get the distribution centres
    virtual psi::SharedMatrix centres() const {return centres_;}

    /// Get the origins
    virtual psi::SharedMatrix origins() const {return origins_;}

    /// Get the distributed charges, dipoles, quadrupoles, octupoles and hexadecapoles
    virtual std::vector<psi::SharedMatrix> charges() const {return charges_;}
    virtual std::vector<psi::SharedMatrix> dipoles() const {return dipoles_;}
    virtual std::vector<psi::SharedMatrix> quadrupoles() const {return quadrupoles_;}
    virtual std::vector<psi::SharedMatrix> octupoles() const {return octupoles_;}
    virtual std::vector<psi::SharedMatrix> hexadecapoles() const {return hexadecapoles_;}
    virtual psi::SharedMatrix charges(int n) const {return charges_.at(n);}
    virtual psi::SharedMatrix dipoles(int n) const {return dipoles_.at(n);}
    virtual psi::SharedMatrix quadrupoles(int n) const {return quadrupoles_.at(n);}
    virtual psi::SharedMatrix octupoles(int n) const {return octupoles_.at(n);}
    virtual psi::SharedMatrix hexadecapoles(int n) const {return hexadecapoles_.at(n);}

    virtual int n_sites() const {return nSites_;}
    virtual int n_dmtp() const {return nDMTPs_;}


    // <--- Mutators ---> //

    void set_charges(std::vector<psi::SharedMatrix> M) {charges_ = M;}
    void set_dipoles(std::vector<psi::SharedMatrix> M) {dipoles_ = M;}
    void set_quadrupoles(std::vector<psi::SharedMatrix> M) {quadrupoles_ = M;}
    void set_octupoles(std::vector<psi::SharedMatrix> M) {octupoles_ = M;}
    void set_hexadecapoles(std::vector<psi::SharedMatrix> M) {hexadecapoles_ = M;}

    void set_charges(psi::SharedMatrix M, int n) {charges_[n] = std::make_shared<psi::Matrix>(M);}
    void set_dipoles(psi::SharedMatrix M, int n) {dipoles_[n] = std::make_shared<psi::Matrix>(M);}
    void set_quadrupoles(psi::SharedMatrix M, int n) {quadrupoles_[n] = std::make_shared<psi::Matrix>(M);}
    void set_octupoles(psi::SharedMatrix M, int n) {octupoles_[n] = std::make_shared<psi::Matrix>(M);}
    void set_hexadecapoles(psi::SharedMatrix M, int n) {hexadecapoles_[n] = std::make_shared<psi::Matrix>(M);}

    /// Change origins of the distributed multipole moments of ith set
    virtual void recenter(psi::SharedMatrix new_origins, int i);
    /// Change origins of the distributed multipole moments of all sets
    virtual void recenter(psi::SharedMatrix new_origins);



    // <--- Computers ---> //

    /// Compute from the one-particle density matrix
    virtual void compute(psi::SharedMatrix D, bool transition, int i) = 0;
    /// Compute from the set of the one-particle density matrices
    void compute(std::vector<psi::SharedMatrix> D, std::vector<bool> transition);
    /// Compute from the *sum* of the ground-state alpha and beta one-particle density matrices (transition=false, i=0)
    void compute(void);

    /// Evaluate generalized interaction energy
    std::vector<double> energy(std::shared_ptr<DMTPole> other, const std::string& type = "R-5");




  protected:

    /// Compute multipole integrals
    void compute_integrals();
    /// Compute order of the integrals
    void compute_order();


    /// Name
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

    /// 
    bool hasCharges_;
    bool hasDipoles_;
    bool hasQuadrupoles_;
    bool hasOctupoles_;
    bool hasHexadecapoles_;

    /// DMTP centres
    psi::SharedMatrix centres_;

    /// DMTP origins
    psi::SharedMatrix origins_;

    /// DMTP Multipoles
    std::vector<psi::SharedMatrix> charges_;
    std::vector<psi::SharedMatrix> dipoles_;
    std::vector<psi::SharedMatrix> quadrupoles_;
    std::vector<psi::SharedMatrix> octupoles_;
    std::vector<psi::SharedMatrix> hexadecapoles_;

    /// Initialize and allocate memory
    virtual void allocate();
};

/** \brief Cumulative Atomic Multipole Moments.
 *
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
