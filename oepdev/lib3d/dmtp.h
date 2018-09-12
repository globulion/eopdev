#ifndef _oepdev_libutil_dmtp_h
#define _oepdev_libutil_dmtp_h
/** @file dmtp.h */

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"

namespace oepdev{

/** \addtogroup OEPDEV_3DFIELDS
 * @{
 */
 
using namespace std;
using namespace psi;
//using SharedMatrix = std::shared_ptr<psi::Matrix>;

/** \brief Distributed Multipole Analysis Container and Computer. Abstract Base.
 *
 */
class DMTPole
{
  public:

    // <--- Qualifiers ---> // ->probably to be deprecated

    ///**
    // * Type of the DMTP distribution. 
    // *
    // * Based on this qualifier the number and positions of the distribution centres and origins are determined.
    // * `CAMM` - The Cumulative Atomic Multipole Moments 
    // * `MMM`  - The Molecular Multipole Moments
    // * `LMTP` - 
    // * `DMA` - The Distributed Multipole Analysis of Stone
    // * `ESP` - Electrostatic potential charges
    // */
    //enum DistributionType {CAMM, MMM, LMTP, DMA, ESP};

    ///**
    // * Whether it contains nuclear part or not.
    // *
    // * `NucleiIncluded` - Nuclei are included and atomic numbers are added to distributed charges at atomic locations.
    // * `NucleiExcluded` - Nuclei are excluded.
    // */
    //enum NuclearPart {NucleiIncluded, NucleiExcluded};


    // <--- Constructors and Destructor ---> //

    /** \brief Construct an empty DMTP object from the molecule.
     *
     */
    DMTPole(std::shared_ptr<psi::Molecule> mol, int n);

    /** \brief Construct an empty DMTP object from the wavefunction.
     *
     */
    DMTPole(std::shared_ptr<psi::Wavefunction> wfn, int n);

    /** \brief Construct an empty DMTP object from the wavefunction.
     *
     */
    static std::shared_ptr<DMTPole> build(std::shared_ptr<psi::Wavefunction> wfn, 
                                          const std::string& type = "CAMM",
                                          int n = 1);


    /// Destructor
    virtual ~DMTPole();
  
 
    // <--- Accessors ---> //
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


    // <--- Computers ---> //

    /// Compute from the one-particle density matrix
    virtual void compute(psi::SharedMatrix D, int n) = 0;
    /// Compute from the set of the one-particle density matrices
    void compute(std::vector<psi::SharedMatrix> D);

  protected:

    /// Name
    std::string name_;
    /// Molecule associated with this DMTP
    psi::SharedMolecule mol_;
    /// Wavefunction associated with this DMTP
    psi::SharedWavefunction wfn_;
    /// Number of DMTP's
    int nDMTPs_;

    /// Number of distribution centres
    int nCentres_;

    /// Number of origins
    int nOrigins_;

    /// 
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

class CAMM : public DMTPole 
{
 public:
   CAMM(psi::SharedWavefunction wfn, int n);
   virtual ~CAMM();
   virtual void compute(psi::SharedMatrix D, int n);
};

/** @}*/

} // EndNameSpace oepdev
#endif //_oepdev_libutil_dmtp_h
