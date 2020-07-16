#ifndef _oepdev_liboep_oep_h_ 
#define _oepdev_liboep_oep_h_ 
/** @file oep.h */

#include<cstdio>
#include<string>
#include<vector>
#include<map>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/local.h"
#include "../libutil/cis.h"
#include "../libpsi/integral.h"
#include "../libpsi/potential.h"
#include "../lib3d/space3d.h"
#include "../lib3d/dmtp.h"

namespace oepdev{

using namespace psi;
using namespace std;
using SharedWavefunction = std::shared_ptr<Wavefunction>;
using SharedBasisSet     = std::shared_ptr<BasisSet>;
using SharedMatrix       = std::shared_ptr<Matrix>;
using SharedVector       = std::shared_ptr<Vector>;
using SharedDMTPole      = std::shared_ptr<DMTPole>;
using SharedLocalizer    = std::shared_ptr<Localizer>;
using SharedCISData      = std::shared_ptr<CISData>;
/** \addtogroup OEPDEV_OEPS
 * @{
 */


/**
 *  \brief Container to handle the type of One-Electron Potentials.
 */
struct OEPType
{
    /// Name of this type of OEP
    std::string name;
    /// Is this OEP DF-based? 
    bool is_density_fitted;
    /// Number of OEP's within a type
    int n;
    /// All OEP's of this type gathered in a matrix form
    SharedMatrix matrix;
    /// Distributed Multipole Object
    SharedDMTPole dmtp;
    /// CIS data
    SharedCISData cis_data;

    /// Initializer
    OEPType() = default;
    /// Initializer from list
    OEPType(std::string, bool, int, SharedMatrix, SharedDMTPole, SharedCISData);
    /// Copy constructor
    OEPType(const OEPType*);
};


/** \brief Generalized One-Electron Potential: Abstract base.
 * 
 *  Manages OEP's in matrix and 3D forms.
 *  Available OEP categories:
 *   - `ELECTROSTATIC ENERGY`
 *   - `REPULSION ENERGY`
 *   - `CHARGE TRANSFER ENERGY`
 *   - `EET COUPLING CONSTANT`
 */
class OEPotential : public std::enable_shared_from_this<OEPotential> 
{

  protected:

    /// Psi4 options
    Options options_;
    /// Wavefunction
    SharedWavefunction wfn_;
    /// Promary Basis set
    SharedBasisSet primary_;
    /// Auxiliary Basis set
    SharedBasisSet auxiliary_;
    /// Intermediate Basis set
    SharedBasisSet intermediate_;
    /// Molecular Orbital Localizer
    SharedLocalizer localizer_;

    /// Name of this OEP;
    std::string name_;
    /// Types of OEP's within the scope of this object
    std::map<std::string, OEPType> oepTypes_;

    /// Integral factory
    std::shared_ptr<psi::IntegralFactory> intsFactory_;
    /// Matrix of potential one-electron integrals
    psi::SharedMatrix potMat_;
    /// One-electron integral shared pointer
    std::shared_ptr<psi::OneBodyAOInt> OEInt_;
    /// One-electron potential shared pointer
    std::shared_ptr<oepdev::PotentialInt> potInt_;
    /// Occupied orbitals: Canonical (CMO)
    psi::SharedMatrix cOcc_;
    /// Virtual orbitals (Canonical or Valence)
    psi::SharedMatrix cVir_;
    /// Occupied orbital energies: Canonical (CMO)
    psi::SharedVector epsOcc_;
    /// Virtual orbital energies (Canonical or Valence)
    psi::SharedVector epsVir_;
    /// Occupied orbitals: Localized (LMO)
    psi::SharedMatrix lOcc_;
    /// Canonical to Occupied orbitals transformation
    psi::SharedMatrix T_;
    /// LMO Centroids
    std::vector<psi::SharedVector> lmoc_;

    /// Is the object initialized? (MOs computed)
    bool initialized_;

  public:

    // <--- Constructors and Destructor ---> //

    /**\brief Fully ESP-based OEP object
     * 
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, Options& options);

    /**\brief General OEP object
     *
     * @param wfn          - wavefunction
     * @param auxiliary    - auxiliary basis set for density fitting of OEP's
     * @param intermediate - intermediate basis set for density fitting of OEP's
     * @param options      - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);

    /// Copy constructor
    OEPotential(const OEPotential*);

    /// Make a deep copy of this object
    virtual std::shared_ptr<OEPotential> clone(void) const = 0;

    /// Destructor
    virtual ~OEPotential();


    // <--- Factories ---> //      

    /**\brief Build fully ESP-based OEP object
     * 
     * @param type    - OEP category
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */

    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, Options& options);

    /**\brief Build general OEP object
     *
     * @param type         - OEP category
     * @param wfn          - wavefunction
     * @param auxiliary    - auxiliary basis set for density fitting of OEP's
     * @param intermediate - intermediate basis set for density fitting of OEP's
     * @param options      - Psi4 options
     */
    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, 
                                              SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);


    // <--- Methods/Computers ---> //

    /// Compute matrix forms of all OEP's within all OEP types
    virtual void compute(void);

    /** \brief Compute matrix forms of all OEP's within a specified OEP type
     */
    virtual void compute(const std::string& oepType) = 0;

    /// Compute value of potential in point x, y, z and save at v
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) = 0;

    /** \brief Create 3D vector field with OEP
     *  @param oepType - type of OEP. ESP-based OEP is assumed.
     *  @return Vector field 3D with the OEP values.
     */
    std::shared_ptr<OEPotential3D<OEPotential>> make_oeps3d(const std::string& oepType);

    /// Write potential to a cube file
    virtual void write_cube(const std::string& oepType, const std::string& fileName);

    /// Localize Occupied MO's
    virtual void localize(void);

    /// Compute MO centroids from LCAO-MO matrix
    virtual std::vector<psi::SharedVector> mo_centroids(psi::SharedMatrix C);

    /// Rotate 
    virtual void rotate(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux);
    /// Translate
    virtual void translate(psi::SharedVector t);
    /// Superimpose
    virtual void superimpose(const Matrix& refGeometry, 
                             const std::vector<int>& supList, 
                             const std::vector<int>& reordList);


    // <--- Accessors ---> //

    /// Retrieve name of this OEP
    std::string name() const { return name_; }

    /// Retrieve the potentials
    OEPType oep(const std::string& oepType) const {return oepTypes_.at(oepType);}

    /// Retrieve the potentials of a particular OEP type in a matrix form
    SharedMatrix matrix(const std::string& oepType) const { return oepTypes_.at(oepType).matrix; }

    /// Retrieve the number of a particular OEP type
    int n(const std::string& oepType) const {return oepTypes_.at(oepType).n;}

    /// Retrieve wavefunction object
    SharedWavefunction wfn() const {return wfn_;}

    /// Retrieve Canonical occupied MOs
    SharedMatrix cOcc() const {return cOcc_;}
    /// Retrieve Canonical virtual MOs
    SharedMatrix cVir() const {return cVir_;}
    /// Retrieve Canonical occupied MO energies
    SharedVector epsOcc() const {return epsOcc_;}
    /// Retrieve Canonical virtual MO energies
    SharedVector epsVir() const {return epsVir_;}

    /// Retrieve Localized occupied MOs
    SharedMatrix lOcc() const {return lOcc_;}
    /// Retrieve Canonical to Localized occupied MO transformation matrix
    SharedMatrix T() const {return T_;}

    /// Retrieve MO Localizer
    SharedLocalizer localizer() const {return localizer_;}

    /// Retrieve LMO Centroids
    std::vector<std::shared_ptr<psi::Vector>> lmoc() const {return lmoc_;}


    // <--- Mutators ---> //

    /// Set the name of this OEP
    void set_name(const std::string& name) {name_ = name;}

    /// Set the localized molecular orbitals in OEP calculation
    void set_localized_orbitals(std::shared_ptr<psi::Localizer> localizer);

    /// Set the localized molecular orbitals in OEP calculation
    void set_localized_orbitals(std::shared_ptr<OEPotential> oep) {if (oep->localizer_) set_localized_orbitals(oep->localizer_);}

    /// Set the occupied canonical orbitals in OEP calculations
    void set_occupied_canonical_orbitals(std::shared_ptr<OEPotential> oep);

    /// Whether to use localized molecular orbitals in OEP calculation; Default: False
    bool use_localized_orbitals;

    /// Whether to use QUAMBO orbitals to construct VVOs; Default: False
    bool use_quambo_orbitals;


    // <--- Printers ---> //

    /// Header information
    virtual void print_header() const; // = 0;

    /// Print the contents (OEP data)
    void print() const;

    /// Initialize the object (expert) 
    virtual void initialize() = 0;

  protected:
   /// Deep-copy the data
   virtual void copy_from(const OEPotential*);
   /// Rotate basic data
   virtual void rotate_basic(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux);
   /// Translate basic data
   virtual void translate_basic(psi::SharedVector t);

   virtual void rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux);
   virtual void translate_oep(psi::SharedVector t);

   /// Compute MOs (used in initialization stage)
   virtual void compute_molecular_orbitals();


  private:

    /// Initialize defaults
    void common_init();

};


/**\brief Generalized One-Electron Potential for Electrostatic Energy.
 * 
 *  Contains the following OEP types:
 *     - `V` 
 */
class ElectrostaticEnergyOEPotential : public OEPotential 
{
  public:
    /// Only ESP-based potential is worth implementing
    ElectrostaticEnergyOEPotential(SharedWavefunction wfn, Options& options);
    ElectrostaticEnergyOEPotential(const ElectrostaticEnergyOEPotential* f);

    virtual ~ElectrostaticEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;
    virtual std::shared_ptr<OEPotential> clone(void) const override;
    virtual void initialize() override;

  protected:
    virtual void rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) override;
    virtual void translate_oep(psi::SharedVector t) override;


  private:
    /// Set defaults
    void common_init();

    /// Auxilary computers
    double compute_3D_V(const double& x, const double& y, const double& z);
};

/**\brief Generalized One-Electron Potential for Pauli Repulsion Energy.
 * 
 *  Contains the following OEP types:
 *   - `Murrell-etal.S1`
 *   - `Otto-Ladik.S2`
 */
class RepulsionEnergyOEPotential : public OEPotential 
{
  public:
    RepulsionEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);
    RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options);
    RepulsionEnergyOEPotential(const RepulsionEnergyOEPotential* f);

    virtual ~RepulsionEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;
    virtual std::shared_ptr<OEPotential> clone(void) const override;

    virtual void initialize() override;


  protected:
    virtual void rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) override;
    virtual void translate_oep(psi::SharedVector t) override;

  private:
    /// Set defaults
    void common_init();

    /// Auxiliary computers
    void compute_murrell_etal_s1();
    void compute_otto_ladik_s2_esp();
    void compute_otto_ladik_s2_camm_a();
    void compute_otto_ladik_s2_camm_A();

    void compute_3D_otto_ladik_s2(const double& x, const double& y, const double& z);
    double* vec_otto_ladik_s2_;
};

/**\brief Generalized One-Electron Potential for Charge-Transfer Interaction Energy.
 * 
 *  Contains the following OEP types:
 *   - `Otto-Ladik.V1.GDF`     - DF-based term (group I)
 *   - `Otto-Ladik.V3.CAMM-nj` - CAMM-based term (group III; truncated on distributed charges)
 *  
 * Group II terms do not require any particular OEP's due to great siplification of this term.
 * Atomic numbers and LMO centroids are sufficient.
 */
class ChargeTransferEnergyOEPotential : public OEPotential 
{
  public:
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options);
    ChargeTransferEnergyOEPotential(const ChargeTransferEnergyOEPotential* f);

    virtual ~ChargeTransferEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;
    virtual std::shared_ptr<OEPotential> clone(void) const override;

    virtual void initialize() override;

  protected:
    virtual void rotate_oep(psi::SharedMatrix, psi::SharedMatrix, psi::SharedMatrix) override;
    virtual void translate_oep(psi::SharedVector) override;

  private:
    /// Set defaults
    void common_init();

    /// Auxiliary computers
    void compute_otto_ladik_v1_gdf();
  //void compute_otto_ladik_v2();
    void compute_otto_ladik_v3_camm_nj();
};


/**\brief Generalized One-Electron Potential for EET coupling calculations.
 * 
 *  Contains the following OEP types:
 *    - `Fujimoto.GDF` - Joint OEP type for ET(L), ET(HL), HT(H) and HT(HL)
 *    - `Fujimoto.CIS` - CIS data
 *    - `Fujimoto.EXCH`- Pure-exchange coupling matrix \f$ G_{\mu\nu} \equiv (\mu\mu | \nu\nu) \f$
 *    - `Fujimoto.CT_M`- \f$ (HH|LL) \f$ integral for the H_34 Hamiltonian matrix elements (CT)
 */
class EETCouplingOEPotential : public OEPotential 
{
  public:
    EETCouplingOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);
    EETCouplingOEPotential(SharedWavefunction wfn, Options& options);
    EETCouplingOEPotential(const EETCouplingOEPotential* f);

    virtual ~EETCouplingOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;
    virtual std::shared_ptr<OEPotential> clone(void) const override; 

    virtual void initialize() override;

  protected:
    virtual void rotate_oep(psi::SharedMatrix, psi::SharedMatrix, psi::SharedMatrix) override;
    virtual void translate_oep(psi::SharedVector) override;

  private:
    /// Set defaults
    void common_init();

    /// Auxiliary computers
    void compute_fujimoto_gdf();
    void compute_fujimoto_cis();
    void compute_fujimoto_exch();
    void compute_fujimoto_ct_m();

};

using SharedOEPotential = std::shared_ptr<OEPotential>;

/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_oep_h_ 
