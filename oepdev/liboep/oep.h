#ifndef _oepdev_liboep_oep_h_ 
#define _oepdev_liboep_oep_h_ 
/** @file oep.h */

#include<cstdio>
#include<string>
#include<vector>
#include<map>

#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libthce/thce.h"
#include "psi4/libcubeprop/csg.h"

#include "../lib3d/space3d.h"

namespace oepdev{

using namespace psi;
using namespace std;
using SharedWavefunction = std::shared_ptr<Wavefunction>;
using SharedBasisSet     = std::shared_ptr<BasisSet>;
using SharedTensor       = std::shared_ptr<Tensor>;
using SharedMatrix       = std::shared_ptr<Matrix>;
using SharedVector       = std::shared_ptr<Vector>;
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

    /// Name of this OEP;
    std::string name_;
    /// Types of OEP's within the scope of this object
    std::map<std::string, OEPType> oepTypes_;

    /// Integral factory
    std::shared_ptr<psi::IntegralFactory> intsFactory_;
    /// Matrix of potential one-electron integrals
    std::shared_ptr<psi::Matrix> potMat_;
    /// One-electron integral shared pointer
    std::shared_ptr<psi::OneBodyAOInt> OEInt_;
    /// One-electron potential shared pointer
    std::shared_ptr<     PotentialInt> potInt_;

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

    /// Write potential to a cube file
    virtual void write_cube(const std::string& oepType, const std::string& fileName);

    /// Rotate 
    virtual void rotate(const Matrix& rotmat);
    /// Translate
    virtual void translate(const Vector& trans);
    /// Superimpose
    virtual void superimpose(const Matrix& refGeometry, 
                             const std::vector<int>& supList, 
                             const std::vector<int>& reordList);


    // <--- Accessors ---> //

    /// Retrieve name of this OEP
    std::string name() const { return name_; }

    /// Retrieve the potentials
    OEPType oep(const std::string& oepType) const {return oepTypes_.at(oepType);}

    /// Retrieve the potentials in a matrix form
    SharedMatrix matrix(const std::string& oepType) const { return oepTypes_.at(oepType).matrix; }

    /// Retrieve wavefunction object
    SharedWavefunction wfn() const {return wfn_;}


    // <--- Mutators ---> //
    void set_name(const std::string& name) {name_ = name;}


    // <--- Printers ---> //
    virtual void print_header() const = 0;


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

    virtual ~ElectrostaticEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;

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

    virtual ~RepulsionEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();

    /// Auxiliary computers
    void compute_murrell_etal_s1();
    void compute_murrell_etal_s2();
};

/**\brief Generalized One-Electron Potential for Charge-Transfer Interaction Energy.
 * 
 *  Contains the following OEP types:
 *   - `Otto-Ladik.V1` - DF-based term
 *   - `Otto-Ladik.V2` - ESP-based term
 *   - `Otto-Ladik.V3` - ESP-based term
 */
class ChargeTransferEnergyOEPotential : public OEPotential 
{
  public:
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~ChargeTransferEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();

    /// Auxiliary computers
    void compute_otto_ladik_v1();
    void compute_otto_ladik_v2();
    void compute_otto_ladik_v3();
};


/**\brief Generalized One-Electron Potential for EET coupling calculations.
 * 
 *  Contains the following OEP types:
 *    - `Fujimoto.ET1` 
 *    - `Fujimoto.ET2` 
 *    - `Fujimoto.HT1` 
 *    - `Fujimoto.HT1`
 *    - `Fujimoto.HT2` 
 *    - `Fujimoto.CT1` 
 *    - `Fujimoto.CT2`
 */
class EETCouplingOEPotential : public OEPotential 
{
  public:
    EETCouplingOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);
    EETCouplingOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~EETCouplingOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();
};

/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_oep_h_ 
