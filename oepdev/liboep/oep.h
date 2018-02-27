#ifndef _oepdev_liboep_liboep_h_ 
#define _oepdev_liboep_liboep_h_ 
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

#include "../libutil/space3d.h"

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

    /// Name of this OEP;
    std::string name_;
    /// Types of OEP's within the scope of this object
    std::vector<std::string> oepTypes_;
    /// OEP's matrix forms for each OEP type
    std::map<std::string, SharedMatrix> oepMatrices_;

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

    /**\brief ESP-based OEP object
     * 
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, Options& options);

    /**\brief DF-based OEP object
     *
     * @param wfn        - wavefunction
     * @param auxiliary  - basis set for density fitting of OEP's
     * @param options    - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);

    /// Destructor
    virtual ~OEPotential();


    // <--- Factories ---> //      

    /**\brief Build ESP-based OEP object
     * 
     * @param type    - OEP category
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */

    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, Options& options);
    /**\brief Build DF-based OEP object
     *
     * @param type       - OEP category
     * @param wfn        - wavefunction
     * @param auxiliary  - basis set for density fitting of OEP's
     * @param options    - Psi4 options
     */
    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, 
                                              SharedBasisSet auxiliary, Options& options);


    // <--- Characterizers ---> //

    /// Is this OEP density-fitted? 
    const bool is_density_fitted;
    /// Is this OEP ESP-based? 
    const bool is_esp_based;


    // <--- Methods/Computers ---> //

    //@{ Compute One Electron Effective Coefficients
    virtual void compute(const std::string& oepType) = 0;
    virtual void compute(void);
    //@}
    //@{ Compute 3D potential
    /** Write potential to a cube file */
    virtual void write_cube(const std::string& oepType, const std::string& fileName);
    /** Compute value of potential in point x, y, z and save at v */
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) = 0;
    //@}
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
    /// Retrieve matrix potential
    SharedMatrix matrix(const std::string& oepType) const { return oepMatrices_.at(oepType); }
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
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
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
    RepulsionEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~RepulsionEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
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
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~ChargeTransferEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
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
    EETCouplingOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    EETCouplingOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~EETCouplingOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();
};

/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_liboep_h_ 
