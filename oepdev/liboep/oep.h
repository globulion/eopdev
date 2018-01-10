#ifndef _oepdev_liboep_liboep_h_ 
#define _oepdev_liboep_liboep_h_ 

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

namespace oepdev{

using namespace psi;
using namespace std;
using SharedWavefunction = std::shared_ptr<Wavefunction>;
using SharedBasisSet     = std::shared_ptr<BasisSet>;
using SharedTensor       = std::shared_ptr<Tensor>;
using SharedMatrix       = std::shared_ptr<Matrix>;
using SharedVector       = std::shared_ptr<Vector>;

/* \brief Generalized One-Electron Potential: Abstract base.
 * 
 *  Contains OEP's in matrix and 3D forms.
 */
class OEPotential {
  private:
    /// Initialize defaults
    void common_init();

  protected:
    /// Psi4 options
    Options options_;
    /// Wavefunction
    const SharedWavefunction wfn_;
    /// Promary Basis set
    const SharedBasisSet primary_;
    /// Auxiliary Basis set
    const SharedBasisSet auxiliary_;

    /// Name of this OEP;
    std::string name_;
    /// Types of OEP's within the scope of this object
    std::vector<std::string> oepTypes_;
    /// OEP's matrix forms for each OEP type
    std::map<std::string, SharedTensor> oepMatrices_;

  public:

    // <--- Constructors and Destructor ---> //

    /* \brief ESP-based OEP object
     * 
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, Options& options);

    /* \brief DF-based OEP object
     *
     * @param wfn        - wavefunction
     * @param auxiliary  - basis set for density fitting of OEP's
     * @param options    - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);

    /// Destructor
    virtual ~OEPotential();


    // <--- Factories ---> //      

    /* \brief ESP-based OEP object
     * 
     * @param type    - OEP category
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */
    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, Options& options);
    /* \brief DF-based OEP object
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
    /// Compute 3D potential
    virtual void compute_3D(const std::string& oepType, const std::string& fileName) = 0;
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
    SharedTensor matrix(const std::string& oepType) const { return oepMatrices_.at(oepType); }


    // <--- Mutators ---> //
    void set_name(const std::string& name) {name_ = name;}


    // <--- Printers ---> //
    virtual void print_header() const = 0;


};


/* \brief Generalized One-Electron Potential for Pauli repulsion energy calculations.
 * 
 *  Contains the following OEP types:
 */
class RepulsionEnergyOEPotential : public OEPotential 
{
  private:
    /// Set defaults
    void common_init();
  protected:
  public:
    RepulsionEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~RepulsionEnergyOEPotential();

    virtual void compute(const std::string& oepType);
    virtual void compute_3D(const std::string& oepType, const std::string& fileName);
    virtual void print_header() const;

};

/* \brief Generalized One-Electron Potential for EET coupling calculations.
 * 
 *  Contains the following OEP types:
 *    "ET1" "ET2" "HT1" "HT1" "HT2" "CT1" "CT2"
 */
class EETCouplingOEPotential : public OEPotential 
{
  private:
    /// Set defaults
    void common_init();
  protected:
  public:
    EETCouplingOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    EETCouplingOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~EETCouplingOEPotential();

    virtual void compute(const std::string& oepType);
    virtual void compute_3D(const std::string& oepType, const std::string& fileName);
    virtual void print_header() const;

};


} // EndNameSpace oepdev

#endif // _oepdev_liboep_liboep_h_ 
