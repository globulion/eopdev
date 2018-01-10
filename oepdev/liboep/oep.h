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
 * 
 */
class OEPotential {
  private:
  protected:
    /// Psi4 options
    Options _options;
    /// Wavefunction
    const SharedWavefunction _wfn;
    /// Promary Basis set
    const SharedBasisSet _primary;
    /// Auxiliary Basis set
    const SharedBasisSet _auxiliary;

    /// Name of this OEP;
    std::string _name;
    /// Types of OEP's within the scope of this object
    std::vector<std::string> _oepTypes;
    /// OEP's matrix forms for each OEP type
    std::map<std::string, SharedTensor> _oepMatrices;

  public:
    // <--- Constructors and Destructor ---> //

    /// ESP-based OEP object
    OEPotential(const SharedWavefunction& wfn, Options& options);
    /// DF-based OEP object
    OEPotential(const SharedWavefunction& wfn, const SharedBasisSet& auxiliary, Options& options);

    /// Destructor
    virtual ~OEPotential();

    // <--- Characterizers ---> //

    /// Is this OEP density-fitted? 
    const bool is_density_fitted;
    /// Is this OEP ESP-based? 
    const bool is_esp_based;

    // <--- Methods/Computers ---> //

    /// Compute One Electron Effective Coefficients
    virtual void compute(std::string oepType) = 0;
    virtual void compute(void);
    /// Compute 3D potential
    virtual void compute_3D(std::string oepType, std::string fileName) = 0;
    /// Rotate 
    virtual void rotate(const Matrix& rotmat) = 0;
    /// Translate
    virtual void translate(const Vector& trans) = 0;
    /// Superimpose
    virtual void superimpose(const Matrix& refGeometry, 
                             const std::vector<int>& supList, 
                             const std::vector<int>& reordList) = 0;

    // <--- Getters ---> //

    /// Retrieve name of this OEP
    std::string name() const { return _name; }
    /// Retrieve matrix potential
    SharedTensor matrix(std::string oepType) const { return _oepMatrices.at(oepType); }
   


};


} // EndNameSpace oepdev

#endif // _oepdev_liboep_liboep_h_ 
