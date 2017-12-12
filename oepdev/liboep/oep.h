#ifndef OEPDEV_LIBOEP_H
#define OEPDEV_LIBOEP_H

#include<cstdio>
#include<string>
#include<vector>

#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/options.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libthce/thce.h"

namespace oepdev_liboep{

using namespace psi;
using namespace std;
using SharedWavefunction = std::shared_ptr<Wavefunction>;
using SharedBasisSet     = std::shared_ptr<BasisSet>;
using SharedTensor       = std::shared_ptr<Tensor>;

/* \brief Generalized One-Electron Potential
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

  public:
    // <--- Constructors and Destructor ---> //

    /// ESP-based OEP object
    OEPotential(SharedWavefunction wfn, Options& options);
    /// DF-based OEP object
    OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);

    /// Destructor
    virtual ~OEPotential();
};


} // EndNameSpace oepdev_liboep
