#ifndef _oepdev_libutil_util_h
#define _oepdev_libutil_util_h
/** @file util.h */

#include<cstdio>
#include<string>
#include<map>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/local.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"

#include "psi4/libscf_solver/rhf.h"
#include "psi4/libdpd/dpd.h"

namespace oepdev{

using namespace psi;
using namespace std;

using SharedMolecule           = std::shared_ptr<Molecule>;        
using SharedSuperFunctional    = std::shared_ptr<SuperFunctional>;
using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedVector             = std::shared_ptr<Vector>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedBasisSet           = std::shared_ptr<BasisSet>;
using SharedMOSpace            = std::shared_ptr<MOSpace>;
using SharedMOSpaceVector      = std::vector<std::shared_ptr<MOSpace>>;
using SharedIntegralTransform  = std::shared_ptr<IntegralTransform>;
using SharedLocalizer          = std::shared_ptr<Localizer>;
/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief Print preambule for module OEPDEV
 */
extern "C" void preambule(void);


/** \brief Set up DFT functional.
 *
 *  Now it accepts only pure HF functional.
 *  @param name name of the functional ("HF" is now only available)
 *  @param options psi::Options object
 *  @return psi::SharedSuperFunctional object with functional.
 */
extern "C" std::shared_ptr<SuperFunctional> 
create_superfunctional(std::string name, Options& options);


/** \brief Extract molecule from dimer.
 *
 *  @param molecule_dimer psi::SharedMolecule object with dimer
 *  @param id index of a molecule (starts from 1)
 *  @return psi::SharedMolecule object with indicated monomer
 */
extern "C" std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id);


/** \brief Solve RHF-SCF equations for a given molecule in a given basis set.
 *
 *  @param molecule psi::SharedMolecule object with molecule
 *  @param primary shared primary basis set 
 *  @param functional DFT functional
 *  @param options psi::Options object
 *  @param psio psi::PSIO object
 *  @return psi::SharedWavefunction SCF wavefunction of the molecule
 */
extern "C" std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
          std::shared_ptr<SuperFunctional> functional,
          Options& options,
          std::shared_ptr<PSIO> psio);

/** @}*/

} // EndNameSpace oepdev_libutil


#endif // _oepdev_libutil_util_h

