#ifndef _oepdev_libutil_cphf_h
#define _oepdev_libutil_cphf_h

#include<cstdio>
#include<string>
#include<vector>

#include "../../include/oepdev_files.h"

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/factory.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/jk_independent.h"
#include "psi4/libfock/link.h"
#include "psi4/libfock/direct_screening.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/local.h"

#if OEPDEV_USE_PSI4_DIIS_MANAGER == 0
  #include "diis.h"
#else
  #include "psi4/libdiis/diismanager.h"
  #include "psi4/libdiis/diisentry.h"
#endif

namespace oepdev{

using namespace std;
using namespace psi;

template< typename... Args >
std::string string_sprintf( const char* format, Args... args ) {
  int length = std::snprintf( nullptr, 0, format, args... );
  assert( length >= 0 );

  char* buf = new char[length + 1];
  std::snprintf( buf, length + 1, format, args... );

  std::string str( buf );
  delete[] buf;
  return std::move(str);
}

/** \addtogroup OEPDEV_UTILITIES
 * @{
 */

/**\brief CPHF solver class.
 *
 * Solves CPHF equations (now only for RHF wavefunction).
 * Computes molecular and polarizabilities associated with the localized
 * molecular orbitals (LMO). 
 * \note Useful options:
 *       - `CPHF_CONVER`    - convergence of CPHF. Default: `1e-8` (au)
 *       - `CPHF_CONVER`    - maximum numberof iterations. Default: `50`
 *       - `CPHF_DIIS`      - wheather use DIIS or not. Default: `True`
 *       - `CPHF_DIIS_DIM`  - dimension of iterative subspace. Default: `3`
 *       - `CPHF_LOCALIZER` - set orbital localization method. Available: `BOYS` and `PIPEK_MEZEY`. Default: `BOYS`
 */
class CPHF {
   protected:
     /// Wavefunction object
     std::shared_ptr<psi::Wavefunction> _wfn;
     /// Orbital localizer
     std::shared_ptr<Localizer> _localizer;
     /// Number of occupied orbitals
     const int _no;
     /// Number of virtual orbitals
     const int _nv;
     /// Number of basis functions
     const int _nn;
     /// Memory
     long int _memory;
     /// Maximum number of iterations
     int _maxiter;
     /// CPHF convergence threshold
     double _conv;
     /// whether use DIIS or not
     bool _with_diis;
     /// Size of subspace 
     const int _diis_dim;
     /// Primary Basis Set
     std::shared_ptr<BasisSet> _primary;
     /// Occupied orbitals
     std::shared_ptr<Matrix> _cocc;
     /// Virtual orbitals
     std::shared_ptr<Matrix> _cvir;
     /// Occupied orbital energies
     std::shared_ptr<Vector> _eps_occ;
     /// Virtual orbital energies
     std::shared_ptr<Vector> _eps_vir;
     /// the DIIS managers for each perturbation operator x, y and z
     #if OEPDEV_USE_PSI4_DIIS_MANAGER == 0
       std::vector<std::shared_ptr<oepdev::DIISManager>> _diis;
     #else
       std::vector<std::shared_ptr<psi::DIISManager>> _diis;
     #endif
     /// Options
     Options _options;

     // <--- target quantities ---> //

     /// Total (molecular) polarizability tensor
     std::shared_ptr<Matrix> _molecularPolarizability;

     /// orbital-associated polarizabilities tensors
     std::vector<std::shared_ptr<Matrix>> _orbitalPolarizabilities;

     /// LMO centroids
     std::vector<std::shared_ptr<Vector>> _orbitalCentroids;
 
   public:
     /// \brief Constructor
     ///
     /// @param ref_wfn reference HF wavefunction
     /// @param options set of Psi4 options 
     ///
     CPHF(SharedWavefunction ref_wfn, Options& options);

     /// Desctructor
    ~CPHF();

     /// run the calculations
     void compute(void);

     /// print to output file
     void print(void) const;

     // <--- getters ---> //

     /// get the number of occupied orbitals
     int nocc(void) const {return _no;}

     /// retrieve the molecular (total) polarizability
     std::shared_ptr<Matrix> polarizability(void) const {return _molecularPolarizability;}

     /// retrieve the *i*-th orbital-associated polarizability
     std::shared_ptr<Matrix> polarizability(int i) const {return _orbitalPolarizabilities[i];}

     /// retrieve the *i*-th orbital (LMO) centroid 
     std::shared_ptr<Vector> lmo_centroid(int i) const {return _orbitalCentroids[i];}

     /// retrieve the orbital localizer
     std::shared_ptr<Localizer> localizer(void) const {return _localizer;}
};

/** @}*/
}// EndNameSpace oepdev

/** \example example_cphf.cc
 *  Shows how to use the oepdev::CPHF solver to compute molecular and LMO-distributed polarizabilities
 *  at RHF level of theory.
 */

#endif // _oepdev_libutil_cphf_h

