#ifndef _oepdev_libutil_cphf_h
#define _oepdev_libutil_cphf_h

#include<cstdio>
#include<string>
#include<vector>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
//#include "psi4/libmints/tensor.h"
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
#ifdef DIIS_BB
  #include "diis.h"
#else
  #include "psi4/libdiis/diismanager.h"
  #include "psi4/libdiis/diisentry.h"
#endif

namespace oepdev_libutil{

using namespace std;
using namespace psi;

/* \brief CPHF solver class.
 *
 * Solves CPHF equations (now only for RHF wavefunction).
 * Computes molecular and orbital-associated polarizabilities.
 *
 * Suggested usage:
 *
 * std::shared_ptr<CPHF> cphf(new CPHF(ref_wfn, options));
 * cphf->compute();
 * std::shared_ptr<Matrix> polarizability = cphf->get_molecular_polarizability();
 * std::shared_ptr<Tensor> orbital_polars = cphf->get_orbital_polarizabilities();
 *
 */
class CPHF {
   protected:
     /// Number of occupied orbitals
     int _no;
     /// Number of virtual orbitals
     int _nv;
     /// Number of basis functions
     int _nn;
     /// Memory
     long int _memory;
     /// Maximum number of iterations
     int _maxiter;
     /// CPHF convergence threshold
     double _conv;
     /// whether use DIIS or not
     bool _with_diis;
     /// Size of subspace 
     int _diis_dim;
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
     #ifdef DIIS_BB
       std::vector<std::shared_ptr<DIIS>> _diis;
     #else
       std::vector<std::shared_ptr<DIISManager>> _diis;
     #endif
     /// Options
     Options _options;

     // <--- target quantities ---> //

     /// Total (molecular) polarizability tensor
     std::shared_ptr<Matrix> _molecular_polarizability;

     /// orbital-associated polarizabilities tensors
     //std::shared_ptr<Tensor> _orbital_polarizabilities;
 
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

     /// retrieve the molecular (total) polarizability
     std::shared_ptr<Matrix> get_molecular_polarizability(void) const {return _molecular_polarizability;};

     /// retrieve the orbital-associated polarizabilities
     //std::shared_ptr<Tensor> get_orbital_polarizabilities(void) const {return _orbital_polarizabilities;};

};


}// EndNameSpace oepdev_libutil

#endif // _oepdev_libutil_cphf_h

