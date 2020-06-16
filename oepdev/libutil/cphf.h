#ifndef _oepdev_libutil_cphf_h
#define _oepdev_libutil_cphf_h

#include<cstdio>
#include<string>
#include<vector>

#include "../../include/oepdev_files.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/factory.h"
#include "psi4/libciomr/libciomr.h"
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
 *       - `CPHF_DIIS`      - wheather use DIIS or not. Default: `true`
 *       - `CPHF_DIIS_DIM`  - dimension of iterative subspace. Default: `3`
 *       - `CPHF_LOCALIZE`  - localize the molecular orbitals? Default: `true`
 *       - `CPHF_LOCALIZER` - set orbital localization method. Available: `BOYS` and `PIPEK_MEZEY`. Default: `BOYS`
 */
class CPHF {
   protected:
     /** \name Basic Data */
     //@{
     /// Wavefunction object
     std::shared_ptr<psi::Wavefunction> _wfn;
     /// Options
     Options& _options;
     /// Primary Basis Set
     std::shared_ptr<BasisSet> _primary;
     /// Orbital localizer
     std::shared_ptr<Localizer> _localizer;
     //@}

     /** \name Sizing Information */
     //@{
     /// Number of occupied orbitals
     const int _no;
     /// Number of virtual orbitals
     const int _nv;
     /// Number of basis functions
     const int _nn;
     /// Memory
     long int _memory;
     //@}

     /** \name Parameters of CPHF Calculations */
     //@{
     /// Maximum number of iterations
     int _maxiter;
     /// CPHF convergence threshold
     double _conv;
     /// whether use DIIS or not
     bool _with_diis;
     /// Size of subspace 
     const int _diis_dim;
     //@}

     /** \name Molecular Orbitals */
     //@{
     /// Occupied orbitals
     std::shared_ptr<Matrix> _cocc;
     /// Virtual orbitals
     std::shared_ptr<Matrix> _cvir;
     /// Occupied orbital energies
     std::shared_ptr<Vector> _eps_occ;
     /// Virtual orbital energies
     std::shared_ptr<Vector> _eps_vir;
     /// Transformation from old to new MO's
     std::shared_ptr<psi::Matrix> _T;
     //@}

     /** \name DIIS Manager */
     //@{
     /// the DIIS managers for each perturbation operator x, y and z
     #if OEPDEV_USE_PSI4_DIIS_MANAGER == 0
       std::vector<std::shared_ptr<oepdev::DIISManager>> _diis;
     #else
       std::vector<std::shared_ptr<psi::DIISManager>> _diis;
     #endif
     //@}

     /** \name Response Properties */
     //@{
     /// Total (molecular) polarizability tensor
     std::shared_ptr<Matrix> _molecularPolarizability;

     /// LMO centroids
     std::vector<std::shared_ptr<Vector>> _orbitalCentroids;

     /// orbital-associated polarizability tensors
     std::vector<std::shared_ptr<Matrix>> _orbitalPolarizabilities;

     /// orbital-orbital charge-transfer polarizability tensors
     std::vector<std::vector<std::shared_ptr<Matrix>>> _orbitalChargeTransferPolarizabilities;

     /// Perturbation X Operator O->V matrices in AO basis
     std::vector<std::shared_ptr<Matrix>> _X_OV_ao_matrices;
     /// Perturbation X Operator O->V matrices in MO basis
     std::vector<std::shared_ptr<Matrix>> _X_OV_mo_matrices;
     /// Electric Field Operator O->V matrices in MO basis
     std::vector<std::shared_ptr<Matrix>> _F_OV_mo_matrices;
     //@}


   public:

     /** \name Constructor and Destructor */
     //@{

     /// \brief Constructor
     ///
     /// @param ref_wfn reference HF wavefunction
     /// @param options set of Psi4 options 
     ///
     CPHF(SharedWavefunction ref_wfn, Options& options);

     /// Desctructor
    ~CPHF();
     //@}

     /** \name Executor */
     //@{
     /// run the calculations
     void compute(void);
     //@}

     /** \name Printer */
     //@{
     /// print to output file
     void print(void) const;
     //@}

     /** \name Accessors */
     //@{
     /// get the number of occupied orbitals
     int nocc(void) const {return _no;}

     /// grab the wavefunction
     std::shared_ptr<Wavefunction> wfn(void) const {return _wfn;}

     /// grab the Psi4 options
     Options& options(void) const {return _options;}

     /// retrieve the molecular (total) polarizability
     std::shared_ptr<Matrix> polarizability(void) const {return _molecularPolarizability;}

     /// retrieve the *i*-th orbital-associated polarizability
     std::shared_ptr<Matrix> polarizability(int i) const {return _orbitalPolarizabilities[i];}

     /// retrieve the charge-transfer polarizability associated with orbitals *i* and *j*
     std::shared_ptr<Matrix> polarizability(int i, int j) const {return _orbitalChargeTransferPolarizabilities[i][j];}

     /// retrieve the X operator O-V perturbation matrix in AO basis for *x*-th component
     std::shared_ptr<Matrix> X(int x) const {return _X_OV_ao_matrices[x];}

     /// retrieve the X operator O-V perturbation matrix in AO basis for all three Cartesian components
     std::vector<std::shared_ptr<Matrix>> X(void) const {return _X_OV_ao_matrices;}

     /// retrieve the X operator O-V perturbation matrix in MO basis for *x*-th component
     std::shared_ptr<Matrix> X_mo(int x) const {return _X_OV_mo_matrices[x];}

     /// retrieve the X operator O-V perturbation matrix in MO basis for all three Cartesian components
     std::vector<std::shared_ptr<Matrix>> X_mo(void) const {return _X_OV_mo_matrices;}

     /// retrieve the F operator O-V perturbation matrix in MO basis for *x*-th component
     std::shared_ptr<Matrix> F_mo(int x) const {return _F_OV_mo_matrices[x];}

     /// retrieve the F operator O-V perturbation matrix in MO basis for all three Cartesian components
     std::vector<std::shared_ptr<Matrix>> F_mo(void) const {return _F_OV_mo_matrices;}

     /// retrieve the transformation from old to new MO's
     std::shared_ptr<Matrix> T(void) const {return _T;}

     /// retrieve the Cocc (always Canonical)
     std::shared_ptr<Matrix> Cocc(void) const {return _cocc;}

     /// retrieve the Cvir
     std::shared_ptr<Matrix> Cvir(void) const {return _cvir;}

     /// retrieve the epsocc (always Canonical)
     std::shared_ptr<Vector> epsocc(void) const {return _eps_occ;}

     /// retrieve the epsvir
     std::shared_ptr<Vector> epsvir(void) const {return _eps_vir;}


     /// retrieve the *i*-th orbital (LMO) centroid 
     std::shared_ptr<Vector> lmo_centroid(int i) const {return _orbitalCentroids[i];}

     /// retrieve the orbital localizer
     std::shared_ptr<Localizer> localizer(void) const {return _localizer;}
     //@}
};

/// CPHF object
using SharedCPHF = std::shared_ptr<CPHF>;

/** @}*/
}// EndNameSpace oepdev

/** \example example_cphf.cc
 *  Shows how to use the oepdev::CPHF solver to compute molecular and LMO-distributed polarizabilities
 *  at RHF level of theory.
 */

#endif // _oepdev_libutil_cphf_h

