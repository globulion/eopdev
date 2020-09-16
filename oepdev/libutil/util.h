#ifndef _oepdev_libutil_util_h
#define _oepdev_libutil_util_h
/** @file util.h */

#include<cstdio>
#include<string>
#include<cmath>
#include<map>
#include<cassert>

#include "psi4/psi4-dec.h"
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
#include "psi4/libmints/mintshelper.h"
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
using SharedBasisSet           = std::shared_ptr<BasisSet>;       
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
extern "C" PSI_API
void preambule(void);

/** \brief Format string output.
 *  Example: std::string text = oepdev::string_sprinff("Test %3d, %13.5f", 5, -10.5425);
 */
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

/** \brief Set up DFT functional.
 *
 *  Now it accepts only pure HF functional.
 *  @param name name of the functional ("HF" is now only available)
 *  @param options psi::Options object
 *  @return psi::SharedSuperFunctional object with functional.
 */
extern "C" PSI_API
std::shared_ptr<SuperFunctional> 
create_superfunctional(std::string name, Options& options);

/** \brief Build BasisSet by Copy.
 *
 *  @param basis_ref - reference basis set
 *  @param molecule_target - target molecule
 *  @return psi::SharedBasisSet object.
 */
extern "C" PSI_API
SharedBasisSet
create_basisset_by_copy(SharedBasisSet basis_ref, SharedMolecule molecule_target);

/** \brief Build BasisSet by Copy for a Particular Atom.
 *
 *  @param basis_ref - reference basis set
 *  @param molecule_target - target molecule (atom in this case)
 *  @param idx_atom - index of an atom in basis_ref->molecule()
 *  @return psi::SharedBasisSet object.
 */
extern "C" PSI_API
SharedBasisSet
create_atom_basisset_by_copy(SharedBasisSet basis_ref, SharedMolecule molecule_target, int idx_atom);

/** \brief Extract molecule from dimer.
 *
 *  @param molecule_dimer psi::SharedMolecule object with dimer
 *  @param id index of a molecule (starts from 1)
 *  @return psi::SharedMolecule object with indicated monomer
 */
extern "C" PSI_API
std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id);

/** \brief Compute distance between two points in nD space.
 *
 *  @param v1 - vector 1
 *  @param v2 - vector 2
 *  @return distance
 *  The vectors have to have the same length.
 */
extern "C" PSI_API
double compute_distance(psi::SharedVector v1,
                        psi::SharedVector v2);


/** \brief Solve RHF-SCF equations for a given molecule in a given basis set.
 *
 *  @param molecule psi::SharedMolecule object with molecule
 *  @param primary basis set 
 *  @param auxiliary basis set 
 *  @param guess basis set 
 *  @param functional DFT functional
 *  @param options psi::Options object
 *  @param psio psi::PSIO object
 *  @param compute_mints Compute integrals (write IWL TOC entry - necessary when transforming integrals)
 *  @return psi::SharedWavefunction SCF wavefunction of the molecule
 */
extern "C" PSI_API
std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
	  std::shared_ptr<BasisSet> auxiliary,
          std::shared_ptr<BasisSet> guess,
          std::shared_ptr<SuperFunctional> functional,
          Options& options,
          std::shared_ptr<PSIO> psio,
	  bool compute_mints = false);

/** \brief Solve RHF-SCF equations for a given molecule in a given basis set.
 *
 *  @param molecule psi::SharedMolecule object with molecule
 *  @param primary shared primary basis set 
 *  @param auxiliary shared auxiliary basis set 
 *  @param sad SAD basis set list
 *  @param sad_fit SAD DF fitting basis set list
 *  @param functional DFT functional
 *  @param options psi::Options object
 *  @param psio psi::PSIO object
 *  @param compute_mints Compute integrals (write IWL TOC entry - necessary when transforming integrals)
 *  @return psi::SharedWavefunction SCF wavefunction of the molecule
 */
extern "C" PSI_API
std::shared_ptr<Wavefunction>
solve_scf_sad(std::shared_ptr<Molecule> molecule, 
              std::shared_ptr<BasisSet> primary,               
	      std::shared_ptr<BasisSet> auxiliary,
              std::vector<std::shared_ptr<BasisSet>> sad,
              std::vector<std::shared_ptr<BasisSet>> sad_fit,
              std::shared_ptr<SuperFunctional> functional,
              Options& options,
              std::shared_ptr<PSIO> psio,
	      bool compute_mints = false);


/** \brief Compute the scalar magnitude of multipole moment.
 *
 *  @param moment - multipole moment vector with unique matrix elements. Now supported only for dipole and quadrupole.
 *  @return       - the average multipole moment value.
 *
 *  The magnitudes of multipole moments are defined here as follows:
 *   - The dipole moment magnitude is just a norm
 *  \f[
 *    \lvert \mu \rvert \equiv \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
 *  \f]
 *   - The quadrupole moment magnitude refers to the traceless moment in Buckingham convention
 *  \f[
 *    \lvert \Theta \rvert \equiv \sqrt{\Theta_{zz}^2 + \frac{1}{3}\left(\Theta_{xx} - \Theta_{yy}\right)^2
 *              + \frac{4}{3} \left( \Theta_{xy}^2 + \Theta_{xz}^2 + \Theta_{yz}^2 \right)}
 *  \f]
 *     In the above equation, the quadrupole moment elements refer to its traceless form.
 */
extern "C" PSI_API
double average_moment(std::shared_ptr<psi::Vector> moment);

/** \brief Compute the Coulomb and exchange integral matrices in MO basis.
 *
 *  Transforms the AO ERI's based on provided C matrix.
 *  @param wfn    - Wavefunction object
 *  @param C      - molecular orbital coefficients (AO x MO)
 *  @return       - vector with J_ij and K_ij matrix
 *
 */
extern "C" PSI_API
std::vector<std::shared_ptr<psi::Matrix>> calculate_JK(std::shared_ptr<psi::Wavefunction> wfn, std::shared_ptr<psi::Matrix> C);

/** \brief Compute the IDF exchange-correlation energy.
 *
 *  TODO
 *  Provide matrices in AO basis!
 */
extern "C" PSI_API
double calculate_idf_xc_energy(
  std::shared_ptr<psi::Wavefunction> wfn, 
  std::shared_ptr<psi::Matrix> D,
  std::vector<std::shared_ptr<psi::Matrix>> f,
  std::vector<std::shared_ptr<psi::Matrix>> g,
  std::shared_ptr<psi::Vector> w, // weights from Gauss-Legendre quadrature
  std::shared_ptr<psi::Vector> o, // omega values from Gauss-Legendre quadrature
  double N, double aN, double xiN, double AN, double wnorm);


/** \brief Compute the Coulomb and exchange integral matrices in MO basis. 
 *
 *  Reads the existing MO ERI's.
 *  @param wfn    - Wavefunction object
 *  @param tr     - IntegralTransform object
 *  @param D      - density matrix in MO basis
 *  @return       - vector with J_ij and K_ij matrix
 *
 */
extern "C" PSI_API
std::vector<std::shared_ptr<psi::Matrix>> calculate_JK_r(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, std::shared_ptr<psi::Matrix> Dij);


/** \brief Compute the derivative of exchange-correlation energy wrt the density matrix in MO-A basis. 
 *
 *  Reads the existing MO ERI's.
 *  @param wfn    - Wavefunction object
 *  @param tr     - IntegralTransform object
 *  @param C      - Transformation matrix MO-B::MO-A (columns are MO-A basis)
 *  @param A      - Vector of matrices A^(n)_{bd}
 *  @return       - derivative matrix in MO-A basis
 *
 */
extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_der_D(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, 
		std::shared_ptr<psi::Matrix> C,
		std::vector<std::shared_ptr<psi::Matrix>> A);

/** \brief Compute the exchange-correlation energy from ERI in MO-SCF basis. 
 *
 *  Reads the existing MO ERI's.
 *  @param wfn    - Wavefunction object
 *  @param tr     - IntegralTransform object
 *  @param f      - f_ij matrix in MO-NEW basis
 *  @param C      - Transformation matrix MO-SCF::MO-NEW (columns are MO-A basis)
 *  @return       - Exchange-correlation energy
 *
 */
extern "C" PSI_API
double calculate_e_xc(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, 
		std::shared_ptr<psi::Matrix> f,
		std::shared_ptr<psi::Matrix> C);

/** \brief Compute the contracted derivative of power of a square and symmetric matrix.
 *
 *  The contracted matrix derivative is defined here as
 *  \f[
 *    {\bf D} = \frac{d {\bf A}^\gamma}{\bf A} : {\mathbb{I}}
 *  \f]
 *  where \f$ {\mathbb{I}} \f$ is the identity matrix.
 *  The derivative, which is the fourth-rank tensor, is computed by 
 *  the forward 2-centre finite difference formula,
 *  \f[
 *    f' = \left( f(h) - f(0) \right) / h
 *  \f]
 *  \notes 
 *    * if \f$ \gamma \f$ is non-integer, input matrix has to be positive-definite.
 *
 *  @param A      - Matrix
 *  @param g      - Power
 *  @param step   - Differentiation step \f$ h \f$
 *  @return       - Contracted derivative (matrix)
 *
 */
extern "C" PSI_API
std::shared_ptr<psi::Matrix>
matrix_power_derivative(std::shared_ptr<psi::Matrix> A,
              		double g, double step);


/** \brief Compute the Effective DFI Potential Matrix Due To Electrons. 
 *
 *  Potential is felt by molecule A and induced by electrons in molecule B.
 *  @param f_aabb - IntegralFactory of type (AA|BB)
 *  @param f_abab - IntegralFactory of type (AB|AB)
 *  @param d_b    - one-particle density matrix in AO basis of B
 *  @return       - V_el(B) matrix in AO basis set of A
 *
 *  If f_abab is nullptr, then only Coulomb matrix is computed.
 *  Otherwise, also exchange contribution is computed.
 */
extern "C"
std::shared_ptr<psi::Matrix> _calculate_DFI_Vel(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::IntegralFactory> f_abab,
                std::shared_ptr<psi::Matrix> d_b);


/** \brief Compute the Effective DFI Coulomb+Exchange Potential Matrix Due To Electrons. 
 *
 *  Potential is felt by molecule A and induced by electrons in molecule B.
 *  @param f_aabb - IntegralFactory of type (AA|BB)
 *  @param f_abab - IntegralFactory of type (AB|AB)
 *  @param d_b    - one-particle density matrix in AO basis of B
 *  @return       - V_el(B) matrix in AO basis set of A
 *
 */
extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_DFI_Vel_JK(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::IntegralFactory> f_abab,
                std::shared_ptr<psi::Matrix> d_b);

/** \brief Compute the Effective DFI Coulomb Potential Matrix Due To Electrons. 
 *
 *  Potential is felt by molecule A and induced by electrons in molecule B.
 *  @param f_aabb - IntegralFactory of type (AA|BB)
 *  @param d_b    - one-particle density matrix in AO basis of B
 *  @return       - V_el(B) matrix in AO basis set of A
 *
 */
extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_DFI_Vel_J(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::Matrix> d_b);


/** \brief Compute the 2-Electron Part of the Effective OEP Matrix for Auxiliary Basis Set Optimization. 
 *
 *  @param nt     - number of test basis functions
 *  @param f_pppt - IntegralFactory of type (PP|PT)
 *  @param ca     - target MOs
 *  @param da     - one-particle density matrix in AO basis
 *  @return       - V matrix
 *
 */
extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_OEP_basisopt_V(const int& nt,
                std::shared_ptr<psi::IntegralFactory> f_pppt,
                std::shared_ptr<psi::Matrix> ca, std::shared_ptr<psi::Matrix> da);


/** \brief Compute the objective function value for auxiliary basis set optimization of OEPs
 *
 *  @param ti - Ti matrix
 *  @param mints - integral helper (instantiated with bsf_i)
 *  @param bsf_m - auxiliary AO basis to optimize
 *  @param bsf_i - intermediate AO basis
 *  @returns value of objective function equal to negative trace of overlap matrix
 */
extern "C" PSI_API
double bs_optimize_projection(std::shared_ptr<psi::Matrix> ti,
                std::shared_ptr<psi::MintsHelper> mints,
                std::shared_ptr<psi::BasisSet> bsf_m, std::shared_ptr<psi::BasisSet> bsf_i);




/** @}*/

} // EndNameSpace oepdev_libutil


#endif // _oepdev_libutil_util_h

