#ifndef _oepdev_libutil_cis_h
#define _oepdev_libutil_cis_h
/** @file cis.h */

#include <string>
#include <utility>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libfock/jk.h"

#include "../lib3d/dmtp.h"
#include "davidson_liu.h"


namespace oepdev{

using SharedMolecule           = std::shared_ptr<psi::Molecule>;
using SharedDMTPole            = std::shared_ptr<oepdev::DMTPole>;
using SharedMatrix             = std::shared_ptr<psi::Matrix>;
using SharedVector             = std::shared_ptr<psi::Vector>;
using SharedMOSpace            = std::shared_ptr<psi::MOSpace>;
using SharedMOSpaceVector      = std::vector<std::shared_ptr<psi::MOSpace>>;
using SharedIntegralTransform  = std::shared_ptr<psi::IntegralTransform>;

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief CIS wavefunction parameters. Container structure.
 *
 */
struct CISData
{
    /// Excitation energy
    double E_ex;
    /// CIS HOMO-LUMO amplitude
    double t_homo_lumo;
    /// Excited state density matrix (sum of alpha and beta)
    SharedMatrix Pe;
    /// Transition ground-to-excited state density matrix (sum of alpha and beta)
    SharedMatrix Peg;
    /// TrCAMM
    SharedDMTPole trcamm;
    /// CAMM for HOMO orbital
    SharedDMTPole camm_homo;
    /// CAMM for LUMO orbital
    SharedDMTPole camm_lumo;

    /// Null Constructor
    CISData(void) = default;
    /// Copy Constructor
    CISData(const CISData*);
};


/** \brief CISComputer
 *
 */
class CISComputer : public DavidsonLiu {

  // --> public interface <-- //
  public:

   /** \brief Build CIS Computer.
    *
    * @param type          - Type of computer
    * @param wfn           - Psi4 wavefunction
    * @param opt           - Psi4 options
    * @param reference     - Reference Slater determinant (`RHF`, `UHF` available).
    *
    * Available computer types:
    *  - `RESTRICTED` or `RCIS` - RHF wavefunction is used as reference state
    *  - `UNRESTRICTED` or `UCIS` - UHF wavefunction is used as reference state
    *
    *
    * # Implementation
    * 
    * The CIS Hamiltonian in the basis space of singly-excited Slater determinants
    * is constructed from canonical molecular orbitals (CMO's)
    *
    * \f{align*}{
    *   \big< \Phi_0 \big| \mathscr{H} \big| \Phi_i^a \big> &= 0 \\
    *   \big< \Phi_j^b \big| \mathscr{H} \big| \Phi_i^a \big> &= \delta_{ij}\delta_{ab}
    *    \left( \varepsilon_a - \varepsilon_i \right) + 
    *    \big< aj \big| ib \big> - \big< aj \big| bi \big>
    * \f}
    *
    * where *i* labels the occupied CMO's whereas *a* labels the virtual CMO's.
    * In the above equation, \f$ \big< aj \big| ib \big> \f$ is the 2-electron 4-centre integral 
    * in physicist's notation.
    * After integrating out the spin coordinate, four blocks of Hamiltonian are explicitly given as
    *
    * \f{align*}{
    *     \big< \Phi_j^b \big| \mathscr{H} \big| \Phi_i^a \big> &= \delta_{ij}\delta_{ab}
    *            \left( \varepsilon_a - \varepsilon_i \right) +
    *   \big[ ia \big| jb \big] - \big[ ab \big| ij \big] \\
    *   \big< \Phi_{\overline{j}}^{\overline{b}} \big| \mathscr{H} \big| \Phi_{\overline{i}}^{\overline{a}} \big> &=
    *      \delta_{\overline{i}\overline{j}}\delta_{\overline{a}\overline{b}}
    *     \left( \varepsilon_{\overline{a}} - \varepsilon_{\overline{i}} \right) +
    *     \big[ {\overline{i}}{\overline{a}} \big| {\overline{j}}{\overline{b}} \big] 
    *   - \big[ {\overline{a}}{\overline{b}} \big| {\overline{i}}{\overline{j}} \big] \\
    *   \big< \Phi_{\overline{j}}^{\overline{b}} \big| \mathscr{H} \big| \Phi_i^a \big> &= 
    *   \big[ ia \big| {\overline{j}}{\overline{b}} \big] \\
    *   \big< \Phi_j^b \big| \mathscr{H} \big| \Phi_{\overline{i}}^{\overline{a}} \big> &=
    *   \big[ {\overline{i}}{\overline{a}} \big| jb \big]
    * \f}
    * where the \f$ \big[ ia \big| jb \big] \f$ is the 2-electron 4-centre integral in the chemist's (Coulomb)
    * notation.
    *
    * Such matrix is diagonalized yelding the excitation energies (wrt HF ground state)
    * as well as the CIS coefficients
    *
    * \f[
    *  \sum_{ij}\sum_{ab} t_{i,I}^a H_{ij}^{ab} t_{j,J}^b = E_I \delta_{IJ}
    * \f]
    * where the summations above extend over alpha and beta electron spin labels
    * and \f$ t_{i,I}^a \f$ is the CIS amplitude for the *I*th excited state,
    * associated with the \f$ i\rightarrow a\f$ excitation with respect to the HF reference determinant.
    * Note that \f$ E_I \f$ is *not* the excited state energy, but the energy relative the the HF reference
    * energy.
    *
    * \see For Davidson-Liu solution to CIS problem, see oepdev::R_CISComputer_DL and oepdev::U_CISComputer_DL.
    *
    * ## Transition density matrix
    *
    * AO basis transition density from ground (HF) to excited (CIS) state is given by
    *
    * \f[
    *   P_{\mu\nu}^{(g\rightarrow e)} = 
    *      \sum_i^{\rm Occ} \sum_a^{\rm Vir} t_{i,e}^a C_{\nu i} C_{\mu a} +
    *      \sum_{\overline{i}}^{\rm Occ} \sum_{\overline{a}}^{\rm Vir} 
    *      t_{{\overline{i}},e}^{\overline{a}} C_{\nu \overline{i}} C_{\mu \overline{a}}
    * \f]
    *
    * ## Excited state density matrix
    *
    * CMO basis excited state density matrix for alpha spin is given by
    *
    * Analogous expression is given for the beta spin. 
    * 
    * AO representation of the CMO excited state density matrix is
    * \f[
    *  P_{\mu\nu}^{(e)} = \sum_{pq} C_{\mu p} P_{pq}^{(e)} C_{\nu q}
    *     + \sum_{{\overline p}{\overline q}} C_{\mu {\overline p}} P_{{\overline p}{\overline q}}^{(e)} C_{\nu {\overline q}}
    * \f]
    * which
    * is the sum of alpha and beta density matrices in CMO basis transformed to AO basis.
    * 
    * The CMO excited state density matrix for spin alpha is given by
    * \f[
    * P_{pq}^{(e)} = \left\{\begin{matrix}
    * \delta_{pq} - \sum_{a}^{\rm Vir} t_{p,e}^a t_{q,e}^a &\text{for p,q } \in {\rm Occ}\\ 
    * \sum_{i}^{\rm Occ} t_{i,e}^p t_{i,e}^q &\text{for p,q } \in {\rm Vir}\\ 
    * 0 &\text{otherwise}
    * \end{matrix}\right.
    * \f]
    * The beta spin density matrix is generated analogously as above.
    *
    * The cumulative atomic multipole moments (CAMM) are computed from
    * the excited state density matrices in AO basis. The nuclear contribution is included.
    *
    * ## Transition multipole moments
    *
    * The transition dipole moment is computed from the AO transition density matrix and
    * the dipole integrals in AO basis, i.e.,
    * \f[
    *  \big< \Phi_0 \big| \hat{\upmu}_u \big| \Psi_e \big> = {\rm Tr}\left[ {\bf d}^{(u)} \cdot {\bf P}^{g\rightarrow e}\right]
    * \f]
    *
    * Oscillator strength is computed from the transition dipole moment via
    * \f[
    *   f^{g\rightarrow e} = \frac{2}{3} E_{e} \Big| \big< \Phi_0 \big| \hat{\boldsymbol{\upmu}} \big| \Psi_e \big> \Big|^2
    * \f]
    *
    * Transition cumulative atomic multipole moments (TrCAMM) are computed from
    * the transition density matrices in AO basis. The nuclear contribution is not included.
    *
    * ## Spin angular momentum
    *
    * The expectation value of the \f$ \hat{S}^2\f$ operator is calculated from the CIS amplitudes
    * and MOs of the reference wavefunction according to D. Maurice and M. Head-Gordon, *Int. J. Quant. Chem.*, 
    * **1995**, 95, 010361-10:
    *
    * \f{multline*}{
    *   \left< \hat{S}^2 \right>_{\rm UCIS} = \left< \hat{S}^2 \right>_{\rm UHF}
    *   - {\rm Tr}\left[ {\bf Q}^{(\alpha)}_{\rm Occ} \cdot \left\{ {\bf P}^{(e,\alpha)}_{\rm Occ} - {\bf 1}\right\} \right]
    *   - {\rm Tr}\left[ {\bf Q}^{(\beta )}_{\rm Occ} \cdot \left\{ {\bf P}^{(e,\beta )}_{\rm Occ} - {\bf 1}\right\} \right]\\
    *   - {\rm Tr}\left[ {\bf Q}^{(\alpha)}_{\rm Vir} \cdot {\bf P}^{(e,\alpha)}_{\rm Vir} \right]
    *   - {\rm Tr}\left[ {\bf Q}^{(\beta )}_{\rm Vir} \cdot {\bf P}^{(e,\beta )}_{\rm Vir} \right] 
    *   - 2\sum_{i}^{\rm Occ}\sum_{a}^{\rm Vir}
    *      \sum_{\overline{j}}^{\rm Occ}\sum_{\overline{b}}^{\rm Vir} 
    *      \Delta^*_{i{\overline{j}}} \Delta_{a{\overline{b}}}
    *      t_{i,e}^a t_{{\overline{j}},e}^{\overline{b}}
    * \f}
    * where
    * \f{align*}{
    *    [{\bf Q}^{(\alpha)}_{\rm Occ}]_{ij}                       &= \sum_{\overline{k}}^{\rm Occ} \Delta^*_{\overline{k}i}\Delta_{\overline{k}j} \\
    *    [{\bf Q}^{(\beta )}_{\rm Occ}]_{\overline{i}\overline{j}} &= \sum_{k}^{\rm Occ} \Delta^*_{k\overline{i}}\Delta_{k\overline{j}} \\
    *    [{\bf Q}^{(\alpha)}_{\rm Vir}]_{ab}                       &= \sum_{\overline{k}}^{\rm Occ} \Delta^*_{\overline{k}a}\Delta_{\overline{k}b} \\
    *    [{\bf Q}^{(\beta )}_{\rm Vir}]_{\overline{a}\overline{b}} &= \sum_{k}^{\rm Occ} \Delta^*_{k\overline{a}}\Delta_{k\overline{b}} 
    * \f}
    * and
    * \f[
    *   \Delta_{pq} = \sum_{\mu\nu} C_{\mu p} S_{\mu\nu} C_{\nu p}
    * \f]
    * The diagnostic for UHF spin contamination is given by
    * \f[
    * \left< \hat{S}^2 \right>_{\rm UHF} = \left< \hat{S}^2 \right>_{\rm exact}
    *  + N_{\beta} - \sum_i^{\rm Occ}\sum_{\overline{j}}^{\rm Occ} \vert \Delta_{i\overline{j}} \vert^2
    * \f]
    * with
    * \f[
    * \left< \hat{S}^2 \right>_{\rm exact} = \frac{N_\alpha-N_\beta}{2} 
    *  \left( \frac{N_\alpha-N_\beta+2}{2} \right)
    * \f]
    * and is also printed out to the output file.
    *
    * \note Useful options:
    *   - `CIS_TYPE`    - Algorithm of CIS. Available: `DAVIDSON_LIU` (Default), `DIRECT_EXPLICIT` (only RHF reference), `EXPLICIT`.
    *   - `CIS_SCHWARTZ_CUTOFF`  - Cutoff for Schwartz ERI screening. Default: 0.0. Relevant if `DAVIDSON_LIU` or `DIRECT_EXPLICIT` are chosen as CIS type.
    *   - `CIS_STANDARDIZE_AMPLITUDES`  - If true, CIS amplitudes of each excited state are rephased so that the leading amplitude is positive. Default: true.
    *   - `OEPDEV_AMPLITUDE_PRINT_THRESHOLD`  - Control threshold how many CIS amplitudes to print to the output. Default: 0.1.
    *   - For UHF references, SAD guess might lead to triplet instabilities. It is then better to set `CORE` as the UHF guess
    */
   static std::shared_ptr<CISComputer> build(const std::string& type, 
                                             std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt,
                                             const std::string& reference = "");
    //->add this option for convenience*   - `CIS_NSTATES` - Number of lowest-energy excited states to include. Default: `-1` (means all states are saved).

   /// Destructor
   virtual ~CISComputer();

   /// Solve the CIS problem
   virtual void compute(void);

   /// Clear DPD instance
   virtual void clear_dpd(void);

   /// Get the total number of excited states
   int nstates(void) const {return nstates_;}

   /// Get the CIS eigenvalues
   psi::SharedVector eigenvalues() const {return E_;}
   psi::SharedVector E() const {return E_;}

   /// Get the CIS eigenvectors
   psi::SharedMatrix eigenvectors() const {return U_;}
   psi::SharedMatrix U() const {return U_;}

   /// Get the HOMO+*h*->LUMO+*l* CIS coefficient for a given excited state *I* for spin alpha and beta
   std::pair<double,double> U_homo_lumo(int I, int h=0, int l=0) const;

   /// Compute MO one-particle alpha density matrix for state *i*
   SharedMatrix Da_mo(int i) const;

   /// Compute MO one-particle beta density matrix for state *i*
   SharedMatrix Db_mo(int i) const;

   /// Compute AO one-particle alpha density matrix for state *i*
   SharedMatrix Da_ao(int i) const;

   /// Compute AO one-particle beta density matrix for state *i*
   SharedMatrix Db_ao(int i) const;

   /// Compute CAMM for *j* excited state
   SharedDMTPole camm(int j, bool symmetrize=false) const;

   /// Compute MO one-particle alpha 0->*j* transition density matrix
   SharedMatrix Ta_ao(int j) const;

   /// Compute MO one-particle beta 0->*j* transition density matrix
   SharedMatrix Tb_ao(int j) const;

   /// Compute MO one-particle alpha *i*->*j* transition density matrix
   SharedMatrix Ta_ao(int i, int j) const;

   /// Compute MO one-particle beta *i*->*j* transition density matrix
   SharedMatrix Tb_ao(int i, int j) const;

   /// Compute TrCAMM for 0->*j* transition
   SharedDMTPole trcamm(int j, bool symmetrize=true) const;

   /// Compute TrCAMM for *i*->*j* transition
   SharedDMTPole trcamm(int i, int j, bool symmetrize=true) const;

   /// Compute transition dipole moment for 0->*j* transition
   SharedVector transition_dipole(int j) const;

   /// Compute transition dipole moment for *i*->*j* transition
   SharedVector transition_dipole(int i, int j) const;

   /// Compute oscillator strength for 0->*j* transition
   double oscillator_strength(int j) const;

   /// Compute oscillator strength for *i*->*j* transition
   double oscillator_strength(int i, int j) const;

   /// Compute <S2> expectation value for the *i*th state
   double s2(int i) const;

   /// Determine electronic state
   void determine_electronic_state(int& I);

   /// Return CIS data structure for a given excited state *I*
   std::shared_ptr<CISData> data(int I, int h, int l, bool symmetrize_trcamm=false);

   /// Slater determinant possible references, that are implemented
   static const std::vector<std::string> reference_types;

  protected:
   /// Reference wavefunction
   std::shared_ptr<psi::Wavefunction> ref_wfn_;
   /// Psi4 Options
 //psi::Options& options_; --> moved to Davidson-Liu

   /// Number of MO's
   const int nmo_;
   /// Number of alpha occupied MO's
   const int naocc_;
   /// Number of beta occupied MO's
   const int nbocc_;
   /// Number of alpha virtual MO's
   const int navir_;
   /// Number of beta virtual MO's
   const int nbvir_;
   /// Number of excited determinants
   int ndets_;
   /// Number of excited states
   int nstates_;

   /// CIS Excited State Hamiltonian in Slater determinantal basis
   SharedMatrix H_;
   // CIS Coefficients \f$ U_{uI} \f$ for each excited state *I* and basis Slater determinant *u*
   SharedMatrix U_;
   // Electronic excitation energies \f$ E_{I} \f$ wrt ground state
   SharedVector E_;

   // Fock matrices: OO, oo, VV and vv blocks
   //SharedMatrix Fa_oo_, Fb_oo_, Fa_vv_, Fb_vv_;

   /// Computer of generalized JK objects
   std::shared_ptr<psi::JK> jk_;

   // Canonical orbital energies
   SharedVector eps_a_o_, eps_a_v_, eps_b_o_, eps_b_v_;

   /// MO Integral Transformation Type
   const psi::IntegralTransform::TransformationType transformation_type_;
   std::shared_ptr<psi::IntegralTransform> inttrans_;

  // --> protected interface <-- //
  protected:

   CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt, psi::IntegralTransform::TransformationType trans_type);

   virtual void print_header_(void);
   virtual void set_nstates_(void);
   virtual void allocate_memory(void);
   virtual void allocate_hamiltonian_(void);
   virtual void prepare_for_cis_(void);
   virtual void build_hamiltonian_(void) = 0;
   virtual void diagonalize_hamiltonian_(void);
   virtual void standardize_amplitudes_(void);
   virtual void print_excited_states_(void);
   virtual void print_excited_state_character_(int I) = 0;

   virtual void set_beta_(void) = 0;
   virtual void transform_integrals_(void);

   virtual void davidson_liu_compute_diagonal_hamiltonian(void);
   virtual void davidson_liu_compute_sigma(void);

  // --> private interface <-- //
  private:
   void common_init(void);
   double compute_s2_reference(void);
   double compute_s2_exact(void);

};

class R_CISComputer: public CISComputer {
  public:
   R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~R_CISComputer(); 
  protected:
   virtual void print_excited_state_character_(int I);
};

class U_CISComputer: public CISComputer {
  public:
   U_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~U_CISComputer(); 
  protected:
   virtual void print_excited_state_character_(int I);
};



class R_CISComputer_Explicit: public R_CISComputer {
  public:
   R_CISComputer_Explicit(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~R_CISComputer_Explicit(); 
  protected:
   virtual void set_beta_(void);
   virtual void build_hamiltonian_(void);
};

/** \brief CIS Computer with RHF reference: Davidson-Liu Solver.
 *
 * Associated options:
 *  - `CIS_TYPE`             - must be set to `DAVIDSON_LIU` (Default).
 *  - `CIS_SCHWARTZ_CUTOFF`  - Cutoff for Schwartz ERI screening. Default: 0.0.
 *
 * # Implementation
 * ## Diagonal Hamiltonian elements
 *
 * They are computed by using direct method with Schwartz screening of AO ERI's.
 * The implementation formula is
 * \f[
 *    H_{ii}^{aa} = \varepsilon_a - \varepsilon_i + 
 *          \sum_{\alpha\beta\gamma\delta} 
 *         (\alpha\beta \vert \gamma\delta) C_{\alpha i} C_{\delta a} 
 *     \left( C_{\beta a} C_{\gamma i} - C_{\beta i} C_{\gamma a}\right)
 * \f]
 * The block associated with beta spin is equal to alpha block.
 *
 * ## Sigma vectors
 *
 * The sigma vectors are computed from
 * \f{align*}{
 *   \sigma_{i}^{a,k} &= (\varepsilon_a - \varepsilon_i) b_{i}^{a,k} 
 *        + J_{i}^{a}({\bf T}^{(k)})
 *        + J_{i}^{a}(\overline{{\bf T}^{(k)}})
 *        - K_{i}^{a}({\bf T}^{(k)}) \\
 *   \sigma_{\overline{i}}^{\overline{a},k} &= (\varepsilon_{a} - \varepsilon_i) b_{\overline{i}}^{\overline{a},k} 
 *        + J_{i}^{a}({\bf T}^{(k)})
 *        + J_{i}^{a}(\overline{{\bf T}^{(k)}})
 *        - K_{i}^{a}(\overline{{\bf T}^{(k)}}) 
 * \f}
 * where *k* labels the vectors and
 * where the generalized one-particle density matrices are defined by
 * \f{align*}{
 *  T_{\gamma\delta}^{(k)} &= \sum_{jb} C_{\delta b} b_{j}^{b,k} C_{\gamma j} \\
 *  \overline{T}_{\gamma\delta}^{(k)}
 *                   &= \sum_{\overline{j}\overline{b}} C_{\delta \overline{b}} 
 *                   b_{\overline{j}}^{\overline{b},k}  C_{\gamma \overline{j}} 
 * \f}
 * The **J** and **K** matrices in AO basis are computed by using the `psi::JK` object, and subsequently 
 * transformed to CMO's.
 *
 */
class R_CISComputer_DL: public R_CISComputer_Explicit {
  public:
   R_CISComputer_DL(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~R_CISComputer_DL(); 
  protected:
   virtual void set_nstates_(void);
   virtual void transform_integrals_(void);
   virtual void allocate_hamiltonian_(void);
   virtual void build_hamiltonian_(void);
   virtual void diagonalize_hamiltonian_(void);

   virtual void davidson_liu_compute_diagonal_hamiltonian(void);
   virtual void davidson_liu_compute_sigma(void);
  private:
   psi::SharedMatrix Ca_occ__, Ca_vir__;
};

class R_CISComputer_Direct: public R_CISComputer_Explicit {
  public:
   R_CISComputer_Direct(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~R_CISComputer_Direct(); 
  protected:
   virtual void build_hamiltonian_(void);
   virtual void transform_integrals_(void);
};

class U_CISComputer_Explicit: public U_CISComputer {
  public:
   U_CISComputer_Explicit(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~U_CISComputer_Explicit(); 
  protected:
   virtual void set_beta_(void);
   virtual void build_hamiltonian_(void);
};

/** \brief CIS Computer with UHF reference: Davidson-Liu Solver.
 *
 * Associated options:
 *  - `CIS_TYPE`             - must be set to `DAVIDSON_LIU` (Default).
 *  - `CIS_SCHWARTZ_CUTOFF`  - Cutoff for Schwartz ERI screening. Default: 0.0.
 *
 * # Implementation
 * ## Diagonal Hamiltonian elements
 *
 * They are computed by using direct method with Schwartz screening of AO ERI's.
 * The implementation formula is
 * \f{align*}{
 *    H_{ii}^{aa} &= \varepsilon_a - \varepsilon_i + 
 *          \sum_{\alpha\beta\gamma\delta} 
 *         (\alpha\beta \vert \gamma\delta) C_{\alpha i} C_{\delta a} 
 *     \left( C_{\beta a} C_{\gamma i} - C_{\beta i} C_{\gamma a}\right)  \\
 *    H_{\overline{i}\overline{i}}^{\overline{a}\overline{a}} &= \varepsilon_{\overline{a}} - \varepsilon_{\overline{i}} + 
 *          \sum_{\alpha\beta\gamma\delta} 
 *         (\alpha\beta \vert \gamma\delta) C_{\alpha \overline{i}} C_{\delta \overline{a}} 
 *     \left( C_{\beta \overline{a}} C_{\gamma \overline{i}} - C_{\beta \overline{i}} C_{\gamma \overline{a}}\right)  
 * \f}
 *
 * ## Sigma vectors
 *
 * The sigma vectors are computed from
 * \f{align*}{
 *   \sigma_{i}^{a,k} &= (\varepsilon_a - \varepsilon_i) b_{i}^{a,k} 
 *        + J_{i}^{a}({\bf T}^{(k)})
 *        + J_{i}^{a}(\overline{{\bf T}^{(k)}})
 *        - K_{i}^{a}({\bf T}^{(k)}) \\
 *   \sigma_{\overline{i}}^{\overline{a},k} &= (\varepsilon_{\overline{a}} - \varepsilon_{\overline{i}}) 
 *                b_{\overline{i}}^{\overline{a},k} 
 *        + J_{\overline{i}}^{\overline{i}}({\bf T}^{(k)})
 *        + J_{\overline{i}}^{\overline{i}}(\overline{{\bf T}^{(k)}})
 *        - K_{\overline{i}}^{\overline{i}}(\overline{{\bf T}^{(k)}}) 
 * \f}
 * where *k* labels the vectors and
 * where the generalized one-particle density matrices are defined by
 * \f{align*}{
 *  T_{\gamma\delta}^{(k)} &= \sum_{jb} C_{\delta b} b_{j}^{b,k} C_{\gamma j} \\
 *  \overline{T}_{\gamma\delta}^{(k)}
 *                   &= \sum_{\overline{j}\overline{b}} C_{\delta \overline{b}} 
 *                   b_{\overline{j}}^{\overline{b},k}  C_{\gamma \overline{j}} 
 * \f}
 * The **J** and **K** matrices in AO basis are computed by using the `psi::JK` object, and subsequently 
 * transformed to CMO's.
 *
 */

class U_CISComputer_DL: public U_CISComputer_Explicit {
  public:
   U_CISComputer_DL(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt);
   virtual ~U_CISComputer_DL(); 
  protected:
   virtual void set_nstates_(void);
   virtual void transform_integrals_(void);
   virtual void allocate_hamiltonian_(void);
   virtual void build_hamiltonian_(void);
   virtual void diagonalize_hamiltonian_(void);

   virtual void davidson_liu_compute_diagonal_hamiltonian(void);
   virtual void davidson_liu_compute_sigma(void);
  private:
   psi::SharedMatrix Ca_occ__, Ca_vir__, Cb_occ__, Cb_vir__;
};



/** @}*/

} // EndNameSpace oepdev


#endif // _oepdev_libutil_cis_h
