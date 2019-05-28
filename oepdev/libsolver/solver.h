#ifndef _oepdev_libutil_solver_h
#define _oepdev_libutil_solver_h
/** @file solver.h */

#include<cstdio>
#include<string>
#include<map>

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "../libutil/wavefunction_union.h"
#include "../libutil/integrals_iter.h"
#include "../libpsi/integral.h"
#include "../liboep/oep.h"

namespace oepdev{

using namespace std;
using namespace psi;

using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedWavefunctionUnion  = std::shared_ptr<WavefunctionUnion>;
using SharedOEPotential        = std::shared_ptr<OEPotential>;
/** \addtogroup OEPDEV_SOLVERS
 * @{
 */

/**\brief Solver of properties of molecular aggregates. Abstract base.
 *
 * Uses only a wavefunction union object to initialize.
 * Available solvers:
 *  - `ELECTROSTATIC ENERGY`
 *  - `REPULSION ENERGY`
 *  - `CHARGE TRANSFER ENERGY`
 *
 * Options controlling the generalized density fitting (GDF):
 *  - `OEPDEV_DF_TYPE` - type of the GDF. Default: `DOUBLE`.
 *  - `DF_BASIS_OEP` - auxiliary basis set. Default: `sto-3g`.
 *  - `DF_BASIS_INT` - intermediate basis set. Relevant only if double GDF is used. Default: `aug-cc-pVDZ-jkfit`. 
 *     Note that intermediate basis set should be nearly complete.
 */
class OEPDevSolver : public std::enable_shared_from_this<OEPDevSolver>
{
 public:

   // <--- Constructor and Destructor ---> //
   
   /**\brief Take wavefunction union and initialize the Solver.
    * 
    * @param wfn_union   - wavefunction union of isolated molecular wavefunctions
    */
   OEPDevSolver(SharedWavefunctionUnion wfn_union);

   /// Destroctor
   virtual ~OEPDevSolver();


   // <--- Factory Methods ---> //

   /**\brief Build a solver of a particular property for given molecular cluster.
    *
    * @param target         - target property
    * @param wfn_union      - wavefunction union of isolated molecular wavefunctions
    *
    * Implemented target properties:
    *  - `ELECTROSTATIC_ENERGY` - Coulombic interaction energy between unperturbed wavefunctions.
    *  - `REPULSION_ENERGY`     - Pauli repulsion interaction energy between unperturbed wavefunctions.
    *
    * \see ElectrostaticEnergySolver
    */
   static std::shared_ptr<OEPDevSolver> build(const std::string& target, SharedWavefunctionUnion wfn_union);

  
   // <--- Computers ---> //

   /**\brief Compute property by using OEP's
    * 
    *  Each solver object has one `DEFAULT` OEP-based method. 
    *  @param method - flavour of OEP model
    */
   virtual double compute_oep_based(const std::string& method = "DEFAULT") = 0;

   /**\brief Compute property by using benchmark method
    *
    *  Each solver object has one `DEFAULT` benchmark method
    *  @param method - benchmark method
    */
   virtual double compute_benchmark(const std::string& method = "DEFAULT") = 0;


 protected:

   /// Wavefunction union
   SharedWavefunctionUnion wfn_union_;

   /// Options
   psi::Options& options_;

   /// Names of all OEP-based methods implemented for a solver
   std::vector<std::string> methods_oepBased_;

   /// Names of all benchmark methods implemented for a solver
   std::vector<std::string> methods_benchmark_;

};

/**\brief Compute the Coulombic interaction energy between unperturbed wavefunctions.
 *
 * The implemented methods are shown in below
 * <table>
 * <caption id="Tab.1">Methods available in the Solver</caption>
 * <tr><th> Keyword  <th>Method Description  
 * <tr><td colspan=2> <center><strong>Benchmark Methods</strong></center>
 * <tr><td> `AO_EXPANDED`  <td>*Default*. Exact Coulombic energy from atomic orbital expansions.
 * <tr><td> `MO_EXPANDED`  <td>Exact Coulombic energy from molecular orbital expansions
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `ESP_SYMMETRIZED` <td>*Default*. Coulombic energy from ESP charges interacting with nuclei and electronic density.
 *                                Symmetrized with respect to monomers.
 * <tr><td> `CAMM` <td>Coulombic energy from CAMM distributions.
 * </table>
 *
 * Below the detailed description of the above methods is given.
 *
 * # Benchmark Methods
 * ## Exact Coulombic energy from atomic orbital expansions.
 *    
 * The Coulombic interaction energy is given by
 * \f[
 *    E^{\rm Coul} = E^{\rm Nuc-Nuc} + E^{\rm Nuc-El} + E^{\rm El-El}
 * \f]
 * where the nuclear-nuclear repulsion energy is
 * \f[
 *     E^{\rm Nuc-Nuc} = \sum_{x\in A}\sum_{y\in B} \frac{Z_xZ_y}{\lvert {\bf r}_x - {\bf r}_y\rvert}
 * \f]
 * the nuclear-electronic attraction energy is
 * \f[
 *     E^{\rm Nuc-El } = \sum_{x\in A}\sum_{\lambda\sigma\in B} Z_x V_{\lambda\sigma}^{(x)} 
 *                       \left(D_{\lambda\sigma}^{(\alpha)} + D_{\lambda\sigma}^{(\beta)}\right)
 *                     + \sum_{y\in B}\sum_{\mu    \nu   \in A} Z_y V_{\mu    \nu   }^{(y)} 
 *                       \left(D_{\mu    \nu   }^{(\alpha)} + D_{\mu    \nu   }^{(\beta)}\right)
 * \f]
 * and the electron-electron repulsion energy is
 * \f[
 *     E^{\rm El-El  } = \sum_{\mu    \nu   \in A} \sum_{\lambda\sigma\in B}
 *                       \left\{ D_{\mu    \nu   }^{(\alpha)} + D_{\mu    \nu   }^{(\beta)} \right\}
 *                       \left\{ D_{\lambda\sigma}^{(\alpha)} + D_{\lambda\sigma}^{(\beta)} \right\}
 *                       \left( \mu\nu \vert \lambda\sigma \right)
 * \f]
 * In the above equations, 
 * \f[
 *    V_{\lambda\sigma}^{(x)} \equiv \int \frac{\varphi_\lambda^{*}({\bf r}) \varphi_\sigma({\bf r})}
 *                                             {\lvert {\bf r} - {\bf r}_x\rvert} d{\bf r} 
 * \f]
 *
 * ## Exact Coulombic energy from molecular orbital expansion.
 * 
 * This approach is fully equivalent to the atomic orbital expansion shown above.
 * For the closed shell case, the Coulombic interaction energy is given by
 * \f[
 *    E^{\rm Coul} = E^{\rm Nuc-Nuc} + E^{\rm Nuc-El} + E^{\rm El-El}
 * \f]
 * where the nuclear-nuclear repulsion energy is
 * \f[
 *     E^{\rm Nuc-Nuc} = \sum_{x\in A}\sum_{y\in B} \frac{Z_xZ_y}{\lvert {\bf r}_x - {\bf r}_y\rvert}
 * \f]
 * the nuclear-electronic attraction energy is
 * \f[
 *     E^{\rm Nuc-El } = 2\sum_{i\in A} \sum_{y\in B} V_{ii}^{(y)} +
 *                       2\sum_{j\in B} \sum_{x\in A} V_{jj}^{(x)}
 * \f]
 * and the electron-electron repulsion energy is
 * \f[
 *     E^{\rm El-El  } = 4\sum_{i\in A}\sum_{j\in B} (ii \vert jj)
 * \f]
 *
 * # OEP-Based Methods
 * ## Coulombic energy from ESP charges interacting with nuclei and electronic density.
 * 
 * In this approach, nuclear and electronic density of either species is approximated
 * by ESP charges. In order to achieve symmetric expression, the interaction
 * is computed twice (ESP of A interacting with density matrix and nuclear charges of B and vice versa)
 * and then divided by 2. Thus,
 * \f[
 *     E^{\rm Coul} \approx \frac{1}{2} 
 *                          \left[
 *                          \sum_{x\in A}\sum_{y\in B} \frac{Z_xq_y}{\lvert {\bf r}_x - {\bf r}_y\rvert}
 *                        + \sum_{y\in B}\sum_{\mu    \nu   \in A} q_y V_{\mu    \nu   }^{(y)} 
 *                          \left(D_{\mu    \nu   }^{(\alpha)} + D_{\mu    \nu   }^{(\beta)}\right)
 *                        + \sum_{y\in B}\sum_{x\in A} \frac{q_xZ_y}{\lvert {\bf r}_x - {\bf r}_y\rvert}
 *                        + \sum_{x\in A}\sum_{\lambda\sigma\in B} q_x V_{\lambda\sigma}^{(x)} 
 *                          \left(D_{\lambda\sigma}^{(\alpha)} + D_{\lambda\sigma}^{(\beta)}\right)
 *                        \right]
 * \f]
 * If the basis set is large and the number of ESP centres \f$ q_{x(y)} \f$ is sufficient,
 * the sum of first two contributions equals the sum of the latter two contributions.
 *
 * *Notes:* 
 *   - This solver also computes and prints the ESP-ESP point charge interaction energy,
 *     \f[
 *       E^{\rm Coul, ESP} \approx \sum_{x\in A}\sum_{y\in B} \frac{q_xq_y}{\lvert {\bf r}_x - {\bf r}_y\rvert}
 *     \f]
 *     for reference purposes.
 *   - In order to construct this solver, **always** use the `OEPDevSolver::build` static factory method.
 */
class ElectrostaticEnergySolver : public OEPDevSolver
{
  public:
    ElectrostaticEnergySolver(SharedWavefunctionUnion wfn_union);
    virtual ~ElectrostaticEnergySolver();

    virtual double compute_oep_based(const std::string& method = "DEFAULT");
    virtual double compute_benchmark(const std::string& method = "DEFAULT");

  private:
    double compute_oep_based_esp_symmetrized();
    double compute_oep_based_camm();
    double compute_benchmark_ao_expanded(); 
    double compute_benchmark_mo_expanded();
};

/**\brief Compute the Pauli-Repulsion interaction energy between unperturbed wavefunctions.
 *
 * The implemented methods are shown below
 * <table>
 * <caption id="Tab.1">Methods available in the Solver</caption>
 * <tr><th> Keyword  <th>Method Description  
 * <tr><td colspan=2> <center><strong>Benchmark Methods</strong></center>
 * <tr><td> `HAYES_STONE`       <td>*Default*. Pauli Repulsion energy at HF level from Hayes and Stone (1984). 
 * <tr><td> `DENSITY_BASED`     <td>Pauli Repulsion energy at HF level from Mandado and Hermida-Ramon (2012).
 * <tr><td> `MURRELL_ETAL`      <td>Approximate Pauli Repulsion energy at HF level from Murrell et al (1967).
 * <tr><td> `OTTO_LADIK`        <td>Approximate Pauli Repulsion energy at HF level from Otto and Ladik (1975).
 * <tr><td> `EFP2`              <td>Approximate Pauli Repulsion energy at HF level from EFP2 model.
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `MURRELL_ETAL_GDF_ESP`  <td>*Default*. OEP-Murrell et al's: S1 term via DF-OEP, S2 term via ESP-OEP.
 * <tr><td> `MURRELL_ETAL_GDF_CAMM` <td>OEP-Murrell et al's: S1 term via DF-OEP, S2 term via CAMM-OEP.
 * <tr><td> `MURRELL_ETAL_ESP`  <td>OEP-Murrell et al's: S1 and S2 via ESP-OEP (not implemented)
 * </table>
 * *Note:*
 *   - This solver also computes and prints the exchange energy at HF level (formulae are given below)
 *     for reference purposes.
 *   - In order to construct this solver, **always** use the `OEPDevSolver::build` static factory method.
 *
 * Below the detailed description of the implemented equations is given 
 * for each of the above provided methods.
 * In the formulae across, it is assumed that the orbitals are real.
 * The Coulomb notation for 
 * electron repulsion integrals (ERI's) is adopted; i.e,
 * \f[
 *  (ac \vert bd) = \iint d{\bf r}_1 d{\bf r}_2 
 *   \phi_a({\bf r}_1) \phi_c({\bf r}_1) \frac{1}{r_{12}} \phi_b({\bf r}_2) \phi_d({\bf r}_2)
 * \f]
 * Greek subscripts denote basis set orbitals whereas Italic subscripts denote the occupied
 * molecular orbitals.
 *
 * # Benchmark Methods
 * ## Pauli Repulsion energy at HF level by Hayes and Stone (1984).
 *    
 * For a closed-shell system, equation of Hayes and Stone (1984)
 * becomes
 * \f[
 *    E^{\rm Rep} = 2\sum_{kl} 
                    \left( V^A_{kl} + V^B_{kl} + T_{kl} \right) 
                    \left[ [{\bf S}^{-1}]_{lk} - \delta_{lk} \right]
                +   \sum_{klmn} 
                    (kl \vert mn) 
                    \left\{ 
      2[{\bf S}^{-1}]_{kl} [{\bf S}^{-1}]_{mn} - 
       [{\bf S}^{-1}]_{kn} [{\bf S}^{-1}]_{lm} -
      2\delta_{kl} \delta_{mn} +
       \delta_{kn} \delta_{lm}
                    \right\}
 * \f]
 * where \f$ {\bf S} \f$ is the overlap matrix between the doubly-occupied
 * orbitals.
 * The exact, pure exchange energy is for a closed shell case given as
 * \f[
     E^{\rm Ex,pure} = -2\sum_{a\in A} \sum_{b\in B} (ab \vert ba)
 * \f]
 * Similarity transformation of molecular orbitals does not affect the resulting energies.
 * The overall exchange-repulsion interaction energy is then (always net repulsive)
 * \f[ 
 *  E^{\rm Ex-Rep} = E^{\rm Ex,pure} + E^{\rm Rep} 
 * \f]
 *
 * ## Repulsion energy of Mandado and Hermida-Ramon (2011)
 *
 * At the Hartree-Fock level, the exchange-repulsion energy from the density-based scheme 
 * of Mandado and Hermida-Ramon (2011) is fully equivalent to the method by Hayes and Stone (1984).
 * However, density-based method enables to compute exchange-repulsion energy 
 * at any level of theory. It is derived based on the Pauli deformation density matrix,
 * \f[
 *   \Delta {\bf D}^{\rm Pauli} \equiv {\bf D}^{oo} - {\bf D}
 * \f]
 * where \f$ {\bf D}^{oo}\f$ and \f$ {\bf D}\f$ are the density matrix formed from
 * mutually orthogonal sets of molecular orbitals within the entire aggregate (formed
 * by symmetric orthogonalization of MO's) and the density matrix of the unperturbed
 * system (that can be understood as a Hadamard sum \f$ {\bf D} \equiv {\bf D}^A \oplus {\bf D}^B\f$).
 *
 * At HF level, the Pauli deformation density matrix is given by
 * \f[
 *   \Delta {\bf D}^{\rm Pauli} = {\bf C} \left[ {\bf S}^{-1} - {\bf 1} \right] {\bf C}^\dagger
 * \f]
 * whereas the density matrix constructed from mutually orthogonal orbitals is
 * \f[
 *    {\bf D}^{oo} = {\bf C} {\bf S}^{-1} {\bf C}^\dagger
 * \f]
 * In the above equations, \f$ {\bf S} \f$ is the overlap matrix between doubly occupied molecular orbitals
 * of the entire aggregate.
 * 
 * Here, the expressions for the exchange-repulsion energy
 * at any level of theory are shown for the case of open-shell system.
 * The net repulsive energy is given as
 * \f[
     E^{\rm Ex-Rep} = E^{\rm Rep,1} + E^{\rm Rep,2} + E^{\rm Ex}
 * \f]
 * where the one- and two-electron part of the repulsion energy is
 * \f{align*}{
      E^{\rm Rep,1} &= E^{\rm Rep,Kin} + E^{\rm Rep,Nuc} \\
      E^{\rm Rep,2} &= E^{\rm Rep,el-\Delta} + E^{\rm Rep,\Delta-\Delta}
 * \f}
 * The kinetic and nuclear contributions are
 * \f{align*}{
     E^{\rm Rep,Kin} &= 2\sum_{\alpha\beta \in A,B} \Delta D_{\alpha\beta}^{\rm Pauli} T_{\alpha\beta}  \\
     E^{\rm Rep,Nuc} &= 2\sum_{\alpha\beta \in A,B} \Delta D_{\alpha\beta}^{\rm Pauli} \sum_{z\in A,B} V_{\alpha\beta}^{(z)} 
 * \f}
 * whereas the electron-deformation and deformation-deformation interaction contributions are
 * \f{align*}{
     E^{\rm Rep,el-\Delta} &= 4\sum_{\alpha\beta\gamma\delta \in A,B}
                              \Delta D_{\alpha\beta}^{\rm Pauli} D_{\gamma\delta}
                              (\alpha\beta \vert \gamma\delta) \\
     E^{\rm Rep,\Delta-\Delta} &= 2\sum_{\alpha\beta\gamma\delta \in A,B}
                              \Delta D_{\alpha\beta}^{\rm Pauli} \Delta D_{\gamma\delta}^{\rm Pauli}
                              (\alpha\beta \vert \gamma\delta) 
 * \f}
 * The associated exchange energy is given by
 * \f[
 *   E^{\rm Ex} = -\sum_{\alpha\beta\gamma\delta \in A,B} 
 *                \left[
 *                 D_{\alpha\delta}^{oo} D_{\beta\gamma}^{oo}
 *                -D_{\alpha\delta}^A    D_{\beta\gamma}^A
 *                -D_{\alpha\delta}^B    D_{\beta\gamma}^B
 *                \right]
 *                (\alpha\beta \vert \gamma\delta)
 * \f]
 * It is important to emphasise that, although, at HF level, the particular 
 * 'repulsive' and 'exchange' energies
 * computed by using either Hayes and Stone or Mandado and Hermida-Ramon methods
 * are not equal to each other, they sum up to exactly the same exchange-repulsion
 * energy, \f$ E^{\rm Ex-Rep}\f$. Therefore, these methods at HF level are fully equivalent but the 
 * nature of partitioning of repulsive and exchange parts is different. It is also
 * noted that the orbital localization does *not* affect the resulting energies, as 
 * opposed to the few approximate methods described below (Otto-Ladik and EFP2 methods).
 * 
 * ## Approximate Pauli Repulsion energy at HF level from Murrell et al.
 * 
 * By expanding the overlap matrix in a Taylor series one can show that 
 * the Pauli repulsion energy is approximately given as
 * \f[
 *    E^{\rm Rep} = E^{\rm Rep}(\mathcal{O}(S)) + E^{\rm Rep}(\mathcal{O}(S^2))
 * \f]
 * where the first-order term is 
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S)) = -2\sum_{a\in A} \sum_{b\in B}
 *                S_{ab} \left\{
 *            V^A_{ab} + \sum_{c\in A} \left[ 2(ab \vert cc) - (ac \vert bc) \right]
 *          + V^B_{ab} + \sum_{d\in B} \left[ 2(ab \vert dd) - (ad \vert bd) \right]
 *                 \right\}
 * \f]
 * whereas the second-order term is
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S^2)) = 2\sum_{a\in A} \sum_{b\in B} S_{ab} \left\{
 *                  \sum_{c\in A} S_{bc}
 *                \left[ V_{ac}^B + 2\sum_{d\in B} (ac \vert dd) \right]
 *           +      \sum_{d\in B} S_{ad}
 *                \left[ V_{bd}^A + 2\sum_{x\in A} (bd \vert cc) \right]
 *                - \sum_{c\in A} \sum_{d\in B} S_{cd} (ac \vert bd)
 *         \right\}
 * \f] 
 * Thus derived repulsion energy is invariant with respect to transformation of molecular
 * orbitals, similarly as Hayes-Stone's method and density-based method.
 * By using OEP technique, the above theory can be exactly re-cast *without* any further approximations.
 *
 * ## Approximate Pauli Repulsion energy at HF level from Otto and Ladik (1975).
 * 
 * The Pauli repulsion energy is approximately given as
 * \f[
 *    E^{\rm Rep} = E^{\rm Rep}(\mathcal{O}(S)) + E^{\rm Rep}(\mathcal{O}(S^2))
 * \f]
 * where the first-order term is 
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S)) = -2\sum_{a\in A} \sum_{b\in B}
 *                S_{ab} \left\{
 *            V^A_{ab} + 2\sum_{c\in A} (ab \vert cc) - (ab \vert aa)
 *          + V^B_{ab} + 2\sum_{d\in B} (ab \vert dd) - (ab \vert bb)
 *                 \right\}
 * \f]
 * whereas the second-order term is
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S^2)) = 2\sum_{a\in A} \sum_{b\in B} S_{ab}^2 
 *         \left\{
 *                  V_{aa}^B + V_{bb}^A 
 *                 + 2\sum_{c\in A} (cc \vert bb) 
 *                 + 2\sum_{d\in B} (aa \vert dd)
 *                 -                (aa \vert bb)
 *         \right\}
 * \f] 
 * Thus derived repulsion energy is *not* invariant with respect to transformation of molecular
 * orbitals, in contrast to Hayes-Stone's method and density-based method.
 * It was shown that good results are obtained when using localized molecular orbitals,
 * whereas using canonical molecular orbitals brings poor results.
 * By using OEP technique, the above theory can be exactly re-cast *without* any further approximations.
 *
 * ## Approximate Pauli Repulsion energy at HF level from Jensen and Gordon (1996).
 * 
 * The Pauli repulsion energy used within the EFP2 approach is approximately given as
 * \f[
 *    E^{\rm Rep} = E^{\rm Rep}(\mathcal{O}(S)) + E^{\rm Rep}(\mathcal{O}(S^2))
 * \f]
 * where the first-order term is 
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S)) = -2\sum_{a\in A} \sum_{b\in B}
 *                S_{ab} \left\{
 *            \sum_{c\in A} F^A_{ac} S_{cb} 
 *          + \sum_{d\in B} F^B_{bd} S_{da}
 *          - 2 T_{ab}
 *                 \right\}
 * \f]
 * whereas the second-order term is
 * \f[
 *    E^{\rm Rep}(\mathcal{O}(S^2)) = 2\sum_{a\in A} \sum_{b\in B} S_{ab}^2 
 *         \left\{
 *                   \sum_{x\in A} \frac{-Z_x}{R_{xb}}
 *                 + \sum_{y\in B} \frac{-Z_y}{R_{ya}}
 *                 + \sum_{c\in A} \frac{2}{R_{bc}}
 *                 + \sum_{d\in B} \frac{2}{R_{ad}}
 *                 -\frac{1}{R_{ab}}
 *         \right\}
 * \f] 
 * Thus derived repulsion energy is *not* invariant with respect to transformation of molecular
 * orbitals, in contrast to Hayes-Stone's method and density-based method.
 * It was shown that good results are obtained when using localized molecular orbitals,
 * whereas using canonical molecular orbitals brings poor results.
 * 
 * In EFP2, exchange energy is approximated by spherical Gaussian approximation (SGO).
 * The result of this is the following formula for the exchange energy:
 * \f[
 *  E^{\rm Ex} \approx -4 \sum_{a\in A} \sum_{b\in B}
            \sqrt{\frac{-2 \ln{\vert S_{ab} \vert} }{\pi}}
            \frac{S_{ab}^2}{R_{ab}}
 * \f]
 * In all the above formulas, \f$ R_{ij} \f$ are distances between position vectors
 * of *i*th and *j*th point. The LMO centroids are defined by
 * \f[
 *   {\bf r}_a = (a \vert {\bf r} \vert a)
 * \f]
 * where *a* denotes the occupied molecular orbital.
 *
 * # OEP-Based Methods
 *
 * The Murrell et al's theory of Pauli repulsion for S-1 term
 * and the Otto-Ladik's theory for S-2 term is here re-cast by introducing OEP's.
 * The S-1 term is expressed via DF-OEP, whereas the S-2 term via ESP-OEP.
 *
 * ## S-1 term (Murrell et al.)
 *
 * The OEP reduction without any approximations leads to the following formula
 * \f[
 *  E^{\rm Rep}(\mathcal{O}(S^{1})) = -2\sum_{a\in A} \sum_{b\in B}
 *                S_{ab}
 *                \left\{
 *                        \sum_{\xi \in A} S_{b \xi } G_{\xi  a}^A
 *                      + \sum_{\eta\in B} S_{a \eta} G_{\eta b}^B
 *                \right\}
 * \f]
 * where the OEP matrices are given as
 * \f[
 *  G_{\xi a}^A = \sum_{\xi' \in A} [{\bf S}^{-1}]_{\xi\xi'} 
 *                \sum_{\alpha\in A} \left\{
 *                  C_{\alpha a} V_{\alpha\xi'}^A 
 *                + \sum_{\mu\nu\in A}
 *                  \left[
 *                   2C_{\alpha a}D_{\mu\nu} - C_{\nu a}D_{\alpha\mu}
 *                  \right] 
 *                  (\alpha\xi'\vert\mu\nu)
 *                                \right\}
 * \f]
 * and analogously for molecule *B*.
 * Here, the nuclear attraction integrals are denoted by \f$ V_{\alpha\xi'}^A\f$.
 *
 * ## S-2 term (Otto-Ladik)
 *
 * After the OEP reduction, this contribution under Otto-Ladik approximation has the following form:
 * \f[
 *  E^{\rm Rep}(\mathcal{O}(S^{2})) = 2\sum_{a\in A} \sum_{b\in B}
 *                S_{ab}^2 
 *                \left\{
 *                        \sum_{x\in A} q_{xa} V^{(x)}_{bb} 
 *                      + \sum_{y\in B} q_{yb} V^{(y)}_{aa} 
 *                \right\}
 * \f]
 * where the ESP charges associated with each occupied molecular orbital 
 * reproduce the *effective potential*
 * of molecule in question, i.e.,
 * \f[ 
 *   \sum_{x\in A} \frac{q_{xa}}{\vert {\rm r} - {\rm r}_x \vert} \cong v_a^A({\bf r})
 * \f]
 * where the potential is given by
 * \f[
 *    v_a^A({\bf r}) = \sum_{x\in A} 
 *                      \frac{-Z_x}{\vert {\rm r} - {\rm r}_x \vert}
 *                   + 2\sum_{c\in A}\int\frac{\phi_c({\rm r}')\phi_c({\rm r}')}{\vert {\rm r} - {\rm r}' \vert}\; d{\bf r}'
 *          -\frac{1}{2}             \int\frac{\phi_a({\rm r}')\phi_a({\rm r}')}{\vert {\rm r} - {\rm r}' \vert}\; d{\bf r}'
 * \f]
 *
 */
class RepulsionEnergySolver : public OEPDevSolver
{
  public:
    RepulsionEnergySolver(SharedWavefunctionUnion wfn_union);
    virtual ~RepulsionEnergySolver();

    virtual double compute_oep_based(const std::string& method = "DEFAULT");
    virtual double compute_benchmark(const std::string& method = "DEFAULT");

  private:
    /// Hayes-Stone (1984) method
    double compute_benchmark_hayes_stone();
    /// Mandado and Hermida-Ramon (2012)
    double compute_benchmark_density_based();
    /// Murrell et al's method (1967)
    double compute_benchmark_murrell_etal();
    /// Otto-Ladik method (1975)
    double compute_benchmark_otto_ladik();
    /// EFP2 method (1996)
    double compute_benchmark_efp2();

    /// Murrell et al's/OEP: S1-DF, S2-ESP (2017-2018)
    double compute_oep_based_murrell_etal_gdf_esp();
    /// Murrell et al's/OEP: S1-DF, S2-CAMM (2017-2018)
    double compute_oep_based_murrell_etal_gdf_camm();
    /// Murrell et al's/OEP: S1-ESP, S2-ESP (2017-2018)
    double compute_oep_based_murrell_etal_esp();

    /// Exchange energy at HF level
    double compute_pure_exchange_energy();
    /// Exchange energy at EFP2 level (SGO-approximated version of the above)
    double compute_efp2_exchange_energy(psi::SharedMatrix, std::vector<psi::SharedVector>, std::vector<psi::SharedVector>);
};

/**\brief Compute the Charge-Transfer interaction energy between unperturbed wavefunctions.
 *
 * The implemented methods are shown below
 * <table>
 * <caption id="Tab.1">Methods available in the Solver</caption>
 * <tr><th> Keyword  <th>Method Description  
 * <tr><td colspan=2> <center><strong>Benchmark Methods</strong></center>
 * <tr><td> `OTTO_LADIK`        <td>*Default*. CT energy at HF level from Otto and Ladik (1975). 
 * <tr><td> `EFP2`              <td>CT energy at HF level from EFP2 model.
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `OTTO_LADIK`       <td>*Default*. OEP-based Otto-Ladik expressions.
 * </table>
 *
 *   > In order to construct this solver, **always** use the `OEPDevSolver::build` static factory method.
 *
 * Below the detailed description of the implemented equations is given 
 * for each of the above provided methods.
 * In the formulae across, it is assumed that the orbitals are real.
 * The Coulomb notation for 
 * electron repulsion integrals (ERI's) is adopted; i.e,
 * \f[
 *  (ac \vert bd) = \iint d{\bf r}_1 d{\bf r}_2 
 *   \phi_a({\bf r}_1) \phi_c({\bf r}_1) \frac{1}{r_{12}} \phi_b({\bf r}_2) \phi_d({\bf r}_2)
 * \f]
 * Greek subscripts denote basis set orbitals whereas Italic subscripts denote the occupied
 * molecular orbitals.
 *
 * The CT energy between molecules *A* and *B* is given by
 * \f[
 *   E^{\rm CT} = E^{\rm A^+B^-} + E^{\rm A^-B^+}
 * \f]
 * 
 * # Benchmark Methods
 * ## CT energy at HF level by Otto and Ladik (1975).
 *    
 * For a closed-shell system, CT energy equation of Otto and Ladik
 * becomes
 * \f[
 *   E^{\rm A^+B^-} \approx 
 *                2\sum_{i\in A}^{\rm Occ_A}\sum_{n\in B}^{\rm Vir_B}
 *     \frac{V^2_{in}}{\varepsilon_i - \varepsilon_n}
 * \f]
 * where
 * \f{multline*}{
 *   V_{in} = V^B_{in} + 2\sum_{j\in B}^{\rm Occ_B} (in \vert jj)
 *           -\sum_{k\in A}^{\rm Occ_A} S_{kn} 
 *             \left\{
 *                    V^B_{ik} + 2\sum_{j\in B}^{\rm Occ_B} (ik \vert jj)
 *             \right\} \\
 *           -\sum_{j\in B}^{\rm Occ_B} \left[ S_{ij}
 *            \left\{
 *               V^A_{nj} + 2\sum_{k\in A}^{\rm Occ_A} (1-\delta_{ik})(nj \vert kk)
 *            \right\}
 *            -(nj \vert ij)
 *            \right]
 *            + \sum_{k\in A}^{\rm Occ_A}\sum_{j\in B}^{\rm Occ_B}
 *            S_{kj}(1-\delta_{ik}) (ik \vert nj)
 * \f}
 * and analogously the twin term.
 *
 * ## CT energy at HF level by EFP2.
 *
 * In EFP2 method, CT energy is given as
 * \f[
 *   E^{\rm A^+B^-} \approx 
 *                2\sum_{i\in A}^{\rm Occ_A}\sum_{n\in B}^{\rm Vir_B}
 *     \frac{V^2_{in}}{F_{ii} - T_{nn}}
 * \f]
 * where 
 * \f[
 *   V^2_{in} = \frac{V^B_{in} - \sum_{m\in A}^{\rm All_A} V_{im}S_{mn}^B}
 *                   {1-\sum_{m\in A}^{\rm All_A} S_{mn}^2}
 *          \left\{
 *             V^B_{in} - \sum_{m\in A}^{\rm All_A} V_{im}^BS_{mn}
 *            + \sum_{j\in B}^{\rm Occ_B} S_{ij} \left( T_{nj} - \sum_{m\in A}^{\rm All_A} S_{nm} T_{mj} \right)
 *          \right\}
 * \f]
 * and analogously the twin term.
 *
 * 
 * # OEP-Based Methods
 * ## OEP-Based Otto-Ladik's theory
 * 
 * After introducing OEP's, the original Otto-Ladik's theory is reformulated *without*
 * approximation as
 * \f[
 *   E^{\rm A^+B^-} \approx 
 *                2\sum_{i\in A}^{\rm Occ_A}\sum_{n\in B}^{\rm Vir_B}
 *     \frac{\left( V_{in}^{\rm DF} + V_{in}^{\rm ESP,A} + V_{in}^{\rm ESP,B}\right)^2}{\varepsilon_i - \varepsilon_n}
 * \f]
 * where
 * \f{align*}{
 *   V_{in}^{\rm DF}    &= \sum_{\eta\in B}^{\rm Aux_B} S_{i\eta} G_{\eta n}^B \\
 *   V_{in}^{\rm ESP,A} &= \sum_{k\in A}^{\rm Occ_A} \sum_{j\in B}^{\rm Occ_B} S_{kj}
 *                         \sum_{x\in A} V_{nj}^{(x)} q_{ik}^{(x)} \\
 *   V_{in}^{\rm ESP,B} &=-\sum_{k\in A}^{\rm Occ_A} S_{kn} V_{ik}^B
 * \f}
 * The OEP matrix for density fitted part is given by
 * \f[
 *   G^B_{\eta n} = \sum_{\eta'\in B}^{\rm Aux_B}
 *                [{\bf S}^{-1}]_{\eta\eta'} 
 *                \left\{
 *                 V^B_{\eta' n} + \sum_{j\in B}^{\rm Occ_B}
 *                    \left[
 *                       2(\eta' n \vert jj) - (\eta' j \vert n j)
 *                    \right]
 *                \right\}
 * \f] 
 * The OEP ESP-A charges are fit to reproduce the OEP potential
 * \f[
 *    v_{ik}^A({\bf r}) \equiv 
 *      (1-\delta_{ik}) \int \frac{\phi_i({\bf r}')\phi_k({\bf r}')}{\vert {\bf r} - {\bf r}' \vert} \; d{\bf r}'
 *       -\delta_{ik} \left(
 *        \sum_{x\in A} \frac{-Z_x}{\vert {\bf r} - {\bf r}_x \vert}
 *        + 2\sum_{k\in A}^{\rm Occ_A} 
 *         \int \frac{\phi_k({\bf r}')\phi_k({\bf r}')}{\vert {\bf r} - {\bf r}' \vert} \; d{\bf r}'
 *        - 2                          
 *         \int \frac{\phi_i({\bf r}')\phi_i({\bf r}')}{\vert {\bf r} - {\bf r}' \vert} \; d{\bf r}'
 *       \right)
 * \f]
 * so that 
 * \f[
 *   v_{ik}^A({\bf r}) \cong \sum_{x\in A} \frac{q_{ik}^{(x)}}{\vert {\bf r} - {\bf r}_x \vert}
 * \f]
 * The OEP ESP-B charges are fit to reproduce the electrostatic potential
 * of molecule *B* (they are standard ESP charges).
 */
class ChargeTransferEnergySolver : public OEPDevSolver
{
  public:
    ChargeTransferEnergySolver(SharedWavefunctionUnion wfn_union);
    virtual ~ChargeTransferEnergySolver();

    virtual double compute_oep_based(const std::string& method = "DEFAULT");
    virtual double compute_benchmark(const std::string& method = "DEFAULT");

  private:
    /// Otto-Ladik method (1975)
    double compute_benchmark_otto_ladik();
    /// EFP2 method (1996)
    double compute_benchmark_efp2();
    /// OEP-based Otto-Ladik (2019)
    double compute_oep_based_murrell_etal();
    /// CT energy component
    double compute_ct_component(std::shared_ptr<psi::Vector> eps_occ_X, 
                                std::shared_ptr<psi::Vector> eps_vir_Y, std::shared_ptr<Matrix> V);
};


/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_solver_h

