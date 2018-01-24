#ifndef _oepdev_libutil_solver_h
#define _oepdev_libutil_solver_h

#include<cstdio>
#include<string>
#include<map>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"

#include "wavefunction_union.h"
#include "integrals_iter.h"
#include "../liboep/oep.h"

namespace oepdev{

using namespace std;
using namespace psi;

using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedWavefunctionUnion  = std::shared_ptr<WavefunctionUnion>;
using SharedOEPotential        = std::shared_ptr<OEPotential>;


/**\brief Solver of properties of molecular aggregates. Abstract base.
 *
 * Uses only a wavefunction union object to initialize.
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
 * <tr><td> `AO_EXPANDED`  <td>Default. Exact Coulombic energy from atomic orbital expansions.
 * <tr><td> `MO_EXPANDED`  <td>Exact Coulombic energy from molecular orbital expansions
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `ESP_SYMMETRIZED` <td>Default. Coulombic energy from ESP charges interacting with nuclei and electronic density.
 *                                Symmetrized with respect to monomers.
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
 *                        + \sum_{x\in A}\sum_{\lambda\sigma\in B} Z_x V_{\lambda\sigma}^{(x)} 
 *                          \left(D_{\lambda\sigma}^{(\alpha)} + D_{\lambda\sigma}^{(\beta)}\right)
 *                        \right]
 * \f]
 * If the basis set is large and the number of ESP centres \f$ q_{x(y)} \f$ is sufficient,
 * the sum of first two contrubutions equals the sum of the latter two contributions.
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
   
};


}      // EndNameSpace oepdev
#endif //_oepdev_libutil_wavefunction_union_h

