/*! \page pimplementedmodels Implemented Models
    \tableofcontents
    \section stargets Target Properties

Detailed list of models which is to be implemented in the OEPDev project is given below:

__Table 1.__ Models subject to be implemented and analyzed within oep-dev.

 | Pauli energy             | Induction energy         | EET Coupling         |
 |--------------------------|--------------------------|----------------------|
 | EFP2-Pauli               | EFP2-Induced Dipoles     | TrCAMM               |
 | Murrel et al.'s theory   | Density Susceptibility   | OEP-ET/HT            |
 | OEP-Murrel et al.'s      |                          | TDFI-TI              |
 |                          |                          | FED                  |
 | Exact (Stone's)          | Exact (incl. CT)         | Exact (ESD)          |


   \section smethods Target, Benchmark and Competing Models

The target models introduced in the Project shall be tested against the
following benchmarks and compared with the following state-of-the-art models:

__Table 2.__ Target models vs benchmarks and competitor models.

| Target Model             | Benchmarks               | Competing Model      |
|--------------------------|--------------------------|----------------------|
| OEP-Murrel et al.'s      | Murrel et al.'s          | EFP2-Pauli           |
|                          | Exact (Stone's)          |                      |
| OEP-ET/HT + TrCAMM       | Exact (ESD)              | TDFI-TI              |
|                          | FED                      | FED                  |
|                          | TDFI-TI                  |                      |
| Density Susceptibility   | Exact (incl. CT)         | EFP2-Induced Dipoles |

*/

/*! \page pprogramming Contributing to oep-dev
    \tableofcontents

  OepDev is a plugin to Psi4. Therefore it should follow the programming etiquette of Psi4. Also,
  oep-dev has additional programming tips to make the code more versatile and easy in further development.
  Here, I emphasise on most important aspects regarding the **programming rules**.

\section smain Main routine and libraries

Oep-dev has only *one* source file in the plugin base directory, i.e., `main.cc`. This is the main
driver routine that handles the functionality of the whole OEP testing platform: specifies options for 
Psi4 input file and implements test routines based on the options. Other sources are stored
in `MODULE/libNAME*` directories where `NAME` is the name of the library with sources and header files, 
whereas `MODULE` is the directory of the oep-dev module.

Things to remember:

  1. **No other sources in base directory.** 
     It is not permitted to place any new source or other files in the plugin base directory 
     (i.e., where `main.cc` resides). 
  2. **Sources in library directories.** 
     Any additional source code has to be placed in `oepdev/libNAME*` directory (either existing one or a new one; in the 
     latter case remember to add the new `*.cc` files to `CMakeLists.txt` in the plugin base directory.
  3. **Miscellanea in special directories.** 
     If you want to add additional documentation, put it in the `doc` directory. 
     If you want to add graphics, put it in the `images` directory.


\section sheader Header files in libraries

Header files are handy in obtaining a quick glimpse of the functionality within certain library. Each library
directory should contain at least one header file in oep-dev. However, header files can be problematic if not managed properly. 

Things to remember:

   1. **Header preprocessor variable**. Define the preprocessor variable specyfying the existence of include 
      of the particular header file. The format of such is
      \code{.cpp}
      #ifndef MODULE_LIBRARY_HEADER_h
      #define MODULE_LIBRARY_HEADER_h
      // rest of your code goes here
      #endif // MODULE_LIBRARY_HEADER_h
      \endcode
      Last line is the **end** of the header file. The preprocessor variables represents
      the directory tree `oepdev/MODULE/LIBRARY/HEADER.h` structure (where `oepdev` is the base plugin directory). 
        * `MODULE` is the plugin module name (e.g. `oepdev`, the 
           name of the module directory)
        * `LIBRARY` is the name of the library (e.g. `libutil`, should be the same as library directory name)
        * `HEADER` is the name of the header in library directory (e.g. `diis` for `diis.h` header file)
   2. **Set module namespace**. To prevent naming clashes with other modules and with Psi4 it is important to operate
      in separate namespace (e.g. for a module). 
      \code{.cpp}
      namespace MODULE {
      // your code goes here
      } // EndNameSpace MODULE
      \endcode
      For instance, all classes and functions in `oepdev` module are implemented within the namespace of the same label.
      Considering addition of other local namespaces within a module can also be useful in certain cases.

\section senviron Environmental variables

Defining the set of intrinsic environmental variables can help in code management 
and conditional compilation. The oep-dev environmental variables are defined in
`include/oepdev_files.h` file. Remember also about psi4 environmental variables
defined in `psi4/psifiles.h` header. As a rule, the oep-dev environmental variable
should have the following format:
\code{.sh}
OEPDEV_XXXX
\endcode
where `XXXX` is the descriptive name of variable.

\section sdocumentation Documenting the code

Code has to be documented (at best at a time it is being created). The place for documentation 
is always in header files. Additional documentation can be also placed in source files. Leaving a chunk of code
for a production run without documentation is unacceptable. 

Use Doxygen style for documentation all the time. Remember that it supports markdown which can make the documentation
even more clear and easy to understand.
Additionally you can create a nice `.rst` documentation file for Sphinx program.
If you are coding equations, always include formulae in the documentation!

Things to remember:

   1. **Descriptions of classes, structures, global functions, etc**. Each programming object should have a description.
   2. **Documentation for function arguments and return object**. 
      Usage of functions and class methods should be explained by providing the description of all arguments 
      (use `\param` and `\return` Doxygen keywords).
   3. **One-line description of class member variables**. Any class member variable should be preceded by 
      a one-liner documentation (starting from `///`).
   4. **Do not be afraid of long names in the code**. Self-documenting code is a bless!

\section snaming Naming conventions

Naming is important because it helps to create more readable and clear self-documented code. 

Some loose suggestions:

   1. **Do not be afraid of long names in the code, but avoid redundancy**. Examples of good and bad names:
      * good name: `get_density_matrix`; bad name: `get_matrix`. Unless there is only one type of matrix
        a particular objects can store, `matrix` is not a good name for a getter method. 
      * good name: `class Wavefunction`, bad name: `class WFN`
      * good name: `int numberOfErrorVectors`, bad name: `int nvec`, bad name: `the_number_of_error_vectors`
      * good name: `class EFPotential`, probably bad name: `class EffectiveFragmentPotential`.
        The latter might be understood by some people as a class that inherits from `EffectiveFragment` class. 
        If it is not the case, compromise between abbreviation and long description is OK.
   2. **Short names are OK in special situations**. In cases meaning of a particular variable is obvious and
      it is frequently used in the code locally, it can be named shortly. Examples are:
      * `i` when iterating
      * `no` number of occupied orbitals, `nv` number of virtual orbitals, etc.
   3. **Clumped names for variables and dashed names for functions**. Try to distinguish between variable name 
      like `sizeOfOEPTypeList` and a method name `get_matrix()` (neither `size_of_OEP_type_list`, nor `getMatrix()`).
      This is little bit cosmetics, but helps in managing the code when it grows.
   4. **Class names start from capital letter**. However, avoid only capital letters in class names, unless it is obvious.
      Avoid also dashes in class names (they are reserved for global functions and class methods). Examples: 
      * good name: `DIISManager`, bad name: `DIIS`.
      * good name: `EETCouplingSolver`, bad name: `EETSolver`, very bad: `EET`.
      
\section soop Use Object-Oriented Programming

Try to organise your creations in objects having special relationships and data structures. Encapsulation
helps in producing self-maintaining code and is much easier to use. Use: 
 - **factory design** for creating objects
 - **container design** for designing data structures
 - **polymorphysm** when dealing with various flavours of one particular feature in the data structure
 
 > *Note:* In Psi4, factories are frequently implemented as static methods of the base classes, 
 > for example `psi::BasisSet::build` static method. It can be followed when building object factories 
 > in oep-dev too.

*/

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
 * <tr><td> `AO_EXPANDED`  <td>*Default*. Exact Coulombic energy from atomic orbital expansions.
 * <tr><td> `MO_EXPANDED`  <td>Exact Coulombic energy from molecular orbital expansions
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `ESP_SYMMETRIZED` <td>*Default*. Coulombic energy from ESP charges interacting with nuclei and electronic density.
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

/**\brief Compute the Pauli-Repulsion interaction energy between unperturbed wavefunctions.
 *
 * The implemented methods are shown below
 * <table>
 * <caption id="Tab.1">Methods available in the Solver</caption>
 * <tr><th> Keyword  <th>Method Description  
 * <tr><td colspan=2> <center><strong>Benchmark Methods</strong></center>
 * <tr><td> `HAYES_STONE`       <td>*Default*. Exact Pauli Repulsion energy at HF level from Hayes and Stone (1984). 
 * <tr><td> `MURRELL_ETAL`      <td>Approximate Pauli Repulsion energy at HF level from Murrell et al.
 * <tr><td> `EFP2`              <td>Approximate Pauli Repulsion energy at HF level from EFP2 model.
 * <tr><td colspan=2> <center><strong>OEP-Based Methods</strong></center>
 * <tr><td> `MURRELL_ETAL_MIX`  <td>*Default*. OEP-Murrell et al's: S1 term via DF-OEP, S2 term via ESP-OEP.
 * <tr><td> `MURRELL_ETAL_ESP`  <td>OEP-Murrell et al's: S1 and S2 via ESP-OEP
 * </table>
 *
 * Below the detailed description of the above methods is given.
 * In the formulae, the Coulomb notation for 
 * electron repulsion integrals (ERI's) in MO basis is adopted; i.e,
 * \f[
 *  (ac \vert bd) = \iint d{\bf r}_1 d{\bf r}_2 
 *   \phi_a({\bf r}_1) \phi_c({\bf r}_1) \frac{1}{r_{12}} \phi_b({\bf r}_2) \phi_d({\bf r}_2)
 * \f]
 * It is also assumed that the orbitals are real.
 *
 * # Benchmark Methods
 * ## Exact Pauli Repulsion energy at HF level.
 *    
 * For a closed-shell system, equation of Hayes and Stone (1984)
 * become
 * \f[
 *    E^{\rm Rep} = 2\sum_{ac\in A} \sum_{bd\in B} 
                    \left( V^A_{ab} + V^B_{ab} + T_{ab} \right) 
                    \left[ [{\bf S}^{-1}]_{ab} - \delta_{ab} \right]
                +   \sum_{ac\in A} \sum_{bd\in B}
                    (ac \vert bd) 
                    \left\{ 
       [{\bf S}^{-1}]_{ab} [{\bf S}^{-1}]_{cd} - 
       [{\bf S}^{-1}]_{ad} [{\bf S}^{-1}]_{bc} -
      2\delta_{ab} \delta_{cd} +
      2\delta_{ad} \delta_{bc}
                    \right\}
 * \f]
 * where \f$ {\bf S} \f$ is the overlap matrix between the doubly-occupied
 * orbitals.
 * The exact exchange energy is for a closed shell case given as
 * \f[
     E^{\rm Ex} = -2\sum_{a\in A} \sum_{b\in B} (ab \vert ab)
 * \f]
 * Similarity transformation of molecular orbitals does not affect the resulting energies.
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
 *                \left[ V^B + 2\sum_{d\in B} (ac \vert dd) \right]
 *           +      \sum_{d\in B} S_{ad}
 *                \left[ V^A + 2\sum_{x\in A} (bd \vert cc) \right]
 *                - \sum_{c\in A} \sum_{d\in B} S_{cd} (ac \vert bd)
 *         \right\}
 * \f] 
 * By using OEP technique, the above theory can be exactly re-cast *without* any further approximations.
 *
 * # OEP-Based Methods
 * The Murrell et al's theory of Pauli repulsion is here re-cast by introducing OEP's.
 * 
 * ## S1 term via DF-OEP, S2 term via ESP-OEP.
 * ## S1 and S2 terms via ESP-OEP.
 *
 * *Notes:* 
 *   - This solver also computes and prints the exchange energy at HF level (formula is given above)
 *     for reference purposes.
 *   - In order to construct this solver, **always** use the `OEPDevSolver::build` static factory method.
 */
class RepulsionEnergySolver : public OEPDevSolver
{
  public:
    RepulsionEnergySolver(SharedWavefunctionUnion wfn_union);
    virtual ~RepulsionEnergySolver();

    virtual double compute_oep_based(const std::string& method = "DEFAULT");
    virtual double compute_benchmark(const std::string& method = "DEFAULT");
};


}      // EndNameSpace oepdev
#endif //_oepdev_libutil_wavefunction_union_h

