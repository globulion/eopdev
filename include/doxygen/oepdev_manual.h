//-----------------------------------------------------------

/*! \page intro Introduction
    \tableofcontents

Exploring biological phenomena at molecular scale is oftentimes indispensable
to develop new drugs and intelligent materials. Most of relevant system 
properties are affected by intermolecular interactions with nearby 
environment such as solvent
or closely bound electronic chromophores.
Studying such molecular aggregates requires rigorous and accurate quantum 
chemistry methods, the cost of which grows very fast with the number 
of electrons. 
Despite many methodologies have been devised to describe 
energetic and dynamical properties of **extended molecular systems** efficiently and 
accurately, there exist particularly difficult cases in which 
modelling is still challenging: 

  - describing electronic transitions in solution or
  - when coupled with other electronic transition via resonance energy transfer, 
  - performing molecular dynamics 
    at very high level of theory including dynamic electron correlation, 
  - vibrational frequency calculations of particular normal mode in 
    condensed phases
 
and so on. The reason behind (sometimes prohibitively) 
high costs of fully *ab initio* calculations in the above areas 
is the complexity of mathematical models often based on wave functions 
rather then (conceptually more straightforward) electronic densities. 
On the other hand, 
it has been pointed out before
that the one-electron density 
distributions are of particular importance in chemistry. It can be 
thus utilized as a means of developing a general model that 
re-expresses the physics of intermolecular interactions in terms 
of effective one-electron functions that are easier to handle 
in practice.

This Project focuses on finding a unified way to simplify various
fragment-based approaches
of Quantum Chemistry of extended molecular systems, i.e., 
molecular aggregates such as interacting chromophores
and molecules solvated by water and other solvents. 
Indeed, one of the important difficulties encountered in Quantum Chemistry
of large systems is the need of evaluation of special kind of numbers
known as *electron repulsion integrals*, or in short, ERI's. In a typical
calculation, the amount of ERI's can be as high as tens or even hundreds of millions (!)
that unfortunately prevents from application of conventional methods when the number of particles
in question is too large. In the Project, the complicated expressions involving ERI's 
shall be greatly simplified to reduce the computational costs as much as possible
while introducing no or minor approximations to the original theories.

\section resprojmeth Research Project Methodology

In this Project the new theoretical protocol based on the 
effective one-electron potentials (EOP's) is developed. 
The main principle is to rewrite arbitrary sum of functions \f$ f \f$ 
of electron repulsion integrals (ERI's) by defining \ref poepdesign "EOP's"
according to the following general prescription:
\f{align*}{
\sum_f f\left[ 
   \left( {\phi_i^A}{\phi_j^A} \vert\vert {\phi_k^B}{\phi_l^B} \right)
 \right] &= \left( {\phi_i^A} \vert {v_{kl}^B} \vert {\phi_j^A} \right) \rightarrow 
 \text{ point charge or density fitting} \\
\sum_f f\left[ 
   \left( {\phi_i^A}{\phi_j^B} \vert\vert {\phi_k^B}{\phi_l^B} \right)
 \right] &= \left( {\phi_i^A} \vert {v_{kl}^B} \vert {\phi_j^B} \right) \rightarrow 
 \text{ density fitting,} 
\f}
where \f$ A \f$ and \f$ B \f$ denote different molecules 
and \f$ \phi_i \f$ is the \f$ i \f$-th molecular orbital
or basis function. Here,
\f$ v_{kl}^B \f$ denotes the \ref poeptypes *ab initio* "EOP matrix element". The technique described above 
will be applied to simplify expressions for 
 - short-range excitation 
   energy transfer couplings between chlorophyll subunits of reaction centres 
   in photosynthesis
 - Pauli interaction repulsion energy 
 - charge-transfer interaction energy
 - electric field-induced charge density polarization of molecules. 

The above developments might be used 
in fragment-based *ab initio* molecular dynamics protocols of new generation. 

\section impact Expected Impact on the Development of Science, Civilization and Society

The proposed EOP's are expected to significantly develop 
the fragment-based methods that are widely used 
in physical chemistry and modelling of biologically 
important systems. Owing to universality of EOP's, 
they could find applications in many branches of chemical science: 
*non-empirical* molecular dynamics, short-range
resonance energy transfer in photosynthesis,
electronic and vibrational solvatochromism, 
multidimensional spectroscopy and so on. In particular:

 - the EOP-based models of Pauli repulsion energy and charge-transfer (CT) energy
   could be used to improve the computational performance of the second generation effective fragment
   potential method (EFP2). Particularly, the EFP2 CT term is relatively time consuming and due to this reason
   it is sometimes omitted in the applications of EFP2 to perform molecular dynamics 
   simulations.\cite OEP-1.2020
 - the EOP-based model of EET couplings could significantly improve modelling of energy transfer
   in the light harvesting complexes. At present, short-range phenomena (Dexter mechanisms of EET)
   are very difficult to efficiently and quantitatively evaluate when performing statistical averaging
   and applying to large molecular aggregates. Such Dexter effects could be computed by using EOP's
   in much more efficient manner without loosing high accuracy of state-of-the-art methods such as
   TDFI-TI method.\cite OEP-2.2020
 - the density matrix polarization (DMS) tensors could be used in new generation fragment-based *ab initio*
   molecular dynamics protocols that rigorously take into consideration electron correlation effects.\cite Blasiak.JCP.2018

Therefore, we believe that the application of
EOP's could have an indirect impact on 
the design of novel drugs and materials for industry.

\section oepdevcode The EOPDev Code

To pursue the above challenges in the field of computational
quantum chemistry of extended molecular aggregates,
the EOPDev platform is developed.
Accurate and efficient *ab initio* \ref pimplementedmodels "models"
based on EOP's are implemented in the EOPDev code, along with
the state-of-the-art benchmark and competiting methods. 
Written in C++ with an extensive Python interface, 
EOPDev is a plugin to Psi4 quantum chemistry package. Therefore,
compilation and running the EOPDev code is straightforward and follows
the API interface similar to the one used in Psi4 with 
just a few \ref pprogramming "specific programing conventions".
The detailed discussion about using the EOPDev code
can be found in \ref advanced "usage section". 
\note The `OEP' abbreviation, rather than the `EOP',
      is used throughout the code. It is because the earlier versions
      of EOPDev utilized the former abbreviation consecutively for the shared libraries,
      modules and class names. The abbreviation was changed to the latter
      in this public release of the code.
      Please treat these two abbreviations as synonyms within the
      project and code, both refering to the *effective one-electron potentials*.
*/

//-----------------------------------------------------------

//-----------------------------------------------------------


/*! \page poepdesign EOP Design.
    \tableofcontents

EOP (One-Electron Potential) is associated with certain quantum one-electron operator 
\f$ \hat{v}^A\f$ that defines the ability of molecule \f$ A \f$ to interact in a particular way with other molecules. 
Technically, EOP can be understood as a __container object__ (associated with the molecule in question)
that stores the information about the above mentioned quantum operator. 
Here, it is assumed that similar EOP
object is also defined for all other molecules in a molecular aggregate. 

In case of interaction between molecules \f$ A \f$ and \f$ B \f$,
EOP object of molecule \f$ A \f$ interacts directly with wavefunction object
of the molecule \f$ B \f$. Defining a 
 * Solver class that handles such interaction 
 * Wavefunction class and
 * EOP class

the universal design of EOP-based approaches can be established and developed.

> **Important:**
>  EOP and Wavefunction classes should not be restricted to Hartree-Fock (HF); in generall any correlated 
>  wavefunction and derived EOP`s should be allowed to work with each other.
>  However, in the current version of the project, only HF wavefunctions are considered.
>

\section soepclasses EOP Classes

There are many types of EOP’s, but the underlying principle is the same and independent of the
type of intermolecular interaction. Therefore, the EOP’s should be implemented by using a multi-level class design.
In turn, this design depends on the way EOP’s enter the mathematical expressions, i.e., on the types
of matrix elements of the one-electron effective operator \f$ \hat{v}^A \f$.

\subsection ssoepunification Structure of possible EOP-based expressions and their unification

Structure of EOP-based mathematical expressions is listed below:

| Type   |                 Matrix Element                     |       Comment            |
|--------|----------------------------------------------------|--------------------------|
| Type 1 | \f$ \left( I \lvert \hat{v}^A \rvert J \right) \f$ | \f$ I\in A,\; J\in B\f$  |
| Type 2 | \f$ \left( J \lvert \hat{v}^A \rvert L \right) \f$ | \f$ J,L\in B\f$          |


In the above table, \f$ I \f$, \f$ J \f$ and \f$ K \f$ indices correspond 
to basis functions or molecular orbitals. 
Basis functions can be primary or auxiliary EOP-specialized density-fitting.
Depending on the type of function and matrix element, there are many subtypes of resulting matrix elements 
that differ in their dimensionality. Examples are given below:

| Matrix Element | DF-based form | DMTP-based form |
|----------------|---------------|----------------|
| \f$ \left(\mu \lvert \hat{v}^{A[\mu]} \rvert \sigma \right) \f$ | \f$ \sum_{\iota\in A} v_{\mu\iota}^AS_{\iota\sigma} \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[\mu]} V_{\mu\sigma}^{(\alpha)}\f$ |
| \f$ \left( i  \lvert \hat{v}^{A[i  ]} \rvert j      \right) \f$ | \f$ \sum_{\iota\in A} v_{i  \iota}^AS_{\iota j    } \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[i  ]} V_{i  j     }^{(\alpha)}\f$ |
| \f$ \left(j   \lvert \hat{v}^{A[i  ]} \rvert l      \right) \f$ | \f$ \sum_{\iota\kappa\in A} S_{j \iota} v_{\iota\kappa}^{A[i]}S_{\kappa l} \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[i]} V_{jl}^{(\alpha)}\f$ |

In the formulae above, the EOP-part (stored by EOP instances) and
the Solver-part (to be computed by the Solver) are separated.
It is apparent that all EOP-parts have the form of 2- or 3-index arrays
with different class of axes (molecular orbitals, primary/auxiliary basis, atomic space). 
Therefore, they can be uniquely defined by a unified *tensor object* 
(storing double precision numbers) and unified *dimension object* storing 
the information of the axes classes.

In Psi4, a perfect candidate for the above is `psi4::Tensor` class declared 
in `psi4/libthce/thce.h`. Except from the numeric content its instances also 
store the information of the dimensions in a form of a vector of `psi4::Dimension` instances.

Another possibility is to use `psi::Matrix` objects, instead of `psi4::Tensor` objects, 
possibly putting them into a `std::vector` container in case there is more than two axes.

\note Currently, the second possibility is used, i.e., matrices. For more complex data structures,
      other types of custom objects are defined.

*/

//----------------------------------------------------------------------------------

/*!\page pdensityfitting Density-fitting Specialized for EOP’s
   \tableofcontents

To get the ab-initio representation of a EOP, one can use a procedure similar to
the typical density fitting or resolution of identity, both of which are nowadays widely used 
to compute electron-repulsion integrals (ERI’s) more efficiently. 

\section sdensfitcompl Fitting in Complete Space

An arbitrary one-electron potential of molecule *A* acting on any state vector 
associated with molecule *A* can be expanded in an *auxiliary space* centered 
on *A* as
\f[
   v\vert i) = \sum_{\xi\eta} v\vert \xi) [{\bf S}^{-1}]_{\xi\eta} ( \eta \vert i)
\f]
under the necessary assumption that the auxiliary basis set is *complete*. 
In a special case when the basis set is orthogonal (e.g., molecular orbitals)
the above relation simplifies to
\f[
   v\vert i) = \sum_{\xi} v\vert \xi) ( \xi \vert i)
\f]
It can be easily shown that the above general and exact expansion can be obtained by 
performing a density fitting in the complete space. We expand the LHS of the first equation 
on this page in a series of the auxiliary basis functions scaled by the undetermined expansion coefficients:
\f[
  v\vert i) = \sum_{\xi} {G_{i\xi}} \vert \xi)
\f]
which we shall refer here as to the matrix form of the EOP operator.
By constructing the least-squares objective function 
\f[
 Z[\{G^{(i)}_\xi\}] = \int d{\bf r}_1
                     \left[    v({\bf r}_1) \phi_i({\bf r}_1) - \sum_\xi G^{(i)}_\xi \varphi_\xi({\bf r}_1) \right]^2
\f]
and requiring that
\f[
 \frac{\partial Z[\{G^{(i)}_\xi\}]}{\partial G^{(i)}_\mu} = 0 \text{ for all $\mu$}
\f]
we find the coefficients \f$ G^{(i)}_\xi \f$ to be
\f[
 {\bf G}^{(i)} = {\bf v}^{(i)} \cdot {\bf S}^{-1}
\f]
where 
\f{align*}{
 v^{(i)}_\eta &= (\eta \vert vi) \\
 S_{\eta\xi}  &= (\eta \vert \xi)
\f}
or explitictly
\f[
  {G_{i\xi}} = \sum_\eta [{\bf S}^{-1}]_{\xi\eta} (\eta\vert v\vert i)
\f]
identical to what we obtained from application of the resolution of identity 
in space spanned by non-orthogonal complete set of basis vectors.

Since matrix elements of an EOP operator in auxiliary space can be computed 
in the same way as the matrix elements with any other basis function, 
one can formally write the following identity
\f[
 (X \vert v\vert i) = \sum_{\xi\eta} S_{X \xi} [{\bf S}^{-1}]_{\xi\eta} (\eta\vert v\vert i)
\f]
where \f$ X \f$ is an arbitrary orbital.
When the other orbital
does not belong to molecule *A* but to the (changing) environment, it is 
straightforward to compute the resulting matrix element, which is simply given as  
\f[
   (j_{\in B} \vert v^A \vert i_{\in A}) = \sum_\xi {S_{j\xi}} {G_{i\xi}}
\f]
where *j* denotes the other (environmental) basis function.

In the above equation, the EOP-part (fragment parameters for molecule *A* only) 
and the Solver-part (subject to be computed by solver 
on the fly) are separated. This then forms a basis for fragment-based 
approach to solve Quantum Chemistry problems related to the extended molecular aggregates.

\section sdensfitincompl Fitting in Incomplete Space

Density fitting scheme from previous section has practical disadvantage of a nearly-complete basis set
being usually very large (spanned by large amount of basis set vectors). Any non-complete basis set
won't work in the previous example. Since most of basis sets used in quantum chemistry do not form a complete
set, it is beneficial to design a modified scheme in which it is possible to obtain the **effective** 
matrix elements of the EOP operator in a **incomplete** auxiliary space. This can be achieved by minimizing 
the following objective function
\f[
 Z[\{G^{(i)}_\xi\}] = \iint d{\bf r}_1 d{\bf r}_2
                      \frac{
                     \left[    v({\bf r}_1) \phi_i({\bf r}_1) - \sum_\xi G^{(i)}_\xi \varphi_\xi({\bf r}_1) \right]
                     \left[    v({\bf r}_2) \phi_i({\bf r}_2) - \sum_\xi G^{(i)}_\eta\varphi_\eta({\bf r}_1) \right]
                           }{\vert {\bf r}_1 - {\bf r}_2 \vert}
\f]
Thus requesting that 
\f[
 \frac{\partial Z[\{G^{(i)}_\xi\}]}{\partial G^{(i)}_\mu} = 0 \text{ for all $\mu$}
\f]
we find the coefficients \f$ G^{(i)}_\xi \f$ to be
\f[
 {\bf G}^{(i)} = {\bf b}^{(i)} \cdot {\bf A}^{-1}
\f]
where 
\f{align*}{
 b^{(i)}_\eta &= (\eta \vert\vert vi) \\
 A_{\eta\xi}  &= (\eta \vert\vert \xi)
\f}
The symbol \f$ \vert\vert \f$ is to denote the operator \f$ r_{12}^{-1}\f$ and double integration over \f$ {\bf r}_1 \f$
and \f$ {\bf r}_2 \f$. Thus, in order to use this generalized density fitting scheme
one must to compute two-centre electron repulsion integrals (implemented in ERI_1_1). 

\section sdensfitincomplalt Fitting in Incomplete Space - Alternative Approach

TODO

*/

//--------------------------------------------------------------------------------------------

/*! \page pimplementedmodels Implemented Models
    \tableofcontents
    \section stargets Fragment-Based Methods

List of most important models implemented in the EOPDev project is given below.
Among the interaction energy models are the
second generation of the effective potential method 
(EFP2) \cite EFP2.2013 \cite Li.Gordon.Jensen.JCP.2006 \cite Xu.Gordon.JCP.2013,
perturbation theories of Murrel et al.\cite Murrell.1965, Otto and Ladik \cite Otto.Ladik.ChemPhys.1975
and Hayes and Stone \cite Hayes.Stone.MolPhys.1984, 
density decomposition scheme (DDS) \cite Mandado.Hermida-Ramon.JCTC.2011, 
reduced variational space (RVS) method \cite Stevens.Fink.CPL.1987.
Among the excitation energy transfer (EET) coupling methods are
the TrCAMM method\cite Blasiak.Maj.Cho.Gora.JCTC.2015
and the transfer integral (TI) method \cite TDFI-TI.2012.


__Table 1.__ Theoretical fragment-based models implemented in EOPDev.

 | Pauli energy            | CT energy                  | EET Coupling         | 
 |-------------------------|----------------------------|----------------------|
 | EFP2                    | EFP2                       | TrCAMM               |
 | Murrel et al.           | Otto-Ladik                 | TI                   |
 | EOP-Murrel et al.       | EOP-Otto-Ladik             | EOP-TI               |
 | Otto-Ladik              |                            |                      |
 | EOP-Otto-Ladik          |                            |                      |
 | DDS                     |                            |                      |
 | Hayes-Stone (exact)     | RVS                        | Exact (ESD)          |

   \section smethods Target, Benchmark and Competing Models

The target models introduced in the Project are tested against the
following benchmarks and compared with the following state-of-the-art models:

__Table 2.__ Target models vs benchmarks and competitor models.

| Target Model             | Benchmarks               | Competing Model      |
|--------------------------|--------------------------|----------------------|
| EOP-Murrel et al. (Pauli)| Murrel et al., DDS, Stone| EFP2 (Pauli)         |
| EOP-Otto-Ladik (CT)      | Otto-Ladik, RVS          | EFP2 (CT)            |
| EOP-TI                   | Exact (ESD), TI          | TI                   |

The target models contain their EOP-based versions, that can be executed in the
`OEPDevSolver::compute_oep_based`, and compared with the corresponding
benchmark models `OEPDevSolver::compute_benchmark`. 

*/

//--------------------------------------------------------------------------------------------

/*! \page pprogramming Contributing to EOPDev
    \tableofcontents

EOPDev is a plugin to Psi4. Therefore it should follow the programming etiquette of Psi4. Also,
EOPDev has additional programming tips to make the code more versatile and easy to develop further.
Here, I emphasise on most important aspects regarding the proposed **programming rules**.

\section scontrmain Main routine and libraries

EOPDev has only *one* source file in the plugin base directory, i.e., `main.cc`. This is the main
driver routine that handles the functionality of the whole EOP testing platform: specifies options for 
Psi4 input file and implements test routines based on the options. 
Include files directly related to `main.cc` are stored in the `include` directory, where only header
files are present. Options are specified in `include/oepdev_options.h` whereas macros and defines
in `include/oepdev_files.h`.
Other sources are stored
in `MODULE/libNAME*` directories where `NAME` is the name of the library with sources and header files, 
whereas `MODULE` is the directory of the EOPDev module.

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
directory should contain at least one header file in EOPDev. However, header files can be problematic if not managed properly. 

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
and conditional compilation. The EOPDev environmental variables are defined in
`include/oepdev_files.h` file. Remember also about psi4 environmental variables
defined in `psi4/psifiles.h` header. As a rule, the EOPDev environmental variable
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
      like `sizeOfEOPTypeList` and a method name `get_matrix()` (neither `size_of_EOP_type_list`, nor `getMatrix()`).
      This is little bit cosmetics, but helps in managing the code when it grows.
   4. **Class names start from capital letter**. However, avoid only capital letters in class names, unless it is obvious.
      Avoid also dashes in class names (they are reserved for global functions and class methods). Examples: 
      * good name: `DIISManager`, bad name: `DIIS`.
      * good name: `EETCouplingSolver`, bad name: `EETSolver`, very bad: `EET`.
     
\section stime Track timing when evaluating the code

It is useful to track time elapsed for performing a particular task
by a computer. For this, use for example `psi::timer_on` and `psi::timer_off` functions defined in `psi4/libqt/qt.h`.
Psi4 always generates the report file `timer.dat` that contains all the defined timings.
For example, 
\code{.cpp}
#include "psi/libqt/qt.h"
psi::timer_on("EOP    E(Paul) Murrell-etal S1  ");
// Your code goes here
psi::timer_off("EOP    E(Paul) Murrell-etal S1  ");
\endcode
To maintain the printout in a neat form, the timing associated with the EOPDev code 
can be generated via `misc/python/timing.py` utility script.

\section memclean Clean memory between independent jobs

If you use scratch disk space to store integrals, clean the scratch in between
independent calculations. From C++ level invoke
\code{.cpp}
#include "psi4/libpsio/psio.hpp"
// ...
psi::PSIOManager::shared_object()->psiclean();
\endcode
whereas from the Python level use
\code{.py}
import psi4
# ..
psi4.core.clean()
\endcode
If the scratch space is not cleaned up before next independent task begins, 
certain computational routines might crash with `PSIOError`
or continue without error, but produce wrong results.

\section soop Use Object-Oriented Programming

Try to organise your creations in objects having special relationships and data structures. Encapsulation
helps in producing self-maintaining code and is much easier to use. Use: 
 - **factory design** for creating objects
 - **container design** for designing data structures
 - **polymorphysm** when dealing with various flavours of one particular feature in the data structure
 
 > *Note:* In Psi4, factories are frequently implemented as static methods of the base classes, 
 > for example `psi::BasisSet::build` static method. It can be followed when building object factories 
 > in EOPDev too.

*/

// --------------------------------------------------------------------------

/*! \page advanced Usage
    \tableofcontents

    EOPDev is addressed for developers. Make sure you have first read \ref intro "the introduction"
    and are familiar with the EOP-based ERI elimination technique
    before proceeding.

\section padvinstall Installation

\subsection spadvpreppsi Preparing Psi4

EOPDev is a Psi4 plugin. It requires 
  - Psi4, version 1.2.1 (git commit `7a58e8c`). Has to be modified (see below).
  - Eigen3, any version.

\note
 Before compiling, make sure EFP is enabled in `CMakeLists.txt` (now it is not
 used in EOPDev but maybe in the future it would).

Recently, Psi4 introduced API visibility management. Only certain Psi4 classes
and functions are *exposed* in the `core.so` library, that is further linked to
Psi4 plugin shared library. Due to this reason, not all Psi4 functionalities
can be directly used from outside Psi4. In order to access local API of Psi4
(also used in the EOPDev code) slight modification of Psi4 code and concomitant
rebuild is necessary.

In order to expose local API used by EOPDev and hidden within Psi4 1.2,
two types of small modifications are necessary:
  - M1: add `PSI_API` macro after required class or function declaration in header file
  - M2: add `#include "psi4/pragma.h"` line at the include section of an appropriate header file 

Modification M1 is obligatory for all affected files whereas modification M2 needs to be done only
in headers that do not have "psi4/pragma.h" included explicitly or implicitly.
The list of some Psi4 header files along with the respective changes 
that need to be done are listed in the table below:

| Psi4 Header File                   | Psi4 Class          | Required Changes |
|------------------------------------|---------------------|------------------|
|`libfunctional/superfunctional.h`   | `Superfunctional`   | M1               |
|`libscf_solver/hf.h`                | `HF`                | M1               |
|`libscf_solver/rhf.h`               | `RHF`               | M1               |
|`libcubeprop/csg.h`                 | `CubicScalarGrid`   | M1               |
|`libmints/onebody.h`                | `OneBodyAOInt`      | M1               | 
|`libmints/potential.h`              | `PotentialInt`      | M1               | 
|`libmints/multipoles.h`             | `MultupoleInt`      | M1               |
|`libmints/multipolesymmetry.h`      | `MultipoleSymmetry` | M1               |
|`libmints/fjt.h`                    | `Taylor_Fjt`        | M1               |
|`libmints/fjt.h`                    | `Fjt`               | M1               |
|`libmints/oeprop.h`                 | `EOProp`            | M1, M2           |
|`libmints/gshell.h`                 | `GaussianShell`     | M1, M2           |

To quickly apply these and other required modifications, use the patch files stored in
`misc/patch` directory.
Please make sure to use a proper patch for a chosen Psi4 version.

\subsection spadvcompile Compiltation

After all the above changes have been done in Psi4 (followed by its rebuild)
compile the EOPDev code by running `compile` script. Make sure Eigen3 path 
is set to environment variable `EIGEN3_INCLUDE_DIR` (instructions will appear on the screen).
After compilation is successful, run `ctest` to check if the code works fine.

\note 
 It may happen that during code development there will be symbol lookup error when
 importing `oepdev.so` (in such case EOPDev compiles without error but Python cannot import the module `oepdev`).
 In such circumstance, probably there some local Psi4 feature that is needed in EOPDev
 is not exposed by `PSI_API` macro. To fix this, run `c++filt [name]` where `[name]`
 is the mangled undefined symbol. This will show you which Psi4 class or function
 is not exposed and requires `PSI_API` (change M1 and perhaps M2 too). Such change requires
 Psi4 rebuild and recompilation of EOPDev code. In any case, please contact me and report
 new undefined symbol (blasiak.bartosz@gmail.com).


\section padvcodestr EOPDev Code Structure

As a plugin to Psi4, EOPDev consists of the `main.cc` file with the plugin main routine,
`include/oepdev_options.h` specifying the options of the plugin, `include/oepdev_files.h`
defining all global macros and environmental variables, as well as the `oepdev` directory.
The latter contains the actual EOPDev code that is divided into several subdirectories called
\ref smod "modules".

\subsection smain Main Routine

Before the actual EOPDev calculations are started, the wavefunction of
the input molecular aggregate is computed by Psi4. See the plugin driver script
`pymodule.py` for more details on how the calculation environment is initialized. 
Subsequently, 
one out of four types of target operations can be performed by the program:

 1. `OEP_BUILD` - Compute the EOP effective parameters for one molecule.
 2. `DMATPOL`   - Compute the generalized density matrix susceptibility tensors (DMS's) for one molecule.
 3. `SOLVER`    - Perform calculations for a molecular aggregate. As for now, only dimers are handled.
 4. `TEST`      - Perform the testing routine.

The first two modes are single molecule calculations. `OEP_BUILD` 
uses the `OEPotential::build` static factory to create EOP objects
whereas `DMATPOL` uses the `GenEffParFactory::build` static factory
to greate generalized effective fragment parameters (GEFP's) for polarization.
\note In the future, `OEP_BUILD` will be handled also by `GenEffParFactory::build`
      since EOP parameters are part of the GEFP's.

`SOLVER` requires at least molecular dimer
and the `WavefunctionUnion` object (being the Hartree product of the unperturbed 
monomer wavefunctions) is constructed at the beginning, which is then passed to 
the `OEPDevSolver::build` static factory.
`TEST` can refer to single- or multiple-molecule calculations, whereby each of the testing routines
is listed in the `cmake/CTestTestfile.cmake.in` file.

\subsection smod Modules

The source code is distributed into directories called modules:
 - `liboep`
 - `libgefp`
 - `libsolver`
 - `libints`
 - `libpsi`
 - `lib3d`
 - `libutil`
 - `libtest`

See \ref Modules "Modules" for a detailed description of each of the modules.

\section padvclasses EOPDev Classes: Overview
\subsection ssclassesmoep EOP Module

The EOP module located in `oepdev/liboep` consists of the following abstract bases:
 - `OEPotential` implementing the EOP,
 - `GeneralizedDensityFit` implementing the GDF technique.

Each of the bases contains static factory method called `build` that creates instances
of chosen subclasses.
The module contains also a structure `EOPType` which is a container storing all the data
associated with a particular EOP: type name, dimensions, EOP coefficients and whether is density-fitted
or not.

\subsubsection OEPotential

It is a container and computer class of EOP. Among others, the most important public method
is `OEPotential::compute` which computes all the EOP's (by iterating over all possible EOP types
within a chosen EOP subclass or category). EOP's can be extracted by `OEPotential::oep` method, for instance.
From protected attributes, each OEPotential instance stores blocks of the LCAO-MO matrices associated with the
occupied (`cOcc_`) and virtual (`cVir_`) MO's. It also contains the pointers to the primary, auxiliary and intermediate
basis sets (`primary_`, `auxiliary_` and `intermediate_`, accordingly). Usage example:
\code{.cpp}
#include "oepdev/liboep/oep.h"
oep = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn, options);
oep->compute();
oep->write_cube("V", "oep_cube_file");
\endcode
So far, four OEPotential subclasses are implemented, from which `ElectrostaticEnergyOEPotential`
and `RepulsionEnergyOEPotential` are fully operative, while the rest is under development.

\subsubsection GeneralizedDensityFit

Implements the density fitting schemes for EOP's.

\subsection ssclassesmgefp GEFP Module

This module deals with the effective fragments consituting an extended molecular aggregate.
It builds the platform to test various generalized effective fragment potentials (GEFP).

\subsubsection GenEffPar

Represents generalized effective fragment parameters.

\subsubsection GenEffParFactory

Implements routines of calculation of effective fragment parameters of various types.

\subsubsection GenEffFrag

Represents one effective fragment.

\subsection ssclassesmsolver EOPDev Solver Module

This module sets up a simple platform of comparing benchmark and EOP-based
fragment-based methods.

\subsubsection OEPDevSolver

This is the main solver which as for now assumes molecular dimers (or bi-fragment systems).
It is based on a union of wavefunctions of unperturbed monomers, `WavefunctionUnion`.

\section sprog Developing EOP's

\note This section is for illustrative purpose. The small details of
      the objects such as `OEPType` and others can change over the years due to
      development of the EOPDev code. However, the overal programing scheme
      remains unchanged and valid.

EOP's are implemented in a suitable subclass of the `OEPotential` base.
Due to the fact that EOP's can be density-based or DMTP-based, the classes 
`GeneralizedDensityFit` as well as `ESPSolver` are usually necessary in the implementations. 
Handling the one-electron integrals (OEI's) and the two-electron integrals (ERI's) 
in AO basis is implemented in `IntegralFactory`. In particular, 
potential integrals evaluated at arbitrary centres can be accessed
by using the `PotentialInt` instances.
Useful iterators for looping over AO ERI's the `ShellCombinationsIterator`
and `AOIntegralsIterator` classes. Transformations of OEI's to MO basis
can be easily achieved by transforming AO integral matrices by `cOcc_` and `cVir_`
members of `OEPotential` instances, e.g., by using the `psi::Matrix::doublet` or `psi::Matrix::triplet`
static methods. Transformations of ERI's to MO basis can be performed by using
the `psi4/libtrans/integraltransform.h` library.

It is recommended that the implementation of all the new EOP's follows the following steps:

 1. **Write the class framework.** 
    This includes choosing a proper name of a OEPotential subclass, 
    sketching the constructors and a destructor, and all the necessary methods.
 2. **Implement EOP types.** Each type of EOP is implemented, including the 3D
    vector field in case DMTP-based EOP's are of use.
 3. **Update base factory method**. Add appropriate entries in the `OEPotential::build`
    static factory method.

Below, we shall go through each of these steps separately and discuss them in detail.

\subsection ssoepsubdraft Drafting an EOP Subclass

This stage is the design of the overall framework of EOP subclass. The name should end with `OEPotential`
to maintain the convention used so far. The template for the header file definition can be depicted as follows:
\code{.cpp}
class SampleOEPotential : public OEPotential 
{
  public:
    // Purely DMTP-based EOP's
    SampleOEPotential(SharedWavefunction wfn, Options& options);

    // GDF-based EOP's
    SampleOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options);

    // Necessary destructor
    virtual ~SampleOEPotential();

    // Necessary computer 
    virtual void compute(const std::string& oepType) override;

    // Necessary computer
    virtual void compute_3D(const std::string& oepType, 
                            const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) override;
    // Necessary printer
    virtual void print_header() const override;

  private:
    // Set defaults - good practice
    void common_init();

    // Auxilary computers - exemplary
    double compute_3D_sample_V(const double& x, const double& y, const double& z);
};
\endcode
The constructors need to call the abstract base constructor and then specialized initializations.
It is a good practice to put the specialized common initializers in a separate private method
`common_init` (which is a convention in Psi4 and is adopted also in EOPDev).
For instance, the exemplary constructor is show below:
\code{.cpp}
SampleOEPotential::SampleOEPotential(SharedWavefunction wfn, 
                                     SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

void SampleOEPotential::common_init() 
{
   int n1 = wfn_->Ca_subset("AO","OCC")->ncol();
   int n2 = auxiliary_->nbf();
   int n3 = wfn_->molecule()->natom();

   psi::SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G(S^{-1})", n2, n1);
   psi::SharedMatrix mat_2 = std::make_shared<psi::Matrix>("G(S^{-2})", n3, n1);

   OEPType type_1 = {"Murrell-etal.S1", true , n1, mat_1};
   OEPType type_2 = {"Otto-Ladik.S2"  , false, n1, mat_2};

   oepTypes_[type_1.name] = type_1; 
   oepTypes_[type_2.name] = type_2;
}
\endcode
Note that the `OEPotential::oepTypes_` attribute, which is a `std::map` of structures `OEPType`,
is initialized here. All the EOP types need to be stated in the constructors.
Destructors usually call nothing, unless dynamically allocated memory is also of use.

It is also a good practice to already sketch the `compute` method here by adding certain private
computers, like in the example below:
\code{.cpp}
void SampleOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.S1") this->compute_murrell_etal_s1();  // calls private method
  else if (oepType ==   "Otto-Ladik.S2") this->compute_otto_ladik_s2();    // calls private method
  else throw psi::PSIEXCEPTION("EOPDEV: Error. Incorrect EOP type specified!\n"); // for safety
}
void SampleOEPotential::compute_murrell_etal_s1() 
{
   psi::timer_on ("EOP    E(Paul) Murrell-etal S1  ");
   /* Your implementation goes here */
   psi::timer_off("EOP    E(Paul) Murrell-etal S1  ");
\endcode

\subsubsection  ssoepsubtypes Implementing EOP Types

Implementation of the inner body of `compute` method requires populating the members of `oepTypes_` with
data. This means, that for each EOP type there has to be a specific implementation of EOP parameters.
GDF-based EOP's need to create the `psi::Matrix` with EOP parameters and put them into `oepTypes_`.
In the case of DMTP-based EOP's `compute_3D` method has to be additionally implemented before `compute`
is fully functional. To implement `compute_3D`, `OEPotential::make_oeps3d` method is of high
relevance: it creates `OEPotential3D<T>` instances, where `T` is the EOP subclass. These instances
are `Field3D` objects that define EOP's in 3D Euclidean space. For example,
\code{.cpp}
void SampleOEPotential::compute_otto_ladik_s2() 
{
      // Switch on timer
      psi::timer_on("EOP    E(Paul) Otto-Ladik S2    ");

      // Create 3D field, automated through `make_oeps3d`. Requires `compute_3D` implementation.
      std::shared_ptr<OEPotential3D<OEPotential>> oeps3d = this->make_oeps3d("Otto-Ladik.S2");
      oeps3d->compute();

      // Perform ESP fit to get EOP effective charges
      ESPSolver esp(oeps3d);
      esp.set_charge_sums(0.5);
      esp.compute();
    
      // Put the EOP coefficients into `oepTypes_` 
      for (int i=0; i<esp.charges()->nrow(); ++i) {
           for (int o=0; o<oepTypes_["Otto-Ladik.S2"].n; ++o) {
                oepTypes_["Otto-Ladik.S2"].matrix->set(i, o, esp.charges()->get(i, o));
           }
      }

      // Switch off timer
      psi::timer_off("EOP    E(Paul) Otto-Ladik S2    ");
}
// Necessary implementation for `make_oeps3d` to work
void SampleOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) 
{
   // Loop over all possibilities for EOP types and exclude illegal names
   if (oepType == "Otto-Ladik.S2") {

       // this computes the actual values of EOP = v(x,y,z) and stores it in `vec_otto_ladik_s2_`
       this->compute_3D_otto_ladik_s2(x, y, z); 

       // Assign final value to the buffer vector
       for (int o = 0; o < oepTypes_["Otto-Ladik.S2"].n; ++o) v->set(o, vec_otto_ladik_s2_[o]);

   }
   else if (oepType == "Murrell-etal.S1" ) {/* Even if it is not DMTP-based EOP, this line is necessary */}

   else {
      throw psi::PSIEXCEPTION("EOPDEV: Error. Incorrect EOP type specified!\n"); // Safety
   }
}
\endcode
Note that `make_oeps3d` is not overridable and is fully defined in the base. Do not call 
`OEPotential3D` constructors in the OEPotential subclass (it can be done only from the level of the abstract base
where all the pointers are dynamically converted to an appropriate data type due to polymorphism)!


*/

