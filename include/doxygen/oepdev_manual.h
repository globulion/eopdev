//-----------------------------------------------------------

/*! \page intro Introduction
This page introduces the user to the topic.
Now you can proceed to the \ref advanced "advanced section".
*/

//-----------------------------------------------------------

//-----------------------------------------------------------


/*! \page poepdesign OEP Design.
    \tableofcontents

OEP (One-Electron Potential) is associated with certain quantum one-electron operator 
\f$ \hat{v}^A\f$ that defines the ability of molecule \f$ A \f$ to interact in a particular way with other molecules. 
Technically, OEP can be understood as a __container object__ (associated with the molecule in question)
that stores the information about the above mentioned quantum operator. 
Here, it is assumed that similar OEP
object is also defined for all other molecules in a molecular aggregate. 

In case of interaction between molecules \f$ A \f$ and \f$ B \f$,
OEP object of molecule \f$ A \f$ interacts directly with wavefunction object
of the molecule \f$ B \f$. Defining a 
 * Solver class that handles such interaction 
 * Wavefunction class and
 * OEP class

the universal design of OEP-based approaches can be established and developed.

> **Important:**
>  OEP and Wavefunction classes should not be restricted to Hartree-Fock; in generall any correlated 
>  wavefunction and derived OEP`s should be allowed to work with each other.
>

\section soepclasses OEP Classes

There are many types of OEP’s, but the underlying principle is the same and independent of the
type of intermolecular interaction. Therefore, the OEP’s should be implemented by using a multi-level class design.
In turn, this design depends on the way OEP’s enter the mathematical expressions, i.e., on the types
of matrix elements of the one-electron effective operator \f$ \hat{v}^A \f$.

\subsection ssoepunification Structure of possible OEP-based expressions and their unification

Structure of OEP-based mathematical expressions is listed below:

| Type   |                 Matrix Element                     |       Comment            |
|--------|----------------------------------------------------|--------------------------|
| Type 1 | \f$ \left( I \lvert \hat{v}^A \rvert J \right) \f$ | \f$ I\in A,\; J\in B\f$  |
| Type 2 | \f$ \left( J \lvert \hat{v}^A \rvert L \right) \f$ | \f$ J,L\in B\f$          |


In the above table, \f$ I \f$, \f$ J \f$ and \f$ K \f$ indices correspond 
to basis functions or molecular orbitals. 
Basis functions can be primary or auxiliary OEP-specialized density-fitting.
Depending on the type of function and matrix element, there are many subtypes of resulting matrix elements 
that differ in their dimensionality. Examples are given below:

| Matrix Element | DF-based form | ESP-based form |
|----------------|---------------|----------------|
| \f$ \left(\mu \lvert \hat{v}^{A[\mu]} \rvert \sigma \right) \f$ | \f$ \sum_{\iota\in A} v_{\mu\iota}^AS_{\iota\sigma} \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[\mu]} V_{\mu\sigma}^{(\alpha)}\f$ |
| \f$ \left( i  \lvert \hat{v}^{A[i  ]} \rvert j      \right) \f$ | \f$ \sum_{\iota\in A} v_{i  \iota}^AS_{\iota j    } \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[i  ]} V_{i  j     }^{(\alpha)}\f$ |
| \f$ \left(j   \lvert \hat{v}^{A[i  ]} \rvert l      \right) \f$ | \f$ \sum_{\iota\kappa\in A} S_{j \iota} v_{\iota\kappa}^{A[i]}S_{\kappa l} \f$ |  \f$ \sum_{\alpha\in A} q_\alpha^{A[i]} V_{jl}^{(\alpha)}\f$ |

In the formulae above, the OEP-part (stored by OEP instances) and
the Solver-part (to be computed by the Solver) are separated.
It is apparent that all OEP-parts have the form of 2nd- or 3rd-rank tensors 
with different class of axes (molecular orbitals, primary/auxiliary basis, atomic space). 
Therefore, they can be uniquely defined by a unified *tensor object* 
(storing double precision numbers) and unified *dimension object* storing 
the information of the axes classes.

In Psi4, a perfect candidate for the above is `psi4::Tensor` class declared 
in `psi4/libthce/thce.h`. Except from the numeric content its instances also 
store the information of the dimensions in a form of a vector of `psi4::Dimension` instances.

Another possibility is to use `psi::Matrix` objects, instead of `psi4::Tensor` objects, 
possibly putting them into a `std::vector` container in case there is more than two axes.

*/

/*! \page poeptypes List of One-Electron Potentals
    \tableofcontents

Here I provide the list of OEP’s that have been already derived within 
the scope of the OEPDev project. 

\section soepcoulomb Electrostatic Energy OEP’s

For electrostatic energy calculations, OEP is simply the electrostatic potential 
due to nuclei and electrons. 

3D form:

\f{align*}{
 v({\bf r}) = \sum_x \frac{Z_x}{\vert {\bf r} - {\bf r}_x \vert}
             +\sum_{\mu\nu\in A} P_{\nu\mu} 
              \int d{\bf r'} \frac{\phi^{*}_\mu({\bf r'}) \phi_\nu({\bf r'})}{\vert {\bf r} - {\bf r'} \vert}
\f}

Matrix form:

\f{align*}{
 v_{ik}  &= \sum_{x\in A} Z_x V^{(x)}_{ik} 
           +\sum_{\mu\nu\in A}  P_{\nu\mu} 
           \left( \mu\nu \vert ik \right)
\f}

\section soeppauli Pauli Repulsion OEP’s

The following potentials are derived for the evaluation of the Pauli 
repulsion energy based on Murrel’s expressions.

\subsection ssoeppauli1 First-order contribution in overlap matrix expansion.

This contribution is simply the electrostatic potential coming from all nuclei and electron density
*except* from electron density from molecular orbital \f$ i \f$ that interacts with the 
generalized overlap density between \f$ i \f$ of molecule \f$ A \f$ and \f$ j \f$ of molecule \f$ B \f$.

3D forms:

\f{align*}{
v({\bf r})^{A[i]}_{S^{-1}} &= 
     -\sum_{x\in A} \frac{Z_x}{\vert {\bf r} - {\bf r}_x\vert} 
     + \sum_{\mu\nu\in A} \left\{ D_{\nu\mu} - C^{*}_{\mu i} C_{\nu i} \right \} 
              \int d{\bf r'} \frac{\phi^{*}_\mu({\bf r'}) \phi_\nu({\bf r'})}{\vert {\bf r} - {\bf r'} \vert}
\f}

Matrix forms:

\f{align*}{
v_{\xi i}({S^{-1}}) &= \sum_{\kappa\in A} C_{i\kappa} 
                      \left\{ -\sum_{x\in A} V^{(x)}_{\kappa\xi} 
                              +\sum_{\mu\nu\in A} \left\{ D_{\nu\mu} - C^{*}_{\mu i} C_{\nu i} \right \} 
                          \left( \mu\nu \vert \xi\kappa \right )\right \}
\f}

\subsection ssoeppauli2 Second-order contribution in overlap matrix expansion.

To be added here!

\section soepeet Excitonic Energy Transfer OEP’s

The following potentials are derived for the evaluation of the short-range EET couplings 
based on Fujimoto’s TDFI-TI method.

\subsection ssoepeet1 ET contributions.

3D forms:

\f{align*}{
v({\bf r})^{A[\mu]}_{1} &= -C^*_{\mu L} \sum_{x\in A} \frac{Z_x}{\vert {\bf r} 
                         - {\bf r}_x\vert} + \sum_{\nu\kappa\in A} 
                         \left\{ C^*_{\mu L} D_{\nu\kappa} - \frac{1}{2} C^{*}_{\nu L} D_{\mu\kappa} \right \} 
                 \int d{\bf r'} \frac{\phi^{*}_\nu({\bf r'}) \phi_\kappa({\bf r'})}{\vert {\bf r} - {\bf r'} \vert} \\
v({\bf r})^{A[\mu]}_{2} &= C_{\kappa H} \sum_{\nu\kappa\in A} 
                          \left\{ 2 C^*_{\nu L} C_{\mu H}^* - C^{*}_{\nu H} C_{\mu L}^* \right \} 
                          \int d{\bf r'} 
                          \frac{\phi^{*}_\nu({\bf r'}) \phi_\kappa({\bf r'})}{\vert {\bf r} - {\bf r'} \vert} \\
v({\bf r})^{A[\mu]}_{3} &= v({\bf r})^{A[\mu]}_{1} + v({\bf r})^{A[\mu]}_{1}
\f}

Matrix forms:

\f{align*}{
v_{\mu\xi}(1) &= -C^*_{\mu L} \sum_{x\in A} V^{x}_{\mu\xi} 
                + \sum_{\nu\kappa\in A} \left\{ C^*_{\mu L} D_{\nu\kappa} 
                - \frac{1}{2} C^{*}_{\nu L} D_{\mu\kappa} \right \} 
                \left( \nu\kappa \vert \mu\xi \right ) \\
v_{\mu\xi}(2) &= C_{\kappa H} \sum_{\nu\kappa\in A} 
                 \left\{ 2 C^*_{\nu L} C_{\mu H}^* - C^{*}_{\nu H} C_{\mu L}^* \right \} 
                 \left( \nu\kappa \vert \mu\xi \right) \\
v_{\mu\xi}(3) &= v_{\mu\xi}(1) + v_{\mu\xi}(2)
\f}

\subsection ssoepeet2 HT contributions.

Do be derived.

\subsection ssoepeet3 CT contributions.

To be derived.

\section soephf Full HF Interaction OEP’s

The following potentials are derived for the evaluation 
of the full Hartree-Fock interaction energy based on the OEPDev equations.

//----------------------------------------------------------------------------------

*/

/*!\page pdensityfitting Density-fitting Specialized for OEP’s
   \tableofcontents

To get the ab-initio representation of a OEP, one can use a procedure similar to
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
which we shall refer here as to the matrix form of the OEP operator.
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

Since matrix elements of an OEP operator in auxiliary space can be computed 
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

In the above equation, the OEP-part (fragment parameters for molecule *A* only) 
and the Solver-part (subject to be computed by solver 
on the fly) are separated. This then forms a basis for fragment-based 
approach to solve Quantum Chemistry problems related to the extended molecular aggregates.

\section sdensfitincompl Fitting in Incomplete Space

Density fitting scheme from previous section has practical disadvantage of a nearly-complete basis set
being usually very large (spanned by large amount of basis set vectors). Any non-complete basis set
won't work in the previous example. Since most of basis sets used in quantum chemistry do not form a complete
set, it is beneficial to design a modified scheme in which it is possible to obtain the **effective** 
matrix elements of the OEP operator in a **incomplete** auxiliary space. This can be achieved by minimizing 
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
and \f$ {\bf r}_2 \f$. Thus, it is clear that in order to use this generalized density fitting scheme
one must to compute two-centre electron repulsion integrals (implemented in oepdev::ERI_1_1) 
as well as four-centre asymmetric electron repulsion integrals of the type \f$ (\alpha\beta\gamma||\eta) \f$
(implemented in oepdev::ERI_3_1).
*/

//--------------------------------------------------------------------------------------------

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

//--------------------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------

/*! \page advanced Advanced Usage
This page is for advanced users.
Make sure you have first read \ref intro "the introduction".
*/

