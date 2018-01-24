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

| Type  | Matrix Element | Comment |
|--------|----|---|
| Type 1 | \f$ \left( I \lvert \hat{v}^A \rvert J \right) \f$ | \f$ I\in A,\; J\in B\f$  |
| Type 2 | \f$ \left( J \lvert \hat{v}^A \rvert L \right) \f$ | \f$ J,L\in B\f$  |


In the above table, \f$ I \f$, \f$ J \f$ and \f$ K \f$ indices correspond 
to basis functions or molecular orbitals. 
Basis functions can be primary or auxiliary OEP-specialized density-fitting.
Depending on the type of function and matrix element, there are many subtypes of resulting matrix elements 
that differ in their dimensionality. Examples are given below:

| Matrix Element | DF-based form | ESP-based form |
|----|---|---|

In the formulae above, the OEP-part (stored by OEP instances) is shown in blue 
whereas the Solver-part (to be computed by the Solver) is shown in brown. 
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



*/

/*!\page pdensityfitting Density-fitting specialized for OEP’s
   \tableofcontents

To get the ab-initio representation of a OEP, one can use a procedure similar to
the typical density fitting or resolution of identity, both of which are nowadays widely used 
to compute electron-repulsion integrals (ERI’s) more efficiently. 

An arbitrary one-electron potential of molecule *A* acting on any state vector 
associated with molecule *A* can be expanded in the *auxiliary space* centered 
on *A* as
\f[
   v\vert i) = \sum_{\xi} v\vert \xi) ( \xi \vert i)
\f]
under the necessary assumption that the auxiliary basis set is *complete*. 
In that case, formally one can write the following identity
\f[
 (\eta\vert v\vert i) = \sum_{\xi} (\eta\vert v\vert \xi) S_{\xi i}
\f]
The matrix elements of the OEP operator in auxiliary space can be computed 
in the same way as the matrix elements with any other basis function. 
In reality, it is almost impossible to reach the completness of the basis set, 
however, but it is possible to obtain the **effective** matrix elements of 
the OEP operator in auxiliary space, rather than compute them as they are 
in the above equation explicitly. We expand the LHS of the first equation 
on this page in a series of the auxiliary basis functions scaled by the 
undetermined expansion coefficients: 
\f[
  v\vert i) = \sum_{\xi} {G_{i\xi}} \vert \xi)
\f]
The expansion coefficients are the effective matrix elements 
of the OEP operator in auxiliary basis set. Now, multiplying both sides 
by another auxiliary basis function and subsequently inverting the equation 
one obtains the expansion coefficients:
\f[
   \boxed{ G_{i\eta} = \sum_\xi (i \vert v \vert \eta) \left[ {\bf S}^{-1} \right]_{\eta\xi} }
\f]
In this way, it is possible to approximately determine the matrix elements 
of the OEP operator with any other basis function in case the auxiliary 
basis set is not complete. **In particular**, when the other basis function 
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
  
*/

#ifndef _oepdev_liboep_liboep_h_ 
#define _oepdev_liboep_liboep_h_ 

#include<cstdio>
#include<string>
#include<vector>
#include<map>

#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libthce/thce.h"
#include "psi4/libcubeprop/csg.h"

#include "../libutil/space3d.h"

namespace oepdev{

using namespace psi;
using namespace std;
using SharedWavefunction = std::shared_ptr<Wavefunction>;
using SharedBasisSet     = std::shared_ptr<BasisSet>;
using SharedTensor       = std::shared_ptr<Tensor>;
using SharedMatrix       = std::shared_ptr<Matrix>;
using SharedVector       = std::shared_ptr<Vector>;


/** \brief Generalized One-Electron Potential: Abstract base.
 * 
 *  Manages OEP's in matrix and 3D forms.
 */
class OEPotential : public std::enable_shared_from_this<OEPotential> 
{

  protected:

    /// Psi4 options
    Options options_;
    /// Wavefunction
    SharedWavefunction wfn_;
    /// Promary Basis set
    SharedBasisSet primary_;
    /// Auxiliary Basis set
    SharedBasisSet auxiliary_;

    /// Name of this OEP;
    std::string name_;
    /// Types of OEP's within the scope of this object
    std::vector<std::string> oepTypes_;
    /// OEP's matrix forms for each OEP type
    std::map<std::string, SharedMatrix> oepMatrices_;

    /// Integral factory
    std::shared_ptr<psi::IntegralFactory> intsFactory_;
    /// Matrix of potential one-electron integrals
    std::shared_ptr<psi::Matrix> potMat_;
    /// One-electron integral shared pointer
    std::shared_ptr<psi::OneBodyAOInt> OEInt_;
    /// One-electron potential shared pointer
    std::shared_ptr<     PotentialInt> potInt_;

  public:

    // <--- Constructors and Destructor ---> //

    /**\brief ESP-based OEP object
     * 
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, Options& options);

    /**\brief DF-based OEP object
     *
     * @param wfn        - wavefunction
     * @param auxiliary  - basis set for density fitting of OEP's
     * @param options    - Psi4 options
     */
    OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);

    /// Destructor
    virtual ~OEPotential();


    // <--- Factories ---> //      

    /**\brief Build ESP-based OEP object
     * 
     * @param type    - OEP category
     * @param wfn     - wavefunction
     * @param options - Psi4 options
     */

    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, Options& options);
    /**\brief Build DF-based OEP object
     *
     * @param type       - OEP category
     * @param wfn        - wavefunction
     * @param auxiliary  - basis set for density fitting of OEP's
     * @param options    - Psi4 options
     */
    static std::shared_ptr<OEPotential> build(const std::string& category, SharedWavefunction wfn, 
                                              SharedBasisSet auxiliary, Options& options);


    // <--- Characterizers ---> //

    /// Is this OEP density-fitted? 
    const bool is_density_fitted;
    /// Is this OEP ESP-based? 
    const bool is_esp_based;


    // <--- Methods/Computers ---> //

    //@{ Compute One Electron Effective Coefficients
    virtual void compute(const std::string& oepType) = 0;
    virtual void compute(void);
    //@}
    //@{ Compute 3D potential
    /** Write potential to a cube file */
    virtual void write_cube(const std::string& oepType, const std::string& fileName);
    /** Compute value of potential in point x, y, z and save at v */
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) = 0;
    //@}
    /// Rotate 
    virtual void rotate(const Matrix& rotmat);
    /// Translate
    virtual void translate(const Vector& trans);
    /// Superimpose
    virtual void superimpose(const Matrix& refGeometry, 
                             const std::vector<int>& supList, 
                             const std::vector<int>& reordList);


    // <--- Accessors ---> //

    /// Retrieve name of this OEP
    std::string name() const { return name_; }
    /// Retrieve matrix potential
    SharedMatrix matrix(const std::string& oepType) const { return oepMatrices_.at(oepType); }
    /// Retrieve wavefunction object
    SharedWavefunction wfn() const {return wfn_;}


    // <--- Mutators ---> //
    void set_name(const std::string& name) {name_ = name;}


    // <--- Printers ---> //
    virtual void print_header() const = 0;


  private:

    /// Initialize defaults
    void common_init();

};


/**\brief Generalized One-Electron Potential for Electrostatic Energy calculations.
 * 
 *  Contains the following OEP types:
 *      "V"
 */
class ElectrostaticEnergyOEPotential : public OEPotential 
{
  public:
    /// Only ESP-based potential is worth implementing
    ElectrostaticEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~ElectrostaticEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();
};

/**\brief Generalized One-Electron Potential for Pauli repulsion energy calculations.
 * 
 *  Contains the following OEP types:
 */
class RepulsionEnergyOEPotential : public OEPotential 
{
  public:
    RepulsionEnergyOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~RepulsionEnergyOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();
};

/**\brief Generalized One-Electron Potential for EET coupling calculations.
 * 
 *  Contains the following OEP types:
 *    "ET1" "ET2" "HT1" "HT1" "HT2" "CT1" "CT2"
 */
class EETCouplingOEPotential : public OEPotential 
{
  public:
    EETCouplingOEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options);
    EETCouplingOEPotential(SharedWavefunction wfn, Options& options);

    virtual ~EETCouplingOEPotential();

    virtual void compute(const std::string& oepType) override;
    virtual void compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) override;
    virtual void print_header() const override;

  private:
    /// Set defaults
    void common_init();
};


} // EndNameSpace oepdev

#endif // _oepdev_liboep_liboep_h_ 
