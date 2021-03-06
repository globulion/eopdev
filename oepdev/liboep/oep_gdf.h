#ifndef _oepdev_liboep_oep_gdf_h_ 
#define _oepdev_liboep_oep_gdf_h_ 
/** @file oep_gdf.h */

#include<cstdio>
#include<string>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "../libpsi/integral.h"

namespace oepdev{

using namespace std;

/** \addtogroup OEPDEV_OEPS
 * @{
 */

/**\brief Generalized Density Fitting Scheme. Abstract Base.
 * 
 * Performs the following map:
 * \f[
 *   \hat{v} \left| i \right) \cong \sum_{\eta} G_{\eta i} \left| \eta \right)
 * \f]
 * where \f$ \hat{v} \f$ is the effective one-electron potential (OEP) operator, 
 * \f$ \left| i \right) \f$ is an arbitrary state vector and
 * \f$ \left| \eta \right) \f$ is an auxiliary basis vector.
 * The coefficients \f$ G_{\eta i} \f$ are stored and define the OEP acting on the
 * state \f$ i \f$.
 * The mapping onto the auxiliary space can be done in two ways:
 *  * **Single Density Fit.** \ref singlegdf "This method" requires the auxiliary basis set to be
 *        nearly complete.
 *  * **Double Density Fit.** \ref doublegdf "This method" can be used to arbitrary auxiliary basis sets.
 */
class GeneralizedDensityFit
{
  public: 
    /**\brief Factory for Single GDF Computer.
     * @param bs_auxiliary - auxiliary basis set
     * @param v_vector     - the matrix with \f$ V_{\xi i} \f$ elements
     * @return Generalized Density Fit Computer.
     */
    static std::shared_ptr<GeneralizedDensityFit> build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                                        std::shared_ptr<psi::Matrix> v_vector);

    /**\brief Factory for Double GDF Computer.
     * @param bs_auxiliary    - auxiliary basis set
     * @param bs_intermediate - intermediate basis set
     * @param v_vector        - the matrix with \f$ V_{\epsilon i} \f$ elements
     * @return Generalized Density Fit Computer.
     */
    static std::shared_ptr<GeneralizedDensityFit> build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                                        std::shared_ptr<psi::BasisSet> bs_intermediate,
                                                        std::shared_ptr<psi::Matrix> v_vector);

    /**\brief Factory for Overlap GDF Computer.
     * @param bs_auxiliary    - auxiliary basis set
     * @param bs_intermediate - intermediate basis set
     * @param v_vector        - the matrix with \f$ V_{\epsilon i} \f$ elements
     * @param dummy           - a dummy variable (not used)
     * @return Generalized Density Fit Computer.
     */
    static std::shared_ptr<GeneralizedDensityFit> build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                                        std::shared_ptr<psi::BasisSet> bs_intermediate,
                                                        std::shared_ptr<psi::Matrix> v_vector, int dummy);


    /// Constructor. Initializes the pointers
    GeneralizedDensityFit();

    /// Destructor
    virtual ~GeneralizedDensityFit();
   
    /**\brief Perform the generalized density fit
     * @return The OEP coefficients \f$ G_{\xi i} \f$
     */
    virtual std::shared_ptr<psi::Matrix> compute(void) = 0;

    /// Extract the \f$ G_{\xi i} \f$ coefficients
    std::shared_ptr<psi::Matrix> G(void) const {return G_;}

  protected: 
    /// The OEP coefficients \f$ G_{\xi i} \f$
    std::shared_ptr<psi::Matrix> G_;
    /// The intermediate DF coefficients for \f$ \hat{v}|i) \f$
    std::shared_ptr<psi::Matrix> H_;
    /// The V matrix \f$ \left( \xi \vert \hat{v} i \right) \f$
    std::shared_ptr<psi::Matrix> V_;

    /// Number of auxiliary basis set functions
    int n_a_;
    /// Number of intermediate basis set functions
    int n_i_;
    /// Number of OEP's
    int n_o_;

    /// Basis set: auxiliary
    std::shared_ptr<psi::BasisSet> bs_a_;
    /// Basis set: intermediate
    std::shared_ptr<psi::BasisSet> bs_i_;

    /// Integral factory: aux - aux
    std::shared_ptr<oepdev::IntegralFactory> ints_aa_;
    /// Integral factory: aux - int
    std::shared_ptr<oepdev::IntegralFactory> ints_ai_;
    /// Integral factory: int - int
    std::shared_ptr<oepdev::IntegralFactory> ints_ii_;

    /// Invert a square matrix and check if the inverse is acceptable
    void invert_matrix(std::shared_ptr<psi::Matrix>& M);
};

/**\brief Generalized Density Fitting Scheme - Single Fit.
 * 
 * The density fitting map projects the OEP onto the auxiliary, nearly complete
 * basis set space through application of the resolution of identity.
 * Refer to \ref sdensfitcompl "density fitting in complete space" for more details.
 *
 * \section singlegdf Determination of the OEP matrix
 *
 * Coefficients \f$ {\bf G} \f$ are computed by using the following relation
 * \f[
 *   {\bf G}^{(i)} = {\bf v}^{(i)} \cdot {\bf S}^{-1}
 * \f]
 * where 
 * \f{align*}{
 *   S_{\xi \eta} &= \left( \xi \vert \eta \right) \\
 *   v^{(i)}_\xi  &= \left( \xi \vert \hat{v} i \right)
 * \f}
 * In the above, \f$ \vert \f$ denotes the single integration over electron coordinate, i.e.,
 * \f[
 *    \left( a \vert b \right) \equiv \int d{\bf r} \phi_a^*({\bf r}) \phi_b({\bf r})
 * \f]
 * whereas the spatial form of the potential operator \f$ \hat{v} \f$ can be expressed by
 * \f[
 *    v({\bf r}) \equiv  \int d{\bf r}' \frac{\rho({\bf r}')}{\vert {\bf r}' - {\bf r} \vert}
 * \f]
 * with \f$ \rho({\bf r}) \f$ being the effective one-electron density associated with \f$ \hat{v} \f$.
 */
class SingleGeneralizedDensityFit : public GeneralizedDensityFit
{
  public: 
    SingleGeneralizedDensityFit(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                std::shared_ptr<psi::Matrix> v_vector);
    virtual ~SingleGeneralizedDensityFit();
    std::shared_ptr<psi::Matrix> compute(void);
};

/**\brief Generalized Density Fitting Scheme - Double Fit.
 * 
 * The density fitting map projects the OEP onto an arbitrary (not necessarily complete) auxiliary
 * basis set space through application of the self energy minimization technique.
 * The resulting three-electron repulsion integrals are computed by utilizing the resolution of identity
 * in an intermediate, nearly-complete basis set space, hence performing an internal density fitting
 * in nearly complete basis.
 * Refer to \ref pdensityfitting "density fitting specialized for OEP's" for more details.
 *
 * \section doublegdf Determination of the OEP matrix
 *
 * Coefficients \f$ {\bf G} \f$ are computed by using the following relation
 * \f[
 *   {\bf G} = {\bf A}^{-1} \cdot {\bf R} \cdot {\bf H}
 * \f]
 * where the intermediate projection matrix is given by
 * \f[
 *  {\bf H} = {\bf S}^{-1} \cdot {\bf V}
 * \f]
 * In the above equations,
 * \f{align*}{
 *   A_{\xi \xi'}           &= \left( \xi \vert\vert \xi' \right) \\
 *   R_{\xi \epsilon}       &= \left( \xi \vert\vert \epsilon \right) \\
 *   S_{\epsilon \epsilon'} &= \left( \epsilon \vert \epsilon' \right) \\
 *   V^{\epsilon i}         &= \left( \epsilon \vert \hat{v} i \right) 
 * \f}
 * The following labeling convention is used here:
 *  * \f$ i        \f$ denotes the arbitrary state vector
 *  * \f$ \xi      \f$ denotes the auxiliary basis set element
 *  * \f$ \epsilon \f$ denotes the intermediate (nearly complete) basis set element
 *
 * In the above, \f$ \vert \f$ denotes the single integration over electron coordinate, i.e.,
 * \f[
 *    \left( a \vert b \right) \equiv \int d{\bf r} \phi_a^*({\bf r}) \phi_b({\bf r})
 * \f]
 * whereas \f$ \vert\vert \f$ acts as shown below:
 * \f[
 *    \left( a \vert\vert b \right) \equiv \iint d{\bf r}' d{\bf r}'' 
 *    \frac{\phi_a^*({\bf r}') \phi_b({\bf r}'')}{\vert {\bf r}' - {\bf r}''\vert}
 * \f]
 * The spatial form of the potential operator \f$ \hat{v} \f$ can be expressed by
 * \f[
 *    v({\bf r}) \equiv  \int d{\bf r}' \frac{\rho({\bf r}')}{\vert {\bf r}' - {\bf r} \vert}
 * \f]
 * with \f$ \rho({\bf r}) \f$ being the effective one-electron density associated with \f$ \hat{v} \f$.
 * 
 * \subsection doublegdfdetails Theory behind the double GDF scheme
 *
 * In order to perform the generalized density fitting in an incomplete auxiliary basis set, 
 * one must apply the following formula:
 * \f[
 *   {\bf G} = {\bf A}^{-1} \cdot {\bf B}
 * \f]
 * where one encounters the need of evaluation of the following *three-electron integrals*
 * \f[
 *  B_{\xi i} = \left( \xi \vert\vert \hat{v} i \right) 
 *  \equiv \iiint d{\bf r}' d{\bf r}'' d{\bf r}'''
 *  \phi_\xi^*({\bf r}') 
 *   \frac{1}{\vert {\bf r}' - {\bf r}'' \vert}
 *  \rho({\bf r}''') 
 *   \frac{1}{\vert {\bf r}''' - {\bf r}'' \vert}
 *  \phi_i({\bf r}'')
 * \f]
 * Computation of all the necessery integrals of this kind is very costly and impractical for larger molecules. 
 * However, one can use the same trick that is a kernel of the OEP technique introduced in the OEPDev project, i.e.,
 * introduce the effective potential 
 * in order to get rid of one integration. 
 * This can be done
 * by performing the generalized density fitting in the nearly complete intermediate basis
 * \f[
 *   \hat{v} \left| i \right) \cong \sum_\epsilon H_{\epsilon i} \left| \epsilon \right)
 * \f]
 * Note that this is done just for the sake of factorizing the triple integral and computing the 
 * OEP matrix for the incomplete auxiliary basis. Therefore, the intermediate basis set
 * is used just for a while during density fitting and is no longer necessary later on.
 * By inserting the above identity to the triple integral one can transform it into
 * a sum of the two-electron integrals that are much easier to evaluate. This leads to equations given in the beginning
 * of this section.
 */
class DoubleGeneralizedDensityFit : public GeneralizedDensityFit
{
  public: 
    DoubleGeneralizedDensityFit(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                std::shared_ptr<psi::BasisSet> bs_intermediate,
                                std::shared_ptr<psi::Matrix> v_vector);
    virtual ~DoubleGeneralizedDensityFit();
    std::shared_ptr<psi::Matrix> compute(void);
};


/**\brief Generalized Density Fitting Scheme - Single Fit Based on Minimal Overlap in MO Basis.
 * 
 * The density fitting map projects the OEP onto an arbitrary (not necessarily complete) auxiliary
 * basis set space through application of the basis projection technique.
 * Refer to \ref pdensityfitting "density fitting specialized for OEP's" for more details.
 *
 * \section overlapgdf Determination of the OEP matrix
 *
 * Coefficients \f$ {\bf G} \f$ are computed by using the following relation
 * \f[
 *   {\bf G} = {\bf T}_{{\rm m}\tilde{\rm B}}  \cdot
 *             {\bf S}_{\tilde{\rm B}\tilde{\rm B}}^{-1} \cdot {\bf T}_{{\rm m}\tilde{\rm B}}^\dagger 
 *             \cdot {\bf S}_{\rm mi}
 *             \cdot {\bf G}_{\rm i}
 * \f]
 * where the intermediate projection matrix is given by
 * \f[
 *  {\bf G}_{\rm i} = {\bf S}_{\rm ii}^{-1} \cdot {\bf V}_{\rm i}
 * \f]
 * See \ref sdensfitincomplalt "density fitting of OEPs" for more details regarding the matrices involved
 * in the above expressions.
 * 
 */
class OverlapGeneralizedDensityFit : public GeneralizedDensityFit
{
  public: 
    OverlapGeneralizedDensityFit(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                                 std::shared_ptr<psi::BasisSet> bs_intermediate,
                                 std::shared_ptr<psi::Matrix> v_vector);
    virtual ~OverlapGeneralizedDensityFit();
    std::shared_ptr<psi::Matrix> compute(void);
};

/** @}*/
} // EndNameSpace oepdev

#endif // _oepdev_liboep_oep_gdf_h_ 
