#ifndef _oepdev_libutil_esp_h
#define _oepdev_libutil_esp_h
/** @file esp.h */

#include "space3d.h"

namespace oepdev{
/** \addtogroup OEPDEV_ESP
 * @{
 */
 
using namespace std;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

/** \brief Charges from Electrostatic Potential (ESP). A solver-type class.
 *
 *  Solves the least-squares problem to fit the generalized charges \f$ q_{m;p} \f$, that reproduce
 *  the reference generalized potential \f$ v_p^{\rm ref}({\bf r}) \f$ supplied by the `Field3D` object:
 *  \f[
 *     \int d{\bf r}' \left[ v_p^{\rm ref}({\bf r}') - \sum_m \frac{q_{m;p}}{\left| {\bf r}' - {\bf r}_m \right|} 
 *                    \right]^2  \rightarrow \text{minimize}
 *  \f]
 *  The charges are subject to the following constraint:
 *  \f[
 *      \sum_m q_{m;p} = Q_p \text{ for all $p$}
 *  \f]
 *  ### Method description.
 *  \f$ M \f$ generalized charges is found by solving the matrix equation
 *  \f[
 *   \begin{pmatrix}
 *    {\bf A} & 1 \\
 *     1      & 0
 *    \end{pmatrix}^{-1}
 *    \cdot
 *    \begin{pmatrix}
 *    {\bf b}_p\\ Q_p
 *    \end{pmatrix}
 *    =
 *    \begin{pmatrix}
 *    {\bf q}_p\\ \lambda
 *    \end{pmatrix} 
 *  \f]
 *  where the \f$ {\bf A} \f$ matrix of dimension \f$ (M+1)\times (M+1) \f$ and \f$ {\bf b}_p \f$ vector or length \f$ M+1 \f$
 *  are given as 
 *  \f{align*}{
 *    A_{mn} &= \sum_i \frac{1}{r_{im} r_{in}} \\
 *    b_{m;p}&= \sum_i \frac{v_p^{\rm ref}({\bf r}_m)}{r_{im}}
 *  \f}
 *  In the above equation, summations run over all sample points, at which reference potential
 *  is known. The solution is stored in the \f$ M \times N \f$ matrix, where \f$ N \f$ is the dimensionality
 *  of the 3D vector field (i.e., the number of potentials supplied, \f$ p_{\rm max} \f$).
 *  As a default, \f$ Q_p = 0 \f$ for all potentials. This can be set by `oepdev::ESPSolver::set_charge_sums` method.
 * 
 * \note Useful options:
 *  - `ESP_PAD_SPHERE`        - Padding spherical radius for random points selection. Default: `10.0` [A.U.]
 *  - `ESP_NPOINTS_PER_ATOM`  - Number of random points per atom in a molecule. Detault: `1500`
 *  - `ESP_VDW_RADIUS_C`      - The vdW radius for carbon atom. Default: `3.0` [A.U.]
 *  - `ESP_VDW_RADIUS_H`      - The vdW radius for hydrogen atom. Default: `4.0` [A.U.]
 *  - `ESP_VDW_RADIUS_N`      - The vdW radius for nitrogen atom. Default: `2.4` [A.U.]
 *  - `ESP_VDW_RADIUS_O`      - The vdW radius for oxygen atom. Default: `5.6` [A.U.]
 *  - `ESP_VDW_RADIUS_F`      - The vdW radius for fluorium atom. Default: `2.3` [A.U.]
 *  - `ESP_VDW_RADIUS_CL`     - The vdW radius for chlorium atom. Default: `2.9` [A.U.]
 */
class ESPSolver
{
  public:

    // <--- Constructors and Destructor ---> //

    /** \brief Construct from 3D vector field.
     *
     *  Assume that the centres are on atoms associated with the 3D vector field.
     *  @param field    - oepdev 3D vector field object
     */
    ESPSolver(SharedField3D field);

    /** \brief Construct from 3D vector field.
     *
     *  Solve ESP equations for a custom set of charge distribution centres.
     *  @param field    - oepdev 3D vector field object
     *  @param centres  - matrix with coordinates of charge distribution centres
     */
    ESPSolver(SharedField3D field, psi::SharedMatrix centres);

    /// Destructor
    virtual ~ESPSolver();

   
    // <--- Accessors ---> //
 
    /// Get the (fit) charges
    virtual psi::SharedMatrix charges() const {return charges_;}

    /// Get the charge distribution centres
    virtual psi::SharedMatrix centres() const {return centres_;}


    // <--- Mutators ---> //

    /// Set the charge sums \f$ Q_p \f$
    virtual void set_charge_sums(psi::SharedVector s);

    /// Set the charge sums \f$ Q_p \f$ (equal to all fields)
    virtual void set_charge_sums(const double& s);


    // <--- Computers ---> //

    /// Perform fitting of effective charges
    virtual void compute();


  protected:

    /// Number of fit centres
    const int nCentres_;

    /// Number of fields to fit
    const int nFields_;

    /// Scalar field
    SharedField3D field_;

    /// Charges to be fit
    psi::SharedMatrix charges_;

    /// Centres, at which fit charges will reside
    psi::SharedMatrix centres_;

    /// Vector of sums of partial charges
    psi::SharedVector charge_sums_;

  private:
 
    /// ESP b vector
    psi::SharedMatrix bvec__;

    /// ESP A matrix 
    psi::SharedMatrix amat__;

    /// Initialize defaults
    void common_init();

    /// Compute all the working matrices
    void compute_matrices(void);

    /// Perform fitting
    void fit(void);
};

/** @}*/

} // EndNameSpace oepdev
#endif //_oepdev_libutil_esp_h
