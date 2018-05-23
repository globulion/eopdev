#ifndef _oepdev_libutil_esp_h
#define _oepdev_libutil_esp_h
/** @file esp.h */

#include "psi4/libmints/vector.h"
#include "space3d.h"

namespace psi{
   using namespace std;                                  
   using SharedVetor          = std::shared_ptr<Vector>;
}


namespace oepdev{
/** \addtogroup OEPDEV_ESP
 * @{
 */
 
using namespace std;

using SharedScalarField3D = std::shared_ptr<ScalarField3D>;

/** \brief Charges from Electrostatic Potential (ESP). A solver-type class.
 *
 *  Solves the least-squares problem to fit the generalized charges \f$ q_m \f$, that reproduce
 *  the reference generalized potential \f$ v^{\rm ref}({\bf r}) \f$ supplied by the `ScalarField3D` object:
 *  \f[
 *     \int d{\bf r}' \left[ v^{\rm ref}({\bf r}') - \sum_m \frac{q_m}{\left| {\bf r}' - {\bf r}_m \right|} 
 *                    \right]^2  \rightarrow \text{minimize}
 *  \f]
 *  The charges are subject to the following constraint:
 *  \f[
 *      \sum_m q_m = 0
 *  \f]
 *  ### Method description.
 *  \f$ M \f$ generalized charges is found by solving the matrix equation
 *  \f[
 *   \begin{pmatrix}
 *    {\bf A} & 0 \\
 *     0      & 1
 *    \end{pmatrix}^{-1}
 *    \cdot
 *    \begin{pmatrix}
 *    {\bf b}\\ 0
 *    \end{pmatrix}
 *    =
 *    \begin{pmatrix}
 *    {\bf q}\\ \lambda
 *    \end{pmatrix} 
 *  \f]
 *  where the \f$ {\bf A} \f$ matrix of dimension \f$ M\times M \f$ and \bf b} vector or length \f$ M \f$
 *  are given as 
 *  \f{align*}{
 *    A_{mn} &= \sum_i \frac{1}{r_{im} r_{in}} \\
 *    b_m    &= \sum_i \frac{v^{\rm ref}({\bf r}_m)}{r_{im}}
 *  \f}
 *  In the above equation, summations run over all sample points, at which reference potential
 *  is known.
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

    /** \brief Construct from scalar field.
     *
     *  Assume that the centres are on atoms associated with the scalar field.
     *  @param field    - oepdev scalar field object
     */
    ESPSolver(SharedScalarField3D field);

    /** \brief Construct from scalar field.
     *
     *  Solve ESP equations for a custom set of charge distribution centres.
     *  @param field    - oepdev scalar field object
     *  @param centres  - matrix with coordinates of charge distribution centres
     */
    ESPSolver(SharedScalarField3D field, psi::SharedMatrix centres);

    /// Destructor
    virtual ~ESPSolver();

   
    // <--- Accessors ---> //
 
    /// Get the (fit) charges
    virtual psi::SharedVector charges() const {return charges_;}

    /// Get the charge distribution centres
    virtual psi::SharedMatrix centres() const {return centres_;}


    // <--- Computers ---> //

    /// Perform fitting of effective charges
    virtual void compute();


  protected:

    /// Number of fit centres
    const int nCentres_;

    /// Scalar field
    SharedScalarField3D field_;

    /// Charges to be fit
    psi::SharedVector charges_;

    /// Centres, at which fit charges will reside
    psi::SharedMatrix centres_;

  private:
 
    /// ESP b vector
    psi::SharedVector bvec__;

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
