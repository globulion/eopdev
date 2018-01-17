#ifndef _oepdev_libutil_esp_h
#define _oepdev_libutil_esp_h

#include "psi4/libmints/vector.h"
#include "space3d.h"

namespace psi{
   using namespace std;                                  
   using SharedVetor          = std::shared_ptr<Vector>;
}


namespace oepdev{
 
using namespace std;

using SharedScalarField3D = std::shared_ptr<ScalarField3D>;

/** \brief Charges from Electrostatic Potential (ESP). A solver-type class.
 *
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

} // EndNameSpace oepdev
#endif //_oepdev_libutil_esp_h
