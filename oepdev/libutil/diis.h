#ifndef _oepdev_libutil_diis_h
#define _oepdev_libutil_diis_h

#include<cstdio>
#include<string>
#include<vector>

#include "psi4/libparallel/parallel.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"


namespace oepdev{

using namespace psi;
using namespace std;


/** \class DIISManager
 *  \brief DIIS manager.
 *
 *  Instance can interact directly with the process of solving 
 *  vector quantities in iterative manner. One needs to pass
 *  the dimensions of solution vector as well as the DIIS subspace size.
 *  The iterative procedure requires providing the current vector
 *  and also an estimate of the error vector. The updated DIIS vector
 *  can be copied to an old vector through the Instance.
 */
class DIISManager {

  // --- Integers ---

  /// Size of DIIS subspace 
  const int _dim;
  /// Order of DIIS matrix equation
  const int _dn;
  /// Number of rows of a solution
  const int _na;
  /// Number of columns of a solution
  const int _nb;
  /// Controller of the current size of vector list
  int _ndiis;

  // --- Vectors and Matrices ---

  /// List of current and previous solution vectors
  std::vector<std::shared_ptr<Matrix>> _vectors;
  /// List of current and previous error vectors
  std::vector<std::shared_ptr<Matrix>> _errors;
  /// Work space
  std::shared_ptr<Matrix> _u;

public:

  /** Constructor.
   *  @param dim Size of DIIS subspace
   *  @param na  Number of solution rows
   *  @param nb  Number of solution columns
   */
  DIISManager(int dim, int na, int nb);

  /// Destructor
 ~DIISManager();
  
  //  --- Methods ---

  /**
   *  Put the current solution to the DIIS manager.
   *  @param error  Shared matrix with current solution error
   *  @param vector Shared matrix with current solution vector
   */
  void put(const std::shared_ptr<const Matrix>& error, 
           const std::shared_ptr<const Matrix>& vector);

  /**
   *  Perform DIIS interpolation.
   */
  void compute(void);

  /**
   *  Update solution vector.
   *  Pass the Shared pointer to current solution. 
   *  Then it will be overriden by the updated DIIS solution.
   */
  void update(std::shared_ptr<Matrix>& other);
};


} // EndNameSpace oepdev

#endif // _oepdev_libutil_diis_h
