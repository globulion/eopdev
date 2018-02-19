#ifndef _oepdev_libints_eri_symm_h
#define _oepdev_libints_eri_symm_h

#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"
#include "recurr.h"

namespace oepdev{
using namespace std;

/*! \ingroup OEPDEV_LIBINTS */

/**\brief General Two Electron Integral. 
 *
 * The integral can be defined for any number of Gaussian centres,
 * thus it is not limited to 2-by-2 four-centre ERI. 
 */
class TwoElectronInt : public psi::TwoBodyAOInt
{
 protected:
    /// Maximum angular momentum
    const int max_am_; 

    /// Maximum number of angular momentum functions
    const int n_max_am_; 

    /// Computes the fundamental
    psi::Fjt *fjt_;

    /// Should we use shell pair information?
    bool use_shell_pairs_;

    /// Map of Cartesian components per each am
    const double cartMap_[60];

    /// Get the angular momentum per Cartesian component
    inline int get_cart_am(int am, int n, int x);

    /// Get the (N,L,M)th McMurchie-Davidson coefficient
    inline double get_R(int N, int L, int M) {return mdh_buffer_R_[R_INDEX(N,L,M,0)];}

    /// Computes the ERI's between four shells.
    virtual size_t compute_quartet(int, int, int, int);

    /// Computes the ERI's between three shells.
    virtual size_t compute_triplet(int, int, int);

    /// Computes the ERI's between three shells.
    virtual size_t compute_doublet(int, int);

    /// Buffer for McMurchie-Davidson-Hermite R coefficents
    double* mdh_buffer_R_;


 public:
   TwoElectronInt(const psi::IntegralFactory* integral, int deriv, bool use_shell_pairs);

   virtual ~TwoElectronInt();

   /// Compute ERI's between 4 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int, int, int);

   /// Compute ERI's between 3 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int, int);

   /// Compute ERI's between 2 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int);

   /// Compute ERIs between 4 shells. Result is stored in buffer.
   virtual size_t compute_shell(const psi::AOShellCombinationsIterator&) {};

   virtual size_t compute_shell_deriv1(int, int, int, int) {};
   virtual size_t compute_shell_deriv2(int, int, int, int) {};

};

/**\brief 4-centre ERI of the form (ab|cd).
  *
  */
class ERI_2_2 : public TwoElectronInt
{
  protected:
   /// Compute ERI's between 4 shells
   size_t compute_quartet(int, int, int, int);

   /// Buffer for McMurchie-Davidson-Hermite coefficents for binomial expansion (shells 1 and 2)
   double* mdh_buffer_12_;

   /// Buffer for McMurchie-Davidson-Hermite coefficents for binomial expansion (shells 3 and 4)
   double* mdh_buffer_34_;

  public:
   /// Constructor. Use oepdev::IntegralFactory to generate this object
   ERI_2_2(const psi::IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
   /// Destructor
  ~ERI_2_2();


 private:
   double get_D(int, int, int);
   void put_int(double, double*);

};

} // EndNameSpace oepdev
#endif //_oepdev_libints_eri_symm_h
