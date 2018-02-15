#ifndef _oepdev_libints_eri_symm_h
#define _oepdev_libints_eri_symm_h

#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"

namespace oepdev{
using namespace std;

/*! \ingroup OEPDEV_LIBINTS */

/**\brief General Two Electron Integral.
  *
  */
class TwoElectronInt : public psi::TwoBodyAOInt
{
 protected:
    /// Maximum cartesian class size.
    ///int max_cart_; 

    /// Computes the fundamental
    psi::Fjt *fjt_;

    /// Should we use shell pair information?
    bool use_shell_pairs_;

    /// Map of Cartesian components per each am
    const double cartMap_[60];

    /// Get the angular momentum per Cartesian component
    inline int get_cart_am(int am, int n, int x);

    /// Computes the ERIs between four shells.
    virtual size_t compute_quartet(int, int, int, int);


 public:
   TwoElectronInt(const psi::IntegralFactory* integral, int deriv, bool use_shell_pairs);

   virtual ~TwoElectronInt();

   /// Compute ERIs between 4 shells. Result is stored in buffer.
   size_t compute_shell(int, int, int, int);

};

/**\brief 4-centre ERI of the form (ab|cd).
  *
  */
class ERI_2_2 : public TwoElectronInt
{
  protected:
   size_t compute_quartet(int, int, int, int);
  public:
   ERI_2_2(const psi::IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
  ~ERI_2_2();
};

} // EndNameSpace oepdev
#endif //_oepdev_libints_eri_symm_h
