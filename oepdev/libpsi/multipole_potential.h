/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Improvement of psi::EFPMultipolePotentialInt 
 * from original version from Psi4-1.2.1
 * Modification log:
 *   20.08.2020     - Adding two new constructors to EFPMultipolePotentialInt
 *                    enabling calculations of potential integrals
 *                    for a custom set of probe charges.
 * @END LICENSE
 */

#ifndef _oepdev_libpsi_multipole_potential_h_
#define _oepdev_libpsi_multipole_potential_h_
/** @file multipole_potential.h */

#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/efpmultipolepotential.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/gshell.h"
#include "osrecur.h"

namespace psi{
class SphericalTransform;
}
namespace oepdev{

using namespace std;
/** \addtogroup OEPDEV_LIBINTS
 * @{
 */

/** 
 *  \class EFPMultipolePotentialInt
 *  \brief Computes potential integrals.
 */

class EFPMultipolePotentialInt : public psi::OneBodyAOInt {
   protected:
    // OS Recursion for this type of potential integral
    oepdev::ObaraSaikaTwoCenterEFPRecursion_New mvi_recur_;

    // maximum multipole order to compute
    int max_k_;
    //
    bool do_octupoles_;
    int nchunk_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const psi::GaussianShell&, const psi::GaussianShell&) override;

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    EFPMultipolePotentialInt(std::vector<psi::SphericalTransform>&, std::shared_ptr<psi::BasisSet>, std::shared_ptr<psi::BasisSet>,
                          int max_k = 3, int deriv = 0);
    //! Virtual destructor
    ~EFPMultipolePotentialInt() override;
};



/** @}*/
} // EndNameSpace oepdev

#endif
