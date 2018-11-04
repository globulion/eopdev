/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Improvement of psi::PotentialInt 
 * from original version from Psi4-1.1.
 * Modification log:
 *   11.01.2018     - Adding two new constructors to PotentialInt
 *                    enabling calculations of potential integrals
 *                    for a custom set of probe charges.
 *   17.01.2018     - Adding mutator setting field of point charges
 *                    as a single point. This is convenient when 
 *                    computing potential integrals wrt many points
 *                    without building shared matrix of Zxyz_ by the
 *                    user.
 *
 * @END LICENSE
 */

#ifndef _oepdev_libpsi_potential_h_
#define _oepdev_libpsi_potential_h_
/** @file potential.h */

#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/osrecur.h"

namespace oepdev{

using namespace std;
/** \addtogroup OEPDEV_LIBINTS
 * @{
 */

/** 
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 */
class PotentialInt : public psi::PotentialInt
{
  public:
    /** \brief Constructor. Initialize identically like in psi::PotentilInt
     *
     *   @param st     - Spherical transform object
     *   @param bs1    - basis set for first space
     *   @param bs2    - basis set for second space
     *   @param deriv  - derivative level
     */
    PotentialInt(std::vector<psi::SphericalTransform>& st, 
                 std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2,
                 int deriv = 0);

    /** \brief Constructor. Takes an arbitrary collection of charges.
     *
     *   @param st     - Spherical transform object
     *   @param bs1    - basis set for first space
     *   @param bs2    - basis set for second space
     *   @param Qxyz   - matrix with charges and their positions
     *   @param deriv  - derivative level
     */
    PotentialInt(std::vector<psi::SphericalTransform>& st, 
                 std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2, 
                 std::shared_ptr<psi::Matrix> Qxyz, int deriv = 0);

    /** \brief Constructor. Computes potential for one point x, y, z for a test particle of charge q.
     *
     *   @param st     - Spherical transform object
     *   @param bs1    - basis set for first space
     *   @param bs2    - basis set for second space
     *   @param x      - x coordinate of q
     *   @param y      - y coordinate of q
     *   @param z      - z coordinate of q
     *   @param q      - value of the probe charge
     *   @param deriv  - derivative level
     */
    PotentialInt(std::vector<psi::SphericalTransform>&, std::shared_ptr<psi::BasisSet>, std::shared_ptr<psi::BasisSet>, 
                 const double& x, const double&y, const double&z, const double& q = 1.0, int deriv = 0);


    /** \brief Mutator. Set the charge field to be a x, y, z point of charge q.
     *
     *   @param x      - x coordinate of q
     *   @param y      - y coordinate of q
     *   @param z      - z coordinate of q
     *   @param q      - value of the probe charge
     */
    void set_charge_field(const double& x, const double&y, const double&z, const double& q = 1.0);
};

/** @}*/
} // EndNameSpace oepdev

#endif
