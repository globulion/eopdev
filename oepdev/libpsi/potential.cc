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
 *
 * @END LICENSE
 */

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/cdsalclist.h"
#include "./potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/physconst.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define VDEBUG 1

namespace oepdev{

PotentialInt::PotentialInt(std::vector<psi::SphericalTransform>& st, 
                           std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2, 
                           int deriv)
 : psi::PotentialInt(st, bs1, bs2, deriv)
{

}

PotentialInt::PotentialInt(std::vector<psi::SphericalTransform>& st, 
                           std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2, 
                           std::shared_ptr<psi::Matrix> Qxyz, int deriv)
 : psi::PotentialInt(st, bs1, bs2, deriv)
{
    // Setup the initial field of partial charges
    Zxyz_ = psi::SharedMatrix(new psi::Matrix(Qxyz));
}

PotentialInt::PotentialInt(std::vector<psi::SphericalTransform>& st, 
                           std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2, 
                           const double& x, const double&y, const double&z, const double& q, int deriv)
 : psi::PotentialInt(st, bs1, bs2, deriv) 
{
    // Setup the initial field of partial charges
    std::shared_ptr<psi::Matrix> Qxyz = std::make_shared<psi::Matrix>("Probe charge (q,x,y,z)", 1, 4);
    Qxyz->set(0, 0, q);
    Qxyz->set(0, 1, x);
    Qxyz->set(0, 2, y);
    Qxyz->set(0, 3, z);
    Zxyz_ = psi::SharedMatrix(new psi::Matrix(Qxyz));
}

void PotentialInt::set_charge_field(const double& x, const double&y, const double&z, const double& q)
{
    std::shared_ptr<psi::Matrix> Qxyz = std::make_shared<psi::Matrix>("Probe charge (q,x,y,z)", 1, 4);
    Qxyz->set(0, 0, q);
    Qxyz->set(0, 1, x);
    Qxyz->set(0, 2, y);
    Qxyz->set(0, 3, z);
    Zxyz_ = psi::SharedMatrix(new psi::Matrix(Qxyz));
}



} // EndNameSpace oepdev
