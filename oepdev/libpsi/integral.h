/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Extension of psi::IntegralFactory
 * from original version from Psi4-1.1.
 * Modification log:
 *   19.02.2018     - Creation of ERI_2_2, ERI_3_1, ERI_2_1 and ERI_1_1 objects
 *                    was added. The constructors are the same as in original
 *                    psi::IntegralFactory.
 *
 * @END LICENSE
 */

#ifndef _oepdev_libpsi_integral_h_
#define _oepdev_libpsi_integral_h_
/** @file integral.h */

#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"

namespace oepdev{

using namespace std;
/** \addtogroup OEPDEV_LIBINTS
 * @{
 */

/** 
 *  \class IntegralFactory
 *  \brief Extended IntegralFactory for computing integrals.
 *
 *  In addition to integrals available in Psi4, oepdev::IntegralFactory 
 *  enables to compute additionally:
 *   - `ERI_2_2` integrals
 *   - `ERI_3_1` integrals
 *   - `ERI_2_1` integrals
 *   - `ERI_1_1` integrals
 */
class IntegralFactory : public psi::IntegralFactory
{
  public:

  // ---> Constructors and Descructor <--- //

    /** \brief Initialize integral factory given a BasisSet for each center.
     *  Becomes (bs1 bs2 | bs3 bs4).
     */
  IntegralFactory(std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2,
                  std::shared_ptr<psi::BasisSet> bs3, std::shared_ptr<psi::BasisSet> bs4);
    /** \brief Initialize integral factory given a BasisSet for two centres. 
     *  Becomes (bs1 bs2 | bs1 bs2).
     */
  IntegralFactory(std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2);
    /** \brief Initialize integral factory given a BasisSet for two centres. 
     *  Becomes (bs1 bs1 | bs1 bs1).
     */
  IntegralFactory(std::shared_ptr<psi::BasisSet> bs1);

  /// Destructor
  virtual ~IntegralFactory();


  // ---> Computers <--- //

  /// Returns an ERI_2_2 integral object
  virtual psi::TwoBodyAOInt* eri_2_2(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_3_1 integral object
  virtual psi::TwoBodyAOInt* eri_3_1(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_2_1 integral object
  virtual psi::TwoBodyAOInt* eri_2_1(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_1_1 integral object
  virtual psi::TwoBodyAOInt* eri_1_1(int deriv=0, bool use_shell_pairs=false);

};

/** @}*/
} // EndNameSpace oepdev

#endif
