/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Extension of psi::IntegralFactory
 * from original version from Psi4-1.2.1.
 * Modification log:
 *   19.02.2018     - Creation of ERI_2_2, ERI_3_1, ERI_2_1 and ERI_1_1 objects
 *                    was added. The constructors are the same as in original
 *                    psi::IntegralFactory.
 *
 *   20.08.2020     - Adding calculation of improved EFP multipole potential
 *                    integrals with specifying maximum multipole order.
 *
 * @END LICENSE
 */

#ifndef _oepdev_libpsi_integral_h_
#define _oepdev_libpsi_integral_h_
/** @file integral.h */

#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "multipole_potential.h"

namespace oepdev{

class IntegralFactory;

using namespace std;
/** \addtogroup OEPDEV_LIBINTS
 * @{
 */

class TwoBodyAOInt : public psi::TwoBodyAOInt
{
  protected:
   TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv=0);
   TwoBodyAOInt(const TwoBodyAOInt & rhs);

  public:
   virtual ~TwoBodyAOInt();
   /** \brief Compute two-body two-centre integral and put it into matrix
    *  @param result - matrix where to store (i||j) two-body integrals
    *  @param ibs1   - first basis set axis
    *  @param ibs2   - second basis set axis
    */
   virtual void compute(std::shared_ptr<psi::Matrix>& result, int ibs1 = 0, int ibs2 = 2);
   /**
    * \overload
    */
   virtual void compute(psi::Matrix& result, int ibs1 = 0, int ibs2 = 2);

   virtual size_t compute_shell(int, int, int, int) = 0;  
   virtual size_t compute_shell(int, int, int) = 0;  
   virtual size_t compute_shell(int, int) = 0;  
   virtual size_t compute_shell_deriv1(int, int, int, int) = 0;
   virtual size_t compute_shell_deriv2(int, int, int, int) = 0;
   virtual size_t compute_shell_deriv1(int, int, int) = 0;
   virtual size_t compute_shell_deriv2(int, int, int) = 0;
   virtual size_t compute_shell_deriv1(int, int) = 0;
   virtual size_t compute_shell_deriv2(int, int) = 0;
};

/** 
 *  \class IntegralFactory
 *  \brief Extended IntegralFactory for computing integrals.
 *
 *  In addition to integrals available in Psi4, oepdev::IntegralFactory 
 *  enables to compute also:
 *  - OEI's:
 *    - none at that moment
 *  - ERI's:
 *    - integrals of type (a|b)   - `oepdev::ERI_1_1` 
 *    - integrals of type (ab|c)  - `oepdev::ERI_2_1` 
 *    - integrals of type (abc|d) - `oepdev::ERI_3_1` 
 *    - integrals of type (ab|cd) - `oepdev::ERI_2_2` (also in Psi4 as `psi::ERI`)
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

  /// Returns an improved EFPMultipolePotentialInt
  virtual psi::OneBodyAOInt* ao_efp_multipole_potential_new(int max_k=3, int deriv=0);

  /// Returns an ERI_1_1 integral object
  virtual oepdev::TwoBodyAOInt* eri_1_1(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_2_1 integral object
  virtual oepdev::TwoBodyAOInt* eri_2_1(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_2_2 integral object
  virtual oepdev::TwoBodyAOInt* eri_2_2(int deriv=0, bool use_shell_pairs=false);

  /// Returns an ERI_3_1 integral object
  virtual oepdev::TwoBodyAOInt* eri_3_1(int deriv=0, bool use_shell_pairs=false);

};

/** @}*/
} // EndNameSpace oepdev

#endif
