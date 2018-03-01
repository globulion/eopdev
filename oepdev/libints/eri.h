#ifndef _oepdev_libints_eri_h
#define _oepdev_libints_eri_h
/** @file eri.h */

#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"
#include "../libpsi/integral.h"
#include "recurr.h"

namespace oepdev{
using namespace std;

/** \addtogroup OEPDEV_LIBINTS 
 * @{
 */

/**\brief General Two Electron Integral. 
 *
 * Implements the McMurchie-Davidson recursive scheme for all integral
 * types. The integral can be defined for any number of Gaussian centres,
 * thus it is not limited to 2-by-2 four-centre ERI. Currently implemented
 * subtypes are:
 *  - oepdev::ERI_1_1 - 2-centre electron-repulsion integral (i|j)
 *  - oepdev::ERI_2_2 - 4-centre electron-repulsion integral (ij|kl)
 *  - oepdev::ERI_3_1 - 4-centre electron-repulsion integral (ijk|l)
 *
 * \see OEPDEV_LIBINTS
 */
class TwoElectronInt : public TwoBodyAOInt
{
 protected:
    /// Maximum angular momentum
    const int max_am_; 

    /// Maximum number of angular momentum functions
    const int n_max_am_; 

    /// Computes the fundamental: Boys function value at *T* for degree *v*.
    psi::Fjt *fjt_;

    /// Should we use shell pair information?
    bool use_shell_pairs_;

    /// Map of Cartesian components per each am
    const double cartMap_[60];

    /// Double factorial array
    const double df_[8];

    /// Get the angular momentum per Cartesian component
    inline int get_cart_am(int am, int n, int x);

    /// Get the (N,L,M)th McMurchie-Davidson coefficient
    inline double get_R(int N, int L, int M) {return mdh_buffer_R_[R_INDEX(N,L,M,0)];}

    /// Computes the ERI's between three shells.
    virtual size_t compute_doublet(int, int);

    /// Computes the ERI's between three shells.
    virtual size_t compute_triplet(int, int, int);

    /// Computes the ERI's between four shells.
    virtual size_t compute_quartet(int, int, int, int);

    /// Buffer for the McMurchie-Davidson-Hermite R coefficents
    double* mdh_buffer_R_;


 public:
   TwoElectronInt(const IntegralFactory* integral, int deriv, bool use_shell_pairs);

   virtual ~TwoElectronInt();

   /// Compute ERI's between 2 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int);

   /// Compute ERI's between 3 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int, int);

   /// Compute ERI's between 4 shells. Result is stored in buffer.
   virtual size_t compute_shell(int, int, int, int);

   /// Compute ERIs between 4 shells. Result is stored in buffer. 
   /// Only for use with ERI_2_2 and the same basis sets, otherwise shell pairs won't be compatible.
   virtual size_t compute_shell(const psi::AOShellCombinationsIterator&);

   /// Compute first derivatives of ERI's between 2 shells
   virtual size_t compute_shell_deriv1(int, int);

   /// Compute second derivatives of ERI's between 2 shells
   virtual size_t compute_shell_deriv2(int, int);

   /// Compute first derivatives of ERI's between 3 shells
   virtual size_t compute_shell_deriv1(int, int, int);

   /// Compute second derivatives of ERI's between 3 shells
   virtual size_t compute_shell_deriv2(int, int, int);

   /// Compute first derivatives of ERI's between 4 shells
   virtual size_t compute_shell_deriv1(int, int, int, int);

   /// Compute second derivatives of ERI's between 4 shells
   virtual size_t compute_shell_deriv2(int, int, int, int);

};

/**\brief 2-centre ERI of the form (a|O(2)|b) where O(2) = 1/r12.
  *
  * ERI's are computed for a shell doublet (P|Q) and stored in the
  * `target_full_` buffer, accessible through `buffer()` method:
  * \f{align*}{
  *  & \text{For each }  {(n_1,l_1,m_1)\in P}: \\
  *  & \quad\text{For each } {(n_2,l_2,m_2)\in Q}: \\
  *  & \quad\quad{\rm ERI} = (A\vert B)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f}
  * For detailed description of the McMurchie-Davidson scheme, refer to \ref OEPDEV_LIBINTS.
  *
  * \section seri11implementation Implementation
  * 
  * A set of ERI's in a shell is decontracted as
  * \f[
  *  (A\vert B)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = \sum_{ij} c_i(\alpha_1)c_j(\alpha_2)
  *   (i\vert j)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f]
  * where the primitive ERI is given by
  * \f{multline*}{
  *  (i\vert j)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = 
  *     \sum_{N_1=0}^{n_1} 
  *     \sum_{L_1=0}^{l_1} 
  *     \sum_{M_1=0}^{m_1} 
  *     \sum_{N_2=0}^{n_2} 
  *     \sum_{L_2=0}^{l_2} 
  *     \sum_{M_2=0}^{m_2} 
  *      d_{N_1}^{n_1}  
  *      d_{L_1}^{l_1}
  *      d_{M_1}^{m_1}
  *      d_{N_2}^{n_2}
  *      d_{L_2}^{l_2}
  *      d_{M_2}^{m_2}
  *     \left[N_1L_1M_1 \vert N_2L_2M_2\right]
  * \f}
  */
class ERI_1_1 : public TwoElectronInt
{
  protected:
   /// Compute ERI's between 2 shells
   size_t compute_doublet(int, int);

   /// Buffer for McMurchie-Davidson-Hermite coefficents for monomial expansion (shell 1)
   double* mdh_buffer_1_;

   /// Buffer for McMurchie-Davidson-Hermite coefficents for monomial expansion (shell 2)
   double* mdh_buffer_2_;

  public:
   /// Constructor. Use oepdev::IntegralFactory to generate this object
   ERI_1_1(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
   /// Destructor
  ~ERI_1_1();


 private:
   /// Get the D1 coefficient
   double get_D1(int, int, int);
   double get_D2(int, int, int);

};


/**\brief 4-centre ERI of the form (ab|O(2)|cd) where O(2) = 1/r12.
  *
  * ERI's are computed for a shell quartet (PQ|RS) and stored in the
  * `target_full_` buffer, accessible through `buffer()` method:
  * \f{align*}{
  *  & \text{For each }  {(n_1,l_1,m_1)\in P}: \\
  *  & \quad\text{For each } {(n_2,l_2,m_2)\in Q}: \\
  *  & \quad\quad\text{For each } {(n_3,l_3,m_3)\in R}: \\
  *  & \quad\quad\quad\text{For each } {(n_4,l_4,m_4)\in S}: \\
  *  & \quad\quad\quad\quad{\rm ERI} = (AB\vert CD)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f}
  * For detailed description of the McMurchie-Davidson scheme, refer to \ref OEPDEV_LIBINTS.
  *
  * \section seri22implementation Implementation
  * 
  * A set of ERI's in a shell is decontracted as
  * \f[
  *  (AB\vert CD)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = \sum_{ijkl} c_i(\alpha_1)c_j(\alpha_2)c_k(\alpha_3)c_l(\alpha_4)
  *   (ij\vert kl)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f]
  * where the primitive ERI is given by
  * \f{multline*}{
  *  (ij\vert kl)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = E_{ij}(\alpha_1,\alpha_2) E_{kl}(\alpha_3,\alpha_4) \\
  * \times
  *     \sum_{N_1=0}^{n_1+n_2} 
  *     \sum_{L_1=0}^{l_1+l_2} 
  *     \sum_{M_1=0}^{m_1+m_2} 
  *     \sum_{N_2=0}^{n_3+n_4} 
  *     \sum_{L_2=0}^{l_3+l_4} 
  *     \sum_{M_2=0}^{m_3+m_4}
  *      d_{N_1}^{n_1n_2}  
  *      d_{L_1}^{l_1l_2}
  *      d_{M_1}^{m_1m_2}
  *      d_{N_2}^{n_3n_4}
  *      d_{L_2}^{l_3l_4}
  *      d_{M_2}^{m_3m_4}
  *     \left[N_1L_1M_1 \vert N_2L_2M_2\right]
  * \f}
  * In the above equation, the multiplicative constants are given as
  * \f{align*}{
  *  E_{ij}(\alpha_1,\alpha_2) &= \exp{\left[-\frac{\alpha_1\alpha_2}
  *                                       {\alpha_1+\alpha_2}\vert {\bf A}-{\bf B}\vert^2\right]} \\
  *  E_{kl}(\alpha_3,\alpha_4) &= \exp{\left[-\frac{\alpha_3\alpha_4}
  *                                       {\alpha_3+\alpha_4}\vert {\bf C}-{\bf D}\vert^2\right]} \\
  * \f}
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
   ERI_2_2(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
   /// Destructor
  ~ERI_2_2();


 private:
   /// Get the D2 coefficient 
   double get_D12(int, int, int, int);
   /// Get the D2 coefficient
   double get_D34(int, int, int, int);

};

/**\brief 4-centre ERI of the form (abc|O(2)|d) where O(2) = 1/r12.
  *
  * ERI's are computed for a shell quartet (PQR|S) and stored in the
  * `target_full_` buffer, accessible through `buffer()` method:
  * \f{align*}{
  *  & \text{For each }  {(n_1,l_1,m_1)\in P}: \\
  *  & \quad\text{For each } {(n_2,l_2,m_2)\in Q}: \\
  *  & \quad\quad\text{For each } {(n_3,l_3,m_3)\in R}: \\
  *  & \quad\quad\quad\text{For each } {(n_4,l_4,m_4)\in S}: \\
  *  & \quad\quad\quad\quad{\rm ERI} = (ABC\vert D)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f}
  * For detailed description of the McMurchie-Davidson scheme, refer to \ref OEPDEV_LIBINTS.
  *
  * \section seri31implementation Implementation
  * 
  * A set of ERI's in a shell is decontracted as
  * \f[
  *  (ABC\vert D)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = \sum_{ijkl} c_i(\alpha_1)c_j(\alpha_2)c_k(\alpha_3)c_l(\alpha_4)
  *   (ijk\vert l)[\{\alpha\},{\bf n},{\bf l},{\bf m}]
  * \f]
  * where the primitive ERI is given by
  * \f{multline*}{
  *  (ijk\vert l)[\{\alpha\},{\bf n},{\bf l},{\bf m}] = E_{ijk}(\alpha_1,\alpha_2, \alpha_3) \\
  * \times
  *     \sum_{N_1=0}^{n_1+n_2+n_3} 
  *     \sum_{L_1=0}^{l_1+l_2+l_3} 
  *     \sum_{M_1=0}^{m_1+m_2+m_3} 
  *     \sum_{N_2=0}^{n_4} 
  *     \sum_{L_2=0}^{l_4} 
  *     \sum_{M_2=0}^{m_4} 
  *      d_{N_1}^{n_1n_2n_3}  
  *      d_{L_1}^{l_1l_2l_3}
  *      d_{M_1}^{m_1m_2m_3}
  *      d_{N_2}^{n_4}
  *      d_{L_2}^{l_4}
  *      d_{M_2}^{m_4}
  *     \left[N_1L_1M_1 \vert N_2L_2M_2\right]
  * \f}
  * In the above equation, the multiplicative constants are given as
  * \f[
  *  E_{ijk}(\alpha_1,\alpha_2,\alpha_3)  = \exp{\left[-\frac{\alpha_1\alpha_2}
  *                                        {\alpha_1+\alpha_2}\vert {\bf A}-{\bf B}\vert^2\right]} 
  *                                         \exp{\left[-\frac{(\alpha_1+\alpha_2)\alpha_3}
  *                                        {\alpha_1+\alpha_2+\alpha_3}
  *                                         \vert {\bf P}-{\bf C}\vert^2\right]} 
  * \f]
  */
class ERI_3_1 : public TwoElectronInt
{
  protected:
   /// Compute ERI's between 4 shells
   size_t compute_quartet(int, int, int, int);

   /// Buffer for McMurchie-Davidson-Hermite coefficents for trinomial expansion (shells 1, 2 and 3)
   double* mdh_buffer_123_;

   /// Buffer for McMurchie-Davidson-Hermite coefficents for monomial expansion (shell 4)
   double* mdh_buffer_4_;

  public:
   /// Constructor. Use oepdev::IntegralFactory to generate this object
   ERI_3_1(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
   /// Destructor
  ~ERI_3_1();


 private:
   /// Get the D3 coefficient 
   double get_D123(int, int, int, int, int);
   /// Get the D1 coefficient
   double get_D4(int, int, int);

};


/** @}*/ 

} // EndNameSpace oepdev
#endif //_oepdev_libints_eri_h
