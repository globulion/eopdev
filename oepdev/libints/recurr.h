#ifndef _oepdev_libints_recurr_h
#define _oepdev_libints_recurr_h
/** @file recurr.h */


namespace oepdev{
using namespace std;

/** \addtogroup OEPDEV_LIBINTS 
 * @{
 */

/*! \def D2_INDEX(x,i,j,n)
    Get the index of McMurchie-Davidson-Hermite D2 coefficient stored in the `mdh_buffer_`, 
    that is attributed to the *x* Cartesian coordinate from angular momenta *i*, *j* 
    of function 1 and 2, and the Hermite index *n*.
*/
#define D2_INDEX(x,i,j,n) ((1377*(x))+(153*(i))+(17*(j))+(n))

/*! \def D3_INDEX(x,i,j,k,n)
    Get the index of McMurchie-Davidson-Hermite D3 coefficient stored in the `mdh_buffer_`, 
    that is attributed to the *x* Cartesian coordinate from angular momenta *i*, *j* and *k* 
    of function 1, 2 and 3, and the Hermite index *n*.
*/
#define D3_INDEX(x,i,j,k,n) ((12393*(x))+(1377*(i))+(153*(j))+(17*(k))+(n))

/*! \def R_INDEX(n,l,m,j)
    Get the index of McMurchie-Davidson R coefficient stored in the `mdh_buffer_R_`
    from angular momenta *n*, *l* and *m* and the Boys index *j*.
*/
#define R_INDEX(n,l,m,j) ((14739*(n))+(867*(l))+(51*(m))+(j))

/**\brief Compute McMurchie-Davidson-Hermite (MDH) coefficient for binomial expansion.
 *
 * @param N  - increment in the summation of MDH series
 * @param n1 - angular momentum of first function
 * @param n2 - angular momentum of second function
 * @param PA - cartesian component of P-A distance 
 * @param PB - cartesian component of P-B distance
 * @param aP - free parameter of MDH expansion
 * @return the McMurchie-Davidson-Hermite coefficient
 */
double d_N_n1_n2(int N, int n1, int n2, double PA, double PB, double aP);

/**\brief Compute the McMurchie-Davidson-Hermite coefficients for binomial expansion.
 * @param n1     - angular momentum of first function
 * @param n2     - angular momentum of second function
 * @param PA     - cartesian components of P-A distance 
 * @param PB     - cartesian components of P-B distance
 * @param buffer - the McMurchie-Davidson-Hermite 4-dimensional array (raveled to vector):
 *         - axis 0: dimension 3 (x, y or z Cartesian component)
 *         - axis 1: dimension n1+1 (0 to n1)
 *         - axis 2: dimension n2+1 (0 to n2)
 *         - axis 3: dimension n1+n2+1 (0 to n1+n2)
 *  \see D2_INDEX
 */
void make_mdh_D2_coeff(int n1, int n2, double aPd, double* PA, double* PB, double* buffer);

/**\brief Compute the McMurchie-Davidson-Hermite coefficients for binomial expansion by explicit recursion.
 * This function makes the same changes to buffers as oepdev::make_mdh_D2_coeff, but
 * implements it through explicit recursion by calls to oepdev::d_N_n1_n2. Therefore,
 * it is slightly slower. Here for debugging purposes.
 * @param n1     - angular momentum of first function
 * @param n2     - angular momentum of second function
 * @param aPd    - parameter equal to 0.500/Pa where Pa is exponent
 * @param PA     - cartesian components of P-A distance 
 * @param PB     - cartesian components of P-B distance
 * @param buffer - the McMurchie-Davidson-Hermite 4-dimensional array (raveled to vector):
 *         - axis 0: dimension 3 (x, y or z Cartesian component)
 *         - axis 1: dimension n1+1 (0 to n1)
 *         - axis 2: dimension n2+1 (0 to n2)
 *         - axis 3: dimension n1+n2+1 (0 to n1+n2)
 *  \see D2_INDEX
 */
void make_mdh_D2_coeff_explicit_recursion(int n1, int n2, double aP, double* PA, double* PB, double* buffer);

/**\brief Compute the McMurchie-Davidson R coefficients.
 * @param N      - increment in the summation of MDH series along *x* direction
 * @param L      - increment in the summation of MDH series along *y* direction
 * @param M      - increment in the summation of MDH series along *z* direction
 * @param alpha  - alpha parameter of R coefficient
 * @param a      - *x* component of PQ vector of R coefficient
 * @param b      - *y* component of PQ vector of R coefficient
 * @param c      - *z* component of PQ vector of R coefficient
 * @param F      - array of Boys function values for given alpha and PQ
 * @param buffer - the McMurchie-Davidson 4-dimensional array (raveled to vector):
 *         - axis 0: dimension N+1
 *         - axis 1: dimension L+1
 *         - axis 2: dimension M+1 
 *         - axis 3: dimension N+L+M+1 (*j*-th element)
 */
void make_mdh_R_coeff(int N, int L, int M, double alpha, double a, double b, double c, double* F, double* buffer);

/** @}*/

} // EndNameSpace oepdev
#endif //_oepdev_libints_recurr_h
