#ifndef _oepdev_libints_recurr_h
#define _oepdev_libints_recurr_h

/*! \def N1_N2_N_TO_D(i,j,n)
    Get the index of McMurchie-Davidson-Hermite coefficient stored in the `mdh_buffer_`
    from angular momenta *i*, *j* of function 1 and 2, and the Hermite index *n*.
*/
#define N1_N2_N_TO_D(i,j,n) (((n1+n2+1)*(n2+1)*(i))+((n1+n2+1)*(j))+(n))

namespace oepdev{
using namespace std;



/**\brief Compute Boys function value at *x* for degree *v*.
  *
  * @param v - degree of the Boys function
  * @param x - argument of the Boys function
  * @return the value of the Boys function
  */
double boys(int v, double x);

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
inline double d_N_n1_n2(int N, int n1, int n2, double PA, double PB, double aP = 1.0);

/**\brief Compute the McMurchie-Davidson-Hermite coefficients for binomial expansion.
 * @param ix     - Cartesian component
 * @param N      - increment in the summation of MDH series 
 * @param n1     - angular momentum of first function
 * @param n2     - angular momentum of second function
 * @param PA     - cartesian component of P-A distance 
 * @param PB     - cartesian component of P-B distance
 * @param buffer - the McMurchie-Davidson-Hermite 4-dimensional array (raveled to vector):
 *         - axis 0: dimension 3 (x, y or z Cartesian component)
 *         - axis 1: dimension n1+1 (0 to n1)
 *         - axis 2: dimension n2+1 (0 to n2)
 *         - axis 3: dimension n1+n2+1 (0 to n1+n2)
 *  \see N1_N2_N_TO_D
 */
void make_mdh_coeff(int ix, int N, int n1, int n2, double PA, double PB, double* buffer);

} // EndNameSpace oepdev
#endif //_oepdev_libints_recurr_h
