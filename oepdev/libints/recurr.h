#ifndef _oepdev_libints_recurr_h
#define _oepdev_libints_recurr_h

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
double d_N_n1_n2(int N, int n1, int n2, double PA, double PB, double aP = 1.0);

} // EndNameSpace oepdev
#endif //_oepdev_libints_recurr_h
