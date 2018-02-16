#include "recurr.h"

namespace oepdev{

double boys(int v, double x) 
{
 /// TODO
 return 0.0;
}
double d_N_n1_n2(int N, int n1, int n2, double PA, double PB, double aP)
{
  if      (N==n1==n2==0)   return 1.0;
  else if (N<0)            return 0.0;
  else if((n1+n2)<N)       return 0.0;
  else {
    if (n1>=n2) {
       return double(N+1) * d_N_n1_n2(N+1, n1-1, n2, PA, PB, aP) +
                     PA   * d_N_n1_n2(N  , n1-1, n2, PA, PB, aP) +
                    0.5/aP* d_N_n1_n2(N-1, n1-1, n2, PA, PB, aP);
    } 
    else {
       return double(N+1) * d_N_n1_n2(N+1, n1, n2-1, PA, PB, aP) +
                     PB   * d_N_n1_n2(N  , n1, n2-1, PA, PB, aP) +
                    0.5/aP* d_N_n1_n2(N-1, n1, n2-1, PA, PB, aP);
    }
  }
}

void make_mdh_coeff(int ix, int N, int n1, int n2, double PA, double PB, double* buffer)
{
   int offs = (n1+1)*(n2+1)*(n1+n2+1);
   double xPA = PA[x];
   double xPB = PB[x];
   for (int i = 0; i < n1+1; ++i) {
        for (int j = 0; j < n2+1; ++j) {
             for (int n = 0; n < n1+n2+1; ++n) {
                  buffer[offs * ix + N1_N2_N_TO_D(i,j,n)] = d_N_n1_n2(N, i, j, xPA, xPB);
             }
        }
   }
}



} //EndNameSpace oepdev
