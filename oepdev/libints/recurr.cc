#include <cmath>
#include <iostream>
#include "recurr.h"

namespace oepdev{

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

void make_mdh_D_coeff(int N, int n1, int n2, double* PA, double* PB, double* buffer)
{
   int offs = (n1+1)*(n2+1)*(n1+n2+1);
   for (unsigned int x=0; x<3; ++x) {
        double xPA = PA[x];                                                                    
        double xPB = PB[x];
        for (int i = 0; i < n1+1; ++i) {
             for (int j = 0; j < n2+1; ++j) {
                  for (int n = 0; n < n1+n2+1; ++n) {
                       buffer[D2_INDEX(x,i,j,n)] = d_N_n1_n2(N, i, j, xPA, xPB);
                  }
             }
        }
   }
}

void make_mdh_R_coeff(int N, int L, int M, double alpha, double a, double b, double c, double* F, double* buffer)
{
  for (int j=0; j<N+L+M+1; ++j) {
       buffer[R_INDEX(0,0,0,j)] = pow((-2.0*alpha),(double)j) * F[j];
       //cout << " f:" << F[j] << endl;
       //cout << " i:" << R_INDEX(0,0,0,j) << endl;
       //cout << " b:" << buffer[R_INDEX(0,0,0,j)] << endl;
  }

  for (int m=0; m<M+1; ++m) {
       for (int j=0; j<N+L+M+1; ++j) {
            double   val  =         c * buffer[R_INDEX(0,0,m  ,j+1)];
            if (m>0) val += (double)m * buffer[R_INDEX(0,0,m-1,j+1)];
            buffer[R_INDEX(0,0,m+1,j)] = val;
       }
  }

  for (int l=0; l<L+1; ++l) {
       for (int m=0; m<M+2; ++m) {
            for (int j=0; j<N+L+M+1; ++j) {
                 double   val  =         b * buffer[R_INDEX(0,l  ,m,j+1)];
                 if (l>0) val += (double)l * buffer[R_INDEX(0,l-1,m,j+1)];
                 buffer[R_INDEX(0,l+1,m,j)] = val;
            }
       }
  }
 
  for (int n=0; n<N+1; ++n) {
       for (int l=0; l<L+2; ++l) {
            for (int m=0; m<M+2; ++m) {
                 for (int j=0; j<N+L+M+1; ++j) {
                      double   val  =         a * buffer[R_INDEX(n  ,l,m,j+1)];
                      if (n>0) val += (double)n * buffer[R_INDEX(n-1,l,m,j+1)];
                      buffer[R_INDEX(n+1,l,m,j)] = val;
                 }
            }
       }
  }

}




} //EndNameSpace oepdev
