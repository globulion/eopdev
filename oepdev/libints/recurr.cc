#include <cmath>
#include <iostream>
#include "recurr.h"

namespace oepdev{

double d_N_n1_n2(int N, int n1, int n2, double PA, double PB, double aP)
{
  if ((N==n1)&&(N==n2)&&(N==0))   return 1.0;
  else if (N<0)                   return 0.0; 
  else if((n1+n2)<N)              return 0.0;
  else {
    if (n1>n2) {
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
//// This is for debug purposes
void make_mdh_D2_coeff_explicit_recursion(int n1, int n2, double aP, double* PA, double* PB, double* buffer)
{
   for (unsigned int x=0; x<3; ++x) {
        double xPA = PA[x];                                                                    
        double xPB = PB[x];
        for (int i = 0; i < n1+1; ++i) {
             for (int j = 0; j < n2+1; ++j) {
                  for (int n = 0; n < n1+n2+1; ++n) {
                       buffer[D2_INDEX(x,i,j,n)] = d_N_n1_n2(n, i, j, xPA, xPB, aP);
                  }
             }
        }
   }
}
void make_mdh_D2_coeff(int n1, int n2, double aPd, double* PA, double* PB, double* buffer)
{
   double xPA = PA[0];
   double yPA = PA[1];
   double zPA = PA[2];
   double xPB = PB[0];
   double yPB = PB[1];
   double zPB = PB[2];
   for (int j=1; j < n2+1; ++j) {
        for (int n = 0; n < n2+1; ++n) {
             //int I = 17*(j-1) + n;
             //double   valx= double(n+1)*buffer[I+1   ] 
             //             + xPB        *buffer[I     ];
             //double   valy= double(n+1)*buffer[I+1378] 
             //             + yPB        *buffer[I+1377];
             //double   valz= double(n+1)*buffer[I+2755] 
             //             + zPB        *buffer[I+2754];
             //if (n>0) {
             //    valx += aPd     *buffer[I-1   ];
             //    valy += aPd     *buffer[I+1377];
             //    valz += aPd     *buffer[I+2754];
             //}
             //buffer[I+17  ] = valx;
             //buffer[I+1394] = valy;
             //buffer[I+2771] = valz;

             double   valx= double(n+1)*buffer[D2_INDEX(0,0,j-1,n+1)] 
                          + xPB        *buffer[D2_INDEX(0,0,j-1,n  )];
             double   valy= double(n+1)*buffer[D2_INDEX(1,0,j-1,n+1)] 
                          + yPB        *buffer[D2_INDEX(1,0,j-1,n  )];
             double   valz= double(n+1)*buffer[D2_INDEX(2,0,j-1,n+1)] 
                          + zPB        *buffer[D2_INDEX(2,0,j-1,n  )];
             if (n>0) {
                 valx += aPd     *buffer[D2_INDEX(0,0,j-1,n-1)];
                 valy += aPd     *buffer[D2_INDEX(1,0,j-1,n-1)];
                 valz += aPd     *buffer[D2_INDEX(2,0,j-1,n-1)];
             }
             buffer[D2_INDEX(0,0,j,n)] = valx;
             buffer[D2_INDEX(1,0,j,n)] = valy;
             buffer[D2_INDEX(2,0,j,n)] = valz;
        }
   }
   for (int i = 1; i < n1+1; ++i) {
        for (int j = 0; j < n2+1; ++j) {
             for (int n = 0; n < n1+n2+1; ++n) {
                  double   valx= double(n+1)*buffer[D2_INDEX(0,i-1,j,n+1)] 
                               + xPA        *buffer[D2_INDEX(0,i-1,j,n  )];
                  double   valy= double(n+1)*buffer[D2_INDEX(1,i-1,j,n+1)] 
                               + yPA        *buffer[D2_INDEX(1,i-1,j,n  )];
                  double   valz= double(n+1)*buffer[D2_INDEX(2,i-1,j,n+1)] 
                               + zPA        *buffer[D2_INDEX(2,i-1,j,n  )];
                  if (n>0) {
                      valx += aPd     *buffer[D2_INDEX(0,i-1,j,n-1)];
                      valy += aPd     *buffer[D2_INDEX(1,i-1,j,n-1)];
                      valz += aPd     *buffer[D2_INDEX(2,i-1,j,n-1)];
                  }
                  buffer[D2_INDEX(0,i,j,n)] = valx;
                  buffer[D2_INDEX(1,i,j,n)] = valy;
                  buffer[D2_INDEX(2,i,j,n)] = valz;
             }
        }
   }
}

void make_mdh_R_coeff(int N, int L, int M, double alpha, double a, double b, double c, double* F, double* buffer)
{
  for (int j=0; j<N+L+M+1; ++j) {
       buffer[R_INDEX(0,0,0,j)] = pow((-2.0*alpha),(double)j) * F[j];
  }
  for (int m=0; m<M; ++m) {
       for (int j=0; j<N+L+M; ++j) {
            double   val  =         c * buffer[R_INDEX(0,0,m  ,j+1)];
            if (m>0) val += (double)m * buffer[R_INDEX(0,0,m-1,j+1)];
            buffer[R_INDEX(0,0,m+1,j)] = val;
       }
  }
  for (int l=0; l<L; ++l) {
       for (int m=0; m<M+1; ++m) {
            for (int j=0; j<N+L+M; ++j) {
                 double   val  =         b * buffer[R_INDEX(0,l  ,m,j+1)];
                 if (l>0) val += (double)l * buffer[R_INDEX(0,l-1,m,j+1)];
                 buffer[R_INDEX(0,l+1,m,j)] = val;
            }
       }
  }
  for (int n=0; n<N; ++n) {
       for (int l=0; l<L+1; ++l) {
            for (int m=0; m<M+1; ++m) {
                 for (int j=0; j<N+L+M; ++j) {
                      double   val  =         a * buffer[R_INDEX(n  ,l,m,j+1)];
                      if (n>0) val += (double)n * buffer[R_INDEX(n-1,l,m,j+1)];
                      buffer[R_INDEX(n+1,l,m,j)] = val;
                 }
            }
       }
  }
}


} //EndNameSpace oepdev
