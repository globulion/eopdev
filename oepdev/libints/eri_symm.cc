#include "psi4/libmints/gshell.h"
#include "recurr.h"
#include "eri_symm.h"

namespace oepdev{
using namespace std;

#define LOCAL_2PI52 34.986836655249724969962699105963

TwoElectronInt::TwoElectronInt(const psi::IntegralFactory* integral, int deriv, bool use_shell_pairs) :
   psi::TwoBodyAOInt(integral, deriv), 
   use_shell_pairs_(use_shell_pairs),
   cartMap_{0,0,0,
            1,0,0, 0,1,0, 0,0,1,
            2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1,
            3,0,0, 0,3,0, 0,0,3, 1,2,0, 2,1,0, 2,0,1, 1,0,2, 0,1,2, 0,2,1, 1,1,1},
   max_am_(std::max({basis1()->max_am(),
                     basis2()->max_am(),
                     basis3()->max_am(),
                     basis4()->max_am()}))

{
    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                  INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    try {
        target_full_ = new double[size];
        target_ = target_full_;
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double) * size);

    try {
        source_full_ = new double[size];
        source_ = source_full_;
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating source_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(source_, 0, sizeof(double) * size);

}

TwoElectronInt::~TwoElectronInt()
{
    delete[] target_full_;
    delete[] source_full_;
}

int TwoElectronInt::get_cart_am(int am, int n, int x)
{
  int n_prev;
  if      (am==0) n_prev = 0;
  else if (am==1) n_prev = 3;
  else if (am==2) n_prev = 12;
  else if (am==3) n_prev = 30;
  else if (am==4) throw psi::PSIEXCEPTION("oepdev::TwoElectronInt: AM 4 not implemented!");
  else throw psi::PSIEXCEPTION("oepdev::TwoElectronInt: wrong AM chosen!");
  return n_prev + 3*n + x;
}

size_t TwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
 return compute_quartet(sh1, sh2, sh3, sh4);
}

size_t TwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
 return 0;
}


//////
ERI_2_2::ERI_2_2(const psi::IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    if (deriv_>0) throw psi::PSIEXCEPTION("oepdev::ERI_2_2: Derivatives are not implemented yet!");

    // The +1 is needed for derivatives to work.
    fjt_ = new psi::Taylor_Fjt(basis1()->max_am() +
                               basis2()->max_am() +
                               basis3()->max_am() +
                               basis4()->max_am() +
                               deriv_+1, 1e-15);

    // Allocate the buffer for McMurchie-Davidson-Hermite coefficients
    size_t size = max_am_ * max_am_;
    try {
        mdh_buffer_12_ = new double[size];
        mdh_buffer_34_ = new double[size];
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating mdh_buffer_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(mdh_buffer_12_, 0, sizeof(double) * size);
    memset(mdh_buffer_34_, 0, sizeof(double) * size);
}

ERI_2_2::~ERI_2_2()
{
    delete fjt_;
    delete[] mdh_buffer_;
}


size_t ERI_2_2::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
    // Shells
    const psi::GaussianShell &s1 = bs1_->shell(sh1);
    const psi::GaussianShell &s2 = bs2_->shell(sh2);
    const psi::GaussianShell &s3 = bs3_->shell(sh3);
    const psi::GaussianShell &s4 = bs4_->shell(sh4);

    // Angular momenta
    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();
    int am  = am1 + am2 + am3 + am4; // total am

    // Number of Cartesian functions
    int nam1= INT_CART(am1);
    int nam2= INT_CART(am2);
    int nam3= INT_CART(am3);
    int nam4= INT_CART(am4);

    // Number of primitives
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    int nprim3 = s3.nprimitive();
    int nprim4 = s4.nprimitive();

    // Exponents and contraction coefficients
    const double* a1s = s1.exps();
    const double* a2s = s2.exps();
    const double* a3s = s3.exps();
    const double* a4s = s4.exps();
    const double* c1s = s1.coefs();
    const double* c2s = s2.coefs();
    const double* c3s = s3.coefs();
    const double* c4s = s4.coefs();

    // Coordinates of atomic centres
    double A[3], B[3], C[3], D[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];
    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];
    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];

    double xAB = A[0]-B[0];
    double yAB = A[1]-B[1];
    double zAB = A[2]-B[2];
    double rAB2= xAB*xAB+yAB*yAB+zAB*zAB;
    double xCD = C[0]-D[0];
    double yCD = C[1]-D[1];
    double zCD = C[2]-D[2];
    double rCD2= xCD*xCD+yCD*yCD+zCD*zCD;

    // How many integrals to compute?
    size_t size = nam1 * nam2 * nam3 * nam4;

    // Iterate over primitives
    size_t nprim = 0;

    double PA[3], PB[3], QC[3], QD[3];
    double P[3], Q[3];

    for (int p1 = 0; p1 < nprim1; ++p1) {
         double a1 = a1s[p1];
         double c1 = c1s[p1];
         for (int p2 = 0; p2 < nprim2; ++p2) {
              double a2 = a2s[p2];
              double c2 = c2s[p2];

              double a12 = a1 + a2;
              double ooz = 1.0 / a12;

              P[0] = (a1 * A[0] + a2 * B[0]) * ooz;
              P[1] = (a1 * A[1] + a2 * B[1]) * ooz;
              P[2] = (a1 * A[2] + a2 * B[2]) * ooz;
              PA[0] = P[0] - A[0];
              PA[1] = P[1] - A[1];
              PA[2] = P[2] - A[2];
              PB[0] = P[0] - B[0];
              PB[1] = P[1] - B[1];
              PB[2] = P[2] - B[2];

              double E12 = exp(-(a1*a2*ooz)*rAB2);

              for (int p3 = 0; p3 < nprim3; ++p3) {
                   double a3 = a3s[p3];
                   double c3 = c3s[p3];
                   for (int p4 = 0; p4 < nprim4; ++p4) { 
                        double a4 = a4s[p4];
                        double c4 = c4s[p4];

                        double a34 = a3 + a4;                 
                        double oox = 1.0 / a34;
                                                              
                        Q[0] = (a3 * C[0] + a4 * D[0]) * oox;
                        Q[1] = (a3 * C[1] + a4 * D[1]) * oox;
                        Q[2] = (a3 * C[2] + a4 * D[2]) * oox;
                        QC[0] = Q[0] - C[0];
                        QC[1] = Q[1] - C[1];
                        QC[2] = Q[2] - C[2];
                        QD[0] = Q[0] - D[0];
                        QD[1] = Q[1] - D[1];
                        QD[2] = Q[2] - D[2];

                        double xPQ = P[0]-Q[0];
                        double yPQ = P[1]-Q[1];
                        double zPQ = P[2]-Q[2];
                        double rPQ2= xPQ*xPQ+yPQ*yPQ+zPQ*zPQ;

                        double E34 = exp(-(a3*a4*oox)*rCD2);

                        double apq = a12 + a34;
                        double ooy = 1.0/apq;
                        double pref = c1*c2*c3*c4*E12*E34;
                        double lambda = LOCAL_2PI52 * (1.0/(a12*a34)) * sqrt(ooy);
                        double alpha = a12*a34*ooy;
                        double T = alpha*rPQ2;

                        // Compute McMurchie-Davidson-Hermite coefficients        
                        make_mdh_coeff(am1+am2, am1, am2, PA, PB, mdh_buffer_12_);
                        make_mdh_coeff(am3+am4, am3, am4, PC, PD, mdh_buffer_34_);

                        // Compute McMurchie-Davidson R-coefficients
                        make_mdh_R(...)

                        // Compute the intermediate ERI's
                        for (int ni = 0; ni < nam1; ++ni) {
                             int nx1 = get_cart_am(am1, ni, 0);
                             int ny1 = get_cart_am(am1, ni, 1);
                             int nx1 = get_cart_am(am1, ni, 2);
                             for (int nj = 0; nj < nam2; ++nj) {
                                  int nx2 = get_cart_am(am2, nj, 0); 
                                  int ny2 = get_cart_am(am2, nj, 1);
                                  int nx2 = get_cart_am(am2, nj, 2);
                                  for (int nk = 0; nk < nam3; ++nk) {
                                       int nx3 = get_cart_am(am3, nk, 0);  
                                       int ny3 = get_cart_am(am3, nk, 1);
                                       int nx3 = get_cart_am(am3, nk, 2);
                                       for (int nl = 0; nl < nam4; ++nl) {
                                            int nx4 = get_cart_am(am4, n4, 0);   
                                            int ny4 = get_cart_am(am4, n4, 1);
                                            int nx4 = get_cart_am(am4, n4, 2);

                                            double integral = 0.0;

                                            // Iterate over Hermite functions
                                            for (int N1 = 0; N1 < (nx1+nx2+1); ++N1) {
                                            for (int N2 = 0; N2 < (nx3+nx4+1); ++N2) {
                                                 for (int L1 = 0; L1 < (ny1+ny2+1); ++L1) {
                                                 for (int L2 = 0; L2 < (ny3+ny4+1); ++L2) {
                                                      for (int M1 = 0; M1 < (nz1+nz1+1); ++M1) {
                                                      for (int M2 = 0; M2 < (nz3+nz4+1); ++M2) {
                                                           double i1 = ((N2+L2+M2)%2) ? -1.0 : 1.0;
                                                           integral += get_D(nx1,nx2,N1) * get_D(nx3,nx4,N2)
                                                                     * get_D(ny1,ny2,L1) * get_D(ny3,ny4,N2)
                                                                     * get_D(nz1,nz2,M1) * get_D(nz3,nz4,M2)
                                                                     * get_R(N1+N2,L1+L2,M1+M2)
                                                                     * lambda * i1;
                                                      }
                                                      }
                                                 }
                                                 }
                                            }
                                            }
                                            put_int(integral*pref, target_full_);
                                       }
                                  }
                             }
                        }

                        //
                        ++nprim;
                   }
              }
         }
    }

    // Finish
    return size;
}

} // EndNameSpace oepdev
