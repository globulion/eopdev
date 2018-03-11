#include "psi4/libmints/gshell.h"
#include "../../include/oepdev_files.h"
#include "eri.h"

namespace oepdev{
using namespace std;

#define LOCAL_2PI52 34.986836655249724969962699105963

TwoElectronInt::TwoElectronInt(const IntegralFactory* integral, int deriv, bool use_shell_pairs) :
   TwoBodyAOInt(integral, deriv), 
   use_shell_pairs_(use_shell_pairs),
   cartMap_{0,0,0,
            1,0,0, 0,1,0, 0,0,1,
            2,0,0, 1,1,0, 1,0,1, 0,2,0, 0,1,1, 0,0,2,
            3,0,0, 2,1,0, 2,0,1, 1,2,0, 1,1,1, 1,0,2, 0,3,0, 0,2,1, 0,1,2, 0,0,3},
   df_{1.0,1.0,3.0,15.0,105.0,945.0,10395.0,135135.0},
   max_am_(std::max({basis1()->max_am(),
                     basis2()->max_am(),
                     basis3()->max_am(),
                     basis4()->max_am()})),
   n_max_am_(2*max_am_+1)
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

    // Allocate the buffer for McMurchie-Davidson R coefficients
    size_t nt = OEPDEV_N_MAX_AM;
    size = nt*nt*nt*(3*nt);
    try {
        mdh_buffer_R_ = new double[size];
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating mdh_buffer_R_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(mdh_buffer_R_, 0, sizeof(double) * size);

}

TwoElectronInt::~TwoElectronInt()
{
    delete[] target_full_;
    delete[] source_full_;
    delete[] mdh_buffer_R_;
}

int TwoElectronInt::get_cart_am(int am, int n, int x)
{
  int n_prev;
  if      (am==0) n_prev = 0;
  else if (am==1) n_prev = 3;
  else if (am==2) n_prev = 12;
  else if (am==3) n_prev = 30;
  else if (am==4) throw psi::PSIEXCEPTION("oepdev::TwoElectronInt: AM 4 not implemented!");
  else            throw psi::PSIEXCEPTION("oepdev::TwoElectronInt: wrong AM chosen!");
  return cartMap_[n_prev + 3*n + x];
}
size_t TwoElectronInt::compute_shell(int sh1, int sh2)
{
 return compute_doublet(sh1, sh2);
}
size_t TwoElectronInt::compute_shell(int sh1, int sh2, int sh3)
{
 return compute_triplet(sh1, sh2, sh3);
}
size_t TwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
 return compute_quartet(sh1, sh2, sh3, sh4);
}
size_t TwoElectronInt::compute_doublet(int sh1, int sh2)
{
 return 0;
}
size_t TwoElectronInt::compute_triplet(int sh1, int sh2, int sh3)
{
 return 0;
}
size_t TwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
 return 0;
}
size_t TwoElectronInt::compute_shell(const psi::AOShellCombinationsIterator &shellIter) {
   return compute_quartet(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}
size_t TwoElectronInt::compute_shell_deriv1(int sh1, int sh2)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: First derivatives are not implemented yet!");
}
size_t TwoElectronInt::compute_shell_deriv1(int sh1, int sh2, int sh3)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: First derivatives are not implemented yet!");
}
size_t TwoElectronInt::compute_shell_deriv1(int sh1, int sh2, int sh3, int sh4)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: First derivatives are not implemented yet!");
}
size_t TwoElectronInt::compute_shell_deriv2(int sh1, int sh2)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: Second derivatives are not implemented yet!");
}
size_t TwoElectronInt::compute_shell_deriv2(int sh1, int sh2, int sh3)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: Second derivatives are not implemented yet!");
}
size_t TwoElectronInt::compute_shell_deriv2(int sh1, int sh2, int sh3, int sh4)
{
   throw psi::PSIEXCEPTION(" oepdev::TwoElectronInt: Second derivatives are not implemented yet!");
}
//--//////////////////////////////////////////////////////////////////////////////
ERI_2_2::ERI_2_2(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    if (deriv_>0) throw psi::PSIEXCEPTION("oepdev::ERI_2_2: Derivatives are not implemented yet!");

    // Boys functor. The +1 is needed for derivatives to work.
    fjt_ = new psi::Taylor_Fjt(basis1()->max_am() +
                               basis2()->max_am() +
                               basis3()->max_am() +
                               basis4()->max_am() +
                               deriv_+1, 1e-15);

    // Allocate the buffer for McMurchie-Davidson-Hermite coefficients (raveled 4-dimensional array)
    //     DIM:   3  Lmax+1            Lmax+1           2*Lmax+1        
    size_t size = 3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1)*OEPDEV_N_MAX_AM;
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

    // Initialize the first (constant) elements of the buffers
    mdh_buffer_12_[D2_INDEX(0,0,0,0)] = 1.0;
    mdh_buffer_12_[D2_INDEX(1,0,0,0)] = 1.0;
    mdh_buffer_12_[D2_INDEX(2,0,0,0)] = 1.0;
    mdh_buffer_34_[D2_INDEX(0,0,0,0)] = 1.0;
    mdh_buffer_34_[D2_INDEX(1,0,0,0)] = 1.0;
    mdh_buffer_34_[D2_INDEX(2,0,0,0)] = 1.0;
}
ERI_2_2::~ERI_2_2()
{
    delete fjt_;
    delete[] mdh_buffer_12_;
    delete[] mdh_buffer_34_;
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

    // Offsets for Cartesian angular momenta buffer
    int nn1 = am1*(am1+1)*(am1+2)/2;
    int nn2 = am2*(am2+1)*(am2+2)/2;
    int nn3 = am3*(am3+1)*(am3+2)/2;
    int nn4 = am4*(am4+1)*(am4+2)/2;

    // Number of Cartesian functions
    int nam1= INT_NCART(am1);
    int nam2= INT_NCART(am2);
    int nam3= INT_NCART(am3);
    int nam4= INT_NCART(am4);

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

    // Clean ERI buffer after previous quartet
    #pragma omp parallel for
    for (int i=0; i<size; ++i) target_full_[i] = 0.0;

    // Iterate over primitives
    double PA[3], PB[3], QC[3], QD[3], PQ[3];
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


              double E12 = exp(-a1*a2*ooz*rAB2);
              if (E12>OEPDEV_CRIT_ERI) {

              // Compute McMurchie-Davidson-Hermite coefficients for first doublet
              make_mdh_D2_coeff(am1, am2, 0.5/a12, PA, PB, mdh_buffer_12_);

              for (int p3 = 0; p3 < nprim3; ++p3) {
                   double a3 = a3s[p3];
                   double c3 = c3s[p3] ;
                   for (int p4 = 0; p4 < nprim4; ++p4) { 
                        double a4 = a4s[p4];
                        double c4 = c4s[p4] ;

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

                        double E34 = exp(-a3*a4*oox*rCD2);

                        double xPQ = P[0]-Q[0];
                        double yPQ = P[1]-Q[1];
                        double zPQ = P[2]-Q[2];
                        double rPQ2= xPQ*xPQ+yPQ*yPQ+zPQ*zPQ;

                        double apq = a12 + a34;
                        double ooy = 1.0/apq;
                        double lambda = LOCAL_2PI52 * (1.0/(a12*a34)) * sqrt(ooy);
                        if (E34>OEPDEV_CRIT_ERI) {
                        double pref = c1*c2*c3*c4*E12*E34;
                        double alpha = a12*a34*ooy;
                        double T = alpha*rPQ2;

                        // Compute McMurchie-Davidson-Hermite coefficients for second doublet
                        make_mdh_D2_coeff(am3, am4, 0.5/a34, QC, QD, mdh_buffer_34_);

                        // Compute McMurchie-Davidson R-coefficients
                        fjt_->set_rho(alpha);
                        double* F = fjt_->values(am, T);
                        make_mdh_R_coeff(am, am, am, alpha, xPQ, yPQ, zPQ, F, mdh_buffer_R_);

                        // Compute the intermediate ERI's
                        int iint = 0;
                        for (int ni = 0; ni < nam1; ++ni) {
                          int inx1 = nn1 + 3*ni;
                          int nx1 = cartMap_[inx1  ];
                          int ny1 = cartMap_[inx1+1];
                          int nz1 = am1 - nx1 - ny1;
                          for (int nj = 0; nj < nam2; ++nj) {
                            int inx2 = nn2 + 3*nj; 
                            int nx2 = cartMap_[inx2  ];
                            int ny2 = cartMap_[inx2+1];
                            int nz2 = am2 - nx2 - ny2;
                            for (int nk = 0; nk < nam3; ++nk) {
                              int inx3 = nn3 + 3*nk;
                              int nx3 = cartMap_[inx3  ];
                              int ny3 = cartMap_[inx3+1];
                              int nz3 = am3 - nx3 - ny3;
                              for (int nl = 0; nl < nam4; ++nl) {
                                int inx4 = nn4 + 3*nl;
                                int nx4 = cartMap_[inx4  ];
                                int ny4 = cartMap_[inx4+1];
                                int nz4 = am4 - nx4 - ny4;

                                double integral = 0.0;

                                // Iterate over Hermite functions
                                for (int N1 = 0; N1 < (nx1+nx2+1); ++N1) {
                                  double Dnx1nx2N1 = get_D12(0,nx1,nx2,N1);
                                  for (int N2 = 0; N2 < (nx3+nx4+1); ++N2) {
                                    double Dnx3nx4N2 = get_D34(0,nx3,nx4,N2);
                                    for (int L1 = 0; L1 < (ny1+ny2+1); ++L1) {
                                      double Dny1ny2L1 = get_D12(1,ny1,ny2,L1);
                                      for (int L2 = 0; L2 < (ny3+ny4+1); ++L2) {
                                        double Dny3ny4L2 = get_D34(1,ny3,ny4,L2);
                                        for (int M1 = 0; M1 < (nz1+nz2+1); ++M1) {
                                          for (int M2 = 0; M2 < (nz3+nz4+1); ++M2) {                       
                                            double i1 = ((N2+L2+M2)%2) ? -1.0 : 1.0;
                                            //integral += get_D12(0,nx1,nx2,N1) * get_D34(0,nx3,nx4,N2)
                                            //          * get_D12(1,ny1,ny2,L1) * get_D34(1,ny3,ny4,L2)
                                            //          * get_D12(2,nz1,nz2,M1) * get_D34(2,nz3,nz4,M2)
                                            //          * get_R(N1+N2,L1+L2,M1+M2)
                                            //          * lambda * i1;
                                            integral += Dnx1nx2N1 * Dnx3nx4N2
                                                      * Dny1ny2L1 * Dny3ny4L2
                                                      * get_D12(2,nz1,nz2,M1) * get_D34(2,nz3,nz4,M2)
                                                      * get_R(N1+N2,L1+L2,M1+M2)
                                                      * lambda * i1;
                                            //integral += mdh_buffer_12_[D2_INDEX(0,nx1,nx2,N1)] 
                                            //          * mdh_buffer_34_[D2_INDEX(0,nx3,nx4,N2)]
                                            //          * mdh_buffer_12_[D2_INDEX(1,ny1,ny2,L1)] 
                                            //          * mdh_buffer_34_[D2_INDEX(1,ny3,ny4,L2)]
                                            //          * mdh_buffer_12_[D2_INDEX(2,nz1,nz2,M1)] 
                                            //          * mdh_buffer_34_[D2_INDEX(2,nz3,nz4,M2)]
                                            //          * mdh_buffer_R_[R_INDEX(N1+N2,L1+L2,M1+M2,0)]
                                            //          * lambda * i1;
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                                target_full_[iint]+= integral*pref;
                                ++iint;
                              }
                            }
                          }
                        }}
                        //
                   }
              }} 
         }
    }

    // Finish
    return size;
}
double ERI_2_2::get_D12(int x, int n1, int n2, int N){
 return mdh_buffer_12_[D2_INDEX(x,n1,n2,N)];
}
double ERI_2_2::get_D34(int x, int n1, int n2, int N){
 return mdh_buffer_34_[D2_INDEX(x,n1,n2,N)];
}
//--//////////////////////////////////////////////////////////////////////////////
ERI_3_1::ERI_3_1(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    if (deriv_>0) throw psi::PSIEXCEPTION("oepdev::ERI_2_2: Derivatives are not implemented yet!");

    // Boys functor. The +1 is needed for derivatives to work.
    fjt_ = new psi::Taylor_Fjt(basis1()->max_am() +
                               basis2()->max_am() +
                               basis3()->max_am() +
                               basis4()->max_am() +
                               deriv_+1, 1e-15);

    // Allocate the buffer for McMurchie-Davidson-Hermite D3 coefficients (raveled 5-dimensional array)
    //     DIM:   3    Lmax+1            Lmax+1             Lmax+1           3*Lmax+1        
    size_t size_1 = 3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1)*(3*OEPDEV_MAX_AM+1);
    // Allocate the buffer for McMurchie-Davidson-Hermite D1 coefficients (raveled 3-dimensional array)
    //     DIM:   3    Lmax+1            Lmax+1
    size_t size_2 = 3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1);
    try {
        mdh_buffer_123_ = new double[size_1];
        mdh_buffer_4_   = new double[size_2];
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating mdh_buffer_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(mdh_buffer_123_, 0, sizeof(double) * size_1);
    memset(mdh_buffer_4_  , 0, sizeof(double) * size_2);

    // Initialize the first (constant) elements of the buffers
    mdh_buffer_123_[D3_INDEX(0,0,0,0,0)] = 1.0;
    mdh_buffer_123_[D3_INDEX(1,0,0,0,0)] = 1.0;
    mdh_buffer_123_[D3_INDEX(2,0,0,0,0)] = 1.0;
    mdh_buffer_4_[D1_INDEX(0,0,0)] = 1.0;
    mdh_buffer_4_[D1_INDEX(1,0,0)] = 1.0;
    mdh_buffer_4_[D1_INDEX(2,0,0)] = 1.0;

}
ERI_3_1::~ERI_3_1()
{
    delete fjt_;
    delete[] mdh_buffer_123_;
    delete[] mdh_buffer_4_;
}
size_t ERI_3_1::compute_quartet(int sh1, int sh2, int sh3, int sh4)
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

    // Offsets for Cartesian angular momenta buffer
    int nn1 = am1*(am1+1)*(am1+2)/2;
    int nn2 = am2*(am2+1)*(am2+2)/2;
    int nn3 = am3*(am3+1)*(am3+2)/2;
    int nn4 = am4*(am4+1)*(am4+2)/2;

    // Number of Cartesian functions
    int nam1= INT_NCART(am1);
    int nam2= INT_NCART(am2);
    int nam3= INT_NCART(am3);
    int nam4= INT_NCART(am4);

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

    // How many integrals to compute?
    size_t size = nam1 * nam2 * nam3 * nam4;

    // Clean ERI buffer after previous quartet
    for (int i=0; i<size; ++i) target_full_[i] = 0.0;

    // Iterate over primitives
    double PA[3], PB[3], RA[3], RB[3], RC[3];
    double P[3], R[3];

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

              double E12 = exp(-a1*a2*ooz*rAB2);
              if (E12>OEPDEV_CRIT_ERI) {

              for (int p3 = 0; p3 < nprim3; ++p3) {
                   double a3 = a3s[p3];
                   double c3 = c3s[p3] ;
                   double a123 = a12 + a3;
                   double ooo  = 1.0/a123;

                   R[0] = (a3 * C[0] + a12 * P[0]) * ooo;
                   R[1] = (a3 * C[1] + a12 * P[1]) * ooo;
                   R[2] = (a3 * C[2] + a12 * P[2]) * ooo;

                   RA[0] = R[0] - A[0];
                   RA[1] = R[1] - A[1];
                   RA[2] = R[2] - A[2];
                   RB[0] = R[0] - B[0];
                   RB[1] = R[1] - B[1];
                   RB[2] = R[2] - B[2];
                   RC[0] = R[0] - C[0];
                   RC[1] = R[1] - C[1];
                   RC[2] = R[2] - C[2];

                   double xPC = P[0] - C[0];
                   double yPC = P[1] - C[1];
                   double zPC = P[2] - C[2];
                   double rPC2 = xPC*xPC + yPC*yPC + zPC*zPC;
                   double E123 = E12*exp(-a12*a3*ooo*rPC2);

                   if (E123 > OEPDEV_CRIT_ERI) {

                   // Compute McMurchie-Davidson-Hermite coefficients for triplet
                   make_mdh_D3_coeff(am1, am2, am3, 0.5*ooo, RA, RB, RC, mdh_buffer_123_);

                   for (int p4 = 0; p4 < nprim4; ++p4) { 
                        double a4 = a4s[p4];
                        double c4 = c4s[p4] ;

                        double xRD = R[0] - D[0];
                        double yRD = R[1] - D[1];
                        double zRD = R[2] - D[2];

                        double rRD2= xRD*xRD+yRD*yRD+zRD*zRD;
                        double oxo = 1.0/(a123+a4);

                        double lambda = LOCAL_2PI52 * (1.0/(a123*a4)) * sqrt(oxo);

                        double pref = c1*c2*c3*c4*E123;
                        double alpha = a123*a4*oxo;
                        double T = alpha*rRD2;

                        // Compute McMurchie-Davidson-Hermite coefficients for 4-th shell
                        make_mdh_D1_coeff(am4, 0.5/a4, mdh_buffer_4_);

                        // Compute McMurchie-Davidson R-coefficients
                        fjt_->set_rho(alpha);
                        double* F = fjt_->values(am, T);
                        make_mdh_R_coeff(am, am, am, alpha, xRD, yRD, zRD, F, mdh_buffer_R_);

                        // Compute the intermediate ERI's
                        int iint = 0;
                        for (int ni = 0; ni < nam1; ++ni) {
                          int inx1 = nn1 + 3*ni;
                          int nx1 = cartMap_[inx1  ];
                          int ny1 = cartMap_[inx1+1];
                          int nz1 = am1 - nx1 - ny1;
                          for (int nj = 0; nj < nam2; ++nj) {
                            int inx2 = nn2 + 3*nj; 
                            int nx2 = cartMap_[inx2  ];
                            int ny2 = cartMap_[inx2+1];
                            int nz2 = am2 - nx2 - ny2;
                            for (int nk = 0; nk < nam3; ++nk) {
                              int inx3 = nn3 + 3*nk;
                              int nx3 = cartMap_[inx3  ];
                              int ny3 = cartMap_[inx3+1];
                              int nz3 = am3 - nx3 - ny3;
                              for (int nl = 0; nl < nam4; ++nl) {
                                int inx4 = nn4 + 3*nl;
                                int nx4 = cartMap_[inx4  ];
                                int ny4 = cartMap_[inx4+1];
                                int nz4 = am4 - nx4 - ny4;

                                double integral = 0.0;

                                // Iterate over Hermite functions
                                for (int N1 = 0; N1 < (nx1+nx2+nx3+1); ++N1) {
                                  double Dnx1nx2nx3N1 = get_D123(0,nx1,nx2,nx3,N1);
                                  for (int N2 = 0; N2 < (nx4+1); ++N2) {
                                    double Dnx4N2 = get_D4(0,nx4,N2);
                                    for (int L1 = 0; L1 < (ny1+ny2+ny3+1); ++L1) {
                                      double Dny1ny2ny3L1 = get_D123(1,ny1,ny2,ny3,L1);
                                      for (int L2 = 0; L2 < (ny4+1); ++L2) {
                                        double Dny4L2 = get_D4(1,ny4,L2);
                                        for (int M1 = 0; M1 < (nz1+nz2+nz3+1); ++M1) {
                                          for (int M2 = 0; M2 < (nz4+1); ++M2) {                       
                                            double i1 = ((N2+L2+M2)%2) ? -1.0 : 1.0;
                                            integral += Dnx1nx2nx3N1 * Dnx4N2
                                                      * Dny1ny2ny3L1 * Dny4L2
                                                      * get_D123(2,nz1,nz2,nz3,M1) * get_D4(2,nz4,M2)
                                                      * get_R(N1+N2,L1+L2,M1+M2)
                                                      * lambda * i1;
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                                target_full_[iint]+= integral*pref;
                                ++iint;
                              }
                            }
                          }
                        }
                        //
                   }} 
              }}
         }
    }

    // Finish
    return size;
}
double ERI_3_1::get_D123(int x, int n1, int n2, int n3, int N){
 return mdh_buffer_123_[D3_INDEX(x,n1,n2,n3,N)];
}
double ERI_3_1::get_D4(int x, int n1, int N){
 return mdh_buffer_4_[D1_INDEX(x,n1,N)];
}
//--//////////////////////////////////////////////////////////////////////////////
ERI_1_1::ERI_1_1(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    if (deriv_>0) throw psi::PSIEXCEPTION("oepdev::ERI_1_1: Derivatives are not implemented yet!");

    // Boys functor. The +1 is needed for derivatives to work.
    fjt_ = new psi::Taylor_Fjt(basis1()->max_am() +
                               basis2()->max_am() +
                               basis3()->max_am() +
                               basis4()->max_am() +
                               deriv_+1, 1e-15);

    // Allocate the buffer for McMurchie-Davidson-Hermite D1 coefficients (raveled 3-dimensional array)
    //     DIM:   3    Lmax+1            Lmax+1
    size_t size = 3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1);
    try {
        mdh_buffer_1_ = new double[size];
        mdh_buffer_2_ = new double[size];
    }
    catch (std::bad_alloc &e) {
        psi::outfile->Printf("Error allocating mdh_buffer_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(mdh_buffer_1_, 0, sizeof(double) * size);
    memset(mdh_buffer_2_, 0, sizeof(double) * size);

    // Initialize the first (constant) elements of the buffers
    mdh_buffer_1_[D1_INDEX(0,0,0)] = 1.0;
    mdh_buffer_1_[D1_INDEX(1,0,0)] = 1.0;
    mdh_buffer_1_[D1_INDEX(2,0,0)] = 1.0;
    mdh_buffer_2_[D1_INDEX(0,0,0)] = 1.0;
    mdh_buffer_2_[D1_INDEX(1,0,0)] = 1.0;
    mdh_buffer_2_[D1_INDEX(2,0,0)] = 1.0;

}
ERI_1_1::~ERI_1_1()
{
    delete fjt_;
    delete[] mdh_buffer_1_;
    delete[] mdh_buffer_2_;
}
size_t ERI_1_1::compute_doublet(int sh1, int sh2)
{
    // Shells
    const psi::GaussianShell &s1 = bs1_->shell(sh1);
    const psi::GaussianShell &s2 = bs2_->shell(sh2);

    // Angular momenta
    int am1 = s1.am();
    int am2 = s2.am();
    int am  = am1 + am2; // total am

    // Offsets for Cartesian angular momenta buffer
    int nn1 = am1*(am1+1)*(am1+2)/2;
    int nn2 = am2*(am2+1)*(am2+2)/2;

    // Number of Cartesian functions
    int nam1= INT_NCART(am1);
    int nam2= INT_NCART(am2);

    // Number of primitives
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();

    // Exponents and contraction coefficients
    const double* a1s = s1.exps();
    const double* a2s = s2.exps();
    const double* c1s = s1.coefs();
    const double* c2s = s2.coefs();

    // Coordinates of atomic centres
    double A[3], B[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    double xAB = A[0]-B[0];
    double yAB = A[1]-B[1];
    double zAB = A[2]-B[2];
    double rAB2= xAB*xAB+yAB*yAB+zAB*zAB;

    // How many integrals to compute?
    size_t size = nam1 * nam2;

    // Clean ERI buffer after previous quartet
    for (int i=0; i<size; ++i) target_full_[i] = 0.0;

    // Iterate over primitives
    for (int p1 = 0; p1 < nprim1; ++p1) {
         double a1 = a1s[p1];
         double c1 = c1s[p1];

         // Compute McMurchie-Davidson-Hermite D1 coefficients for first shell
         make_mdh_D1_coeff(am1, 0.5/a1, mdh_buffer_1_);

         for (int p2 = 0; p2 < nprim2; ++p2) {
              double a2 = a2s[p2];
              double c2 = c2s[p2];

              double a12 = a1 + a2;
              double oxo = 1.0 / a12;

              double lambda = LOCAL_2PI52 * (1.0/(a1*a2)) * sqrt(oxo);

              double pref = c1*c2;
              double alpha = a1*a2*oxo;
              double T = alpha*rAB2;

              // Compute McMurchie-Davidson-Hermite D1 coefficients for second shell
              make_mdh_D1_coeff(am2, 0.5/a2, mdh_buffer_2_);

              // Compute McMurchie-Davidson R-coefficients
              fjt_->set_rho(alpha);
              double* F = fjt_->values(am, T);
              make_mdh_R_coeff(am, am, am, alpha, xAB, yAB, zAB, F, mdh_buffer_R_);

              // Compute the intermediate ERI's
              int iint = 0;
              for (int ni = 0; ni < nam1; ++ni) {
                int inx1 = nn1 + 3*ni;
                int nx1 = cartMap_[inx1  ];
                int ny1 = cartMap_[inx1+1];
                int nz1 = am1 - nx1 - ny1;
                for (int nj = 0; nj < nam2; ++nj) {
                  int inx2 = nn2 + 3*nj; 
                  int nx2 = cartMap_[inx2  ];
                  int ny2 = cartMap_[inx2+1];
                  int nz2 = am2 - nx2 - ny2;

                  double integral = 0.0;

                  // Iterate over Hermite functions
                  for (int N1 = 0; N1 < (nx1+1); ++N1) {
                    double Dnx1N1 = get_D1(0,nx1,N1);
                    for (int N2 = 0; N2 < (nx2+1); ++N2) {
                      double Dnx2N2 = get_D2(0,nx2,N2);
                      for (int L1 = 0; L1 < (ny1+1); ++L1) {
                        double Dny1L1 = get_D1(1,ny1,L1);
                        for (int L2 = 0; L2 < (ny2+1); ++L2) {
                          double Dny2L2 = get_D2(1,ny2,L2);
                          for (int M1 = 0; M1 < (nz1+1); ++M1) {
                            for (int M2 = 0; M2 < (nz2+1); ++M2) {                       
                              double i1 = ((N2+L2+M2)%2) ? -1.0 : 1.0;
                              integral += Dnx1N1 * Dnx2N2
                                        * Dny1L1 * Dny2L2
                                        * get_D1(2,nz1,M1) * get_D2(2,nz2,M2)
                                        * get_R(N1+N2,L1+L2,M1+M2)
                                        * lambda * i1;
                            }
                          }
                        }
                      }
                    }
                  }
                  target_full_[iint]+= integral*pref;
                  ++iint;
                }
              }
         }
    }

    // Finish
    return size;
}
double ERI_1_1::get_D1(int x, int n1, int N){
 return mdh_buffer_1_[D1_INDEX(x,n1,N)];
}
double ERI_1_1::get_D2(int x, int n2, int N){
 return mdh_buffer_2_[D1_INDEX(x,n2,N)];
}


} // EndNameSpace oepdev
