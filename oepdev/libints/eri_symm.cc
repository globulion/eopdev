#include "psi4/libmints/gshell.h"
#include "recurr.h"
#include "eri_symm.h"

namespace oepdev{
using namespace std;


TwoElectronInt::TwoElectronInt(const psi::IntegralFactory* integral, int deriv, bool use_shell_pairs) :
   psi::TwoBodyAOInt(integral, deriv), 
   use_shell_pairs_(use_shell_pairs),
   cartMap_{0,0,0,
            1,0,0, 0,1,0, 0,0,1,
            2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1,
            3,0,0, 0,3,0, 0,0,3, 1,2,0, 2,1,0, 2,0,1, 1,0,2, 0,1,2, 0,2,1, 1,1,1}
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

}

ERI_2_2::~ERI_2_2()
{
    delete fjt_;
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
    int am = am1 + am2 + am3 + am4; // total am

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

    // How many integrals to compute?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

    // Iterate over primitives
    size_t nprim = 0;

    //
    for (int p1 = 0; p1 < nprim1; ++p1) {
         double a1 = a1s[p1];
         double c1 = c1s[p1];
         for (int p2 = 0; p2 < nprim2; ++p2) {
              double a2 = a2s[p2];
              double c2 = c2s[p2];
              for (int p3 = 0; p3 < nprim3; ++p3) {
                   double a3 = a3s[p3];
                   double c3 = c3s[p3];
                   for (int p4 = 0; p4 < nprim4; ++p4) { 
                        double a4 = a4s[p4];
                        double c4 = c4s[p4];

                        // Compute 
                        ++nprim;
                   }
              }
         }
    }

    // Finish
    return size;
}

} // EndNameSpace oepdev
