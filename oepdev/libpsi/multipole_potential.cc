/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Improvement of psi::EFPMultipolePotentialInt 
 * from original version from Psi4-1.2.1
 * Modification log:
 *   20.08.2020     - Adding two new constructors to EFPMultipolePotentialInt
 *                    enabling calculations of potential integrals
 *                    for a custom set of probe charges.
 *
 * @END LICENSE
 */

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/physconst.h"
#include "multipole_potential.h"
#include <iostream>

namespace oepdev{

#define have_moment(k, max_k) (k <= max_k);


EFPMultipolePotentialInt::EFPMultipolePotentialInt(std::vector<psi::SphericalTransform> &spherical_transforms,
                                             std::shared_ptr<psi::BasisSet> bs1, std::shared_ptr<psi::BasisSet> bs2, int max_k,
                                             int nderiv)
 // : psi::EFPMultipolePotentialInt(spherical_transforms, bs1, bs2, deriv), max_k_(max_k) {
    : psi::OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv),
      mvi_recur_(bs1->max_am(), bs2->max_am(), max_k),
      max_k_(max_k), do_octupoles_(false), nchunk_(10)
{

    if (max_k == 3) { do_octupoles_ = true; nchunk_ = 20; }

    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (nderiv == 0) {
        buffer_ = new double[nchunk_*maxnao1*maxnao2];
        set_chunks(nchunk_);
    }
    else
      //throw psi::FeatureNotImplemented("LibMints", "MultipolePotentialInts called with deriv > 0",  __FILE__, __LINE__);
        throw psi::PSIEXCEPTION("MultipolePotentialInts called with deriv > 0 not implemented!");

}

EFPMultipolePotentialInt::~EFPMultipolePotentialInt() { delete[] buffer_; }

void EFPMultipolePotentialInt::compute_pair(const psi::GaussianShell &s1, const psi::GaussianShell &s2) {

    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // Not sure if these are needed.
    int size =  INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, nchunk_ * size * sizeof(double));

    double ***q   = mvi_recur_.q();
    double ***x   = mvi_recur_.x();
    double ***y   = mvi_recur_.y();
    double ***z   = mvi_recur_.z();
    double ***xx  = mvi_recur_.xx();
    double ***yy  = mvi_recur_.yy();
    double ***zz  = mvi_recur_.zz();
    double ***xy  = mvi_recur_.xy();
    double ***xz  = mvi_recur_.xz();
    double ***yz  = mvi_recur_.yz();

    double ***xxx = nullptr;
    double ***yyy = nullptr;
    double ***zzz = nullptr;
    double ***xxy = nullptr;
    double ***xxz = nullptr;
    double ***xyy = nullptr;
    double ***yyz = nullptr;
    double ***xzz = nullptr;
    double ***yzz = nullptr;
    double ***xyz = nullptr;

    if (do_octupoles_) {
    xxx = mvi_recur_.xxx();
    yyy = mvi_recur_.yyy();
    zzz = mvi_recur_.zzz();
    xxy = mvi_recur_.xxy();
    xxz = mvi_recur_.xxz();
    xyy = mvi_recur_.xyy();
    yyz = mvi_recur_.yyz();
    xzz = mvi_recur_.xzz();
    yzz = mvi_recur_.yzz();
    xyz = mvi_recur_.xyz();
    }

    double Cx = origin_[0];
    double Cy = origin_[1];
    double Cz = origin_[2];

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0 / gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;
            double PC[3];

            PC[0] = P[0] - Cx;
            PC[1] = P[1] - Cy;
            PC[2] = P[2] - Cz;

            // Get recursive
            mvi_recur_.compute(PA, PB, PC, gamma, am1, am2);

            // Gather contributions.
            ao12 = 0;
            for (int ii = 0; ii <= am1; ++ii) {
                int l1 = am1 - ii;
                for (int jj = 0; jj <= ii; ++jj) {
                    int m1 = ii - jj;
                    int n1 = jj;

                    for (int kk = 0; kk <= am2; ++kk) {
                        int l2 = am2 - kk;
                        for (int ll = 0; ll <= kk; ++ll) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // Compute location in the recursion
                            int iind = l1 * ixm + m1 * iym + n1 * izm;
                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

                            buffer_[ao12 + 0 * size]  += q[iind][jind][0] * over_pf;
                            buffer_[ao12 + 1 * size]  += x[iind][jind][0] * over_pf;
                            buffer_[ao12 + 2 * size]  += y[iind][jind][0] * over_pf;
                            buffer_[ao12 + 3 * size]  += z[iind][jind][0] * over_pf;
                            buffer_[ao12 + 4 * size]  += xx[iind][jind][0] * over_pf;
                            buffer_[ao12 + 5 * size]  += yy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 6 * size]  += zz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 7 * size]  += xy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 8 * size]  += xz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 9 * size]  += yz[iind][jind][0] * over_pf;
                            if (do_octupoles_) {
                            buffer_[ao12 + 10 * size] += xxx[iind][jind][0] * over_pf;
                            buffer_[ao12 + 11 * size] += yyy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 12 * size] += zzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 13 * size] += xxy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 14 * size] += xxz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 15 * size] += xyy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 16 * size] += yyz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 17 * size] += xzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 18 * size] += yzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 19 * size] += xyz[iind][jind][0] * over_pf;
                            }
                            ao12++;
                        }
                    }
                }
            }
        }
    }

}

} // EndNameSpace oepdev
