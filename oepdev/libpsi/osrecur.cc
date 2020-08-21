/*
 * @BEGIN LICENSE
 *
 * Addon to Psi4: an open-source quantum chemistry software package
 *
 * BARTOSZ BŁASIAK (blasiak.bartosz@gmail.com)
 * Improvement of osrecur.h 
 * from original version from Psi4-1.2.1.
 * Modification log:
 *   20.08.2020     - Adding Obara-Saika recursion 
 *                    for improved oepdev::EFPMultipolePotentialInt
 *
 * @END LICENSE
 */

#include <cmath>
#include <stdexcept>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/exception.h"
//#include "psi4/libmints/osrecur.h"
#include "osrecur.h"
#include <iostream>
#include <ctime>


namespace oepdev{

double dfxxx[MAX_DF];


double ***init_box(int a, int b, int c) {
    int i, j;
    double ***box;

    box = (double ***)malloc(sizeof(double **) * a);
    for (i = 0; i < a; i++) box[i] = (double **)malloc(sizeof(double *) * b);
    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            box[i][j] = (double *)malloc(sizeof(double) * c);
            memset((void *)box[i][j], '\0', sizeof(double) * c);
        }
    }
    return box;
}

void zero_box(double ***box, int a, int b, int c) {
    int i, j;
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            memset((void *)box[i][j], 0, sizeof(double) * c);
        }
    }
}

void free_box(double ***box, int a, int b) {
    int i, j;

    for (i = 0; i < a; i++)
        for (j = 0; j < b; j++) free(box[i][j]);

    for (i = 0; i < a; i++) free(box[i]);

    free(box);
}


ObaraSaikaTwoCenterEFPRecursion_New::ObaraSaikaTwoCenterEFPRecursion_New(int max_am1, int max_am2,
                                                                                               int max_k)
    : max_am1_(max_am1), max_am2_(max_am2), do_octupoles_(false) {
    if (max_am1 < 0)
      //throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMVIRecursion -- max_am1 must be nonnegative", __FILE__,
      //                       __LINE__);
       throw psi::PSIEXCEPTION("ObaraSaikaTwoCenterEFPRecursion -- max_am1 must be nonnegative!");
    if (max_am2 < 0)
      //throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMVIRecursion -- max_am2 must be nonnegative", __FILE__,
      //                       __LINE__);
       throw psi::PSIEXCEPTION("ObaraSaikaTwoCenterEFPRecursion -- max_am2 must be nonnegative!");

    if (max_k == 3) do_octupoles_ = true;
// this is just ridiculous! 
//clock_t t_time = -clock();
dfxxx[0] = 1.0;
dfxxx[1] = 1.0;
dfxxx[2] = 1.0;
for (int i = 3; i < MAX_DF; ++i) {
    dfxxx[i] = (i - 1) * dfxxx[i - 2];
}
//t_time += clock(); // Clock END
//std::cout << " o TIME FORMING OF DFXXX : " << ((double)t_time/CLOCKS_PER_SEC) << std::endl;

    size_ = max_am1 > max_am2 ? max_am1 : max_am2;
    size_ += 1;
    size_ = (size_-1)*size_*(size_+1)+1;
    q_   = init_box(size_, size_, max_am1_ + max_am2_ + 4);
    x_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    y_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    z_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    xx_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    yy_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    zz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    xy_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    xz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    yz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    if (do_octupoles_) {
    xxx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    zzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xxy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xxz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    }
}

ObaraSaikaTwoCenterEFPRecursion_New::~ObaraSaikaTwoCenterEFPRecursion_New() {
    free_box(q_  , size_, size_);
    free_box(x_  , size_, size_);
    free_box(y_  , size_, size_);
    free_box(z_  , size_, size_);
    free_box(xx_ , size_, size_);
    free_box(yy_ , size_, size_);
    free_box(zz_ , size_, size_);
    free_box(xy_ , size_, size_);
    free_box(xz_ , size_, size_);
    free_box(yz_ , size_, size_);
    if (do_octupoles_) {
    free_box(xxx_, size_, size_);
    free_box(yyy_, size_, size_);
    free_box(zzz_, size_, size_);
    free_box(xxy_, size_, size_);
    free_box(xxz_, size_, size_);
    free_box(xyy_, size_, size_);
    free_box(yyz_, size_, size_);
    free_box(xzz_, size_, size_);
    free_box(yzz_, size_, size_);
    free_box(xyz_, size_, size_);
    }
}

#define EPS 1.0e-17

void ObaraSaikaTwoCenterEFPRecursion_New::calculate_f(double *F, int n, double t) {
    int i, m;
    int m2;
    double t2;
    double num;
    double sum;
    double term1;
    static double K = 1.0/M_2_SQRTPI;
    double et;


    if (t>20.0){
        t2 = 2*t;
        et = exp(-t);
        t = sqrt(t);
        F[0] = K*erf(t)/t;
        for(m=0; m<=n-1; m++){
            F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
        }
    }
    else {
        et = exp(-t);
        t2 = 2*t;
        m2 = 2*n;
        num = dfxxx[m2];
        i=0;
        sum = 1.0/(m2+1);
        do{
            i++;
            num = num*t2;
            term1 = num/dfxxx[m2+2*i+2];
            sum += term1;
        } while (std::fabs(term1) > EPS && i < MAX_FAC);
        F[n] = sum*et;
        for(m=n-1;m>=0;m--){
            F[m] = (t2*F[m+1] + et)/(2*m+1);
        }
    }

}

void ObaraSaikaTwoCenterEFPRecursion_New::compute(double PA[3], double PB[3], double PC[3], double zeta,
                                                             int am1, int am2) {
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2 + 3;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1]; // TODO: Move this allocation into constructor

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    double z2 = 2.0*zeta;
    double z2x= z2 * PC[0];
    double z2y= z2 * PC[1];
    double z2z= z2 * PC[2];
    double zz2x = 2.0*zeta*z2x;
    double zz2y = 2.0*zeta*z2y;
    double zz2z = 2.0*zeta*z2z;

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        q_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        double q = q_[0][0][m+1];
        x_[0][0][m] = z2x*q;
        y_[0][0][m] = z2y*q;
        z_[0][0][m] = z2z*q;
    }
    for (m=0; m<=mmax-2; ++m) {
        double q1= q_[0][0][m+1];
        double q2= q_[0][0][m+2];
        double z2q1 = z2*q1;
        xx_[0][0][m] = z2x * z2x * q2 - z2q1;
        yy_[0][0][m] = z2y * z2y * q2 - z2q1;
        zz_[0][0][m] = z2z * z2z * q2 - z2q1;
        xy_[0][0][m] = z2x * z2y * q2;
        xz_[0][0][m] = z2x * z2z * q2;
        yz_[0][0][m] = z2y * z2z * q2;
    }
    if (do_octupoles_) {
    for (m=0; m<=mmax-3; ++m) {
        double q2= q_[0][0][m+2];
        double q3= q_[0][0][m+3];

        xxx_[0][0][m] = z2x * z2x * z2x *q3 - 3.0*zz2x*q2;
        yyy_[0][0][m] = z2y * z2y * z2y *q3 - 3.0*zz2y*q2;
        zzz_[0][0][m] = z2z * z2z * z2z *q3 - 3.0*zz2z*q2;
        xxy_[0][0][m] = z2x * z2x * z2y *q3 -     zz2y*q2;  
        xxz_[0][0][m] = z2x * z2x * z2z *q3 -     zz2z*q2;
        xyy_[0][0][m] = z2x * z2y * z2y *q3 -     zz2x*q2;
        yyz_[0][0][m] = z2y * z2y * z2z *q3 -     zz2z*q2;
        xzz_[0][0][m] = z2x * z2z * z2z *q3 -     zz2x*q2;
        yzz_[0][0][m] = z2y * z2z * z2z *q3 -     zz2y*q2;
        xyz_[0][0][m] = z2x * z2y * z2z *q3;
    }
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) { /* Electrostatic potential integrals */
                        q_[0][bind][m] = PB[2] * q_[0][bind-bzm][m] - PC[2] * q_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) { /* Electric field integrals */
                        x_[0][bind][m] = PB[2] * x_[0][bind-bzm][m] - PC[2] * x_[0][bind-bzm][m+1];
                        y_[0][bind][m] = PB[2] * y_[0][bind-bzm][m] - PC[2] * y_[0][bind-bzm][m+1];
                        z_[0][bind][m] = PB[2] * z_[0][bind-bzm][m] - PC[2] * z_[0][bind-bzm][m+1] + q_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-2; ++m) { /* Gradients of the electric field */
                        xx_[0][bind][m] = PB[2]*xx_[0][bind-bzm][m] - PC[2]*xx_[0][bind-bzm][m+1];
                        yy_[0][bind][m] = PB[2]*yy_[0][bind-bzm][m] - PC[2]*yy_[0][bind-bzm][m+1];
                        zz_[0][bind][m] = PB[2]*zz_[0][bind-bzm][m] - PC[2]*zz_[0][bind-bzm][m+1] + 2*z_[0][bind-bzm][m+1];
                        xy_[0][bind][m] = PB[2]*xy_[0][bind-bzm][m] - PC[2]*xy_[0][bind-bzm][m+1];
                        xz_[0][bind][m] = PB[2]*xz_[0][bind-bzm][m] - PC[2]*xz_[0][bind-bzm][m+1] + x_[0][bind-bzm][m+1];
                        yz_[0][bind][m] = PB[2]*yz_[0][bind-bzm][m] - PC[2]*yz_[0][bind-bzm][m+1] + y_[0][bind-bzm][m+1];
                    }
                    if (do_octupoles_) {
                    for (m=0; m <=mmax-b-3; ++m) { /* Hessians of the electric field */
                        xxx_[0][bind][m] = PB[2]*xxx_[0][bind-bzm][m] - PC[2] * xxx_[0][bind-bzm][m+1];
                        yyy_[0][bind][m] = PB[2]*yyy_[0][bind-bzm][m] - PC[2] * yyy_[0][bind-bzm][m+1];
                        zzz_[0][bind][m] = PB[2]*zzz_[0][bind-bzm][m] - PC[2] * zzz_[0][bind-bzm][m+1] + 3 * zz_[0][bind-bzm][m+1];
                        xxy_[0][bind][m] = PB[2]*xxy_[0][bind-bzm][m] - PC[2] * xxy_[0][bind-bzm][m+1];
                        xxz_[0][bind][m] = PB[2]*xxz_[0][bind-bzm][m] - PC[2] * xxz_[0][bind-bzm][m+1] + 1 * xx_[0][bind-bzm][m+1];
                        xyy_[0][bind][m] = PB[2]*xyy_[0][bind-bzm][m] - PC[2] * xyy_[0][bind-bzm][m+1];
                        yyz_[0][bind][m] = PB[2]*yyz_[0][bind-bzm][m] - PC[2] * yyz_[0][bind-bzm][m+1] + 1 * yy_[0][bind-bzm][m+1];
                        xzz_[0][bind][m] = PB[2]*xzz_[0][bind-bzm][m] - PC[2] * xzz_[0][bind-bzm][m+1] + 2 * xz_[0][bind-bzm][m+1];
                        yzz_[0][bind][m] = PB[2]*yzz_[0][bind-bzm][m] - PC[2] * yzz_[0][bind-bzm][m+1] + 2 * yz_[0][bind-bzm][m+1];
                        xyz_[0][bind][m] = PB[2]*xyz_[0][bind-bzm][m] - PC[2] * xyz_[0][bind-bzm][m+1] + 1 * xy_[0][bind-bzm][m+1];
                    }
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bz-1) * (q_[0][bind-2*bzm][m] - q_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bz-1) * (x_[0][bind-2*bzm][m] - x_[0][bind-2*bzm][m+1]);
                            y_[0][bind][m] += ooz * (bz-1) * (y_[0][bind-2*bzm][m] - y_[0][bind-2*bzm][m+1]);
                            z_[0][bind][m] += ooz * (bz-1) * (z_[0][bind-2*bzm][m] - z_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-2; m++) {
                            xx_[0][bind][m] += ooz*(bz-1)*(xx_[0][bind-2*bzm][m] - xx_[0][bind-2*bzm][m+1]);
                            yy_[0][bind][m] += ooz*(bz-1)*(yy_[0][bind-2*bzm][m] - yy_[0][bind-2*bzm][m+1]);
                            zz_[0][bind][m] += ooz*(bz-1)*(zz_[0][bind-2*bzm][m] - zz_[0][bind-2*bzm][m+1]);
                            xy_[0][bind][m] += ooz*(bz-1)*(xy_[0][bind-2*bzm][m] - xy_[0][bind-2*bzm][m+1]);
                            xz_[0][bind][m] += ooz*(bz-1)*(xz_[0][bind-2*bzm][m] - xz_[0][bind-2*bzm][m+1]);
                            yz_[0][bind][m] += ooz*(bz-1)*(yz_[0][bind-2*bzm][m] - yz_[0][bind-2*bzm][m+1]);
                        }
                        if (do_octupoles_) {
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(bz-1)*(xxx_[0][bind-2*bzm][m] - xxx_[0][bind-2*bzm][m+1]);
                            yyy_[0][bind][m] += ooz*(bz-1)*(yyy_[0][bind-2*bzm][m] - yyy_[0][bind-2*bzm][m+1]);
                            zzz_[0][bind][m] += ooz*(bz-1)*(zzz_[0][bind-2*bzm][m] - zzz_[0][bind-2*bzm][m+1]);
                            xxy_[0][bind][m] += ooz*(bz-1)*(xxy_[0][bind-2*bzm][m] - xxy_[0][bind-2*bzm][m+1]);
                            xxz_[0][bind][m] += ooz*(bz-1)*(xxz_[0][bind-2*bzm][m] - xxz_[0][bind-2*bzm][m+1]);
                            xyy_[0][bind][m] += ooz*(bz-1)*(xyy_[0][bind-2*bzm][m] - xyy_[0][bind-2*bzm][m+1]);
                            yyz_[0][bind][m] += ooz*(bz-1)*(yyz_[0][bind-2*bzm][m] - yyz_[0][bind-2*bzm][m+1]);
                            xzz_[0][bind][m] += ooz*(bz-1)*(xzz_[0][bind-2*bzm][m] - xzz_[0][bind-2*bzm][m+1]);
                            yzz_[0][bind][m] += ooz*(bz-1)*(yzz_[0][bind-2*bzm][m] - yzz_[0][bind-2*bzm][m+1]);
                            xyz_[0][bind][m] += ooz*(bz-1)*(xyz_[0][bind-2*bzm][m] - xyz_[0][bind-2*bzm][m+1]);
                        }
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[1] * q_[0][bind-bym][m] - PC[1] * q_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[1] * x_[0][bind-bym][m] - PC[1] * x_[0][bind-bym][m+1];
                        y_[0][bind][m] = PB[1] * y_[0][bind-bym][m] - PC[1] * y_[0][bind-bym][m+1] + q_[0][bind-bym][m+1];
                        z_[0][bind][m] = PB[1] * z_[0][bind-bym][m] - PC[1] * z_[0][bind-bym][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        xx_[0][bind][m] = PB[1]*xx_[0][bind-bym][m] - PC[1]*xx_[0][bind-bym][m+1];
                        yy_[0][bind][m] = PB[1]*yy_[0][bind-bym][m] - PC[1]*yy_[0][bind-bym][m+1] + 2*y_[0][bind-bym][m+1];
                        zz_[0][bind][m] = PB[1]*zz_[0][bind-bym][m] - PC[1]*zz_[0][bind-bym][m+1];
                        xy_[0][bind][m] = PB[1]*xy_[0][bind-bym][m] - PC[1]*xy_[0][bind-bym][m+1] + x_[0][bind-bym][m+1];
                        xz_[0][bind][m] = PB[1]*xz_[0][bind-bym][m] - PC[1]*xz_[0][bind-bym][m+1];
                        yz_[0][bind][m] = PB[1]*yz_[0][bind-bym][m] - PC[1]*yz_[0][bind-bym][m+1] + z_[0][bind-bym][m+1];
                    }
                    if (do_octupoles_) {
                    for (m=0; m <=mmax-b-3; ++m) {
                        xxx_[0][bind][m] = PB[1]*xxx_[0][bind-bym][m] - PC[1] * xxx_[0][bind-bym][m+1];
                        yyy_[0][bind][m] = PB[1]*yyy_[0][bind-bym][m] - PC[1] * yyy_[0][bind-bym][m+1] + 3 * yy_[0][bind-bym][m+1];
                        zzz_[0][bind][m] = PB[1]*zzz_[0][bind-bym][m] - PC[1] * zzz_[0][bind-bym][m+1];
                        xxy_[0][bind][m] = PB[1]*xxy_[0][bind-bym][m] - PC[1] * xxy_[0][bind-bym][m+1] + 1 * xx_[0][bind-bym][m+1];
                        xxz_[0][bind][m] = PB[1]*xxz_[0][bind-bym][m] - PC[1] * xxz_[0][bind-bym][m+1];
                        xyy_[0][bind][m] = PB[1]*xyy_[0][bind-bym][m] - PC[1] * xyy_[0][bind-bym][m+1] + 2 * xy_[0][bind-bym][m+1];
                        yyz_[0][bind][m] = PB[1]*yyz_[0][bind-bym][m] - PC[1] * yyz_[0][bind-bym][m+1] + 2 * yz_[0][bind-bym][m+1];
                        xzz_[0][bind][m] = PB[1]*xzz_[0][bind-bym][m] - PC[1] * xzz_[0][bind-bym][m+1];
                        yzz_[0][bind][m] = PB[1]*yzz_[0][bind-bym][m] - PC[1] * yzz_[0][bind-bym][m+1] + 1 * zz_[0][bind-bym][m+1];
                        xyz_[0][bind][m] = PB[1]*xyz_[0][bind-bym][m] - PC[1] * xyz_[0][bind-bym][m+1] + 1 * xz_[0][bind-bym][m+1];
                    }
                    }

                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (by-1) * (q_[0][bind-2*bym][m] - q_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (by-1) * (x_[0][bind-2*bym][m] - x_[0][bind-2*bym][m+1]);
                            y_[0][bind][m] += ooz * (by-1) * (y_[0][bind-2*bym][m] - y_[0][bind-2*bym][m+1]);
                            z_[0][bind][m] += ooz * (by-1) * (z_[0][bind-2*bym][m] - z_[0][bind-2*bym][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            xx_[0][bind][m] += ooz*(by-1)*(xx_[0][bind-2*bym][m] - xx_[0][bind-2*bym][m+1]);
                            yy_[0][bind][m] += ooz*(by-1)*(yy_[0][bind-2*bym][m] - yy_[0][bind-2*bym][m+1]);
                            zz_[0][bind][m] += ooz*(by-1)*(zz_[0][bind-2*bym][m] - zz_[0][bind-2*bym][m+1]);
                            xy_[0][bind][m] += ooz*(by-1)*(xy_[0][bind-2*bym][m] - xy_[0][bind-2*bym][m+1]);
                            xz_[0][bind][m] += ooz*(by-1)*(xz_[0][bind-2*bym][m] - xz_[0][bind-2*bym][m+1]);
                            yz_[0][bind][m] += ooz*(by-1)*(yz_[0][bind-2*bym][m] - yz_[0][bind-2*bym][m+1]);
                        }
                        if (do_octupoles_) {
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(by-1)*(xxx_[0][bind-2*bym][m] - xxx_[0][bind-2*bym][m+1]);
                            yyy_[0][bind][m] += ooz*(by-1)*(yyy_[0][bind-2*bym][m] - yyy_[0][bind-2*bym][m+1]);
                            zzz_[0][bind][m] += ooz*(by-1)*(zzz_[0][bind-2*bym][m] - zzz_[0][bind-2*bym][m+1]);
                            xxy_[0][bind][m] += ooz*(by-1)*(xxy_[0][bind-2*bym][m] - xxy_[0][bind-2*bym][m+1]);
                            xxz_[0][bind][m] += ooz*(by-1)*(xxz_[0][bind-2*bym][m] - xxz_[0][bind-2*bym][m+1]);
                            xyy_[0][bind][m] += ooz*(by-1)*(xyy_[0][bind-2*bym][m] - xyy_[0][bind-2*bym][m+1]);
                            yyz_[0][bind][m] += ooz*(by-1)*(yyz_[0][bind-2*bym][m] - yyz_[0][bind-2*bym][m+1]);
                            xzz_[0][bind][m] += ooz*(by-1)*(xzz_[0][bind-2*bym][m] - xzz_[0][bind-2*bym][m+1]);
                            yzz_[0][bind][m] += ooz*(by-1)*(yzz_[0][bind-2*bym][m] - yzz_[0][bind-2*bym][m+1]);
                            xyz_[0][bind][m] += ooz*(by-1)*(xyz_[0][bind-2*bym][m] - xyz_[0][bind-2*bym][m+1]);
                        }
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[0] * q_[0][bind-bxm][m] - PC[0] * q_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[0] * x_[0][bind-bxm][m] - PC[0] * x_[0][bind-bxm][m+1] + q_[0][bind-bxm][m+1];
                        y_[0][bind][m] = PB[0] * y_[0][bind-bxm][m] - PC[0] * y_[0][bind-bxm][m+1];
                        z_[0][bind][m] = PB[0] * z_[0][bind-bxm][m] - PC[0] * z_[0][bind-bxm][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        xx_[0][bind][m] = PB[0]*xx_[0][bind-bxm][m] - PC[0]*xx_[0][bind-bxm][m+1] + 2*x_[0][bind-bxm][m+1];
                        yy_[0][bind][m] = PB[0]*yy_[0][bind-bxm][m] - PC[0]*yy_[0][bind-bxm][m+1];
                        zz_[0][bind][m] = PB[0]*zz_[0][bind-bxm][m] - PC[0]*zz_[0][bind-bxm][m+1];
                        xy_[0][bind][m] = PB[0]*xy_[0][bind-bxm][m] - PC[0]*xy_[0][bind-bxm][m+1] + y_[0][bind-bxm][m+1];
                        xz_[0][bind][m] = PB[0]*xz_[0][bind-bxm][m] - PC[0]*xz_[0][bind-bxm][m+1] + z_[0][bind-bxm][m+1];
                        yz_[0][bind][m] = PB[0]*yz_[0][bind-bxm][m] - PC[0]*yz_[0][bind-bxm][m+1];
                    }
                    if (do_octupoles_) {
                    for (m=0; m <=mmax-b-3; ++m) {
                        xxx_[0][bind][m] = PB[0]*xxx_[0][bind-bxm][m] - PC[0] * xxx_[0][bind-bxm][m+1] + 3 * xx_[0][bind-bxm][m+1];
                        yyy_[0][bind][m] = PB[0]*yyy_[0][bind-bxm][m] - PC[0] * yyy_[0][bind-bxm][m+1];
                        zzz_[0][bind][m] = PB[0]*zzz_[0][bind-bxm][m] - PC[0] * zzz_[0][bind-bxm][m+1];
                        xxy_[0][bind][m] = PB[0]*xxy_[0][bind-bxm][m] - PC[0] * xxy_[0][bind-bxm][m+1] + 2 * xy_[0][bind-bxm][m+1];
                        xxz_[0][bind][m] = PB[0]*xxz_[0][bind-bxm][m] - PC[0] * xxz_[0][bind-bxm][m+1] + 2 * xz_[0][bind-bxm][m+1];
                        xyy_[0][bind][m] = PB[0]*xyy_[0][bind-bxm][m] - PC[0] * xyy_[0][bind-bxm][m+1] + 1 * yy_[0][bind-bxm][m+1];
                        yyz_[0][bind][m] = PB[0]*yyz_[0][bind-bxm][m] - PC[0] * yyz_[0][bind-bxm][m+1];
                        xzz_[0][bind][m] = PB[0]*xzz_[0][bind-bxm][m] - PC[0] * xzz_[0][bind-bxm][m+1] + 1 * zz_[0][bind-bxm][m+1];
                        yzz_[0][bind][m] = PB[0]*yzz_[0][bind-bxm][m] - PC[0] * yzz_[0][bind-bxm][m+1];
                        xyz_[0][bind][m] = PB[0]*xyz_[0][bind-bxm][m] - PC[0] * xyz_[0][bind-bxm][m+1] + 1 * yz_[0][bind-bxm][m+1];
                    }
                    }

                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bx-1) * (q_[0][bind-2*bxm][m] - q_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bx-1) * (x_[0][bind-2*bxm][m] - x_[0][bind-2*bxm][m+1]);
                            y_[0][bind][m] += ooz * (bx-1) * (y_[0][bind-2*bxm][m] - y_[0][bind-2*bxm][m+1]);
                            z_[0][bind][m] += ooz * (bx-1) * (z_[0][bind-2*bxm][m] - z_[0][bind-2*bxm][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            xx_[0][bind][m] += ooz*(bx-1)*(xx_[0][bind-2*bxm][m] - xx_[0][bind-2*bxm][m+1]);
                            yy_[0][bind][m] += ooz*(bx-1)*(yy_[0][bind-2*bxm][m] - yy_[0][bind-2*bxm][m+1]);
                            zz_[0][bind][m] += ooz*(bx-1)*(zz_[0][bind-2*bxm][m] - zz_[0][bind-2*bxm][m+1]);
                            xy_[0][bind][m] += ooz*(bx-1)*(xy_[0][bind-2*bxm][m] - xy_[0][bind-2*bxm][m+1]);
                            xz_[0][bind][m] += ooz*(bx-1)*(xz_[0][bind-2*bxm][m] - xz_[0][bind-2*bxm][m+1]);
                            yz_[0][bind][m] += ooz*(bx-1)*(yz_[0][bind-2*bxm][m] - yz_[0][bind-2*bxm][m+1]);
                        }
                        if (do_octupoles_) {
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(bx-1)*(xxx_[0][bind-2*bxm][m] - xxx_[0][bind-2*bxm][m+1]);
                            yyy_[0][bind][m] += ooz*(bx-1)*(yyy_[0][bind-2*bxm][m] - yyy_[0][bind-2*bxm][m+1]);
                            zzz_[0][bind][m] += ooz*(bx-1)*(zzz_[0][bind-2*bxm][m] - zzz_[0][bind-2*bxm][m+1]);
                            xxy_[0][bind][m] += ooz*(bx-1)*(xxy_[0][bind-2*bxm][m] - xxy_[0][bind-2*bxm][m+1]);
                            xxz_[0][bind][m] += ooz*(bx-1)*(xxz_[0][bind-2*bxm][m] - xxz_[0][bind-2*bxm][m+1]);
                            xyy_[0][bind][m] += ooz*(bx-1)*(xyy_[0][bind-2*bxm][m] - xyy_[0][bind-2*bxm][m+1]);
                            yyz_[0][bind][m] += ooz*(bx-1)*(yyz_[0][bind-2*bxm][m] - yyz_[0][bind-2*bxm][m+1]);
                            xzz_[0][bind][m] += ooz*(bx-1)*(xzz_[0][bind-2*bxm][m] - xzz_[0][bind-2*bxm][m+1]);
                            yzz_[0][bind][m] += ooz*(bx-1)*(yzz_[0][bind-2*bxm][m] - yzz_[0][bind-2*bxm][m+1]);
                            xyz_[0][bind][m] += ooz*(bx-1)*(xyz_[0][bind-2*bxm][m] - xyz_[0][bind-2*bxm][m+1]);
                        }
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[2] * q_[aind-azm][bind][m] - PC[2] * q_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[2] * x_[aind-azm][bind][m] - PC[2] * x_[aind-azm][bind][m+1];
                                    y_[aind][bind][m] = PA[2] * y_[aind-azm][bind][m] - PC[2] * y_[aind-azm][bind][m+1];
                                    z_[aind][bind][m] = PA[2] * z_[aind-azm][bind][m] - PC[2] * z_[aind-azm][bind][m+1] + q_[aind-azm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {
                                    xx_[aind][bind][m] = PA[2]*xx_[aind-azm][bind][m] - PC[2]*xx_[aind-azm][bind][m+1];
                                    yy_[aind][bind][m] = PA[2]*yy_[aind-azm][bind][m] - PC[2]*yy_[aind-azm][bind][m+1];
                                    zz_[aind][bind][m] = PA[2]*zz_[aind-azm][bind][m] - PC[2]*zz_[aind-azm][bind][m+1] + 2*z_[aind-azm][bind][m+1];
                                    xy_[aind][bind][m] = PA[2]*xy_[aind-azm][bind][m] - PC[2]*xy_[aind-azm][bind][m+1];
                                    xz_[aind][bind][m] = PA[2]*xz_[aind-azm][bind][m] - PC[2]*xz_[aind-azm][bind][m+1] + x_[aind-azm][bind][m+1];
                                    yz_[aind][bind][m] = PA[2]*yz_[aind-azm][bind][m] - PC[2]*yz_[aind-azm][bind][m+1] + y_[aind-azm][bind][m+1];
                                }
                                if (do_octupoles_) {
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[2]*xxx_[aind-azm][bind][m] - PC[2]*xxx_[aind-azm][bind][m+1];
                                    yyy_[aind][bind][m] = PA[2]*yyy_[aind-azm][bind][m] - PC[2]*yyy_[aind-azm][bind][m+1];
                                    zzz_[aind][bind][m] = PA[2]*zzz_[aind-azm][bind][m] - PC[2]*zzz_[aind-azm][bind][m+1] + 3*zz_[aind-azm][bind][m+1];
                                    xxy_[aind][bind][m] = PA[2]*xxy_[aind-azm][bind][m] - PC[2]*xxy_[aind-azm][bind][m+1];
                                    xxz_[aind][bind][m] = PA[2]*xxz_[aind-azm][bind][m] - PC[2]*xxz_[aind-azm][bind][m+1] + 1*xx_[aind-azm][bind][m+1];
                                    xyy_[aind][bind][m] = PA[2]*xyy_[aind-azm][bind][m] - PC[2]*xyy_[aind-azm][bind][m+1];
                                    yyz_[aind][bind][m] = PA[2]*yyz_[aind-azm][bind][m] - PC[2]*yyz_[aind-azm][bind][m+1] + 1*yy_[aind-azm][bind][m+1];
                                    xzz_[aind][bind][m] = PA[2]*xzz_[aind-azm][bind][m] - PC[2]*xzz_[aind-azm][bind][m+1] + 2*xz_[aind-azm][bind][m+1];
                                    yzz_[aind][bind][m] = PA[2]*yzz_[aind-azm][bind][m] - PC[2]*yzz_[aind-azm][bind][m+1] + 2*yz_[aind-azm][bind][m+1];
                                    xyz_[aind][bind][m] = PA[2]*xyz_[aind-azm][bind][m] - PC[2]*xyz_[aind-azm][bind][m+1] + 1*xy_[aind-azm][bind][m+1];
                                }
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (az-1) * (q_[aind-2*azm][bind][m] - q_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (az-1) * (x_[aind-2*azm][bind][m] - x_[aind-2*azm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (az-1) * (y_[aind-2*azm][bind][m] - y_[aind-2*azm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (az-1) * (z_[aind-2*azm][bind][m] - z_[aind-2*azm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(az-1)*(xx_[aind-2*azm][bind][m] - xx_[aind-2*azm][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(az-1)*(yy_[aind-2*azm][bind][m] - yy_[aind-2*azm][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(az-1)*(zz_[aind-2*azm][bind][m] - zz_[aind-2*azm][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(az-1)*(xy_[aind-2*azm][bind][m] - xy_[aind-2*azm][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(az-1)*(xz_[aind-2*azm][bind][m] - xz_[aind-2*azm][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(az-1)*(yz_[aind-2*azm][bind][m] - yz_[aind-2*azm][bind][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(az-1)*(xxx_[aind-2*azm][bind][m] - xxx_[aind-2*azm][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(az-1)*(yyy_[aind-2*azm][bind][m] - yyy_[aind-2*azm][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(az-1)*(zzz_[aind-2*azm][bind][m] - zzz_[aind-2*azm][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(az-1)*(xxy_[aind-2*azm][bind][m] - xxy_[aind-2*azm][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(az-1)*(xxz_[aind-2*azm][bind][m] - xxz_[aind-2*azm][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(az-1)*(xyy_[aind-2*azm][bind][m] - xyy_[aind-2*azm][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(az-1)*(yyz_[aind-2*azm][bind][m] - yyz_[aind-2*azm][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(az-1)*(xzz_[aind-2*azm][bind][m] - xzz_[aind-2*azm][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(az-1)*(yzz_[aind-2*azm][bind][m] - yzz_[aind-2*azm][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(az-1)*(xyz_[aind-2*azm][bind][m] - xyz_[aind-2*azm][bind][m+1]);
                                    }
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bz * (q_[aind-azm][bind-bzm][m] - q_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bz * (x_[aind-azm][bind-bzm][m] - x_[aind-azm][bind-bzm][m+1]);
                                        y_[aind][bind][m] += ooz * bz * (y_[aind-azm][bind-bzm][m] - y_[aind-azm][bind-bzm][m+1]);
                                        z_[aind][bind][m] += ooz * bz * (z_[aind-azm][bind-bzm][m] - z_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*bz*(xx_[aind-azm][bind-bzm][m] - xx_[aind-azm][bind-bzm][m+1]);
                                        yy_[aind][bind][m] += ooz*bz*(yy_[aind-azm][bind-bzm][m] - yy_[aind-azm][bind-bzm][m+1]);
                                        zz_[aind][bind][m] += ooz*bz*(zz_[aind-azm][bind-bzm][m] - zz_[aind-azm][bind-bzm][m+1]);
                                        xy_[aind][bind][m] += ooz*bz*(xy_[aind-azm][bind-bzm][m] - xy_[aind-azm][bind-bzm][m+1]);
                                        xz_[aind][bind][m] += ooz*bz*(xz_[aind-azm][bind-bzm][m] - xz_[aind-azm][bind-bzm][m+1]);
                                        yz_[aind][bind][m] += ooz*bz*(yz_[aind-azm][bind-bzm][m] - yz_[aind-azm][bind-bzm][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*bz*(xxx_[aind-azm][bind-bzm][m] - xxx_[aind-azm][bind-bzm][m+1]);
                                        yyy_[aind][bind][m] += ooz*bz*(yyy_[aind-azm][bind-bzm][m] - yyy_[aind-azm][bind-bzm][m+1]);
                                        zzz_[aind][bind][m] += ooz*bz*(zzz_[aind-azm][bind-bzm][m] - zzz_[aind-azm][bind-bzm][m+1]);
                                        xxy_[aind][bind][m] += ooz*bz*(xxy_[aind-azm][bind-bzm][m] - xxy_[aind-azm][bind-bzm][m+1]);
                                        xxz_[aind][bind][m] += ooz*bz*(xxz_[aind-azm][bind-bzm][m] - xxz_[aind-azm][bind-bzm][m+1]);
                                        xyy_[aind][bind][m] += ooz*bz*(xyy_[aind-azm][bind-bzm][m] - xyy_[aind-azm][bind-bzm][m+1]);
                                        yyz_[aind][bind][m] += ooz*bz*(yyz_[aind-azm][bind-bzm][m] - yyz_[aind-azm][bind-bzm][m+1]);
                                        xzz_[aind][bind][m] += ooz*bz*(xzz_[aind-azm][bind-bzm][m] - xzz_[aind-azm][bind-bzm][m+1]);
                                        yzz_[aind][bind][m] += ooz*bz*(yzz_[aind-azm][bind-bzm][m] - yzz_[aind-azm][bind-bzm][m+1]);
                                        xyz_[aind][bind][m] += ooz*bz*(xyz_[aind-azm][bind-bzm][m] - xyz_[aind-azm][bind-bzm][m+1]);
                                    }
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[1] * q_[aind-aym][bind][m] - PC[1] * q_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[1] * x_[aind-aym][bind][m] - PC[1] * x_[aind-aym][bind][m+1];
                                    y_[aind][bind][m] = PA[1] * y_[aind-aym][bind][m] - PC[1] * y_[aind-aym][bind][m+1] + q_[aind-aym][bind][m+1];
                                    z_[aind][bind][m] = PA[1] * z_[aind-aym][bind][m] - PC[1] * z_[aind-aym][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {
                                    xx_[aind][bind][m] = PA[1]*xx_[aind-aym][bind][m] - PC[1]*xx_[aind-aym][bind][m+1];
                                    yy_[aind][bind][m] = PA[1]*yy_[aind-aym][bind][m] - PC[1]*yy_[aind-aym][bind][m+1] + 2*y_[aind-aym][bind][m+1];
                                    zz_[aind][bind][m] = PA[1]*zz_[aind-aym][bind][m] - PC[1]*zz_[aind-aym][bind][m+1];
                                    xy_[aind][bind][m] = PA[1]*xy_[aind-aym][bind][m] - PC[1]*xy_[aind-aym][bind][m+1] + x_[aind-aym][bind][m+1];
                                    xz_[aind][bind][m] = PA[1]*xz_[aind-aym][bind][m] - PC[1]*xz_[aind-aym][bind][m+1];
                                    yz_[aind][bind][m] = PA[1]*yz_[aind-aym][bind][m] - PC[1]*yz_[aind-aym][bind][m+1] + z_[aind-aym][bind][m+1];
                                }
                                if (do_octupoles_) {
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[1]*xxx_[aind-aym][bind][m] - PC[1]*xxx_[aind-aym][bind][m+1];
                                    yyy_[aind][bind][m] = PA[1]*yyy_[aind-aym][bind][m] - PC[1]*yyy_[aind-aym][bind][m+1] + 3*yy_[aind-aym][bind][m+1];
                                    zzz_[aind][bind][m] = PA[1]*zzz_[aind-aym][bind][m] - PC[1]*zzz_[aind-aym][bind][m+1];
                                    xxy_[aind][bind][m] = PA[1]*xxy_[aind-aym][bind][m] - PC[1]*xxy_[aind-aym][bind][m+1] + 1*xx_[aind-aym][bind][m+1];
                                    xxz_[aind][bind][m] = PA[1]*xxz_[aind-aym][bind][m] - PC[1]*xxz_[aind-aym][bind][m+1];
                                    xyy_[aind][bind][m] = PA[1]*xyy_[aind-aym][bind][m] - PC[1]*xyy_[aind-aym][bind][m+1] + 2*xy_[aind-aym][bind][m+1];
                                    yyz_[aind][bind][m] = PA[1]*yyz_[aind-aym][bind][m] - PC[1]*yyz_[aind-aym][bind][m+1] + 2*yz_[aind-aym][bind][m+1];
                                    xzz_[aind][bind][m] = PA[1]*xzz_[aind-aym][bind][m] - PC[1]*xzz_[aind-aym][bind][m+1];
                                    yzz_[aind][bind][m] = PA[1]*yzz_[aind-aym][bind][m] - PC[1]*yzz_[aind-aym][bind][m+1] + 1*zz_[aind-aym][bind][m+1];
                                    xyz_[aind][bind][m] = PA[1]*xyz_[aind-aym][bind][m] - PC[1]*xyz_[aind-aym][bind][m+1] + 1*xz_[aind-aym][bind][m+1];
                                }
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ay-1) * (q_[aind-2*aym][bind][m] - q_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ay-1) * (x_[aind-2*aym][bind][m] - x_[aind-2*aym][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ay-1) * (y_[aind-2*aym][bind][m] - y_[aind-2*aym][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ay-1) * (z_[aind-2*aym][bind][m] - z_[aind-2*aym][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(ay-1)*(xx_[aind-2*aym][bind][m] - xx_[aind-2*aym][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(ay-1)*(yy_[aind-2*aym][bind][m] - yy_[aind-2*aym][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(ay-1)*(zz_[aind-2*aym][bind][m] - zz_[aind-2*aym][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(ay-1)*(xy_[aind-2*aym][bind][m] - xy_[aind-2*aym][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(ay-1)*(xz_[aind-2*aym][bind][m] - xz_[aind-2*aym][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(ay-1)*(yz_[aind-2*aym][bind][m] - yz_[aind-2*aym][bind][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(ay-1)*(xxx_[aind-2*aym][bind][m] - xxx_[aind-2*aym][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(ay-1)*(yyy_[aind-2*aym][bind][m] - yyy_[aind-2*aym][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(ay-1)*(zzz_[aind-2*aym][bind][m] - zzz_[aind-2*aym][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(ay-1)*(xxy_[aind-2*aym][bind][m] - xxy_[aind-2*aym][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(ay-1)*(xxz_[aind-2*aym][bind][m] - xxz_[aind-2*aym][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(ay-1)*(xyy_[aind-2*aym][bind][m] - xyy_[aind-2*aym][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(ay-1)*(yyz_[aind-2*aym][bind][m] - yyz_[aind-2*aym][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(ay-1)*(xzz_[aind-2*aym][bind][m] - xzz_[aind-2*aym][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(ay-1)*(yzz_[aind-2*aym][bind][m] - yzz_[aind-2*aym][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(ay-1)*(xyz_[aind-2*aym][bind][m] - xyz_[aind-2*aym][bind][m+1]);
                                    }
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * by * (q_[aind-aym][bind-bym][m] - q_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * by * (x_[aind-aym][bind-bym][m] - x_[aind-aym][bind-bym][m+1]);
                                        y_[aind][bind][m] += ooz * by * (y_[aind-aym][bind-bym][m] - y_[aind-aym][bind-bym][m+1]);
                                        z_[aind][bind][m] += ooz * by * (z_[aind-aym][bind-bym][m] - z_[aind-aym][bind-bym][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*by*(xx_[aind-aym][bind-bym][m] - xx_[aind-aym][bind-bym][m+1]);
                                        yy_[aind][bind][m] += ooz*by*(yy_[aind-aym][bind-bym][m] - yy_[aind-aym][bind-bym][m+1]);
                                        zz_[aind][bind][m] += ooz*by*(zz_[aind-aym][bind-bym][m] - zz_[aind-aym][bind-bym][m+1]);
                                        xy_[aind][bind][m] += ooz*by*(xy_[aind-aym][bind-bym][m] - xy_[aind-aym][bind-bym][m+1]);
                                        xz_[aind][bind][m] += ooz*by*(xz_[aind-aym][bind-bym][m] - xz_[aind-aym][bind-bym][m+1]);
                                        yz_[aind][bind][m] += ooz*by*(yz_[aind-aym][bind-bym][m] - yz_[aind-aym][bind-bym][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*by*(xxx_[aind-aym][bind-bym][m] - xxx_[aind-aym][bind-bym][m+1]);
                                        yyy_[aind][bind][m] += ooz*by*(yyy_[aind-aym][bind-bym][m] - yyy_[aind-aym][bind-bym][m+1]);
                                        zzz_[aind][bind][m] += ooz*by*(zzz_[aind-aym][bind-bym][m] - zzz_[aind-aym][bind-bym][m+1]);
                                        xxy_[aind][bind][m] += ooz*by*(xxy_[aind-aym][bind-bym][m] - xxy_[aind-aym][bind-bym][m+1]);
                                        xxz_[aind][bind][m] += ooz*by*(xxz_[aind-aym][bind-bym][m] - xxz_[aind-aym][bind-bym][m+1]);
                                        xyy_[aind][bind][m] += ooz*by*(xyy_[aind-aym][bind-bym][m] - xyy_[aind-aym][bind-bym][m+1]);
                                        yyz_[aind][bind][m] += ooz*by*(yyz_[aind-aym][bind-bym][m] - yyz_[aind-aym][bind-bym][m+1]);
                                        xzz_[aind][bind][m] += ooz*by*(xzz_[aind-aym][bind-bym][m] - xzz_[aind-aym][bind-bym][m+1]);
                                        yzz_[aind][bind][m] += ooz*by*(yzz_[aind-aym][bind-bym][m] - yzz_[aind-aym][bind-bym][m+1]);
                                        xyz_[aind][bind][m] += ooz*by*(xyz_[aind-aym][bind-bym][m] - xyz_[aind-aym][bind-bym][m+1]);
                                    }
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[0] * q_[aind-axm][bind][m] - PC[0] * q_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[0] * x_[aind-axm][bind][m] - PC[0] * x_[aind-axm][bind][m+1] + q_[aind-axm][bind][m+1];
                                    y_[aind][bind][m] = PA[0] * y_[aind-axm][bind][m] - PC[0] * y_[aind-axm][bind][m+1];
                                    z_[aind][bind][m] = PA[0] * z_[aind-axm][bind][m] - PC[0] * z_[aind-axm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
                                    xx_[aind][bind][m] = PA[0]*xx_[aind-axm][bind][m] - PC[0]*xx_[aind-axm][bind][m+1] + 2*x_[aind-axm][bind][m+1];
                                    yy_[aind][bind][m] = PA[0]*yy_[aind-axm][bind][m] - PC[0]*yy_[aind-axm][bind][m+1];
                                    zz_[aind][bind][m] = PA[0]*zz_[aind-axm][bind][m] - PC[0]*zz_[aind-axm][bind][m+1];
                                    xy_[aind][bind][m] = PA[0]*xy_[aind-axm][bind][m] - PC[0]*xy_[aind-axm][bind][m+1] + y_[aind-axm][bind][m+1];
                                    xz_[aind][bind][m] = PA[0]*xz_[aind-axm][bind][m] - PC[0]*xz_[aind-axm][bind][m+1] + z_[aind-axm][bind][m+1];
                                    yz_[aind][bind][m] = PA[0]*yz_[aind-axm][bind][m] - PC[0]*yz_[aind-axm][bind][m+1];
                                }
                                if (do_octupoles_) {
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[0]*xxx_[aind-axm][bind][m] - PC[0]*xxx_[aind-axm][bind][m+1] + 3*xx_[aind-axm][bind][m+1];
                                    yyy_[aind][bind][m] = PA[0]*yyy_[aind-axm][bind][m] - PC[0]*yyy_[aind-axm][bind][m+1];
                                    zzz_[aind][bind][m] = PA[0]*zzz_[aind-axm][bind][m] - PC[0]*zzz_[aind-axm][bind][m+1];
                                    xxy_[aind][bind][m] = PA[0]*xxy_[aind-axm][bind][m] - PC[0]*xxy_[aind-axm][bind][m+1] + 2*xy_[aind-axm][bind][m+1];
                                    xxz_[aind][bind][m] = PA[0]*xxz_[aind-axm][bind][m] - PC[0]*xxz_[aind-axm][bind][m+1] + 2*xz_[aind-axm][bind][m+1];
                                    xyy_[aind][bind][m] = PA[0]*xyy_[aind-axm][bind][m] - PC[0]*xyy_[aind-axm][bind][m+1] + 1*yy_[aind-axm][bind][m+1];
                                    yyz_[aind][bind][m] = PA[0]*yyz_[aind-axm][bind][m] - PC[0]*yyz_[aind-axm][bind][m+1];
                                    xzz_[aind][bind][m] = PA[0]*xzz_[aind-axm][bind][m] - PC[0]*xzz_[aind-axm][bind][m+1] + 1*zz_[aind-axm][bind][m+1];
                                    yzz_[aind][bind][m] = PA[0]*yzz_[aind-axm][bind][m] - PC[0]*yzz_[aind-axm][bind][m+1];
                                    xyz_[aind][bind][m] = PA[0]*xyz_[aind-axm][bind][m] - PC[0]*xyz_[aind-axm][bind][m+1] + 1*yz_[aind-axm][bind][m+1];
                                }
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ax-1) * (q_[aind-2*axm][bind][m] - q_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ax-1) * (x_[aind-2*axm][bind][m] - x_[aind-2*axm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ax-1) * (y_[aind-2*axm][bind][m] - y_[aind-2*axm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ax-1) * (z_[aind-2*axm][bind][m] - z_[aind-2*axm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(ax-1)*(xx_[aind-2*axm][bind][m] - xx_[aind-2*axm][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(ax-1)*(yy_[aind-2*axm][bind][m] - yy_[aind-2*axm][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(ax-1)*(zz_[aind-2*axm][bind][m] - zz_[aind-2*axm][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(ax-1)*(xy_[aind-2*axm][bind][m] - xy_[aind-2*axm][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(ax-1)*(xz_[aind-2*axm][bind][m] - xz_[aind-2*axm][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(ax-1)*(yz_[aind-2*axm][bind][m] - yz_[aind-2*axm][bind][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(ax-1)*(xxx_[aind-2*axm][bind][m] - xxx_[aind-2*axm][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(ax-1)*(yyy_[aind-2*axm][bind][m] - yyy_[aind-2*axm][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(ax-1)*(zzz_[aind-2*axm][bind][m] - zzz_[aind-2*axm][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(ax-1)*(xxy_[aind-2*axm][bind][m] - xxy_[aind-2*axm][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(ax-1)*(xxz_[aind-2*axm][bind][m] - xxz_[aind-2*axm][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(ax-1)*(xyy_[aind-2*axm][bind][m] - xyy_[aind-2*axm][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(ax-1)*(yyz_[aind-2*axm][bind][m] - yyz_[aind-2*axm][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(ax-1)*(xzz_[aind-2*axm][bind][m] - xzz_[aind-2*axm][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(ax-1)*(yzz_[aind-2*axm][bind][m] - yzz_[aind-2*axm][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(ax-1)*(xyz_[aind-2*axm][bind][m] - xyz_[aind-2*axm][bind][m+1]);
                                    }
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bx * (q_[aind-axm][bind-bxm][m] - q_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bx * (x_[aind-axm][bind-bxm][m] - x_[aind-axm][bind-bxm][m+1]);
                                        y_[aind][bind][m] += ooz * bx * (y_[aind-axm][bind-bxm][m] - y_[aind-axm][bind-bxm][m+1]);
                                        z_[aind][bind][m] += ooz * bx * (z_[aind-axm][bind-bxm][m] - z_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*bx*(xx_[aind-axm][bind-bxm][m] - xx_[aind-axm][bind-bxm][m+1]);
                                        yy_[aind][bind][m] += ooz*bx*(yy_[aind-axm][bind-bxm][m] - yy_[aind-axm][bind-bxm][m+1]);
                                        zz_[aind][bind][m] += ooz*bx*(zz_[aind-axm][bind-bxm][m] - zz_[aind-axm][bind-bxm][m+1]);
                                        xy_[aind][bind][m] += ooz*bx*(xy_[aind-axm][bind-bxm][m] - xy_[aind-axm][bind-bxm][m+1]);
                                        xz_[aind][bind][m] += ooz*bx*(xz_[aind-axm][bind-bxm][m] - xz_[aind-axm][bind-bxm][m+1]);
                                        yz_[aind][bind][m] += ooz*bx*(yz_[aind-axm][bind-bxm][m] - yz_[aind-axm][bind-bxm][m+1]);
                                    }
                                    if (do_octupoles_) {
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*bx*(xxx_[aind-axm][bind-bxm][m] - xxx_[aind-axm][bind-bxm][m+1]);
                                        yyy_[aind][bind][m] += ooz*bx*(yyy_[aind-axm][bind-bxm][m] - yyy_[aind-axm][bind-bxm][m+1]);
                                        zzz_[aind][bind][m] += ooz*bx*(zzz_[aind-axm][bind-bxm][m] - zzz_[aind-axm][bind-bxm][m+1]);
                                        xxy_[aind][bind][m] += ooz*bx*(xxy_[aind-axm][bind-bxm][m] - xxy_[aind-axm][bind-bxm][m+1]);
                                        xxz_[aind][bind][m] += ooz*bx*(xxz_[aind-axm][bind-bxm][m] - xxz_[aind-axm][bind-bxm][m+1]);
                                        xyy_[aind][bind][m] += ooz*bx*(xyy_[aind-axm][bind-bxm][m] - xyy_[aind-axm][bind-bxm][m+1]);
                                        yyz_[aind][bind][m] += ooz*bx*(yyz_[aind-axm][bind-bxm][m] - yyz_[aind-axm][bind-bxm][m+1]);
                                        xzz_[aind][bind][m] += ooz*bx*(xzz_[aind-axm][bind-bxm][m] - xzz_[aind-axm][bind-bxm][m+1]);
                                        yzz_[aind][bind][m] += ooz*bx*(yzz_[aind-axm][bind-bxm][m] - yzz_[aind-axm][bind-bxm][m+1]);
                                        xyz_[aind][bind][m] += ooz*bx*(xyz_[aind-axm][bind-bxm][m] - xyz_[aind-axm][bind-bxm][m+1]);
                                    }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;

}

} // EndNameSpace oepdev
