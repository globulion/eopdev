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

#ifndef _oepdev_libpsi_osrecur_h_
#define _oepdev_libpsi_osrecur_h_
/** @file osrecur.h */

#define MAX_DF 500
#define MAX_FAC 100


namespace oepdev{
//using namespace std;

extern double dfxxx[MAX_DF];


double ***init_box(int a, int b, int c);
void zero_box(double ***box, int a, int b, int c);
void free_box(double ***box, int a, int b);


/** \addtogroup OEPDEV_LIBINTS
 * @{
 */

/** 
 *  \class ObaraSaikaTwoCenterEFPRecursion_New
 *  \brief Obara-Saika recursion formulae for improved EFP multipole potential integrals.
 */
class ObaraSaikaTwoCenterEFPRecursion_New
{
protected:
    int max_am1_;
    int max_am2_;
    int size_;

    bool do_octupoles_;

    double*** q_;
    double*** x_;
    double*** y_;
    double*** z_;
    double*** xx_;
    double*** xy_;
    double*** xz_;
    double*** yy_;
    double*** yz_;
    double*** zz_;
    double*** xxx_;
    double*** xxy_;
    double*** xxz_;
    double*** xyy_;
    double*** xyz_;
    double*** xzz_;
    double*** yyy_;
    double*** yyz_;
    double*** yzz_;
    double*** zzz_;

    // Forms Fm(U) from A20 (OS 1986)
    void calculate_f(double *F, int n, double t);

public:
    // No default constructor
    ObaraSaikaTwoCenterEFPRecursion_New() {};
    // No assignment operator
    ObaraSaikaTwoCenterEFPRecursion_New& operator=(const ObaraSaikaTwoCenterEFPRecursion_New&) {};

public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterEFPRecursion_New(int max_am1, int max_am2, int max_k);
    virtual ~ObaraSaikaTwoCenterEFPRecursion_New();

    /// Returns the potential integral 3D matrix
    double*** q  () const { return q_;   }
    double*** x  () const { return x_;   }
    double*** y  () const { return y_;   }
    double*** z  () const { return z_;   }
    double*** xx () const { return xx_;  }
    double*** yy () const { return yy_;  }
    double*** zz () const { return zz_;  }
    double*** xy () const { return xy_;  }
    double*** xz () const { return xz_;  }
    double*** yz () const { return yz_;  }
    double*** xxx() const { return xxx_; }
    double*** yyy() const { return yyy_; }
    double*** zzz() const { return zzz_; }
    double*** xxy() const { return xxy_; }
    double*** xxz() const { return xxz_; }
    double*** xyy() const { return xyy_; }
    double*** yyz() const { return yyz_; }
    double*** xzz() const { return xzz_; }
    double*** yzz() const { return yzz_; }
    double*** xyz() const { return xyz_; }

    /// Computes the potential integral 3D matrix using the data provided.
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);

};

/** @}*/
} // EndNameSpace oepdev

#endif
