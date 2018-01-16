#ifndef _oepdev_libutil_potential_h
#define _oepdev_libutil_potential_h

#include <iostream>
#include <memory>
#include <cstdio>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <map>

#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/integral.h"
#include "psi4/libcubeprop/csg.h"
#include "psi4/liboptions/liboptions.h"
//#include "../liboep/oep.h"
#include "../libpsi/potential.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace psi{
 using namespace std;
 using SharedBasisSet = std::shared_ptr<BasisSet>;
 using SharedMolecule = std::shared_ptr<Molecule>;
 using SharedMatrix   = std::shared_ptr<Matrix  >;
 using SharedWavefunction = std::shared_ptr<Wavefunction>;
}

namespace oepdev{
 
using namespace std;

// Points3DIterator
class Points3DIterator
{
  public:
    Points3DIterator(const int& np);
    virtual ~Points3DIterator();

    // Cube
    static shared_ptr<Points3DIterator> build(const int& nx, const int& ny, const int& nz,
                                              const double& dx, const double& dy, const double& dz,
                                              const double& ox, const double& oy, const double& oz);
    // Random
    static shared_ptr<Points3DIterator> build(const int& np, const double& cx, const double& cy, const double& cz,
                                              const double& radius);
    // Random
    static shared_ptr<Points3DIterator> build(const int& np, const double& pad, psi::SharedMolecule mol);

    virtual double x() const {return current_.x;}
    virtual double y() const {return current_.y;}
    virtual double z() const {return current_.z;}
    virtual int index() const {return current_.index;}

    virtual bool is_done() {return done_;}


    virtual void first() = 0;
    virtual void next() = 0;


  protected:
    struct Point {
       double x, y, z;
       int index;
    };
    Point current_;
    const int np_;
    bool done_;
    int index_;
};


class CubePoints3DIterator : public Points3DIterator
{
  public:
    CubePoints3DIterator(const int& nx, const int& ny, const int& nz,
                         const double& dx, const double& dy, const double& dz,
                         const double& ox, const double& oy, const double& oz);
    virtual ~CubePoints3DIterator();

    virtual void first();
    virtual void next();

  protected:
    const int nx_, ny_, nz_;
    const double dx_, dy_, dz_;
    const double ox_, oy_, oz_;
    int ii_, jj_, kk_;

  private:
    const int ix_max__, iy_max__, iz_max__;
};

class RandomPoints3DIterator : public Points3DIterator
{
  public:
    RandomPoints3DIterator(const int& np, const double& cx, const double& cy, const double& cz,
                                          const double& radius);
    RandomPoints3DIterator(const int& np, const double& pad, psi::SharedMolecule mol);

    virtual ~RandomPoints3DIterator();

    virtual void first();
    virtual void next();                        

  protected:
    double cx_, cy_, cz_;
    double radius_;
    double r_, phi_, theta_, x_, y_, z_;

    psi::SharedMatrix excludeSpheres_;
    std::map<std::string, double> vdwRadius_;

    std::default_random_engine randomNumberGenerator_;
    std::uniform_real_distribution<double> randomDistribution_;

    virtual double random_double() {return randomDistribution_(randomNumberGenerator_);}
    virtual void draw_random_point();
    virtual bool is_in_vdWsphere(double x, double y, double z) const;
};


// Distribution 3D
class Distribution3D
{
  public:
    enum Distribution {Random, Cube};

    Distribution3D(Distribution distribution, int& np);
    Distribution3D(Distribution distribution, const int& np);

    virtual ~Distribution3D();

    static shared_ptr<Distribution3D> build(const int& npoints, const double& cx, const double& cy, const double& cz,
                                            const double& radius);
    static shared_ptr<Distribution3D> build(const int& npoints, const double& padding, psi::SharedMolecule mol);
    static shared_ptr<Distribution3D> build(const int& nx, const int& ny, const int& nz,
                                            const double& px, const double& py, const double& pz,
                                            psi::SharedBasisSet bs, psi::Options& options);

    virtual int npoints() const {return np_;}
    virtual shared_ptr<Points3DIterator> points_iterator() const {return pointsIterator_;}
    virtual Distribution get_type() const {return type_;}

    virtual void print() const = 0;

  protected:
    const int np_;
    Distribution type_;
    shared_ptr<Points3DIterator> pointsIterator_;
};


class RandomDistribution3D : public Distribution3D
{
  public:
    RandomDistribution3D(Distribution distribution, const int& npoints,
                         const double& cx, const double& cy, const double& cz, const double& radius);
    RandomDistribution3D(Distribution distribution, const int& npoints, const double& padding, psi::SharedMolecule mol);
    virtual ~RandomDistribution3D();

    virtual void print() const;

  protected:
};


// Cube Distribution
class CubeDistribution3D : public Distribution3D, public psi::CubicScalarGrid
{
  public:
    CubeDistribution3D(Distribution distribution, const int& nx, const int& ny, const int& nz,
                                                  const double& px, const double& py, const double& pz,
                                                  psi::SharedBasisSet bs, psi::Options& options);

    virtual ~CubeDistribution3D();
    virtual void print() const;

    virtual void write_cube_file(psi::SharedMatrix v, const std::string& name);
  protected:
};


/** \brief 3-dimensional scalar potential manifold: Abstract base.
 *
 */
class Potential3D 
{
  public:
    Potential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius);
    Potential3D(const int& np, const double& pad, psi::SharedMolecule mol);
    Potential3D(const int& nx, const int& ny, const int& nz,
                const double& px, const double& py, const double& pz,
                std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options);

    virtual ~Potential3D();

    static shared_ptr<Potential3D> build(const std::string& type, const int& np,
                                         const double& cx, const double& cy, const double& cz, const double& radius);
    static shared_ptr<Potential3D> build(const std::string& type, const int& np,
                                         const double& pad, psi::SharedMolecule mol);
    static shared_ptr<Potential3D> build(const std::string& type, const int& nx, const int& ny, const int& nz,
                                         const double& px, const double& py, const double& pz,
                                         psi::SharedWavefunction bs, psi::Options& options);

    virtual int npoints() const {return distribution_->npoints();}
    virtual std::shared_ptr<psi::Matrix> data() const {return data_;}
    virtual std::shared_ptr<psi::Wavefunction> wfn() const {return wfn_;}
    virtual std::shared_ptr<Distribution3D> distribution() const {return distribution_;}

    virtual void compute();
    virtual void write_cube_file(const std::string& name);

    virtual double compute_xyz(const double& x, const double& y, const double& z) = 0;
    virtual void print() const = 0;

  protected:
    std::shared_ptr<Distribution3D> distribution_;
    std::shared_ptr<psi::Matrix> data_;
    std::shared_ptr<psi::Wavefunction> wfn_;

    psi::Matrix geom_;
    std::shared_ptr<psi::IntegralFactory> fact_;
    std::shared_ptr<psi::Matrix> pot_;
    std::shared_ptr<psi::OneBodyAOInt> ints_;
    std::shared_ptr<psi::BasisSet> primary_;
    int nbf_;
    
};

class EPotential3D : public Potential3D
{
  public:
    EPotential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius);
    EPotential3D(const int& np, const double& padding, psi::SharedMolecule mol);
    EPotential3D(const int& nx, const int& ny, const int& nz,
                 const double& px, const double& py, const double& pz,
                 psi::SharedWavefunction wfn, psi::Options& options);

    virtual ~EPotential3D();

    virtual double compute_xyz(const double& x, const double& y, const double& z);
    virtual void print() const;
  protected:
};


/** \brief Class template for OEP scalar fields.
 *
 *  @tparam T the compatible class (e.g. OEPotential)
 */
template<class T>
class OEPotential3D : public Potential3D
{
  public:
    // <--- Constructors and Destructor ---> //

    /** \brief Construct random spherical collection of scalar field of type T.
     *
     *  The points are drawn according to uniform distrinution in 3D space.
     *  @param np         - number of points to draw
     *  @param padding    - spherical padding distance (au)
     *  @param oep        - OEP object of type T
     *  @param oepType    - type of OEP
     */
    OEPotential3D(const int& np, const double& padding, std::shared_ptr<T> oep, const std::string& oepType);

    /** \brief Construct ordered 3D collection of scalar field of type T.
     *
     *  The points are generated according to Gaussian cube file format.
     *  @param nx         - number of points along x direction
     *  @param ny         - number of points along y direction
     *  @param nz         - number of points along z direction
     *  @param px         - padding distance along p direction
     *  @param py         - padding distance along p direction
     *  @param pz         - padding distance along p direction
     *  @param oep        - OEP object of type T
     *  @param oepType    - type of OEP
     *  @param options    - Psi4 options object
     */
    OEPotential3D(const int& nx, const int& ny, const int& nz, 
                  const double& px, const double& py, const double& pz,
                  std::shared_ptr<T> oep, const std::string& oepType, psi::Options& options);

    /// Destructor
    virtual ~OEPotential3D();
    
    virtual double compute_xyz(const double& x, const double& y, const double& z);
    virtual void print() const;

  protected:
    std::shared_ptr<T> oep_;
    std::string oepType_;
    
};

template <class T>
OEPotential3D<T>::OEPotential3D(const int& np, const double& padding, std::shared_ptr<T> oep, const std::string& oepType)
 : Potential3D(np, padding, oep->wfn()->molecule()), oep_(oep), oepType_(oepType)
{

}

template <class T>
OEPotential3D<T>::OEPotential3D(const int& nx, const int& ny, const int& nz,
                                const double& px, const double& py, const double& pz,
                                std::shared_ptr<T> oep, const std::string& oepType, psi::Options& options)
 : Potential3D(nx, ny, nz, px, py, pz, oep->wfn(), options), oep_(oep), oepType_(oepType)
{

}

template <class T>
OEPotential3D<T>::~OEPotential3D() 
{

}

template <class T>
void OEPotential3D<T>::print() const 
{

}

template <class T>
double OEPotential3D<T>::compute_xyz(const double& x, const double& y, const double& z)
{
   double val;
   oep_->compute_3D(oepType_, x, y, z, val);
   return val;
}





//class BarePotential3D : public Potential3D
//{
//  public:
//    BarePotential3D(SharedMatrix pot);
//    virtual void compute();
//  protected:
//  private:
//    void common_init();
//};
//
///** \brief Charges from Electrostatic Potential (ESP): abstract base of ESP method.
// *
// */
////class ESP
////{
////  public:
////    //@{ Use atoms as ESP centres
////    ESP(const Wavefunction& wfn);
////    ESP(const OEPotential& oep, const std::string& oepType);
////
////    static std::shared_ptr<ESP> build(const std::string& method, const Wavefunction& wfn);
////    static std::shared_ptr<ESP> build(const std::string& method, const OEPotential& oep, const std::string& oepType);
////    //@}
////
////    virtual ~ESP();
////
////    virtual void fit();
////
////    virtual void compute_potential() = 0;
////    virtual void print() const = 0;
////
////  protected:
////    SharedPotential3D potRef_;
////    SharedPotential3D pot_;
////    
////  private:
////    void common_init(void);
////};
//
////class OEPotentialESP : public ESP
////{
////  public:
////    OEPotentialESP(SharedOEPotential oep, const std::string& oepType);
////  protected:
////  private:
////};

} // EndNameSpace oepdev
#endif //_oepdev_libutil_potential_h
