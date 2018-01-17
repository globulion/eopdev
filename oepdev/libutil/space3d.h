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
#include "../libpsi/potential.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace psi{
 using namespace std;
 using SharedBasisSet     = std::shared_ptr<BasisSet    >; 
 using SharedMolecule     = std::shared_ptr<Molecule    >;
 using SharedMatrix       = std::shared_ptr<Matrix      >;
 using SharedWavefunction = std::shared_ptr<Wavefunction>;
}

namespace oepdev{
 
using namespace std;

/** \brief Iterator over a collection of points in 3D space. Abstract base.
 *
 *  Points3DIterators are constructed either as iterators over:
 *   - a random collections or
 *   - an ordered (g09 cube-like) collections.
 *  __Note:__ Always create instances by using static factory methods.
 */
class Points3DIterator
{
  public:

    // <--- Constructor and Destructor ---> //

    /** \brief Plain constructor. Initializes the abstract features.
      *  @param np   - number of points this iterator is constructed for
      */
    Points3DIterator(const int& np);

    /// Destructor
    virtual ~Points3DIterator();


    // <--- Factory Methods ---> //

    /** \brief G09 Cube collection iterator.
      *
      *  The points are generated according to Gaussian cube file format. 
      *  @param nx         - number of points along x direction
      *  @param ny         - number of points along y direction
      *  @param nz         - number of points along z direction
      *  @param dx         - spacing distance along x direction
      *  @param dy         - spacing distance along y direction
      *  @param dz         - spacing distance along y direction
      *  @param ox         - coordinate x of cube origin 
      *  @param oy         - coordinate y of cube origin               
      *  @param oz         - coordinate z of cube origin               
      */
    static shared_ptr<Points3DIterator> build(const int& nx, const int& ny, const int& nz,
                                              const double& dx, const double& dy, const double& dz,
                                              const double& ox, const double& oy, const double& oz);

    /** \brief Random collection iterator.
      *
      *  The points are drawn according to uniform distrinution in 3D space.
      *  @param np         - number of points to draw
      *  @param radius     - sphere radius inside which points are to be drawn
      *  @param cx         - coordinate x of sphere's centre
      *  @param cy         - coordinate y of sphere's centre
      *  @param cz         - coordinate z of sphere's centre
      */
    static shared_ptr<Points3DIterator> build(const int& np, const double& radius, 
                                              const double& cx, const double& cy, const double& cz);

    /** \brief Random collection iterator.
      *
      *  The points are drawn according to uniform distrinution in 3D space
      *  enclosing a molecule given. All drawn points lie outside the van der Waals
      *  volume.
      *  @param np         - number of points to draw
      *  @param pad        - radius padding of a minimal sphere enclosing the molecule
      *  @param mol        - Psi4 molecule object
      */
    static shared_ptr<Points3DIterator> build(const int& np, const double& pad, psi::SharedMolecule mol);


    // <--- Accessors ---> //

    //@{ Retrieve current point coordinates and point index
    virtual double x() const {return current_.x;}
    virtual double y() const {return current_.y;}
    virtual double z() const {return current_.z;}
    virtual int index() const {return current_.index;}
    //@}

    /// Check if iteration is finished
    virtual bool is_done() {return done_;}


    // <--- Mutators ---> //

    /// Initialize first iteration
    virtual void first() = 0;

    /// Step to next iteration
    virtual void next() = 0;


  protected:

    //@{ Current point
    struct Point {
       double x, y, z;
       int index;
    };

    Point current_;
    //@}

    /// Number of points
    const int np_;
    /// Status of the iterator
    bool done_;
    /// Current index
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
    RandomPoints3DIterator(const int& np, const double& radius, const double& cx, const double& cy, const double& cz);
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


/** \brief Collection of points in 3D space. Abstract base.
 *
 *  Create random or ordered (g09 cube-like) collections of points in 3D space.
 *  
 *  __Note:__ Always create instances by using static factory methods.
 */
class PointsCollection3D
{
  public:

    // <--- Enumerators ---> //

    /// Public descriptior of collection type
    enum Collection {Random, Cube};


    // <--- Constructors and Destructor ---> //

    /** \brief Initialize abstract features.
      *  @param np   - number of points to be created
      */
    PointsCollection3D(Collection collectionType, int& np);
    PointsCollection3D(Collection collectionType, const int& np);

    /// Destructor
    virtual ~PointsCollection3D();


    // <--- Factories ---> //

    /** \brief Random collection of points.
      *
      *  Points uniformly span a sphere.
      *  @param npoints    - number of points to draw
      *  @param radius     - sphere radius inside which points are to be drawn
      *  @param cx         - coordinate x of sphere's centre
      *  @param cy         - coordinate y of sphere's centre
      *  @param cz         - coordinate z of sphere's centre
      */
    static shared_ptr<PointsCollection3D> build(const int& npoints, const double& radius,
                                                const double& cx = 0.0, 
                                                const double& cy = 0.0, 
                                                const double& cz = 0.0);

    /** \brief Random collection of points.
      *
      *  Points uniformly span space inside a sphere enclosing a molecule. 
      *  exluding the van der Waals volume.
      *  @param np         - number of points to draw
      *  @param padding    - radius padding of a minimal sphere enclosing the molecule
      *  @param mol        - Psi4 molecule object
      */
    static shared_ptr<PointsCollection3D> build(const int& npoints, const double& padding, psi::SharedMolecule mol);

    /** \brief G09 Cube collection of points.
      *
      *  The points span a parallelpiped according to Gaussian cube file format. 
      *  @param nx         - number of points along x direction
      *  @param ny         - number of points along y direction
      *  @param nz         - number of points along z direction
      *  @param px         - padding distance along x direction
      *  @param py         - padding distance along y direction
      *  @param pz         - padding distance along z direction
      *  @param bs         - Psi4 basis set object
      *  @param options    - Psi4 options object
      */
    static shared_ptr<PointsCollection3D> build(const int& nx, const int& ny, const int& nz,
                                                const double& px, const double& py, const double& pz,
                                                psi::SharedBasisSet bs, psi::Options& options);


    // <--- Accessors ---> //

    /// Get the number of points
    virtual int npoints() const {return np_;}

    /// Get the iterator over this collection of points
    virtual shared_ptr<Points3DIterator> points_iterator() const {return pointsIterator_;}

    /// Get the collection type
    virtual Collection get_type() const {return collectionType_;}

    virtual void print() const = 0;


  protected:

    /// Number of points
    const int np_;

    /// Collection type
    Collection collectionType_;

    /// iterator over points collection
    shared_ptr<Points3DIterator> pointsIterator_;
};


class RandomPointsCollection3D : public PointsCollection3D
{
  public:
    RandomPointsCollection3D(Collection collectionType, const int& npoints, const double& radius,
                         const double& cx, const double& cy, const double& cz);
    RandomPointsCollection3D(Collection collectionType, const int& npoints, const double& padding, psi::SharedMolecule mol);
    virtual ~RandomPointsCollection3D();

    virtual void print() const;

  protected:
};


// Cube Distribution
class CubePointsCollection3D : public PointsCollection3D, public psi::CubicScalarGrid
{
  public:
    CubePointsCollection3D(Collection collectionType, const int& nx, const int& ny, const int& nz,
                                                  const double& px, const double& py, const double& pz,
                                                  psi::SharedBasisSet bs, psi::Options& options);

    virtual ~CubePointsCollection3D();
    virtual void print() const;

    virtual void write_cube_file(psi::SharedMatrix v, const std::string& name);
  protected:
};


/** \brief Scalar field in 3D space. Abstract base.
 *
 *  Create various scalar fields with points distributed randomly or 
 *  as an ordered g09 cube-like collection.
 *  
 *  __Note:__ Always create instances by using static factory methods.
 */

class ScalarField3D
{
  public:
    
    // <--- Constructors and Destructor ---> //

    // TODO!
    ScalarField3D(const int& np, const double& radius, const double& cx, const double& cy, const double& cz);
    ScalarField3D(const int& np, const double& pad, psi::SharedWavefunction wfn, psi::Options& options);
    ScalarField3D(const int& nx, const int& ny, const int& nz,
                  const double& px, const double& py, const double& pz,
                  std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options);

    /// Destructor
    virtual ~ScalarField3D();


    // <--- Factory Methods ---> //

    // TODO!
    static shared_ptr<ScalarField3D> build(const std::string& type, const int& np, const double& radius,
                                           const double& cx, const double& cy, const double& cz);
    static shared_ptr<ScalarField3D> build(const std::string& type, const int& np,
                                           const double& pad, psi::SharedWavefunction wfn, psi::Options& options);
    static shared_ptr<ScalarField3D> build(const std::string& type, const int& nx, const int& ny, const int& nz,
                                           const double& px, const double& py, const double& pz,
                                           psi::SharedWavefunction wfn, psi::Options& options);


    // <--- Accessors ---> //

    /// Get the number of points at which the scalar field is defined
    virtual int npoints() const {return pointsCollection_->npoints();}

    /// Get the collection of points
    virtual std::shared_ptr<PointsCollection3D> points_collection() const {return pointsCollection_;}

    /// Get the data matrix in a form { [x, y, z, f(x, y, z)] }
    virtual std::shared_ptr<psi::Matrix> data() const {return data_;}

    /// Get the wavefunction
    virtual std::shared_ptr<psi::Wavefunction> wfn() const {return wfn_;}

    /// Get the information if data is already computed or not
    virtual bool is_computed() const {return isComputed_;}


    // <--- Computers ---> //

    /// Compute the scalar field in each point from the point collection
    virtual void compute();

    /// Compute a value of scalar field at point (x, y, z)
    virtual double compute_xyz(const double& x, const double& y, const double& z) = 0;


    // <--- Printers ---> //

    /// Write the cube file (only for Cube collections, otherwise does nothing)
    virtual void write_cube_file(const std::string& name);

    /// Print information of the object to Psi4 output
    virtual void print() const = 0;


  protected:

    /// Collection of points at which the scalar field is to be computed
    std::shared_ptr<PointsCollection3D> pointsCollection_;

    /// The data matrix in a form { [x, y, z, f(x, y, z)] }
    std::shared_ptr<psi::Matrix> data_;

    /// Wavefunction
    std::shared_ptr<psi::Wavefunction> wfn_;

    /// Geometry of a molecule
    psi::Matrix geom_;
 
    /// Integral factory
    std::shared_ptr<psi::IntegralFactory> fact_;

    /// Matrix of potential one-electron integrals
    std::shared_ptr<psi::Matrix> pot_;

    /// One-electron integral shared pointer
    std::shared_ptr<psi::OneBodyAOInt> oneInt_;

    /// One-electron potential shared pointer
    std::shared_ptr<     PotentialInt> potInt_;

    /// Basis set
    std::shared_ptr<psi::BasisSet> primary_;

    /// Number of basis functions
    int nbf_;

    /// Has data already computed?
    bool isComputed_;


  private:

    /// Initialize common defaults
    void common_init(void);

}; 

class ElectrostaticPotential3D : public ScalarField3D
{
  public:
    ElectrostaticPotential3D(const int& np, const double& radius, const double& cx, const double& cy, const double& cz);
    ElectrostaticPotential3D(const int& np, const double& padding, psi::SharedWavefunction wfn, psi::Options& options);
    ElectrostaticPotential3D(const int& nx, const int& ny, const int& nz,
                             const double& px, const double& py, const double& pz,
                             psi::SharedWavefunction wfn, psi::Options& options);

    virtual ~ElectrostaticPotential3D();

    virtual double compute_xyz(const double& x, const double& y, const double& z);
    virtual void print() const;
};


/** \brief Class template for OEP scalar fields.
 *
 *  @tparam T the compatible class (e.g. OEPotential)
 */
template<class T>
class OEPotential3D : public ScalarField3D
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
    //TODO!
    OEPotential3D(const int& np, const double& padding, std::shared_ptr<T> oep, const std::string& oepType);

    /** \brief Construct ordered 3D collection of scalar field of type T.
     *
     *  The points are generated according to Gaussian cube file format.
     *  @param nx         - number of points along x direction
     *  @param ny         - number of points along y direction
     *  @param nz         - number of points along z direction
     *  @param px         - padding distance along x direction
     *  @param py         - padding distance along y direction
     *  @param pz         - padding distance along z direction
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
 : ScalarField3D(np, padding, oep->wfn()->molecule()), oep_(oep), oepType_(oepType)
{

}

template <class T>
OEPotential3D<T>::OEPotential3D(const int& nx, const int& ny, const int& nz,
                                const double& px, const double& py, const double& pz,
                                std::shared_ptr<T> oep, const std::string& oepType, psi::Options& options)
 : ScalarField3D(nx, ny, nz, px, py, pz, oep->wfn(), options), oep_(oep), oepType_(oepType)
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

} // EndNameSpace oepdev
#endif //_oepdev_libutil_potential_h
