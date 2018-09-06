#ifndef _oepdev_libutil_potential_h
#define _oepdev_libutil_potential_h
/** @file spce3d.h */

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
#include "psi4/libmints/vector.h"
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
/** \addtogroup OEPDEV_3DFIELDS
 * @{
 */


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
     *   @param np   - number of points this iterator is constructed for
     */
    Points3DIterator(const int& np);

    /// Destructor
    virtual ~Points3DIterator();


    // <--- Factory Methods ---> //

    /** \brief Build G09 Cube collection iterator.
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

    /** \brief Build random collection iterator.
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

    /** \brief Build random collection iterator.
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


/** \brief Iterator over a collection of points in 3D space. g09 Cube-like order.
 *
 *  __Note:__ Always create instances by using static factory method from Points3DIterator.
 *            Do not use constructor of this class.
 */
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


/** \brief Iterator over a collection of points in 3D space. Random collection.
 *
 *  __Note:__ Always create instances by using static factory method from Points3DIterator.
 *            Do not use constructors of this class.
 */
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

    /** @brief Build random collection of points.
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

    /** @brief Build random collection of points.
     *
     *  Points uniformly span space inside a sphere enclosing a molecule. 
     *  exluding the van der Waals volume.
     *  @param np         - number of points to draw
     *  @param padding    - radius padding of a minimal sphere enclosing the molecule
     *  @param mol        - Psi4 molecule object
     */
    static shared_ptr<PointsCollection3D> build(const int& npoints, const double& padding, psi::SharedMolecule mol);

    /** @brief Build G09 Cube collection of points.
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

    /// Print the information to Psi4 output file
    virtual void print() const = 0;


  protected:

    /// Number of points
    const int np_;

    /// Collection type
    Collection collectionType_;

    /// iterator over points collection
    shared_ptr<Points3DIterator> pointsIterator_;
};


/** \brief Collection of random points in 3D space.
 *
 * __Note:__ Do not use constructors of this class explicitly. Instead,
 * use static factory methods of the superclass to create instances.
 */
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


/** \brief G09 cube-like ordered collection of points in 3D space.
 *
 * __Note:__ Do not use constructors of this class explicitly. Instead,
 * use static factory methods of the superclass to create instances.
 */
class CubePointsCollection3D : public PointsCollection3D, public psi::CubicScalarGrid
{
  public:
    CubePointsCollection3D(Collection collectionType, const int& nx, const int& ny, const int& nz,
                                                      const double& px, const double& py, const double& pz,
                                                      psi::SharedBasisSet bs, psi::Options& options);

    virtual ~CubePointsCollection3D();

    virtual void print() const;

    virtual void write_cube_file(psi::SharedMatrix v, const std::string& name, const int& col = 0);
  protected:
};

/** \brief General Vector Dield in 3D Space. Abstract base.
 *
 *  Create vector field defined at points distributed randomly or as an ordered 
 *  g09 cube-like collection. Currently implemented fields are:
 *
 *   - Electrostatic potential     - computes electrostatic potential (requires wavefunction)
 *   - Template of generic classes - compute custom vector fields (requires generic object 
 *                                   that is able to compute the field in 3D space)
 *
 *  __Note:__ Always create instances by using static factory methods `build`.
 *  The following types of 3D vector fields are currently implemented:
 *   - `ELECTROSTATIC POTENTIAL`
 */
class Field3D
{
  public:
    
    // <--- Constructors and Destructor ---> //

    /// Construct potential on random grid by providing wavefunction. Excludes space within vdW volume
    Field3D(const int& ndim, const int& np, const double& pad, psi::SharedWavefunction wfn, psi::Options& options);

    /// Construct potential on random grid by providing molecule. Excludes space within vdW volume
    /// Field3D(const int& ndim, const int& np, const double& pad, psi::SharedMolecule mol, psi::Options& options);

    /// Construct potential on cube grid by providing wavefunction
    Field3D(const int& ndim, 
            const int& nx, const int& ny, const int& nz,
            const double& px, const double& py, const double& pz,
            std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options);

    /// Destructor
    virtual ~Field3D();


    // <--- Factory Methods ---> //

    /**\brief Build 3D field of random points. vdW volume is excluded.
     *
     *  @param ndim    - dimensionality of 3D field (1: scalar field, >2: vector field)
     *  @param type    - type of 3D field                                      
     *  @param np      - number of points 
     *  @param pad     - radius padding of a minimal sphere enclosing the molecule
     *  @param wfn     - Psi4 Wavefunction containing the molecule
     *  @param options - Psi4 options
     */
    static shared_ptr<Field3D> build(const std::string& type, const int& np,
                                     const double& pad, psi::SharedWavefunction wfn, psi::Options& options,
                                     const int& ndim = 1);

    /**\brief Build 3D field of points on a g09-cube grid.
     *
     *  @param ndim    - dimensionality of 3D field (1: scalar field, >2: vector field)
     *  @param type    - type of 3D field                                      
     *  @param nx      - number of points along x direction
     *  @param ny      - number of points along y direction
     *  @param nz      - number of points along z direction
     *  @param px      - padding distance along x direction
     *  @param py      - padding distance along y direction
     *  @param pz      - padding distance along z direction
     *  @param wfn     - Psi4 Wavefunction containing the molecule
     *  @param options - Psi4 options
     */
    static shared_ptr<Field3D> build(const std::string& type, const int& nx, const int& ny, const int& nz,
                                     const double& px, const double& py, const double& pz,
                                     psi::SharedWavefunction wfn, psi::Options& options,
                                     const int& ndim = 1);


    // <--- Accessors ---> //

    /// Get the number of points at which the 3D field is defined
    virtual int npoints() const {return pointsCollection_->npoints();}

    /// Get the collection of points
    virtual std::shared_ptr<PointsCollection3D> points_collection() const {return pointsCollection_;}

    /// Get the data matrix in a form { [x, y, z, f_1(x, y, z), f_2(x, y, z), ... , f_n(x, y, z) ] } where n = ndim
    virtual std::shared_ptr<psi::Matrix> data() const {return data_;}

    /// Get the wavefunction
    virtual std::shared_ptr<psi::Wavefunction> wfn() const {return wfn_;}

    /// Get the information if data is already computed or not
    virtual bool is_computed() const {return isComputed_;}

    /// Get the number of fields
    int dimension() const {return nDim_;}


    // <--- Computers ---> //

    /// Compute the 3D field in each point from the point collection
    virtual void compute();

    /// Compute a value of 3D field at point (x, y, z)
    virtual std::shared_ptr<psi::Vector> compute_xyz(const double& x, const double& y, const double& z) = 0;


    // <--- Printers ---> //

    /// Write the cube file (only for Cube collections, otherwise does nothing)
    virtual void write_cube_file(const std::string& name);

    /// Print information of the object to Psi4 output
    virtual void print() const = 0;


  protected:

    /// Collection of points at which the 3D field is to be computed
    std::shared_ptr<PointsCollection3D> pointsCollection_;

    /// The data matrix in a form { [x, y, z, f_1(x, y, z), f_2(x, y, z), ..., f_n(x, y, z) ] } where n = nDim_
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

    /// Dimensionality of the 3D field (1: scalar field, 2>: vector field)
    int nDim_;

    /// Has data already computed?
    bool isComputed_;


  private:

    /// Initialize common defaults
    void common_init(void);

}; 

/** \brief Electrostatic potential of a molecule.
 *
 *  Computes the electrostatic potential of a molecule directly from the wavefunction.
 *  The electrostatic potential \f$ v({\bf r})\f$ at point \f${\bf r}\f$ is computed 
 *  from the following formula:
 *  \f[
 *     v({\bf r}) = v_{\rm nuc}({\bf r}) + v_{\rm el}({\bf r})
 *  \f]
 *  where the nuclear and electronic contributions are defined accordingly as
 *  \f{align*}{
 *    v_{\rm nuc} ({\bf r})&= \sum_x  \frac{Z_x}{\left| {\bf r} - {\bf r}_x \right|} \\
 *    v_{\rm el}  ({\bf r})&= \sum_{\mu\nu} \left\{ D_{\mu\nu}^{(\alpha)} + D_{\mu\nu}^{(\beta)}\right\} 
 *                              V_{\nu\mu}({\bf r})
 *  \f}
 *  In the above equations, \f$ Z_x \f$ denotes the charge of \f$ x \f$th nucleus, 
 *  \f$ D_{\mu\nu}^{(\omega)} \f$ is the one-particle (relaxed) density matrix element
 *  in AO basis associated with the \f$ \omega \f$ electron spin, and \f$ V_{\mu\nu}({\bf r}) \f$
 *  is the potential one-electron integral defined by
 *  \f[
 *    V_{\nu\mu}({\bf r}) \equiv \int d{\bf r}' \varphi^{*}_\nu({\bf r}') 
 *                                              \frac{1}{\left| {\bf r} - {\bf r}' \right|} 
 *                                               \varphi_\mu({\bf r}')
 *  \f]
 *  
 */
class ElectrostaticPotential3D : public Field3D
{
  public:
    ElectrostaticPotential3D(const int& np, const double& padding, psi::SharedWavefunction wfn, psi::Options& options);
    ElectrostaticPotential3D(const int& nx, const int& ny, const int& nz,
                             const double& px, const double& py, const double& pz,
                             psi::SharedWavefunction wfn, psi::Options& options);

    virtual ~ElectrostaticPotential3D();

    virtual std::shared_ptr<psi::Vector> compute_xyz(const double& x, const double& y, const double& z);
    virtual void print() const;
};

/** \brief Class template for OEP 3D fields.
 *
 *  Used for special type of classes T that contain following public member functions:
 *
 *  \code{.cpp}
 *   class T : public std::enable_shared_from_this<T> {
 * 
 *     public:
 *        void compute_3D(const std::string& descriptor, 
 *                        const double& x, const double& y, const double& z,
 *                        double& v);
 *                                                                              
 *        shared_ptr<psi::Wavefunction> wfn() const {return wfn_;}
 * 
 *     /// The rest of the declaration of T
 *   };
 *  \endcode
 *
 *  with the `descriptor` of a certain 3D field type,
 *  `x`, `y`, `z` the points in 3D space in which the scalar or vector field
 *  has to be computed and stored at `v`. Instances of `T` should store shared
 *  pointer to wavefunction object. List of classes `T` that are
 *  compatible with this class template and are currently implemented in oepdev 
 *  is given below:
 *    - `oepdev::OEPotential` abstract base (do not use derived classes as `T`)
 *
 *  Template parameters:
 *  @tparam T the compatible class (e.g. `oepdev::OEPotential`)
 */
template<class T>
class OEPotential3D : public Field3D
{
  public:

    // <--- Constructors and Destructor ---> //

    /** \brief Construct random spherical collection of 3D field of type T.
     *
     *  The points are drawn according to uniform distrinution in 3D space.
     *  @param ndim       - dimensionality of 3D field (1: scalar field, >2: vector field)
     *  @param np         - number of points to draw
     *  @param padding    - spherical padding distance (au)
     *  @param oep        - OEP object of type T
     *  @param oepType    - type of OEP
     */
    //TODO!
    OEPotential3D(const int& ndim, const int& np, const double& padding, std::shared_ptr<T> oep, const std::string& oepType);

    /** \brief Construct ordered 3D collection of 3D field of type T.
     *
     *  The points are generated according to Gaussian cube file format.
     *  @param ndim       - dimensionality of 3D field (1: scalar field, >2: vector field)
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
    OEPotential3D(const int& ndim, 
                  const int& nx, const int& ny, const int& nz, 
                  const double& px, const double& py, const double& pz,
                  std::shared_ptr<T> oep, const std::string& oepType, psi::Options& options);

    /// Destructor
    virtual ~OEPotential3D();
    
    virtual std::shared_ptr<psi::Vector> compute_xyz(const double& x, const double& y, const double& z);
    virtual void print() const;

  protected:
    /// Shared pointer to the instance of class `T`
    std::shared_ptr<T> oep_;
    /// Descriptor of the 3D field type stored in instance of `T`
    std::string oepType_;
    
};

template <class T>
OEPotential3D<T>::OEPotential3D(const int& ndim, const int& np, const double& padding, std::shared_ptr<T> oep, const std::string& oepType)
 : Field3D(ndim, np, padding, oep->wfn(), oep->wfn()->options()), oep_(oep), oepType_(oepType)
{

}

template <class T>
OEPotential3D<T>::OEPotential3D(const int& ndim, const int& nx, const int& ny, const int& nz,
                                const double& px, const double& py, const double& pz,
                                std::shared_ptr<T> oep, const std::string& oepType, psi::Options& options)
 : Field3D(ndim, nx, ny, nz, px, py, pz, oep->wfn(), options), oep_(oep), oepType_(oepType)
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
std::shared_ptr<psi::Vector> OEPotential3D<T>::compute_xyz(const double& x, const double& y, const double& z)
{
   std::shared_ptr<psi::Vector> v = std::make_shared<psi::Vector>("", nDim_);
   oep_->compute_3D(oepType_, x, y, z, v);
   return v;
}


/** @}*/ 

} // EndNameSpace oepdev


#endif //_oepdev_libutil_potential_h
