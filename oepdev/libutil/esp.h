#ifndef _oepdev_libutil_esp_h
#define _oepdev_libutil_esp_h

#include<string>

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "../liboep/oep.h"

#include <iostream>
#include <memory>
#include <cstdio>
#include <stdio.h>
#include <random>
#include <cmath>

#define _USE_MATH_DEFINES
 
using namespace std;

class CubeGrid
{
  public:
   CubeGrid(int, int, int, double, double, double, double, double, double);
   void print_cube() const;
  protected:
   const int nx_, ny_, nz_;
   const double dx_, dy_, dz_, ox_, oy_, oz_;
};

CubeGrid::CubeGrid(int nx, int ny, int nz, double dx, double dy, double dz, double ox, double oy, double oz) 
 : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), ox_(ox), oy_(oy), oz_(oz) {}
void CubeGrid::print_cube() const {cout << " Cube Props: nx, ny, nz = " << nx_ << ny_ << nz_ << "\n";}

// Points3DIterator
class Points3DIterator
{
  public:
    Points3DIterator(const int& np);
    virtual ~Points3DIterator() {};

    // Cube
    static shared_ptr<Points3DIterator> build(const int& nx, const int& ny, const int& nz,
                                             const double& dx, const double& dy, const double& dz,
                                             const double& ox = 0.0, const double& oy = 0.0, const double& oz = 0.0);
    // Random
    static shared_ptr<Points3DIterator> build(const int& np, const double& cx, const double& cy, const double& cz,
                                              const double& radius);

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

Points3DIterator::Points3DIterator(const int& np)
 : np_(np), done_(false), 
   index_(0)
{ 
    current_.index = 0;
}


class CubePoints3DIterator : public Points3DIterator
{
  public:
    CubePoints3DIterator(const int& nx, const int& ny, const int& nz,
                         const double& dx, const double& dy, const double& dz,
                         const double& ox, const double& oy, const double& oz);

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


CubePoints3DIterator::CubePoints3DIterator(
                        const int& nx, const int& ny, const int& nz,
                        const double& dx, const double& dy, const double& dz,
                        const double& ox, const double& oy, const double& oz)
 : Points3DIterator(nx*ny*nz), nx_(nx), ny_(ny), nz_(nz),
                              dx_(dx), dy_(dy), dz_(dz),
                              ox_(ox), oy_(oy), oz_(oz),
                              ii_(0), jj_(0), kk_(0),
                              ix_max__(nx-1), iy_max__(ny-1), iz_max__(nz-1)
{
}

void CubePoints3DIterator::first() 
{
  current_.x = ox_;
  current_.y = oy_;
  current_.z = oz_;
}
void CubePoints3DIterator::next() 
{
  ++kk_; ++index_;
  if(kk_ > iz_max__) {
     kk_ = 0;
     ++jj_;
     if(jj_ > iy_max__) {
        jj_ = 0;
        ++ii_;
        if(ii_ > ix_max__) {
           done_ = true;
        }
     }
  }
  current_.x = ox_ + dx_ * (double)ii_;
  current_.y = oy_ + dy_ * (double)jj_;
  current_.z = oz_ + dz_ * (double)kk_;
  current_.index = index_;
}

// Random Points Iterator
class RandomPoints3DIterator : public Points3DIterator
{
  public:
    RandomPoints3DIterator(const int& np, const double& cx, const double& cy, const double& cz,
                                          const double& radius);

    virtual void first();
    virtual void next();                        

  protected:
    const double cx_, cy_, cz_;
    const double radius_;
    double r_, phi_, theta_;

    virtual double random_double() {return random_distribution_(generator_);}
    virtual void draw_random_point();
    std::default_random_engine generator_;
    std::uniform_real_distribution<double> random_distribution_;
};

RandomPoints3DIterator::RandomPoints3DIterator(const int& np, const double& cx, const double& cy, const double& cz,
                                                              const double& radius)
 : Points3DIterator(np), cx_(cx), cy_(cy), cz_(cz), radius_(radius),
   random_distribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   
}

void RandomPoints3DIterator::draw_random_point() 
{
  double val;
  theta_     = acos(random_double());
  phi_       = 2.0 * M_PI * random_double();
  r_         = radius_  * cbrt( random_double() );
  current_.x = cx_ + r_ * sin(theta_) * cos(phi_);
  current_.y = cy_ + r_ * sin(theta_) * sin(phi_);
  current_.z = cz_ + r_ * cos(theta_);
}

void RandomPoints3DIterator::first()
{
   draw_random_point();
}

void RandomPoints3DIterator::next()
{
   ++index_;
   if (index_ == np_) done_ = true;
   draw_random_point();
}


// Distribution 3D
class Distribution3D
{
  public:
    enum Distribution {Random, Cube};
    Distribution3D(Distribution distribution);
    virtual ~Distribution3D() {};
    static shared_ptr<Distribution3D> build(const int& npoints, const double& cx, const double& cy, const double& cz,
                                            const double& radius);
    static shared_ptr<Distribution3D> build(const int& nx, const int& ny, const int& nz,
                                            const double& dx, const double& dy, const double& dz,
                                            const double& ox = 0.0, const double& oy = 0.0, const double& oz = 0.0);

    virtual shared_ptr<Points3DIterator> points_iterator() const {return pointsIterator_;}

    virtual void print_header() const = 0;
  protected:
    Distribution type_;
    shared_ptr<Points3DIterator> pointsIterator_;
};


Distribution3D::Distribution3D(Distribution distribution) : type_(distribution) {}
void Distribution3D::print_header() const {}

// Random Distribution
class RandomDistribution3D : public Distribution3D
{
  public:
    RandomDistribution3D(Distribution distribution, const int& npoints,
                         const double& cx, const double& cy, const double& cz, const double& radius);
    virtual void print_header() const;
  protected:
    const int npoints_;
};

RandomDistribution3D::RandomDistribution3D(Distribution distribution, const int& npoints,
                const double& cx, const double& cy, const double& cz,
                const double& radius)
 : Distribution3D(distribution), npoints_(npoints)
{
  pointsIterator_ = make_shared<RandomPoints3DIterator>(npoints, cx, cy, cz, radius);
}

void RandomDistribution3D::print_header() const {cout << "I am a RANDOM distribution with np = " << npoints_ << "\n";}

// Cube Distribution
class CubeDistribution3D : public Distribution3D, public CubeGrid
{
  public:
    CubeDistribution3D(Distribution distribution, const int& nx, const int& ny, const int& nz,
                                                  const double& dx, const double& dy, const double& dz,
                                                  const double& ox, const double& oy, const double& oz);
    virtual void print_header() const;
  protected:
};

CubeDistribution3D::CubeDistribution3D(Distribution distribution, const int& nx, const int& ny, const int& nz,
                                                                  const double& dx, const double& dy, const double& dz,
                                                                  const double& ox, const double& oy, const double& oz)
 : Distribution3D(distribution), CubeGrid(nx, ny, nz, dx, dy, dz, ox, oy, oz)
{
   pointsIterator_ = make_shared<CubePoints3DIterator>(nx, ny, nz, dx, dy, dz, ox, oy, oz);
  // common_init(); 
}

void CubeDistribution3D::print_header() const {
  cout << "I am a CUBE distribution with nx, ny, nz:\n";
  print_cube();
}

// 3D Potential!
class Potential3D 
{
  public:
    Potential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius);
    Potential3D(const int& nx, const int& ny, const int& nz, 
                const double& dx, const double& dy, const double& dz,
                const double& ox, const double& oy, const double& oz);
    virtual ~Potential3D() {};
    static shared_ptr<Potential3D> build(const std::string& type, const int& np,
                                         const double& cx, const double& cy, const double& cz, const double& radius);
    static shared_ptr<Potential3D> build(const std::string& type, const int& nx, const int& ny, const int& nz,
                                         const double& dx, const double& dy, const double& dz,
                                         const double& ox = 0.0, const double& oy = 0.0, const double& oz = 0.0);

    virtual void compute();
    virtual shared_ptr<Distribution3D> distribution() const {return distribution_;}

    virtual double compute_xyz(const double& x, const double& y, const double& z) = 0;
    virtual void print_header() const = 0;

  protected:
    shared_ptr<Distribution3D> distribution_;
};

class EPotential3D : public Potential3D
{
  public:
    EPotential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius);
    EPotential3D(const int& nx, const int& ny, const int& nz,
                 const double& dx, const double& dy, const double& dz,
                 const double& ox, const double& oy, const double& oz);

    virtual double compute_xyz(const double& x, const double& y, const double& z);
    virtual void print_header() const;
  protected:
};

void Potential3D::compute(){
  shared_ptr<Points3DIterator> iter = distribution()->points_iterator();
  if (!iter) cout << "ERROR!!!\n";
  for (iter->first(); iter->is_done()==false; iter->next()) {
       double val = compute_xyz(iter->x(), iter->y(), iter->z());
       cout << iter->x() << "\n";
  ////     distribution_->set(iter.index(), val);
  }
  iter.reset();
}

Potential3D::Potential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius) 
 : distribution_(Distribution3D::build(np, cx, cy, cz, radius))
{
   
}

Potential3D::Potential3D(const int& nx, const int& ny, const int& nz,
                         const double& dx, const double& dy, const double& dz,
                         const double& ox, const double& oy, const double& oz)
  : distribution_(Distribution3D::build(nx, ny, nz, dx, dy, dz, ox, oy, oz))
{
  //distribution_ = Distribution3D::build(nx, ny, nz, dx, dy, dz, ox, oy, oz);
}

EPotential3D::EPotential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius) 
 : Potential3D(np, cx, cy, cz, radius)
{

}

EPotential3D::EPotential3D(const int& nx, const int& ny, const int& nz,
                           const double& dx, const double& dy, const double& dz,
                           const double& ox, const double& oy, const double& oz) 
 : Potential3D(nx, ny, nz, dx, dy, dz, ox, oy, oz)
{

}

void EPotential3D::print_header() const {cout << " I am Electrostatic potential:\n";distribution_->print_header();}

double EPotential3D::compute_xyz(const double& x, const double& y, const double& z)
{
   double val = 1.0;
   return val;
}




// Factory methods
shared_ptr<Distribution3D> Distribution3D::build(const int& npoints, 
           const double& cx, const double& cy, const double& cz, const double& radius)
{
  shared_ptr<Distribution3D> distribution = make_shared<RandomDistribution3D>(Distribution3D::Random,npoints,cx,cy,cz,radius);
  return distribution;
}

shared_ptr<Distribution3D> Distribution3D::build(const int& nx, const int& ny, const int& nz,
                                                 const double& dx, const double& dy, const double& dz,
                                                 const double& ox, const double& oy, const double& oz){
  shared_ptr<Distribution3D> distribution = make_shared<CubeDistribution3D>(Distribution3D::Cube, nx, ny, nz, dx, dy, dz, ox, oy, oz);
  return distribution;
}

shared_ptr<Potential3D> Potential3D::build(const std::string& type, 
                                           const int& np,
                                           const double& cx, const double& cy, const double& cz,
                                           const double& radius)
{
   shared_ptr<Potential3D> potential;
   if          (type == "E") potential = make_shared<EPotential3D>(np, cx, cy, cz, radius);
   else cout << " ERROR! No such potential type <" << type << "> !\n";
   return potential;
}

shared_ptr<Potential3D> Potential3D::build(const std::string& type, 
                                           const int& nx, const int& ny, const int& nz,
                                           const double& dx, const double& dy, const double& dz,
                                           const double& ox, const double& oy, const double& oz)
{
   shared_ptr<Potential3D> potential;
   if          (type == "E") potential = make_shared<EPotential3D>(nx, ny, nz, dx, dy, dz, ox, oy, oz);
   else cout << " ERROR! No such potential type <" << type << "> !\n";
   return potential;
}

shared_ptr<Points3DIterator> Points3DIterator::build(const int& nx, const int& ny, const int& nz,
                                                     const double& dx, const double& dy, const double& dz,
                                                     const double& ox, const double& oy, const double& oz)
{
  shared_ptr<Points3DIterator> iterator = make_shared<CubePoints3DIterator>(nx, ny, nz, dx, dy, dz, ox, oy, oz);
  return iterator;
}

shared_ptr<Points3DIterator> Points3DIterator::build(const int& np, 
                                                     const double& cx, const double& cy, const double& cz,
                                                     const double& radius)
{
 shared_ptr<Points3DIterator> iterator = make_shared<RandomPoints3DIterator>(np, cx, cy, cz, radius);
 return iterator;
}







// -------------------------------------------------- MAIN ------------------------------------ 
int main(int argc, char *argv[])
{
  cout << "Welcome to the Błasiak group!\n";
  // Create bare grid
  unique_ptr<CubeGrid> grid(new CubeGrid(3,4,5,1.0,1.0,1.0,0.0,0.0,0.0));
  grid->print_cube();

  // Create some distributions
  shared_ptr<Distribution3D> distrRand = Distribution3D::build(45, 0.0, 0.0, 0.0, 1.0);
  shared_ptr<Distribution3D> distrCube = Distribution3D::build(1,2,9,1.0,1.0,1.0);
  distrRand->print_header();
  distrCube->print_header();

  // Create some potentials
  shared_ptr<Potential3D> epot = Potential3D::build("E",2,2,2,1.0,1.0,1.0);
  epot->print_header();
  epot->compute();

  // Create some iterators
  //shared_ptr<Points3DIterator> iter = Points3DIterator::build(2,2,2,0.2,0.2,0.2);
  //for (iter->first(); iter->is_done()==false; iter->next()) {
  //     fprintf(stdout, " Point X, Y, Z: %13.5f %13.5f %13.5f\n", iter->x(), iter->y(), iter->z());
  //}
  shared_ptr<Points3DIterator> iter = Points3DIterator::build(10000, 0.0, 0.0, 0.0, 1.0);
  for (iter->first(); iter->is_done()==false; iter->next()) {
       fprintf(stdout, " Point X, Y, Z: %13.5f %13.5f %13.5f\n", iter->x(), iter->y(), iter->z());
  }
  

  
  

  return 0;
}



//--------------------------------------------------------------------------
namespace oepdev{

using namespace psi;
using namespace std;


using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedOEPotential        = std::shared_ptr<OEPotential>;

// Points3DIterator
class Points3DIterator
{
  public:
    Points3DIterator(const int& np);
    virtual ~Points3DIterator() {};

    static shared_ptr<Points3DIterator> build(const int& nx, const int& ny, const int& nz,
                                             const double& dx, const double& dy, const double& dz,
                                             const double& ox = 0.0, const double& oy = 0.0, const double& oz = 0.0);
    static shared_ptr<Points3DIterator> build(SharedMatrix points);

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
    int ii_, jj_, kk_, index_;
};

Points3DIterator::Points3DIterator(const int& np)
 : np_(np), done_(false), 
   index_(0), ii_(0), jj_(0), kk_(0)
{ 
    current_.index = 0;
}


class CubePoints3DIterator : public Points3DIterator
{
  public:
    CubePoints3DIterator(const int& nx, const int& ny, const int& nz,
                         const double& dx, const double& dy, const double& dz,
                         const double& ox, const double& oy, const double& oz);

    virtual void first();
    virtual void next();                        
  protected:
    const int nx_, ny_, nz_;
    const double dx_, dy_, dz_;
    const double ox_, oy_, oz_;
  private:
    const int ix_max__, iy_max__, iz_max__;
};


CubePoints3DIterator::CubePoints3DIterator(
                        const int& nx, const int& ny, const int& nz,
                        const double& dx, const double& dy, const double& dz,
                        const double& ox, const double& oy, const double& oz)
 : Points3DIterator(nx*ny*nz), nx_(nx), ny_(ny), nz_(nz),
                               dx_(dx), dy_(dy), dz_(dz),
                               ox_(ox), oy_(oy), oz_(oz),
                               ix_max__(nx-1), iy_max__(ny-1), iz_max__(nz-1)
{
}

void CubePoints3DIterator::first() 
{
  current_.x = ox_;
  current_.y = oy_;
  current_.z = oz_;
}
void CubePoints3DIterator::next() 
{
  ++kk_; ++index_;
  if(kk_ > iz_max__) {
     kk_ = 0;
     ++jj_;
     if(jj_ > iy_max__) {
        jj_ = 0;
        ++ii_;
        if(ii_ > ix_max__) {
           done_ = true;
        }
     }
  }
  current_.x = ox_ + dx_ * (double)ii_;
  current_.y = oy_ + dy_ * (double)jj_;
  current_.z = oz_ + dz_ * (double)kk_;
  current_.index = index_;
}

class RandomPoints3DIterator : public Points3DIterator
{
  public:
    RandomPoints3DIterator(const int& np);

    virtual void first();
    virtual void next();                        
  protected:
  private:
};


RandomPoints3DIterator::RandomPoints3DIterator(const int& np)
 : Points3DIterator(np)
{
}

void RandomPoints3DIterator::first() 
{
  current_.x = 0.0;
  current_.y = 0.0;
  current_.z = 0.0;
}
void RandomPoints3DIterator::next() 
{
  //current_.x = ox_ + dx_ * (double)ii_;
  //current_.y = oy_ + dy_ * (double)jj_;
  //current_.z = oz_ + dz_ * (double)kk_;
  //current_.index = index_;
}



/** \brief 3-dimensional scalar potential manifold: Abstract base.
 *
 */
class Potential3D
{
  public:
   /**
    * The distribution of points in space.
    *
    * Random = randomized uniform distribution of points excluding vdW spheres
    * Cube   = distribution in a form of a cube grid
    */
    enum Distribution {Random, Cube};

    Potential3D(SharedMatrix pot);
    Potential3D(SharedWavefunction wfn, Distribution distribution,
                const double& padx, const double& pady, const double& padz);
    Potential3D(SharedOEPotential oep, const std::string& oepType, Distribution distribution,
                const double& padx, const double& pady, const double& padz);

    /// Build RandomBare potential
    static std::shared_ptr<Potential3D> build(SharedMatrix pot);
    static std::shared_ptr<Potential3D> build(const std::string& potType, SharedWavefunction wfn, 
                                              const double& padx, 
                                              const double& pady,
                                              const double& padz);
    static std::shared_ptr<Potential3D> build(SharedOEPotential oep, const std::string& oepType, 
                                              Distribution d = Random);

    virtual ~Potential3D();

    virtual SharedMatrix get(void) const {return pot_;}

    virtual void compute() = 0;

  protected:
    SharedMatrix pot_;
  private:
    void common_init();
};

class BarePotential3D : public Potential3D
{
  public:
    BarePotential3D(SharedMatrix pot);
    virtual void compute();
  protected:
  private:
    void common_init();
};

class ElectrostaticPotential3D : public Potential3D
{
  public:
    ElectrostaticPotential3D(SharedWavefunction wfn);
    virtual void compute();
  protected:
  private:
    void common_init();
};

class OEPotential3D : public Potential3D
{
  public:
    OEPotential3D(SharedOEPotential oep, const std::string& oepType);
    virtual void compute();
  protected:
  private:
    void common_init();
};

/** \brief Charges from Electrostatic Potential (ESP): abstract base of ESP method.
 *
 */
//class ESP
//{
//  public:
//    //@{ Use atoms as ESP centres
//    ESP(const Wavefunction& wfn);
//    ESP(const OEPotential& oep, const std::string& oepType);
//
//    static std::shared_ptr<ESP> build(const std::string& method, const Wavefunction& wfn);
//    static std::shared_ptr<ESP> build(const std::string& method, const OEPotential& oep, const std::string& oepType);
//    //@}
//
//    virtual ~ESP();
//
//    virtual void fit();
//
//    virtual void compute_potential() = 0;
//    virtual void print_header() const = 0;
//
//  protected:
//    SharedPotential3D potRef_;
//    SharedPotential3D pot_;
//    
//  private:
//    void common_init(void);
//};

//class OEPotentialESP : public ESP
//{
//  public:
//    OEPotentialESP(SharedOEPotential oep, const std::string& oepType);
//  protected:
//  private:
//};

} // EndNameSpace oepdev
#endif //_oepdev_libutil_esp_h
