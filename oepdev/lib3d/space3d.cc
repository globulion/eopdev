#include "space3d.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "../libutil/util.h"


namespace oepdev{

using namespace psi;
using namespace std;

Points3DIterator::Points3DIterator(const int& np)
 : np_(np), done_(false), 
   index_(0)
{ 
    current_.index = 0;
}

Points3DIterator::~Points3DIterator()
{

}

void Points3DIterator::rewind() {
  done_ = false;
  index_= 0;
  current_.index = 0;
}

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

CubePoints3DIterator::~CubePoints3DIterator()
{

}

void CubePoints3DIterator::first() 
{
  ii_ = 0; jj_ = 0; kk_ = 0;
  current_.x = ox_;
  current_.y = oy_;
  current_.z = oz_;
  rewind();
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

RandomPoints3DIterator::RandomPoints3DIterator(const int& np, const double& radius, 
                                               const double& cx, const double& cy, const double& cz)
 : Points3DIterator(np), cx_(cx), cy_(cy), cz_(cz), radius_(radius),
   randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   
}

RandomPoints3DIterator::RandomPoints3DIterator(const int& np, const double& pad, psi::SharedMolecule mol)
 : Points3DIterator(np), randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   // populate vdwRadius
   vdwRadius_["C"] = psi::Process::environment.options.get_double("ESP_VDW_RADIUS_C" );
   vdwRadius_["H"] = psi::Process::environment.options.get_double("ESP_VDW_RADIUS_H" );
   vdwRadius_["HE"] = psi::Process:environment.options.get_double("ESP_VDW_RADIUS_HE" );
   vdwRadius_["NE"] = psi::Process:environment.options.get_double("ESP_VDW_RADIUS_NE" );
   vdwRadius_["AR"] = psi::Process:environment.options.get_double("ESP_VDW_RADIUS_AR" );
   vdwRadius_["N"] = psi::Process::environment.options.get_double("ESP_VDW_RADIUS_N" );
   vdwRadius_["O"] = psi::Process::environment.options.get_double("ESP_VDW_RADIUS_O" );
   vdwRadius_["F"] = psi::Process::environment.options.get_double("ESP_VDW_RADIUS_F" );
   vdwRadius_["CL"]= psi::Process::environment.options.get_double("ESP_VDW_RADIUS_CL");
   vdwRadius_["LI"]= psi::Process::environment.options.get_double("ESP_VDW_RADIUS_LI");

   // Set the centre of the sphere
   Matrix       geom = mol->geometry();
   Vector3       com = mol->center_of_mass();
   cx_               = com.get(0);
   cy_               = com.get(1);
   cz_               = com.get(2);
  
   // compute radius (including padding)
   double maxval = geom.get(0,0) - com.get(0);
   double minval = geom.get(0,0) - com.get(0);
   double val, rad;
   for (int i = 0; i < geom.nrow(); ++i) {
        for (int j = 0; j < geom.ncol(); ++j) {
             val = geom.get(i,j) - com.get(j);
             if (val > maxval) maxval = val;
             if (val < minval) minval = val;
        }
   }
   rad = std::max(std::abs(maxval), std::abs(minval));
   radius_ = rad + pad;

   // populate excludeSpheres Matrix
   excludeSpheres_ = std::make_shared<Matrix>("van der Waals Exclude Spheres", geom.nrow(), 4);
   for (int i=0; i< geom.nrow(); ++i) {
        for (int j=0; j< geom.ncol(); ++j) {
             excludeSpheres_->set(i, j, geom.get(i,j));
        }
        excludeSpheres_->set(i, 3, vdwRadius_[mol->symbol(i)]);
   }
}

RandomPoints3DIterator::~RandomPoints3DIterator()
{

}

bool RandomPoints3DIterator::is_in_vdWsphere(double x, double y, double z) const
{
   bool isInside = false;
   double rx, ry, rz, r, rw;
   for (int i=0; i<excludeSpheres_->nrow(); ++i) {
        rx = excludeSpheres_->get(i, 0);
        ry = excludeSpheres_->get(i, 1);
        rz = excludeSpheres_->get(i, 2);
        rw = excludeSpheres_->get(i, 3);
        r  = sqrt(std::pow(rx-x, 2.0) + std::pow(ry-y, 2.0) + std::pow(rz-z, 2.0));
        if (r < rw) {
            isInside = true;
            break;
        }
   }
   return isInside;
}

void RandomPoints3DIterator::draw_random_point() 
{
  theta_     = acos(random_double());
  phi_       = 2.0 * M_PI * random_double();
  r_         = radius_  * cbrt( random_double() );
  x_         = cx_ + r_ * sin(theta_) * cos(phi_);
  y_         = cy_ + r_ * sin(theta_) * sin(phi_);
  z_         = cz_ + r_ * cos(theta_);

  if(excludeSpheres_) {
     while (is_in_vdWsphere(x_, y_, z_) == true) {               
            theta_     = acos(random_double());              
            phi_       = 2.0 * M_PI * random_double();
            r_         = radius_  * cbrt( random_double() );
            x_         = cx_ + r_ * sin(theta_) * cos(phi_); 
            y_         = cy_ + r_ * sin(theta_) * sin(phi_);
            z_         = cz_ + r_ * cos(theta_);
     }
  }
  
  current_.x = x_; 
  current_.y = y_;
  current_.z = z_;
}

void RandomPoints3DIterator::first()
{
   rewind();
   draw_random_point();
}

void RandomPoints3DIterator::next()
{
   ++index_;
   if (index_ == np_) done_ = true;
   draw_random_point();
   current_.index = index_;
}

PointsCollection3D::PointsCollection3D(Collection collectionType, int& np) : collectionType_(collectionType), np_(np) 
{

}
PointsCollection3D::PointsCollection3D(Collection collectionType, const int& np) : collectionType_(collectionType), np_(np) 
{

}


PointsCollection3D::~PointsCollection3D()
{

}

void PointsCollection3D::print() const 
{

}

RandomPointsCollection3D::RandomPointsCollection3D(Collection collectionType, 
                                           const int& npoints, const double& radius,                   
                                           const double& cx, const double& cy, const double& cz)
 : PointsCollection3D(collectionType, npoints)
{
  pointsIterator_ = make_shared<RandomPoints3DIterator>(npoints, radius, cx, cy, cz);
}

RandomPointsCollection3D::RandomPointsCollection3D(Collection collectionType, 
                                           const int& npoints,                                   
                                           const double& padding,
                                           psi::SharedMolecule mol)
 : PointsCollection3D(collectionType, npoints)
{
  pointsIterator_ = make_shared<RandomPoints3DIterator>(npoints, padding, mol);
}


RandomPointsCollection3D::~RandomPointsCollection3D()
{

}

void RandomPointsCollection3D::print() const 
{
   cout << "I am a RANDOM distribution with np = " << np_ << "\n";
}

CubePointsCollection3D::CubePointsCollection3D(Collection collectionType, 
                                       const int& nx, const int& ny, const int& nz,
                                       const double& px, const double& py, const double& pz, 
                                       psi::SharedBasisSet bs, psi::Options& options)
 : PointsCollection3D(collectionType, nx*ny*nz), CubicScalarGrid(bs, options)
{
  // Determine XYZmax and XYZmin
  double xmax = mol_->x(0);
  double ymax = mol_->y(0);
  double zmax = mol_->z(0);
  double xmin = xmax;
  double ymin = ymax;
  double zmin = zmax;
  double x, y, z;
  for (int i = 0; i < mol_->natom(); ++i) {
       x = mol_->x(i);
       y = mol_->y(i);
       z = mol_->z(i);
       if (x >= xmax) xmax = x;
       if (y >= ymax) ymax = y;
       if (z >= zmax) zmax = z;
       if (x <= xmin) xmin = x;
       if (y <= ymin) ymin = y;
       if (z <= zmin) zmin = z;
  }
  xmin -= px;  xmax += px;
  ymin -= py;  ymax += py;
  zmin -= pz;  zmax += pz;
  
  // Determine spacings
  double dx = (xmax - xmin) / (double)nx;
  double dy = (ymax - ymin) / (double)ny;
  double dz = (zmax - zmin) / (double)nz;
  
  // Determine cube origin
  double ox = xmin;
  double oy = ymin;
  double oz = zmin;
  
  // Point3D iterator
  pointsIterator_ = make_shared<CubePoints3DIterator>(nx, ny, nz, dx, dy, dz, ox, oy, oz);

  // Initialize the Grid
  int    N[3] = {nx-1, ny-1, nz-1};
  double D[3] = {dx  , dy  , dz  };
  double O[3] = {ox  , oy  , oz  }; 
  build_grid(".", N, D, O);
}

CubePointsCollection3D::~CubePointsCollection3D()
{

}

void CubePointsCollection3D::print() const {
  //psi::outfile->Printf(" ===> Cube 3D Collection <===\n");
}

void CubePointsCollection3D::write_cube_file(psi::SharedMatrix v, const std::string& name, const int& col){
    // => Drop the grid out <= //

    std::stringstream ss;
    ss << filepath_ << "/" << name << ".cube";

    // Is filepath a valid directory?
    std::ifstream infile(filepath_.c_str());
    if (!infile.good()) {
    //if (filesystem::path(filepath_).make_absolute().is_directory() == false) {
        printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath_.c_str());
        outfile->Printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath_.c_str());
        //outfile->Flush(); // -> incompatible with psi4-1.2 and newer
        exit(EXIT_FAILURE);
    }

    FILE* fh = fopen(ss.str().c_str(), "w");

    // Two comment lines
    fprintf(fh, "OepDev Gaussian Cube File.\n");
    fprintf(fh, "Property: %s\n", name.c_str());

    // Number of atoms plus origin of data
    fprintf(fh, "%5d %11.6f %11.6f %11.6f\n", mol_->natom(), O_[0], O_[1], O_[2]);

    //// Number of points along axis, displacement along x,y,z
    fprintf(fh, "%5d %11.6f %11.6f %11.6f\n", N_[0] + 1, D_[0], 0.0, 0.0);
    fprintf(fh, "%5d %11.6f %11.6f %11.6f\n", N_[1] + 1, 0.0, D_[1], 0.0);
    fprintf(fh, "%5d %11.6f %11.6f %11.6f\n", N_[2] + 1, 0.0, 0.0, D_[2]);

    //// Atoms of molecule (Z, Q?, x, y, z)
    for (int A = 0; A < mol_->natom(); A++) {
        fprintf(fh, "%5d %11.6f %11.6f %11.6f %11.6f\n", (int) mol_->Z(A), 0.0, mol_->x(A), mol_->y(A), mol_->z(A));
    }

    //// Data, striped (x, y, z)
    for (int ind = 0; ind < npoints_; ++ind) {
        fprintf(fh, "%13.5E", v->get(ind,3+col));
        if (ind % 6 == 5) fprintf(fh,"\n");
    }

    fclose(fh);
}

void Field3D::compute(){
  shared_ptr<Points3DIterator> iter = points_collection()->points_iterator();
  if (!iter) throw psi::PSIEXCEPTION("ERROR!!! No iterator in Potential3D object!\n");
  for (iter->first(); iter->is_done()==false; iter->next()) {
       std::shared_ptr<psi::Vector> v = compute_xyz(iter->x(), iter->y(), iter->z());
       data_->set(iter->index(), 0  , iter->x());
       data_->set(iter->index(), 1  , iter->y());
       data_->set(iter->index(), 2  , iter->z());
       for (int i=0; i<nDim_; ++i) 
       data_->set(iter->index(), 3+i, v->get(i));
  }
  isComputed_ = true;
  iter->rewind();
}

Field3D::Field3D(const int& ndim, const int& np, const double& pad, psi::SharedWavefunction wfn, psi::Options& opt)
 : pointsCollection_(PointsCollection3D::build(np, pad, wfn->molecule())),
   data_(std::make_shared<Matrix>("XYZ and Potential Values", np, 3+ndim)),
   geom_(psi::Matrix(wfn->molecule()->geometry())), 
   wfn_(wfn),
   nDim_(ndim)
{
   common_init();
}

Field3D::Field3D(const int& ndim, const int& nx, const int& ny, const int& nz,
                             const double& px, const double& py, const double& pz,
                             std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options)
  : pointsCollection_(PointsCollection3D::build(nx, ny, nz, px, py, pz, wfn->basisset(), options)),
    data_(std::make_shared<Matrix>("XYZ and Potential Values", nx*ny*nz, 3+ndim)),
    geom_(psi::Matrix(wfn->molecule()->geometry())),
    wfn_(wfn),
    nDim_(ndim)
{
    common_init();
}

void Field3D::common_init()
{
  primary_= wfn_->basisset();
  nbf_    = wfn_->basisset()->nbf();
  fact_   = std::make_shared<psi::IntegralFactory>(primary_);
  pot_    = std::make_shared<psi::Matrix>("POT", nbf_, nbf_);
  potInt_ = std::make_shared<PotentialInt>(fact_->spherical_transform(), primary_, primary_, 0);
  isComputed_ = false;
}

void Field3D::write_cube_file(const std::string& name)
{
   if (points_collection()->get_type() == PointsCollection3D::Cube) {
       compute();
       std::shared_ptr<CubePointsCollection3D> col = std::dynamic_pointer_cast<CubePointsCollection3D>(points_collection());
       for (int i=0; i<nDim_; ++i) {
            std::string name_i = oepdev::string_sprintf("%s.%d", name.c_str(), i);
	    col->write_cube_file(data(), name_i, i);
       }
   }
}


Field3D::~Field3D()
{

}

ElectrostaticPotential3D::ElectrostaticPotential3D(const int& npoints, const double& padding, psi::SharedWavefunction wfn, psi::Options& options) 
 : Field3D(1, npoints, padding, wfn, options)
{

}


ElectrostaticPotential3D::ElectrostaticPotential3D(const int& nx, const int& ny, const int& nz,
                                                   const double& px, const double& py, const double& pz,
                                                   psi::SharedWavefunction wfn, psi::Options& options)
 : Field3D(1, nx, ny, nz, px, py, pz, wfn, options)
{ 

}

ElectrostaticPotential3D::~ElectrostaticPotential3D()
{

}

void ElectrostaticPotential3D::print() const 
{
   cout << " I am Electrostatic potential:\n";
   pointsCollection_->print();
}

std::shared_ptr<psi::Vector> ElectrostaticPotential3D::compute_xyz(const double& x, const double& y, const double& z)
{
   double val = 0.0;
   // ===> Nuclear contribution <=== //
   for (int i=0; i<wfn_->molecule()->natom(); ++i) {
        val+= (double)wfn_->molecule()->Z(i) /
                  sqrt( pow(geom_.get(i,0)-x, 2.0) 
                      + pow(geom_.get(i,1)-y, 2.0) 
                      + pow(geom_.get(i,2)-z, 2.0) );
   }

   // ===> Electronic contribution <=== //
   {double v;
   potInt_->set_charge_field(x, y, z);
   oneInt_ = potInt_;
   oneInt_->compute(pot_);
   for (int i=0; i<nbf_; ++i) {
        for (int j=0; j<=i; ++j) {
             v= (wfn_->Da()->get(i,j) + wfn_->Db()->get(i,j)) * pot_->get(i,j);
             val += 2.0*v;
             if (i==j) val -= v;
        }
   }
   pot_->zero();}
   
   // save
   std::shared_ptr<psi::Vector> v = std::make_shared<psi::Vector>("", 1);
   v->set(0, val);
   return v;
}



//
// Factory methods
// 

std::shared_ptr<Points3DIterator> Points3DIterator::build(const int& nx, const int& ny, const int& nz,
                                                     const double& dx, const double& dy, const double& dz,
                                                     const double& ox, const double& oy, const double& oz)
{
  std::shared_ptr<Points3DIterator> iterator = std::make_shared<CubePoints3DIterator>(nx, ny, nz, dx, dy, dz, ox, oy, oz);
  return iterator;
}

std::shared_ptr<Points3DIterator> Points3DIterator::build(const int& np, 
                                                     const double& cx, const double& cy, const double& cz,
                                                     const double& radius)
{
 std::shared_ptr<Points3DIterator> iterator = std::make_shared<RandomPoints3DIterator>(np, cx, cy, cz, radius);
 return iterator;
}

shared_ptr<Points3DIterator> Points3DIterator::build(const int& np, const double& pad, psi::SharedMolecule mol)
{
   std::shared_ptr<Points3DIterator> iterator = std::make_shared<RandomPoints3DIterator>(np, pad, mol);
   return iterator;
}


std::shared_ptr<PointsCollection3D> PointsCollection3D::build(const int& npoints, const double& radius,
           const double& cx, const double& cy, const double& cz)
{
  std::shared_ptr<PointsCollection3D> pointsCollection = 
              std::make_shared<RandomPointsCollection3D>(PointsCollection3D::Random, npoints, radius, cx, cy, cz);
  return pointsCollection;
}

std::shared_ptr<PointsCollection3D> PointsCollection3D::build(const int& npoints, const double& padding, psi::SharedMolecule mol)
{
  std::shared_ptr<PointsCollection3D> pointsCollection = std::make_shared<RandomPointsCollection3D>(PointsCollection3D::Random, npoints, padding, mol);
  return pointsCollection;
}


std::shared_ptr<PointsCollection3D> PointsCollection3D::build(const int& nx, const int& ny, const int& nz,
                                                         const double& px, const double& py, const double& pz,
                                                         psi::SharedBasisSet bs, psi::Options& options){
  shared_ptr<PointsCollection3D> pointsCollection = std::make_shared<CubePointsCollection3D>(PointsCollection3D::Cube, nx, ny, nz, px, py, pz, bs, options);
  return pointsCollection;
}

// Build Random with molecule (vdW radii include)
std::shared_ptr<Field3D> Field3D::build(const std::string& type, 
                                               const int& np,
                                               const double& padding,
                                               psi::SharedWavefunction wfn,
                                               psi::Options& options, const int& ndim)
// ndim is not used at that moment
{
   std::shared_ptr<Field3D> field;
   if          (type == "ELECTROSTATIC POTENTIAL") field = std::make_shared<ElectrostaticPotential3D>(np, padding, wfn, options);
   else throw psi::PSIEXCEPTION(" ERROR! Incorrect scalar field type requested!\n");
   return field;
}

// Build CUBE collection with wavefunction
std::shared_ptr<Field3D> Field3D::build(const std::string& type, 
                                        const int& nx, const int& ny, const int& nz,
                                        const double& px, const double& py, const double& pz,
                                        psi::SharedWavefunction wfn, psi::Options& options, const int& ndim)
// ndim is not used at that moment
{
   std::shared_ptr<Field3D> field;
   if          (type == "ELECTROSTATIC POTENTIAL") field = std::make_shared<ElectrostaticPotential3D>(nx, ny, nz, px, py, pz, wfn, options);
   else throw psi::PSIEXCEPTION(" ERROR! Incorrect scalar field type requested!\n");
   return field;
}




///// BarePotential3D
//BarePotential3D::BarePotential3D(SharedMatrix pot)
// : Potential3D(pot)
//{
//  common_init();
//}
//
//void BarePotential3D::common_init()
//{
//  
//}
//
//void BarePotential3D::compute() 
//{
//
//}

} // EndNameSpace oepdev
