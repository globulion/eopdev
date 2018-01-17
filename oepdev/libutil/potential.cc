#include "potential.h"


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

RandomPoints3DIterator::RandomPoints3DIterator(const int& np, const double& cx, const double& cy, const double& cz,
                                                              const double& radius)
 : Points3DIterator(np), cx_(cx), cy_(cy), cz_(cz), radius_(radius),
   randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   
}

RandomPoints3DIterator::RandomPoints3DIterator(const int& np, const double& pad, psi::SharedMolecule mol)
 : Points3DIterator(np), randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   // populate vdwRadius
   vdwRadius_["C"] = 3.000;
   vdwRadius_["H"] = 1.000;
   vdwRadius_["N"] = 2.400;
   vdwRadius_["O"] = 2.600;
   vdwRadius_["F"] = 2.300;
   vdwRadius_["Cl"]= 2.900;

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
   draw_random_point();
}

void RandomPoints3DIterator::next()
{
   ++index_;
   if (index_ == np_) done_ = true;
   draw_random_point();
}

Distribution3D::Distribution3D(Distribution distribution, int& np) : type_(distribution), np_(np) 
{

}
Distribution3D::Distribution3D(Distribution distribution, const int& np) : type_(distribution), np_(np) 
{

}


Distribution3D::~Distribution3D()
{

}

void Distribution3D::print() const 
{

}

RandomDistribution3D::RandomDistribution3D(Distribution distribution, 
                                           const int& npoints,                                   
                                           const double& cx, const double& cy, const double& cz,
                                           const double& radius)
 : Distribution3D(distribution, npoints)
{
  pointsIterator_ = make_shared<RandomPoints3DIterator>(npoints, cx, cy, cz, radius);
}

RandomDistribution3D::RandomDistribution3D(Distribution distribution, 
                                           const int& npoints,                                   
                                           const double& padding,
                                           psi::SharedMolecule mol)
 : Distribution3D(distribution, npoints)
{
  pointsIterator_ = make_shared<RandomPoints3DIterator>(npoints, padding, mol);
}


RandomDistribution3D::~RandomDistribution3D()
{

}

void RandomDistribution3D::print() const 
{
   cout << "I am a RANDOM distribution with np = " << np_ << "\n";
}

CubeDistribution3D::CubeDistribution3D(Distribution distribution, 
                                       const int& nx, const int& ny, const int& nz,
                                       const double& px, const double& py, const double& pz, 
                                       psi::SharedBasisSet bs, psi::Options& options)
 : Distribution3D(distribution, nx*ny*nz), CubicScalarGrid(bs, options)
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
       if (x <  xmin) xmin = x;
       if (y <  ymin) ymin = y;
       if (z <  zmin) zmin = z;
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
  double D[3] = {dx, dy, dz};
  double O[3] = {ox, oy, oz}; 
  build_grid(".", N, D, O);
}

CubeDistribution3D::~CubeDistribution3D()
{

}

void CubeDistribution3D::print() const {
  //psi::outfile->Printf(" ===> Cube 3D Collection <===\n");
}

void CubeDistribution3D::write_cube_file(psi::SharedMatrix v, const std::string& name){
    // => Drop the grid out <= //

    std::stringstream ss;
    ss << filepath_ << "/" << name << ".cube";

    // Is filepath a valid directory?
    std::ifstream infile(filepath_.c_str());
    if (!infile.good()) {
    //if (filesystem::path(filepath_).make_absolute().is_directory() == false) {
        printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath_.c_str());
        outfile->Printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath_.c_str());
        outfile->Flush();
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
        fprintf(fh, "%13.5E", v->get(ind,3));
        if (ind % 6 == 5) fprintf(fh,"\n");
    }

    fclose(fh);
}

void Potential3D::compute(){
  shared_ptr<Points3DIterator> iter = distribution()->points_iterator();
  if (!iter) throw psi::PSIEXCEPTION("ERROR!!! No iterator in Potential3D object!\n");
  for (iter->first(); iter->is_done()==false; iter->next()) {
       //double val = compute_xyz(iter->x(), iter->y(), iter->z());
       //cout << "Point: " << iter->x() << "  " << iter->y() << "  " << iter->z() << "\n";
       data_->set(iter->index(), 0, iter->x());
       data_->set(iter->index(), 1, iter->y());
       data_->set(iter->index(), 2, iter->z());
       data_->set(iter->index(), 3, compute_xyz(iter->x(), iter->y(), iter->z()));
  }
}

Potential3D::Potential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius) 
 : distribution_(Distribution3D::build(np, cx, cy, cz, radius)), 
   data_(std::make_shared<Matrix>("XYZ and Scalar Potential Values", np, 4))
{
   
}

Potential3D::Potential3D(const int& np, const double& pad, psi::SharedMolecule mol)
 : distribution_(Distribution3D::build(np, pad, mol)),
   data_(std::make_shared<Matrix>("XYZ and Scalar Potential Values", np, 4)),
   geom_(psi::Matrix(mol->geometry()))
{

}

Potential3D::Potential3D(const int& nx, const int& ny, const int& nz,
                         const double& px, const double& py, const double& pz,
                         std::shared_ptr<psi::Wavefunction> wfn, psi::Options& options)
  : distribution_(Distribution3D::build(nx, ny, nz, px, py, pz, wfn->basisset(), options)),
    data_(std::make_shared<Matrix>("XYZ and Scalar Potential Values", nx*ny*nz, 4)),
    wfn_(wfn), 
    nbf_(wfn->basisset()->nbf()), 
    geom_(psi::Matrix(wfn->molecule()->geometry())),
    primary_(wfn_->basisset()),
    fact_(std::make_shared<psi::IntegralFactory>(wfn->basisset(), wfn->basisset())),
    pot_(std::make_shared<psi::Matrix>("POT", wfn->basisset()->nbf(), wfn->basisset()->nbf()))
{
    potInt_ = std::make_shared<PotentialInt>(fact_->spherical_transform(), primary_, primary_, 0);
}

void Potential3D::write_cube_file(const std::string& name)
{
   if (distribution()->get_type() == Distribution3D::Cube) {
       compute();
       std::shared_ptr<CubeDistribution3D> distr = std::dynamic_pointer_cast<CubeDistribution3D>(distribution());
       distr->write_cube_file(data(), name);
   }
}


Potential3D::~Potential3D()
{

}

EPotential3D::EPotential3D(const int& np, const double& cx, const double& cy, const double& cz, const double& radius) 
 : Potential3D(np, cx, cy, cz, radius)
{

}

EPotential3D::EPotential3D(const int& npoints, const double& padding, psi::SharedMolecule molecule) 
 : Potential3D(npoints, padding, molecule)
{

}


EPotential3D::EPotential3D(const int& nx, const int& ny, const int& nz,
                           const double& px, const double& py, const double& pz,
                           psi::SharedWavefunction wfn, psi::Options& options)
 : Potential3D(nx, ny, nz, px, py, pz, wfn, options)
{ 

}

EPotential3D::~EPotential3D()
{

}

void EPotential3D::print() const 
{
   cout << " I am Electrostatic potential:\n";
   distribution_->print();
}

double EPotential3D::compute_xyz(const double& x, const double& y, const double& z)
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
   double v;
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
   pot_->zero();
   
   return val;
}



//
// Factory methods
// 

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

shared_ptr<Points3DIterator> Points3DIterator::build(const int& np, const double& pad, psi::SharedMolecule mol)
{
   shared_ptr<Points3DIterator> iterator = make_shared<RandomPoints3DIterator>(np, pad, mol);
   return iterator;
}


shared_ptr<Distribution3D> Distribution3D::build(const int& npoints, 
           const double& cx, const double& cy, const double& cz, const double& radius)
{
  shared_ptr<Distribution3D> distribution = make_shared<RandomDistribution3D>(Distribution3D::Random,npoints,cx,cy,cz,radius);
  return distribution;
}

shared_ptr<Distribution3D> Distribution3D::build(const int& npoints, const double& padding, psi::SharedMolecule mol)
{
  shared_ptr<Distribution3D> distribution = make_shared<RandomDistribution3D>(Distribution3D::Random, npoints, padding, mol);
  return distribution;
}


shared_ptr<Distribution3D> Distribution3D::build(const int& nx, const int& ny, const int& nz,
                                                 const double& px, const double& py, const double& pz,
                                                 psi::SharedBasisSet bs, psi::Options& options){
  shared_ptr<Distribution3D> distribution = make_shared<CubeDistribution3D>(Distribution3D::Cube, nx, ny, nz, px, py, pz, bs, options);
  return distribution;
}

// Build Random without molecule (no vdW radii)
shared_ptr<Potential3D> Potential3D::build(const std::string& type, 
                                           const int& np,
                                           const double& cx, const double& cy, const double& cz,
                                           const double& radius)
{
   shared_ptr<Potential3D> potential;
   if          (type == "ELECTROSTATIC") potential = make_shared<EPotential3D>(np, cx, cy, cz, radius);
   else cout << " ERROR! No such potential type <" << type << "> !\n";
   return potential;
}

// Build Random with molecule (vdW radii include)
shared_ptr<Potential3D> Potential3D::build(const std::string& type, 
                                           const int& np,
                                           const double& padding,
                                           psi::SharedMolecule mol)
{
   shared_ptr<Potential3D> potential;
   if          (type == "ELECTROSTATIC") potential = make_shared<EPotential3D>(np, padding, mol);
   else cout << " ERROR! No such potential type <" << type << "> !\n";
   return potential;
}

// Build CUBE collection with wavefunction
shared_ptr<Potential3D> Potential3D::build(const std::string& type, 
                                           const int& nx, const int& ny, const int& nz,
                                           const double& px, const double& py, const double& pz,
                                           psi::SharedWavefunction wfn, psi::Options& options)
{
   shared_ptr<Potential3D> potential;
   if          (type == "ELECTROSTATIC") potential = make_shared<EPotential3D>(nx, ny, nz, px, py, pz, wfn, options);
   else cout << " ERROR! No such potential type <" << type << "> !\n";
   return potential;
}





////////////////////////////////////////////////////////////////////////////////////////////////////
//Potential3D::Potential3D(SharedMatrix pot) : pot_(std::make_shared<Matrix>(pot)) 
//{
//  common_init();
//}
//
//Potential3D::Potential3D(SharedWavefunction wfn) 
//{
//  common_init();
//}
//
//Potential3D::Potential3D(SharedOEPotential oep, const std::string& oepType)
//{
//  common_init();
//}
//
//void Potential3D::common_init()
//{
//
//}
//
//Potential3D::~Potential3D() 
//{
//
//}
//
//std::shared_ptr<Potential3D> Potential3D::build(SharedMatrix pot)
//{
//   std::shared_ptr<Potential3D> potential = std::make_shared<BarePotential3D>(pot);
//   return potential;
//}
//std::shared_ptr<Potential3D> Potential3D::build(const std::string& type, SharedWavefunction wfn)
//{
//   std::shared_ptr<Potential3D> potential;
//   if     (type == "ELECTROSTATIC") potential = std::make_shared<ElectrostaticPotential3D>(wfn);
//   else throw PSIEXCEPTION("OEPDEV: Potential3D init. Unknown type of potential specified.");
//   return potential;
//}
//std::shared_ptr<Potential3D> Potential3D::build(SharedOEPotential oep, const std::string& oepType)
//{
//   std::shared_ptr<Potential3D> potential = std::make_shared<OEPotential3D>(oep, oepType);
//   return potential;
//}
//
//void Potential3D::compute(bool cube) 
//{
//
//}
//
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
//
//
///// ElectrostaticPotential3D
//ElectrostaticPotential3D::ElectrostaticPotential3D(SharedWavefunction wfn)
// : Potential3D(wfn)
//{
//  common_init();
//}
//
//void ElectrostaticPotential3D::common_init()
//{
//
//}
//
//void ElectrostaticPotential3D::compute() 
//{
//   
//   // ===> Nuclear Contribution <=== //
//}
//
//
///// OEPotential3D
//OEPotential3D::OEPotential3D(SharedOEPotential oep, const std::string& oepType) 
// : Potential3D(oep, oepType)
//{
//  common_init();
//}
//
//void OEPotential3D::common_init()
//{
//
//}
//
//void OEPotential3D::compute() 
//{
//
//}


} // EndNameSpace oepdev
