#include <iostream>
#include <random>
#include "gefp.h"

using namespace std;


void oepdev::GenEffPar::allocate_dipole_polarizability(int nsites, int nbf)
{
   // Here also the distributed centres are allocated.
   std::map<int, char> m;
   m[0] = 'X'; m[1] = 'Y'; m[2] = 'Z';

   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susc;
   std::vector<std::shared_ptr<psi::Vector>> centres;
   for (int n=0; n<nsites; ++n) {
        std::vector<std::shared_ptr<psi::Matrix>> susc_n;
        centres.push_back(std::make_shared<psi::Vector>(oepdev::string_sprintf("Centre (%d)", n+1), 3));
        for (int z=0; z<3; ++z) {
             std::string name = oepdev::string_sprintf("Density Matrix Dipole Polarizability B[%c](%d)", m[z], n+1);
             std::shared_ptr<psi::Matrix> susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
             susc_n.push_back(susc_nz);
        }
        susc.push_back(susc_n);
   }
   set_dipole_polarizability(susc);
   set_centres(centres);
}
void oepdev::GenEffPar::allocate_dipole_dipole_hyperpolarizability(int nsites, int nbf)
{
   std::map<int, char> m;
   m[0] = 'X'; m[1] = 'Y'; m[2] = 'Z';

   //mm[0] = "XX"; mm[1] = "XY"; mm[2] = "XZ";
   //mm[3] = "YX"; mm[4] = "YY"; mm[5] = "YZ";
   //mm[6] = "ZX"; mm[7] = "ZY"; mm[8] = "ZZ";

   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susc;
   for (int n=0; n<nsites; ++n) {
        std::vector<std::shared_ptr<psi::Matrix>> susc_n;
        for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             std::string name = oepdev::string_sprintf("Density Matrix Dipole-Dipole Hyperpolarizability B[%c%c](%d)", m[z1], m[z2], n+1);
             std::shared_ptr<psi::Matrix> susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
             susc_n.push_back(susc_nz);
        }}
        susc.push_back(susc_n);
   }
   set_dipole_dipole_hyperpolarizability(susc);
}
void oepdev::GenEffPar::allocate_quadrupole_polarizability(int nsites, int nbf)
{
   std::map<int, char> m;
   m[0] = 'X'; m[1] = 'Y'; m[2] = 'Z';

   std::vector<std::vector<std::shared_ptr<psi::Matrix>>> susc;
   for (int n=0; n<nsites; ++n) {
        std::vector<std::shared_ptr<psi::Matrix>> susc_n;
        for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             std::string name = oepdev::string_sprintf("Density Matrix Quadrupole Polarizability B[%c%c](%d)", m[z1], m[z2], n+1);
             std::shared_ptr<psi::Matrix> susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
             susc_n.push_back(susc_nz);
        }}
        susc.push_back(susc_n);
   }
   set_quadrupole_polarizability(susc);
}
std::shared_ptr<psi::Matrix> oepdev::GenEffPar::compute_density_matrix(std::shared_ptr<psi::Vector> field)
{
   return oepdev::GenEffPar::compute_density_matrix(field->get(0), field->get(1), field->get(2));
}
std::shared_ptr<psi::Matrix> oepdev::GenEffPar::compute_density_matrix(double fx, double fy, double fz)
{
   int nsites = densityMatrixDipolePolarizability_.size();
   std::vector<std::shared_ptr<psi::Vector>> fields;
   for (int n=0; n<nsites; ++n) {
        std::shared_ptr<psi::Vector> field = std::make_shared<psi::Vector>("", 3); 
        field->set(0, fx);
        field->set(1, fy);
        field->set(2, fz);
        fields.push_back(field);
   }
   return oepdev::GenEffPar::compute_density_matrix(fields);
}
std::shared_ptr<psi::Matrix> oepdev::GenEffPar::compute_density_matrix(std::vector<std::shared_ptr<psi::Vector>> fields) 
{
   if (!hasDensityMatrixDipolePolarizability_) throw psi::PSIEXCEPTION("Density Matrix Dipole Polarizability is not set!");
   int nbf = densityMatrixDipolePolarizability_[0][0]->nrow();
   int nsites = densityMatrixDipolePolarizability_.size();

   std::shared_ptr<psi::Matrix> D = std::make_shared<psi::Matrix>("Density Matrix Change", nbf, nbf);
   for (int n=0; n<nsites; ++n) {
        double fx = fields[n]->get(0);
        double fy = fields[n]->get(1);
        double fz = fields[n]->get(2);
        D->axpy(fx, densityMatrixDipolePolarizability_[n][0]);
        D->axpy(fy, densityMatrixDipolePolarizability_[n][1]);
        D->axpy(fz, densityMatrixDipolePolarizability_[n][2]);
        //cout << densityMatrixDipolePolarizability_[n][0]->get(0,0) << endl;
        //cout << densityMatrixDipolePolarizability_[n][1]->get(0,0) << endl;
        //cout << densityMatrixDipolePolarizability_[n][2]->get(0,0) << endl;
        if (hasDensityMatrixDipoleDipoleHyperpolarizability_) {
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][0]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][4]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][8]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][1]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][3]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][2]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][6]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][5]->get(0,0) << endl;
            //cout << densityMatrixDipoleDipoleHyperpolarizability_[n][7]->get(0,0) << endl;
            D->axpy(fx*fx, densityMatrixDipoleDipoleHyperpolarizability_[n][0]);
            D->axpy(fx*fy, densityMatrixDipoleDipoleHyperpolarizability_[n][1]);
            D->axpy(fx*fz, densityMatrixDipoleDipoleHyperpolarizability_[n][2]);
            D->axpy(fy*fx, densityMatrixDipoleDipoleHyperpolarizability_[n][3]);
            D->axpy(fy*fy, densityMatrixDipoleDipoleHyperpolarizability_[n][4]);
            D->axpy(fy*fz, densityMatrixDipoleDipoleHyperpolarizability_[n][5]);
            D->axpy(fz*fx, densityMatrixDipoleDipoleHyperpolarizability_[n][6]);
            D->axpy(fz*fy, densityMatrixDipoleDipoleHyperpolarizability_[n][7]);
            D->axpy(fz*fz, densityMatrixDipoleDipoleHyperpolarizability_[n][8]);
        }
   }
   return D;
}
std::shared_ptr<psi::Matrix> oepdev::GenEffPar::compute_density_matrix(std::vector<std::shared_ptr<psi::Vector>> fields,
                                                                       std::vector<std::shared_ptr<psi::Matrix>> grads) 
{
   std::shared_ptr<psi::Matrix> D = oepdev::GenEffPar::compute_density_matrix(fields);
   if (hasDensityMatrixQuadrupolePolarizability_) {
       int nsites = densityMatrixDipolePolarizability_.size();
       for (int n=0; n<nsites; ++n) {                                   
            double fxx = grads[n]->get(0, 0);
            double fxy = grads[n]->get(0, 1);
            double fxz = grads[n]->get(0, 2);
            double fyx = grads[n]->get(1, 0);
            double fyy = grads[n]->get(1, 1);
            double fyz = grads[n]->get(1, 2);
            double fzx = grads[n]->get(2, 0);
            double fzy = grads[n]->get(2, 1);
            double fzz = grads[n]->get(2, 2);
            D->axpy(fxx, densityMatrixQuadrupolePolarizability_[n][0]);
            D->axpy(fxy, densityMatrixQuadrupolePolarizability_[n][1]);
            D->axpy(fxz, densityMatrixQuadrupolePolarizability_[n][2]);
            D->axpy(fyx, densityMatrixQuadrupolePolarizability_[n][3]);
            D->axpy(fyy, densityMatrixQuadrupolePolarizability_[n][4]);
            D->axpy(fyz, densityMatrixQuadrupolePolarizability_[n][5]);
            D->axpy(fzx, densityMatrixQuadrupolePolarizability_[n][6]);
            D->axpy(fzy, densityMatrixQuadrupolePolarizability_[n][7]);
            D->axpy(fzz, densityMatrixQuadrupolePolarizability_[n][8]);
       }
   }
   return D;
}
//-- GenEffPar --///////////////////////////////////////////////////////////////////////////////////////
oepdev::GenEffFrag::GenEffFrag(std::string name) : 
  name_(name),
  densityMatrixSusceptibilityGEF_(nullptr),
  electrostaticEnergyGEF_(nullptr),
  repulsionEnergyGEF_(nullptr),
  chargeTransferEnergyGEF_(nullptr),
  EETCouplingConstantGEF_(nullptr)
{
  parameters["POLARIZATION"   ] = densityMatrixSusceptibilityGEF_;
  parameters["COULOMBIC   "   ] = electrostaticEnergyGEF_;
  parameters["REPULSION"      ] = repulsionEnergyGEF_;
  parameters["CHARGE_TRANSFER"] = chargeTransferEnergyGEF_;
  parameters["EET_COUPLING"   ] = EETCouplingConstantGEF_;
}
oepdev::GenEffFrag::~GenEffFrag() {}
void oepdev::GenEffFrag::rotate(std::shared_ptr<psi::Matrix> R)
{
  throw psi::NOT_IMPLEMENTED_EXCEPTION();
}
void oepdev::GenEffFrag::translate(std::shared_ptr<psi::Vector> T)
{
  throw psi::NOT_IMPLEMENTED_EXCEPTION();
}
void oepdev::GenEffFrag::superimpose(std::shared_ptr<psi::Matrix> targetXYZ, std::vector<int> supList)
{
  throw psi::NOT_IMPLEMENTED_EXCEPTION();
}
//-- GenEffParFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::GenEffParFactory::GenEffParFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 wfn_(wfn),
 options_(opt),
 nbf_(wfn->basisset()->nbf()),
 randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0)),
 cphfSolver_(nullptr)
{
   // --> Future: encapsulate the ESP point space with vdW excluded spheres into an object and put to 3D points library <-- //

   // Padding
   double pad = options_.get_double("ESP_PAD_SPHERE");

   // Populate vdwRadius
   vdwRadius_["C"] = options_.get_double("ESP_VDW_RADIUS_C" );
   vdwRadius_["H"] = options_.get_double("ESP_VDW_RADIUS_H" );
   vdwRadius_["N"] = options_.get_double("ESP_VDW_RADIUS_N" );
   vdwRadius_["O"] = options_.get_double("ESP_VDW_RADIUS_O" );
   vdwRadius_["F"] = options_.get_double("ESP_VDW_RADIUS_F" );
   vdwRadius_["Cl"]= options_.get_double("ESP_VDW_RADIUS_CL");

   // Set the centre of the sphere
   psi::Matrix  geom = wfn_->molecule()->geometry();
   psi::Vector3  com = wfn_->molecule()->center_of_mass();
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
   excludeSpheres_ = std::make_shared<psi::Matrix>("van der Waals Exclude Spheres", geom.nrow(), 4);
   for (int i=0; i< geom.nrow(); ++i) {
        for (int j=0; j< geom.ncol(); ++j) {
             excludeSpheres_->set(i, j, geom.get(i,j));
        }
        excludeSpheres_->set(i, 3, vdwRadius_[wfn_->molecule()->symbol(i)]);
   }
}
oepdev::GenEffParFactory::~GenEffParFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::GenEffParFactory::compute()
{

}
bool oepdev::GenEffParFactory::is_in_vdWsphere(double x, double y, double z) const
{ 
   // --> Encapsulate this method in a separate class and put it to 3D points library <-- //
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

std::shared_ptr<psi::Vector> oepdev::GenEffParFactory::draw_random_point() 
{
  double theta     = acos(random_double());
  double phi       = 2.0 * M_PI * random_double();
  double r         = radius_  * cbrt( random_double() );
  double x         = cx_ + r * sin(theta) * cos(phi);
  double y         = cy_ + r * sin(theta) * sin(phi);
  double z         = cz_ + r * cos(theta);

  while (is_in_vdWsphere(x, y, z) == true) {               
         theta     = acos(random_double());              
         phi       = 2.0 * M_PI * random_double();
         r         = radius_  * cbrt( random_double() );
         x         = cx_ + r * sin(theta) * cos(phi); 
         y         = cy_ + r * sin(theta) * sin(phi);
         z         = cz_ + r * cos(theta);
  }
  std::shared_ptr<psi::Vector> point = std::make_shared<psi::Vector>("",3);
  point->set(0, x);
  point->set(1, y);
  point->set(2, z);
  return point;
}
// Static factory method
std::shared_ptr<oepdev::GenEffParFactory> oepdev::GenEffParFactory::build(const std::string& type, 
                                                  std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
{
   if (type == "POLARIZATION") {
       const std::string mode  = opt.get_str("DMATPOL_TRAINING_MODE");
       const int rank_field    = opt.get_int("DMATPOL_FIELD_RANK");
       const int rank_gradient = opt.get_int("DMATPOL_GRADIENT_RANK");

       if (mode == "NONE") {
          {return std::make_shared<oepdev::AbInitioPolarGEFactory>(wfn, opt);}
       } else {
        std::string notsupported = "Susceptibilities for this rank are not supported yet.";                                          
        if        (rank_field == 0) { 
          throw psi::PSIEXCEPTION("Trivially vanishing susceptibility!");
        } else if (rank_field == 1) { // Linear wrt electric field
          if      (rank_gradient == 0) {
            if      (mode == "EFIELD")  {return std::make_shared<oepdev::LinearUniformEFieldPolarGEFactory>(wfn, opt);}
            else if (mode == "CHARGES") {return std::make_shared<oepdev::LinearNonUniformEFieldPolarGEFactory>(wfn, opt);}
            else {throw psi::PSIEXCEPTION(notsupported);}
          }
          else if (rank_gradient == 1) {
            if (mode == "EFIELD") throw psi::PSIEXCEPTION("Options: Gradient rank 1 and uniform EFIELD exclude each other.");
            else if (mode == "CHARGES") return std::make_shared<oepdev::LinearGradientNonUniformEFieldPolarGEFactory>(wfn, opt);
            else {throw psi::PSIEXCEPTION(notsupported);}
          }
          else {throw psi::PSIEXCEPTION(notsupported);}
        } else if (rank_field == 2) {  // Quadratic wrt electric field
          if      (rank_gradient == 0) {
            if      (mode == "EFIELD")  {return std::make_shared<oepdev::QuadraticUniformEFieldPolarGEFactory>(wfn, opt);}
            else if (mode == "CHARGES") {return std::make_shared<oepdev::QuadraticNonUniformEFieldPolarGEFactory>(wfn, opt);}
          }
          else if (rank_gradient == 1) {
            if (mode == "EFIELD") throw psi::PSIEXCEPTION("Options: Gradient rank 1 and uniform EFIELD exclude each other.");
            else if (mode == "CHARGES") return std::make_shared<oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory>(wfn, opt);
            else {throw psi::PSIEXCEPTION(notsupported);} 
          }
          else {throw psi::PSIEXCEPTION(notsupported);}
        } else {
             throw psi::PSIEXCEPTION(notsupported);
        }
       }
   } else {
     throw psi::PSIEXCEPTION("Invalid factory type chosen!");
   }
}
