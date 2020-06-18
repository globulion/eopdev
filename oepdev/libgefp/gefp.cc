#include <iostream>
#include <random>
#include "gefp.h"
#include "../libutil/kabsch_superimposer.h"
#include "../libutil/basis_rotation.h"

using namespace std;

oepdev::GenEffPar::GenEffPar(const GenEffPar* f) {
  name_ = f->name_;
  type_ = f->type_;
  hasDensityMatrixDipolePolarizability_ = f->hasDensityMatrixDipolePolarizability_;
  hasDensityMatrixDipoleDipoleHyperpolarizability_ = f->hasDensityMatrixDipoleDipoleHyperpolarizability_;
  hasDensityMatrixQuadrupolePolarizability_ = f->hasDensityMatrixQuadrupolePolarizability_;
  copy_from(f);
}
void oepdev::GenEffPar::copy_from(const GenEffPar* f) {
  //
  distributedCentres_.clear();
  for (unsigned int i=0; i<f->distributedCentres_.size(); ++i) {
       psi::SharedVector centre = std::make_shared<psi::Vector>(*f->distributedCentres_[i]);
       distributedCentres_.push_back(centre);
  }
  //
  data_matrix_.clear();
  for (auto const& x : f->data_matrix_) {
     std::string key = x.first;
     psi::SharedMatrix mat = std::make_shared<psi::Matrix>(x.second);
     data_matrix_[key] = mat;
  }
  //
  data_vector_.clear();
  for (auto const& x : f->data_vector_) {
     std::string key = x.first;
     psi::SharedVector vec = std::make_shared<psi::Vector>(*(x.second));
     data_vector_[key] = vec;
  }
  //
  data_dpol_.clear();
  for (auto const&x : f->data_dpol_) {
       std::string key = x.first;
       std::vector<psi::SharedMatrix> v;
       for (unsigned int i=0; i<x.second.size(); ++i) {
            psi::SharedMatrix m = std::make_shared<psi::Matrix>(x.second[i]);
            v.push_back(m);
       }
       data_dpol_[key] = v;
  }
  //
  data_dmtp_.clear();
  for (auto const& x : f->data_dmtp_) {
     std::string key = x.first;
     oepdev::SharedDMTPole dmtp = x.second->clone();
     data_dmtp_[key] = dmtp;
  }
  //
  data_oep_.clear();
  for (auto const& x : f->data_oep_) {
     std::string key = x.first;
     oepdev::SharedOEPotential oep = x.second->clone();
     data_oep_[key] = oep;
  }
  //
  densityMatrixDipolePolarizability_.clear();
  for (unsigned int i=0; i<f->densityMatrixDipolePolarizability_.size(); ++i) {
       std::vector<psi::SharedMatrix> v;
       for (unsigned int a=0; a<f->densityMatrixDipolePolarizability_[i].size(); ++a) {
            psi::SharedMatrix mat = std::make_shared<psi::Matrix>(f->densityMatrixDipolePolarizability_[i][a]);
            v.push_back(mat);
       }
       densityMatrixDipolePolarizability_.push_back(v);
  }
  //
  densityMatrixDipoleDipoleHyperpolarizability_.clear();
  for (unsigned int i=0; i<f->densityMatrixDipoleDipoleHyperpolarizability_.size(); ++i) {
       std::vector<psi::SharedMatrix> v;
       for (unsigned int a=0; a<f->densityMatrixDipoleDipoleHyperpolarizability_[i].size(); ++a) {
            psi::SharedMatrix mat = std::make_shared<psi::Matrix>(f->densityMatrixDipoleDipoleHyperpolarizability_[i][a]);
            v.push_back(mat);
       }
       densityMatrixDipoleDipoleHyperpolarizability_.push_back(v);
  }
  //
  densityMatrixQuadrupolePolarizability_.clear();
  for (unsigned int i=0; i<f->densityMatrixQuadrupolePolarizability_.size(); ++i) {
       std::vector<psi::SharedMatrix> v;
       for (unsigned int a=0; a<f->densityMatrixQuadrupolePolarizability_[i].size(); ++a) {
            psi::SharedMatrix mat = std::make_shared<psi::Matrix>(f->densityMatrixQuadrupolePolarizability_[i][a]);
            v.push_back(mat);
       }
       densityMatrixQuadrupolePolarizability_.push_back(v);
  }
}

void oepdev::GenEffPar::allocate_dipole_polarizability(int nsites, int nbf)
{
   // Here also the distributed centres are allocated.
   std::map<int, char> m;
   m[0] = 'X'; m[1] = 'Y'; m[2] = 'Z';

   std::vector<std::vector<psi::SharedMatrix>> susc;
   std::vector<psi::SharedVector> centres;
   for (int n=0; n<nsites; ++n) {
        std::vector<psi::SharedMatrix> susc_n;
        centres.push_back(std::make_shared<psi::Vector>(oepdev::string_sprintf("Centre (%d)", n+1), 3));
        for (int z=0; z<3; ++z) {
             std::string name = oepdev::string_sprintf("Density Matrix Dipole Polarizability B[%c](%d)", m[z], n+1);
             psi::SharedMatrix susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
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

   std::vector<std::vector<psi::SharedMatrix>> susc;
   for (int n=0; n<nsites; ++n) {
        std::vector<psi::SharedMatrix> susc_n;
        for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             std::string name = oepdev::string_sprintf("Density Matrix Dipole-Dipole Hyperpolarizability B[%c%c](%d)", m[z1], m[z2], n+1);
             psi::SharedMatrix susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
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

   std::vector<std::vector<psi::SharedMatrix>> susc;
   for (int n=0; n<nsites; ++n) {
        std::vector<psi::SharedMatrix> susc_n;
        for (int z1=0; z1<3; ++z1) {
        for (int z2=0; z2<3; ++z2) {
             std::string name = oepdev::string_sprintf("Density Matrix Quadrupole Polarizability B[%c%c](%d)", m[z1], m[z2], n+1);
             psi::SharedMatrix susc_nz = std::make_shared<psi::Matrix>(name, nbf, nbf);
             susc_n.push_back(susc_nz);
        }}
        susc.push_back(susc_n);
   }
   set_quadrupole_polarizability(susc);
}
psi::SharedMatrix oepdev::GenEffPar::compute_density_matrix(psi::SharedVector field)
{
   return oepdev::GenEffPar::compute_density_matrix(field->get(0), field->get(1), field->get(2));
}
psi::SharedMatrix oepdev::GenEffPar::compute_density_matrix(double fx, double fy, double fz)
{
   int nsites = densityMatrixDipolePolarizability_.size();
   std::vector<psi::SharedVector> fields;
   for (int n=0; n<nsites; ++n) {
        psi::SharedVector field = std::make_shared<psi::Vector>("", 3); 
        field->set(0, fx);
        field->set(1, fy);
        field->set(2, fz);
        fields.push_back(field);
   }
   return oepdev::GenEffPar::compute_density_matrix(fields);
}
psi::SharedMatrix oepdev::GenEffPar::compute_density_matrix(std::vector<psi::SharedVector> fields) 
{
   if (!hasDensityMatrixDipolePolarizability_) throw psi::PSIEXCEPTION("Density Matrix Dipole Polarizability is not set!");
   int nbf = densityMatrixDipolePolarizability_[0][0]->nrow();
   int nsites = densityMatrixDipolePolarizability_.size();

   psi::SharedMatrix D = std::make_shared<psi::Matrix>("Density Matrix Change", nbf, nbf);
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
        //densityMatrixDipolePolarizability_[n][0]->print();
        //densityMatrixDipolePolarizability_[n][1]->print();
        //densityMatrixDipolePolarizability_[n][2]->print();
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
            //densityMatrixDipoleDipoleHyperpolarizability_[n][0]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][1]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][2]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][3]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][4]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][5]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][6]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][7]->print();
            //densityMatrixDipoleDipoleHyperpolarizability_[n][8]->print();
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
psi::SharedMatrix oepdev::GenEffPar::compute_density_matrix(std::vector<psi::SharedVector> fields,
                                                                       std::vector<psi::SharedMatrix> grads) 
{
   psi::SharedMatrix D = oepdev::GenEffPar::compute_density_matrix(fields);
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
void oepdev::GenEffPar::rotate(psi::SharedMatrix R) {
  // Nothing to do here
  throw psi::PSIEXCEPTION("  OepDEV Error: Rotation feature not implemented in GenEffPar instance!\n");
}
void oepdev::GenEffPar::translate(psi::SharedVector t) {
  // Nothing to do here
  throw psi::PSIEXCEPTION("  OepDEV Error: Translation feature not implemented in GenEffPar instance!\n");
}
void oepdev::GenEffPar::superimpose(psi::SharedMatrix targetXYZ, std::vector<int> supList) {
   // Initialize Kabsch Solver
   KabschSuperimposer sup = KabschSuperimposer();

   // Determine the overlapping structure slices
   psi::SharedMatrix initial_xyz, final_xyz;
   if (supList.empty()) {
       initial_xyz = this->data_matrix_.at("pos");
       final_xyz = targetXYZ;
   } else {
       const int n = supList.size();
       initial_xyz = std::make_shared<psi::Matrix>("", n, 3);
       final_xyz   = std::make_shared<psi::Matrix>("", n, 3);
       for (int i=0; i<n; ++i) {
            initial_xyz->set_row(0, supList[i], this->data_matrix_.at("pos")->get_row(0, supList[i]));
            final_xyz  ->set_row(0, supList[i],                    targetXYZ->get_row(0, supList[i]));
       }
   }

   // Run Kabsch Algorithm
   sup.compute(initial_xyz, final_xyz);
   psi::SharedMatrix r = sup.rotation;
   psi::SharedVector t = sup.translation;
   const double tx = t->get(0);
   const double ty = t->get(1);
   const double tz = t->get(2);
   outfile->Printf("\n Kabsch Superimposition RMS = %14.5f [a.u.]\n", sup.rms());

  // Superimpose POS
  psi::SharedMatrix pos_new = psi::Matrix::doublet(this->data_matrix_.at("pos"), r, false, false);
  for (int i=0; i<pos_new->nrow(); ++i) {
       pos_new->set(i,0,pos_new->get(i,0)+tx);
       pos_new->set(i,1,pos_new->get(i,1)+ty);
       pos_new->set(i,2,pos_new->get(i,2)+tz);
  }
  this->data_matrix_["pos"] = pos_new;

  // Superimpose DMTP
  if (this->data_dmtp_.find("camm") != this->data_dmtp_.end()) {
      outfile->Printf("  Superimposing CAMM in Parameters %s\n", this->name_.c_str());
    //this->data_dmtp_.at("camm")->superimpose(targetXYZ, supList);
      this->data_dmtp_.at("camm")->rotate(r);
      this->data_dmtp_.at("camm")->translate(t);

  }
  // Superimpose LMOC
  if (this->data_matrix_.find("lmoc") != this->data_matrix_.end()) {
      outfile->Printf("  Superimposing LMOC in Parameters %s\n", this->name_.c_str());
      psi::SharedMatrix lmoc = psi::Matrix::doublet(this->data_matrix_.at("lmoc"), r, false, false);
      double** l = lmoc->pointer();
      for (int i=0; i<lmoc->nrow(); ++i) {
           l[i][0] += tx; 
           l[i][1] += ty; 
           l[i][2] += tz; 
      }
      this->data_matrix_["lmoc"] = lmoc;
  }
  // Superimpose LMOC-OEP
  if (this->data_matrix_.find("lmoc-oep") != this->data_matrix_.end()) {
      outfile->Printf("  Superimposing LMOC-OEP in Parameters %s\n", this->name_.c_str());
      psi::SharedMatrix lmoc = psi::Matrix::doublet(this->data_matrix_.at("lmoc-oep"), r, false, false);
      double** l = lmoc->pointer();
      for (int i=0; i<lmoc->nrow(); ++i) {
           l[i][0] += tx; 
           l[i][1] += ty; 
           l[i][2] += tz; 
      }
      this->data_matrix_["lmoc-oep"] = lmoc;
  }
  // Superimpose DPOL
  if (this->data_dpol_.find("0") != this->data_dpol_.end()) {
      outfile->Printf("  Superimposing DPOL in Parameters %s\n", this->name_.c_str());
      std::vector<psi::SharedMatrix> dpol_new;
      for (int i=0; i<this->data_dpol_.at("0").size(); ++i) {
           psi::SharedMatrix A = psi::Matrix::doublet(this->data_dpol_.at("0")[i], r, false, false);
           psi::SharedMatrix a = psi::Matrix::doublet(r, A, true, false);
           dpol_new.push_back(a);
      }
      this->data_dpol_["0"] = dpol_new;
  }

   // Compute AO rotation matrix
   psi::SharedMatrix R_primary, Ri_primary;
   if (this->data_basisset_.find("primary") != this->data_basisset_.end()) {
       R_primary = oepdev::ao_rotation_matrix(r, this->data_basisset_.at("primary"));
       R_primary->set_name("AO Rotation Matrix: Primary Basis");
     //R_primary->print();

       Ri_primary = R_primary->clone(); Ri_primary->invert(); Ri_primary->transpose_this();
   }

   // Superimpose orbitals
   if (this->data_matrix_.find("lmoo") != this->data_matrix_.end()) {
       outfile->Printf("  Superimposing LMOO in Parameters %s\n", this->name_.c_str());
       psi::SharedMatrix lmoo = psi::Matrix::doublet(Ri_primary, this->data_matrix_.at("lmoo"), true, false);
       this->data_matrix_["lmoo"] = lmoo;
   }

   if (this->data_matrix_.find("lmoo-oep") != this->data_matrix_.end()) {
       outfile->Printf("  Superimposing LMOO-OEP in Parameters %s\n", this->name_.c_str());
       psi::SharedMatrix lmoo = psi::Matrix::doublet(Ri_primary, this->data_matrix_.at("lmoo-oep"), true, false);
       this->data_matrix_["lmoo-oep"] = lmoo;
   }

   if (this->data_matrix_.find("cmoo") != this->data_matrix_.end()) {
       outfile->Printf("  Superimposing CMOO in Parameters %s\n", this->name_.c_str());
       psi::SharedMatrix cmoo = psi::Matrix::doublet(Ri_primary, this->data_matrix_.at("cmoo"), true, false);
       this->data_matrix_["cmoo"] = cmoo;
   }

   if (this->data_matrix_.find("cmov") != this->data_matrix_.end()) {
       outfile->Printf("  Superimposing CMOV in Parameters %s\n", this->name_.c_str());
       psi::SharedMatrix cmov = psi::Matrix::doublet(Ri_primary, this->data_matrix_.at("cmov"), true, false);
       this->data_matrix_["cmov"] = cmov;
   }

   if (this->data_matrix_.find("fock_ao") != this->data_matrix_.end()) {
       outfile->Printf("  Superimposing FOCK_AO in Parameters %s\n", this->name_.c_str());
       psi::SharedMatrix fock = psi::Matrix::triplet(R_primary, this->data_matrix_.at("fock_ao"), R_primary, true, false, false);
       this->data_matrix_["fock_ao"] = fock;
   }


   // --> Superimpose OEPs <-- //
   psi::SharedMatrix R_auxiliary;
   if (this->data_basisset_.find("auxiliary") != this->data_basisset_.end()) {
       R_auxiliary = oepdev::ao_rotation_matrix(r, this->data_basisset_.at("auxiliary"));
       R_auxiliary->set_name("AO Rotation Matrix: Auxiliary Basis");
   }

   // Superimpose OEP-REP
   if (this->data_oep_.find("rep") != this->data_oep_.end()) {
       outfile->Printf("  Superimposing OEP-REP in Parameters %s\n", this->name_.c_str());
       this->data_oep_.at("rep")->rotate(r, R_primary, R_auxiliary);
       this->data_oep_.at("rep")->translate(t);
   }

   // Superimpose OEP-CT
   if (this->data_oep_.find("ct") != this->data_oep_.end()) {
       outfile->Printf("  Superimposing OEP-CT in Parameters %s\n", this->name_.c_str());
       this->data_oep_.at("ct")->rotate(r, R_primary, R_auxiliary);
       this->data_oep_.at("ct")->translate(t);
   }


}

//-- GenEffParFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::GenEffParFactory::GenEffParFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 wfn_(wfn),
 options_(opt),
 nbf_(wfn->basisset()->nbf()),
 randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0)),
 cphfSolver_(nullptr),
 abInitioPolarizationSusceptibilitiesFactory_(nullptr)
{
   // --> Future: encapsulate the ESP point space with vdW excluded spheres into an object and put to 3D points library <-- //

   // Padding
   double pad = options_.get_double("ESP_PAD_SPHERE");

   // Populate vdwRadius
   vdwRadius_["C"] = options_.get_double("ESP_VDW_RADIUS_C" );
   vdwRadius_["H"] = options_.get_double("ESP_VDW_RADIUS_H" );
   vdwRadius_["HE"] = options_.get_double("ESP_VDW_RADIUS_HE" );
   vdwRadius_["NE"] = options_.get_double("ESP_VDW_RADIUS_NE" );
   vdwRadius_["AR"] = options_.get_double("ESP_VDW_RADIUS_AR" );
   vdwRadius_["N"] = options_.get_double("ESP_VDW_RADIUS_N" );
   vdwRadius_["O"] = options_.get_double("ESP_VDW_RADIUS_O" );
   vdwRadius_["F"] = options_.get_double("ESP_VDW_RADIUS_F" );
   vdwRadius_["CL"]= options_.get_double("ESP_VDW_RADIUS_CL");
   vdwRadius_["LI"]= options_.get_double("ESP_VDW_RADIUS_LI");

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
oepdev::SharedGenEffPar oepdev::GenEffParFactory::compute()
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

psi::SharedVector oepdev::GenEffParFactory::draw_random_point() 
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
  psi::SharedVector point = std::make_shared<psi::Vector>("",3);
  point->set(0, x);
  point->set(1, y);
  point->set(2, z);
  return point;
}
// Static factory method
oepdev::SharedGenEffParFactory oepdev::GenEffParFactory::build(const std::string& type,
                               psi::SharedWavefunction wfn, psi::Options& opt,
                               psi::SharedBasisSet auxiliary, psi::SharedBasisSet intermediate) {
   oepdev::SharedGenEffParFactory factory;

   if ((!auxiliary) && (!intermediate)) {
       factory = oepdev::GenEffParFactory::build(type, wfn, opt);
   } else {
     // when auxiliary and intermediate are given, OEP-based factory will be of interest
     if (type == "OEP-EFP2") {
         factory = std::make_shared<oepdev::OEP_EFP2_GEFactory>(wfn, opt, auxiliary, intermediate);
     } else {
         factory = oepdev::GenEffParFactory::build(type, wfn, opt);
     }
   }
   return factory;
}
oepdev::SharedGenEffParFactory oepdev::GenEffParFactory::build(const std::string& type, 
                               psi::SharedWavefunction wfn, psi::Options& opt)
{
   if (type == "POLARIZATION") {
       const std::string mode  = opt.get_str("DMATPOL_TRAINING_MODE");
       const int rank_field    = opt.get_int("DMATPOL_FIELD_RANK");
       const int rank_gradient = opt.get_int("DMATPOL_GRADIENT_RANK");

       // Check if gradient models are requested
     //if (rank_gradient > 0) 
     //    throw psi::PSIEXCEPTION("Gradient models are not available now and probably will be deprecated in the future.");

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
   } else if (type=="EFP2") {
     return std::make_shared<oepdev::EFP2_GEFactory>(wfn, opt);
   } else if (type=="OEP-EFP2") {
     return std::make_shared<oepdev::OEP_EFP2_GEFactory>(wfn, opt);
   } else {
     throw psi::PSIEXCEPTION("Invalid factory type chosen!");
   }
}
