#include <iostream>
#include <random>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

//-- FFAbInitioPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::FFAbInitioPolarGEFactory::FFAbInitioPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::PolarGEFactory(wfn, opt)
{
}
oepdev::FFAbInitioPolarGEFactory::~FFAbInitioPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::FFAbInitioPolarGEFactory::compute()
{
  // Sizing
  int nbf = wfn_->basisset()->nbf();

  // Allocate
  std::shared_ptr<MatrixFactory> matrixFactory = std::make_shared<psi::MatrixFactory>();
  matrixFactory->init_with(1, &nbf, &nbf);

  // FF step
  const double s = options_.get_double("DMATPOL_EFIELD_FF_STEP");

  // Prepare finite-difference displacements
  std::vector<std::shared_ptr<psi::Vector>> Fields;
  std::shared_ptr<psi::Vector> field_1 = std::make_shared<psi::Vector>("", 3); field_1->set(0, s);
  std::shared_ptr<psi::Vector> field_2 = std::make_shared<psi::Vector>("", 3); field_2->set(0,-s);
  std::shared_ptr<psi::Vector> field_3 = std::make_shared<psi::Vector>("", 3); field_3->set(1, s);
  std::shared_ptr<psi::Vector> field_4 = std::make_shared<psi::Vector>("", 3); field_4->set(1,-s);
  std::shared_ptr<psi::Vector> field_5 = std::make_shared<psi::Vector>("", 3); field_5->set(2, s);
  std::shared_ptr<psi::Vector> field_6 = std::make_shared<psi::Vector>("", 3); field_6->set(2,-s);

  Fields.push_back(field_1);
  Fields.push_back(field_2);
  Fields.push_back(field_3);
  Fields.push_back(field_4);
  Fields.push_back(field_5);
  Fields.push_back(field_6);

  std::vector<std::shared_ptr<psi::Matrix>> Dmatset;
  for (int n=0; n<6; ++n) {
       cout << " Calculations of FF displacement " << (n+1) << " out of 6 ...\n";
       Dmatset.push_back(perturbed_state(Fields[n])->Da());
  }
 
  // Compute The Density Matrix Dipole Susceptibility Tensors 
  const double h = 1.0/(2.0 * s);

  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipolePolarizability;
  std::vector<std::shared_ptr<psi::Matrix>> B;

  std::shared_ptr<psi::Matrix> Bx = matrixFactory->create_shared_matrix("B[x] Non-Orthogonal AO Basis"); 
  std::shared_ptr<psi::Matrix> By = matrixFactory->create_shared_matrix("B[y] Non-Orthogonal AO Basis"); 
  std::shared_ptr<psi::Matrix> Bz = matrixFactory->create_shared_matrix("B[z] Non-Orthogonal AO Basis"); 

  Bx->copy(Dmatset[0]); Bx->subtract(Dmatset[1]); Bx->scale(h);
  By->copy(Dmatset[2]); By->subtract(Dmatset[3]); By->scale(h);
  Bz->copy(Dmatset[4]); Bz->subtract(Dmatset[5]); Bz->scale(h);

  B.push_back(Bx);
  B.push_back(By);
  B.push_back(Bz);

  densityMatrixDipolePolarizability.push_back(B);
  cout << " Density Matrix Dipole Polarizability computed!\n";

  // Compute The Density Matrix Dipole-Dipole Hyperpolarizability Tensors 
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixDipoleDipoleHyperpolarizability;
  if ((options_.get_int("DMATPOL_FIELD_RANK") == 2) && (options_.get_bool("DMATPOL_FF_AB_INITIO"))) {
     const double t = 1.0 / (2.0*s*s); /* Taylor 1/2 coefficient is already included here */
     std::vector<std::shared_ptr<psi::Matrix>> C;

     std::vector<std::shared_ptr<psi::Vector>> SderFields;                                                                   
     std::shared_ptr<psi::Vector> sfield_1 = std::make_shared<psi::Vector>("", 3); sfield_1->set(0, s); sfield_1->set(1, s);
     std::shared_ptr<psi::Vector> sfield_2 = std::make_shared<psi::Vector>("", 3); sfield_2->set(0,-s); sfield_2->set(1,-s);
     std::shared_ptr<psi::Vector> sfield_3 = std::make_shared<psi::Vector>("", 3); sfield_3->set(0, s); sfield_3->set(2, s);
     std::shared_ptr<psi::Vector> sfield_4 = std::make_shared<psi::Vector>("", 3); sfield_4->set(0,-s); sfield_4->set(2,-s);
     std::shared_ptr<psi::Vector> sfield_5 = std::make_shared<psi::Vector>("", 3); sfield_5->set(1, s); sfield_5->set(2, s);
     std::shared_ptr<psi::Vector> sfield_6 = std::make_shared<psi::Vector>("", 3); sfield_6->set(1,-s); sfield_6->set(2,-s);
     std::shared_ptr<psi::Vector> sfield_7 = std::make_shared<psi::Vector>("", 3);
                                                                                                                             
     SderFields.push_back(sfield_1);
     SderFields.push_back(sfield_2);
     SderFields.push_back(sfield_3);
     SderFields.push_back(sfield_4);
     SderFields.push_back(sfield_5);
     SderFields.push_back(sfield_6);
     SderFields.push_back(sfield_7);
                                                                                                                             
     for (int n=0; n<7; ++n) {
          cout << " Calculations of FF displacement " << (6+n+1) << " out of 13 ...\n";
          Dmatset.push_back(perturbed_state(SderFields[n])->Da());
     }

     std::shared_ptr<psi::Matrix> Cxx = matrixFactory->create_shared_matrix("C[xx] Non-Orthogonal AO Basis");  
     std::shared_ptr<psi::Matrix> Cxy = matrixFactory->create_shared_matrix("C[xy] Non-Orthogonal AO Basis"); 
     std::shared_ptr<psi::Matrix> Cxz = matrixFactory->create_shared_matrix("C[xz] Non-Orthogonal AO Basis"); 
     std::shared_ptr<psi::Matrix> Cyx = matrixFactory->create_shared_matrix("C[yx] Non-Orthogonal AO Basis");  
     std::shared_ptr<psi::Matrix> Cyy = matrixFactory->create_shared_matrix("C[yy] Non-Orthogonal AO Basis"); 
     std::shared_ptr<psi::Matrix> Cyz = matrixFactory->create_shared_matrix("C[yz] Non-Orthogonal AO Basis"); 
     std::shared_ptr<psi::Matrix> Czx = matrixFactory->create_shared_matrix("C[zx] Non-Orthogonal AO Basis");  
     std::shared_ptr<psi::Matrix> Czy = matrixFactory->create_shared_matrix("C[zy] Non-Orthogonal AO Basis"); 
     std::shared_ptr<psi::Matrix> Czz = matrixFactory->create_shared_matrix("C[zz] Non-Orthogonal AO Basis"); 
                                                                                                       
     Cxx->axpy(-2.0, Dmatset[12]); Cxx->add(Dmatset[0]); Cxx->add(Dmatset[1]);  // XX
     Cyy->axpy(-2.0, Dmatset[12]); Cyy->add(Dmatset[2]); Cyy->add(Dmatset[3]);  // YY
     Czz->axpy(-2.0, Dmatset[12]); Czz->add(Dmatset[4]); Czz->add(Dmatset[5]);  // ZZ

     Cxy->axpy( 2.0, Dmatset[12]); Cxy->add(Dmatset[6]); Cxy->add(Dmatset[7]);  // XY
     Cxy->subtract(Dmatset[0]); Cxy->subtract(Dmatset[1]); Cxy->subtract(Dmatset[2]); Cxy->subtract(Dmatset[3]); 

     Cxz->axpy( 2.0, Dmatset[12]); Cxz->add(Dmatset[8]); Cxz->add(Dmatset[9]);  // XZ
     Cxz->subtract(Dmatset[0]); Cxz->subtract(Dmatset[1]); Cxz->subtract(Dmatset[4]); Cxz->subtract(Dmatset[5]); 

     Cyz->axpy( 2.0, Dmatset[12]); Cyz->add(Dmatset[10]); Cyz->add(Dmatset[11]);  // YZ
     Cyz->subtract(Dmatset[2]); Cyz->subtract(Dmatset[3]); Cyz->subtract(Dmatset[4]); Cyz->subtract(Dmatset[5]); 

     Cxx->scale(t); Cyy->scale(t); Czz->scale(t);
     Cxy->scale(t*0.5); Cxz->scale(t*0.5); Cyz->scale(t*0.5);
     Cyx->copy(Cxy); Czx->copy(Cxz); Czy->copy(Cyz);
                                                                                                      
     C.push_back(Cxx);
     C.push_back(Cxy);
     C.push_back(Cxz);
     C.push_back(Cyx);
     C.push_back(Cyy);
     C.push_back(Cyz);
     C.push_back(Czx);
     C.push_back(Czy);
     C.push_back(Czz);
     densityMatrixDipoleDipoleHyperpolarizability.push_back(C);
     cout << " Density Matrix Dipole-Dipole Hyperpolarizability computed!\n";
  }

  // Construct The Effective Fragment Parameters Object
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");

  par->set_dipole_polarizability(densityMatrixDipolePolarizability);
  if ((options_.get_int("DMATPOL_FIELD_RANK") == 2) && (options_.get_bool("DMATPOL_FF_AB_INITIO"))) 
  par->set_dipole_dipole_hyperpolarizability(densityMatrixDipoleDipoleHyperpolarizability);

  std::vector<std::shared_ptr<psi::Vector>> centres;
  std::shared_ptr<psi::Vector> com = std::make_shared<psi::Vector>("", 3);
  com->set(0, wfn_->molecule()->center_of_mass().get(0));
  com->set(1, wfn_->molecule()->center_of_mass().get(1));
  com->set(2, wfn_->molecule()->center_of_mass().get(2));
  centres.push_back(com);
  par->set_centres(centres);

  // Return
  return par;
}
