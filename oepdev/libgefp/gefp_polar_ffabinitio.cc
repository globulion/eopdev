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
 
  // Compute The Density Matrix Susceptibility Tensors 
  const double h = 1.0/(2.0 * s);

  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixSusceptibility;
  std::vector<std::shared_ptr<psi::Matrix>> B;

  std::shared_ptr<psi::Matrix> Bx = matrixFactory->create_shared_matrix("B[x] Orthogonal AO Basis"); 
  std::shared_ptr<psi::Matrix> By = matrixFactory->create_shared_matrix("B[y] Orthogonal AO Basis"); 
  std::shared_ptr<psi::Matrix> Bz = matrixFactory->create_shared_matrix("B[z] Orthogonal AO Basis"); 

  Bx->copy(Dmatset[0]); Bx->subtract(Dmatset[1]); Bx->scale(h);
  By->copy(Dmatset[2]); By->subtract(Dmatset[3]); By->scale(h);
  Bz->copy(Dmatset[4]); Bz->subtract(Dmatset[5]); Bz->scale(h);

  B.push_back(Bx);
  B.push_back(By);
  B.push_back(Bz);

  densityMatrixSusceptibility.push_back(B);
  cout << " Density Matrix Susceptibility computed!\n";

  // Construct The Effective Fragment Parameters Object
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  par->set_dipole_polarizability(densityMatrixSusceptibility);
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
