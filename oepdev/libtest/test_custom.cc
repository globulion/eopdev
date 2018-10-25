#include "test.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "../libutil/cphf.h"
#include "../libutil/unitary_optimizer.h"

using namespace std;

double oepdev::test::Test::test_custom(void)
{
  /* Compute the dipole integrals */
  std::vector<std::shared_ptr<psi::Matrix>> F;
  for (int z=0; z<3; ++z) F.push_back(std::make_shared<psi::Matrix>("",wfn_->basisset()->nbf(),wfn_->basisset()->nbf()));
  std::shared_ptr<psi::IntegralFactory> ints = std::make_shared<psi::IntegralFactory>(wfn_->basisset());
  ints->ao_dipole()->compute(F);

  /* Solve the CPHF equations */
  std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn_, options_);
  solver->compute();

  // Accumulate errors
  double r_sum = 0.0;

  /* Compute G tensor */
  std::vector<std::shared_ptr<psi::Matrix>> G;
  for (int z1=0; z1<3; ++z1) {
       for (int z2=0; z2<3; ++z2) {
            int idx = 3*z1 + z2;
            G.push_back(psi::Matrix::doublet(solver->F_mo(z1), solver->X_mo(z2), false, true));
       }
  }

  // Compute R tensor and minimize Z[X] to get X
  cout << " Computing R tensor...\n";
  int nocc = solver->nocc();
  int n_ = nocc;
  int n2_= n_*n_;
  int n3_= n_*n2_;
  int n4_= n_*n3_;
  int n5_= n_*n4_;
  int n6_= n_*n5_;
  //double* R = new double[n3_];
  //double* P = new double[n6_];
  double R[15625] ;
  double P[125] ;
  for (int i=0; i<n3_; ++i) P[i] = 0.0;
  for (int i=0; i<n6_; ++i) R[i] = 0.0;
  // ...
  for (int i=0; i<nocc; ++i) {
    for (int n=0; n<nocc; ++n) {
      if (i==n) {
          for (int j=0; j<nocc; ++j) {
          for (int k=0; k<nocc; ++k) {
          for (int l=0; l<nocc; ++l) {
          for (int m=0; m<nocc; ++m) {
          //     R[IDX6(i,n,j,k,l,m)] = G[0]->get(j,k) * G[0]->get(l,m) + 
          //                            G[1]->get(j,k) * G[3]->get(l,m) + 
          //                            G[3]->get(j,k) * G[1]->get(l,m) + 
          //                            G[2]->get(j,k) * G[6]->get(l,m) + 
          //                            G[6]->get(j,k) * G[2]->get(l,m) + 
          //                            G[4]->get(j,k) * G[4]->get(l,m) + 
          //                            G[5]->get(j,k) * G[7]->get(l,m) + 
          //                            G[7]->get(j,k) * G[5]->get(l,m) + 
          //                            G[8]->get(j,k) * G[8]->get(l,m) -
          //                             (G[0]->get(j,k) + G[4]->get(j,k) + G[8]->get(j,k)) *
          //                             (G[0]->get(l,m) + G[4]->get(l,m) + G[8]->get(l,m)) ;
              R[IDX6(i,n,j,k,l,m)] = (G[0]->get(j,k) + G[4]->get(j,k) + G[8]->get(j,k)) *
                                     (G[0]->get(l,m) + G[4]->get(l,m) + G[8]->get(l,m)) ;
          }
          }
          }
          }
      }
    }
  }

  // Compute the distributed polarizabilities and XPOL-LMO centroids
  cout << " Minimization of Z[X]...\n";
  std::shared_ptr<oepdev::UnitaryOptimizer_4_2> optimizer = std::make_shared<oepdev::UnitaryOptimizer_4_2>(R, P, nocc, 1.0e-14, 5000, true);
  // Release
  //delete[] P;
  //delete[] R;
  bool success_min = optimizer->minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer->X();
  X_min->print();

  /* Check the final Z value */ 
  double Z=0.0; double** x=X_min->pointer();
  for (int i=0; i<nocc; ++i) {
          for (int j=0; j<nocc; ++j) {
          for (int k=0; k<nocc; ++k) {
          for (int l=0; l<nocc; ++l) {
          for (int m=0; m<nocc; ++m) {
               //Z += X[j][i] * X[k][i] * X[l][i] * X[m][i] * R[IDX6(0,0,j,k,l,m)];
               Z += x[j][i] * x[k][i] * x[l][i] * x[m][i] * 
                                     (G[0]->get(j,k) + G[4]->get(j,k) + G[8]->get(j,k)) *
                                     (G[0]->get(l,m) + G[4]->get(l,m) + G[8]->get(l,m)) ;
          }
          }
          }
          }
  }
  cout << Z << endl;

  /* Compute the transformation matrix */
  std::shared_ptr<psi::Matrix> X = psi::Matrix::doublet(solver->T(), X_min, false, false);
  X->set_name("QUO"); X->print();
  std::shared_ptr<psi::Matrix> U = psi::Matrix::triplet(solver->Cocc(), solver->T(), X_min, false, false, false);
  solver->Cocc()->print();
  U->set_name("Nex LCAO-MO matrix");
  U->print();

  /* Compute the LMO centroids */ 
  std::shared_ptr<psi::Matrix> r_lmo = std::make_shared<psi::Matrix>("LMO centroids (Angstroms)", nocc, 3);
  for (int z=0; z<3; ++z) {
       std::shared_ptr<psi::Matrix> r = psi::Matrix::triplet(U, F[z], U, true, false, false);
       for (int o=0; o<nocc; ++o) r_lmo->set(o, z, r->get(o,o));
  }
  r_lmo->scale(-0.5291772086);
  r_lmo->print();

  /* Compute the distributed polarizabilities */
  std::vector<std::shared_ptr<psi::Matrix>> a_lmo;
  std::shared_ptr<psi::Matrix> a_sum = std::make_shared<psi::Matrix>("Total Sum", 3, 3);
  double z;
  for (int o=0; o<nocc; ++o) {
       std::shared_ptr<psi::Matrix> a = std::make_shared<psi::Matrix>("DPOL", 3, 3);
       for (int z1=0; z1<3; ++z1) {
            for (int z2=0; z2<3; ++z2) {
                 std::shared_ptr<psi::Matrix> A = psi::Matrix::triplet(X_min, G[z1*3+z2], X_min, true, false, false);
                 a->set(z1, z2, A->get(o,o));
            }
       }
       a_lmo.push_back(a);
       a_sum->add(a);
       a->print();
       z += a->trace()*a->trace();
  }
  cout << " Check: z = " << z << endl;
  a_sum->print();
  solver->polarizability()->print();

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}

//  // Compute dipole integrals
//  std::vector<std::shared_ptr<psi::Matrix>> Mao;
//  int nbf = wfn_->basisset()->nbf();
//  for (int z=0;z<3;++z) Mao.push_back(std::make_shared<psi::Matrix>("",nbf,nbf));
//  psi::IntegralFactory ints(wfn_->basisset());
//  std::shared_ptr<psi::OneBodyAOInt> dipInt(ints.ao_dipole());
//  dipInt->compute(Mao);
//
//  // Compute numerical derivatives of dipole integrals
//  const double df = 0.0005; double h = 1./(2.0*df);
//  std::shared_ptr<psi::SuperFunctional> func = oepdev::create_superfunctional("HF", options_);
//  std::shared_ptr<psi::Matrix> Ca_x1, Ca_x2, Ca_y1, Ca_y2, Ca_z1, Ca_z2;
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation( df   , 0.000, 0.000);
//  scf->compute_energy(); Ca_x1 = std::make_shared<psi::Matrix>(scf->Ca());}
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation(-df   , 0.000, 0.000);
//  scf->compute_energy(); Ca_x2 = std::make_shared<psi::Matrix>(scf->Ca());}
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation( 0.000, df   , 0.000); 
//  scf->compute_energy(); Ca_y1 = std::make_shared<psi::Matrix>(scf->Ca());}
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation( 0.000,-df   , 0.000);
//  scf->compute_energy(); Ca_y2 = std::make_shared<psi::Matrix>(scf->Ca());}
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation( 0.000, 0.000, df   );
//  scf->compute_energy(); Ca_z1 = std::make_shared<psi::Matrix>(scf->Ca());}
// {std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, func, options_, wfn_->psio());
//  scf->set_perturbation( 0.000, 0.000,-df   );
//  scf->compute_energy(); Ca_z2 = std::make_shared<psi::Matrix>(scf->Ca());}
//
//  std::shared_ptr<psi::Matrix> M_x_x1 = psi::Matrix::triplet(Ca_x1, Mao[0], Ca_x1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_x_x2 = psi::Matrix::triplet(Ca_x2, Mao[0], Ca_x2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_x1 = psi::Matrix::triplet(Ca_x1, Mao[1], Ca_x1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_x2 = psi::Matrix::triplet(Ca_x2, Mao[1], Ca_x2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_x1 = psi::Matrix::triplet(Ca_x1, Mao[2], Ca_x1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_x2 = psi::Matrix::triplet(Ca_x2, Mao[2], Ca_x2, true, false, false);
//  //
//  std::shared_ptr<psi::MaRrix> M_x_y1 = psi::Matrix::triplet(Ca_y1, Mao[0], Ca_y1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_x_y2 = psi::Matrix::triplet(Ca_y2, Mao[0], Ca_y2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_y1 = psi::Matrix::triplet(Ca_y1, Mao[1], Ca_y1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_y2 = psi::Matrix::triplet(Ca_y2, Mao[1], Ca_y2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_y1 = psi::Matrix::triplet(Ca_y1, Mao[2], Ca_y1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_y2 = psi::Matrix::triplet(Ca_y2, Mao[2], Ca_y2, true, false, false);
//  //
//  std::shared_ptr<psi::MaRrix> M_x_z1 = psi::Matrix::triplet(Ca_z1, Mao[0], Ca_z1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_x_z2 = psi::Matrix::triplet(Ca_z2, Mao[0], Ca_z2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_z1 = psi::Matrix::triplet(Ca_z1, Mao[1], Ca_z1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_y_z2 = psi::Matrix::triplet(Ca_z2, Mao[1], Ca_z2, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_z1 = psi::Matrix::triplet(Ca_z1, Mao[2], Ca_z1, true, false, false);
//  std::shared_ptr<psi::Matrix> M_z_z2 = psi::Matrix::triplet(Ca_z2, Mao[2], Ca_z2, true, false, false);
//  //
//  // derivatives
//  //
//  std::shared_ptr<psi::MaRrix> dM_x_x = std::make_shared<psi::Matrix>(M_x_x1); dM_x_x->subtract(M_x_x2); dM_x_x->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_x_y = std::make_shared<psi::Matrix>(M_x_y1); dM_x_y->subtract(M_x_y2); dM_x_y->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_x_z = std::make_shared<psi::Matrix>(M_x_z1); dM_x_z->subtract(M_x_z2); dM_x_z->scale(h);
//  //
//  std::shared_ptr<psi::MaRrix> dM_y_x = std::make_shared<psi::Matrix>(M_y_x1); dM_y_x->subtract(M_y_x2); dM_y_x->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_y_y = std::make_shared<psi::Matrix>(M_y_y1); dM_y_y->subtract(M_y_y2); dM_y_y->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_y_z = std::make_shared<psi::Matrix>(M_y_z1); dM_y_z->subtract(M_y_z2); dM_y_z->scale(h);
//  //
//  std::shared_ptr<psi::MaRrix> dM_z_x = std::make_shared<psi::Matrix>(M_z_x1); dM_z_x->subtract(M_z_x2); dM_z_x->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_z_y = std::make_shared<psi::Matrix>(M_z_y1); dM_z_y->subtract(M_z_y2); dM_y_y->scale(h);
//  std::shared_ptr<psi::MaRrix> dM_z_z = std::make_shared<psi::Matrix>(M_z_z1); dM_z_z->subtract(M_z_z2); dM_z_z->scale(h);


