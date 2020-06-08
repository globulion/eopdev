#include <iostream>
#include <random>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"
#include "../libutil/unitary_optimizer.h"

using namespace std;

#define IDX_3(i,j,k) (n2*(i)+n1*(j)+(k))
#define IDX_6(i,j,k,l,m,n) (n5*(i)+n4*(j)+n3*(k)+n2*(l)+n1*(m)+(n))

//-- AbInitioPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::AbInitioPolarGEFactory::AbInitioPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::PolarGEFactory(wfn, opt)
{
  cphfSolver_ = std::make_shared<oepdev::CPHF>(wfn_, options_);
  cphfSolver_->compute();
}
oepdev::AbInitioPolarGEFactory::~AbInitioPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::AbInitioPolarGEFactory::compute()
{
  // Sizing
  int nbf = wfn_->basisset()->nbf();
  int nocc= cphfSolver_->nocc();

  // Allocate
  std::map<int, char> mm;
  mm[0] = 'X'; mm[1] = 'Y'; mm[2] = 'Z';
  std::shared_ptr<MatrixFactory> matrixFactory = std::make_shared<psi::MatrixFactory>();
  matrixFactory->init_with(1, &nbf, &nbf);
  psi::IntegralFactory integrals(wfn_->basisset());
  std::shared_ptr<psi::Matrix> Sao = std::make_shared<psi::Matrix>("Overlap AO Ints", nbf, nbf);
  std::shared_ptr<psi::Matrix> X   = std::make_shared<psi::Matrix>("Lowdin Symmetric Orthogonalizer: S^{-1/2}", nbf, nbf);
  std::shared_ptr<psi::Matrix> Y   = std::make_shared<psi::Matrix>("Lowdin Symmetric Deorthogonalizer: S^{+1/2}", nbf, nbf);
  std::vector<std::shared_ptr<psi::Matrix>> Mao, Mao_bar, Mmo_ao_bar, L, G;
  for (int z = 0; z<3; ++z) {
       Mao.push_back(std::make_shared<psi::Matrix>("Dipole AO Ints", nbf, nbf));
       G  .push_back(std::make_shared<psi::Matrix>("Left Inverse L Tensor", nocc, nbf));
  }
  std::shared_ptr<psi::Matrix> eigvec   = matrixFactory->create_shared_matrix("Eigenvectors of Sao");
  std::shared_ptr<psi::Matrix> eigval_m = matrixFactory->create_shared_matrix("Eigenvalues of Sao");
  std::shared_ptr<psi::Matrix> temp     = matrixFactory->create_shared_matrix("Temporary");
  std::shared_ptr<psi::Matrix> proj     = matrixFactory->create_shared_matrix("Projector in Orthogonal AO Space");
  std::shared_ptr<psi::Vector> eigval    (matrixFactory->create_vector());

  // Compute Dipole Integrals
  std::shared_ptr<psi::OneBodyAOInt> dipInt(integrals.ao_dipole());
  dipInt->compute(Mao);

  // Compute Overlap Integrals
  std::shared_ptr<psi::OneBodyAOInt> ovlInt(integrals.ao_overlap());
  ovlInt->compute(Sao);

  // Compute Lowdin Symmetric Orthogonalizer
  Sao->diagonalize(eigvec, eigval);
  double min_S = fabs(eigval->get(0,0));
  for (int i=0; i<nbf; ++i) {
       if (min_S > eigval->get(0, i)) 
           min_S = eigval->get(0, i);
       double v = 1.0 / sqrt(eigval->get(0, i));
       eigval->set(0, i, v);
  }
  outfile->Printf("  Minimum eigenvalue in the overlap matrix is %14.10E.\n", min_S);

  eigval_m->set_diagonal(eigval);
  temp->gemm(false, true, 1.0, eigval_m, eigvec, 0.0);
  X->gemm(false, false, 1.0, eigvec, temp, 0.0);

  // Compute Deorthogonalizer
  Y->copy(X);
  Y->invert();

  // LCAO-LMO Coefficients in Non-Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> U = psi::Matrix::doublet(wfn_->Ca_subset("AO","OCC"), cphfSolver_->localizer()->U(), false, false);

  // LCAO-LMO Coefficients in Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> Ubar = psi::Matrix::doublet(Y, U, false, false);
  Ubar->set_name("LCAO-LMO Coefficients in Orthogonal AO Basis");

  // One-Particle Density Matrix in Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> Dbar = psi::Matrix::doublet(Ubar, Ubar, false, true);
  Dbar->set_name("One-Particle Density Matrix in Orthogonal AO Basis");

  // Transform Dipole Integrals to Orthogonal Basis
  for (int z = 0; z<3; ++z) {
       Mao[z]->scale(-1.0);
       Mao_bar.push_back(psi::Matrix::triplet(X, Mao[z], X, true, false, false));
  }

  // Transform Left Axis of Dipole AO Integrals to LMO Basis
  for (int z = 0; z<3; ++z) {
       Mmo_ao_bar.push_back(psi::Matrix::doublet(Ubar, Mao_bar[z], true, false));
  }

  // Compute L Vector of Matrices
  proj->identity();
  proj->subtract(Dbar);
  for (int z = 0; z<3; ++z) {
       L.push_back(psi::Matrix::doublet(Mmo_ao_bar[z], proj, false, false));
  }

  // Compute Collection of Left L Inverse Matrces
  for (int o = 0; o < nocc; ++o) {
       std::shared_ptr<psi::Matrix> Li = std::make_shared<psi::Matrix>("", nbf, 3);
       for (int z = 0; z < 3; ++z) {
            for (int n = 0; n < nbf; ++n) {
                 double v = L[z]->get(o, n);
                 Li->set(n, z, v);
            }
       }
       std::shared_ptr<psi::Matrix> Q  = psi::Matrix::doublet(Li, Li, true, false);
       Q->invert();
       std::shared_ptr<psi::Matrix> Li_left = psi::Matrix::doublet(Q, Li, false, true);
       for (int z = 0; z < 3; ++z) {
            for (int n = 0; n < nbf; ++n) {
                 double v = Li_left->get(z, n);
                 G[z]->set(o, n, v);
            }
       }
  }

  // Transform G by Distributed Polarizabilities
  for (int o = 0; o < nocc; ++o) {
       std::shared_ptr<psi::Matrix> pol = cphfSolver_->polarizability(o);
       for (int n = 0; n < nbf; ++n) {
            double v_x = pol->get(0,0) * G[0]->get(o, n) + 
                         pol->get(1,0) * G[1]->get(o, n) +
                         pol->get(2,0) * G[2]->get(o, n);
            double v_y = pol->get(0,1) * G[0]->get(o, n) + 
                         pol->get(1,1) * G[1]->get(o, n) +
                         pol->get(2,1) * G[2]->get(o, n);
            double v_z = pol->get(0,2) * G[0]->get(o, n) + 
                         pol->get(1,2) * G[1]->get(o, n) +
                         pol->get(2,2) * G[2]->get(o, n);

            G[0]->set(o, n, v_x);
            G[1]->set(o, n, v_y);
            G[2]->set(o, n, v_z);
       }
  }

  // Compute The Density Matrix Susceptibility Tensors 
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> densityMatrixSusceptibility;
  for (int o = 0; o < nocc; ++o) {
       std::vector<std::shared_ptr<psi::Matrix>> Bi;
       for (int z = 0; z < 3; ++z) {
            std::shared_ptr<psi::Matrix> Biz = matrixFactory->create_shared_matrix("Bi[z] Orthogonal AO Basis"); 
            for (int i = 0; i < nbf; ++i) { 
                 for (int j = 0; j < nbf; ++j) {
                      double v =  G[z]->get(o, i) * Ubar->get(j, o) + 
                                  G[z]->get(o, j) * Ubar->get(i, o) ;
                      for (int k = 0; k < nbf; ++k) {
                           v += G[z]->get(o, k) * (Dbar->get(i, k) * Ubar->get(j, o) +
                                                   Dbar->get(k, j) * Ubar->get(i, o) );
                      }
                      Biz->set(i, j, v);
                 }
            }
            std::shared_ptr<psi::Matrix> Biz_t = psi::Matrix::triplet(X, Biz, X, false, false, false);
            Biz_t->set_name(oepdev::string_sprintf("B[%c](%d) in Non-Orthogonal AO Basis", mm[z], o+1));
            Biz_t->scale(-0.25000000);
            Bi.push_back(Biz_t);
       }
       densityMatrixSusceptibility.push_back(Bi);
  }

  // Construct The Effective Fragment Parameters Object
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  par->set_dipole_polarizability(densityMatrixSusceptibility);
  std::vector<std::shared_ptr<psi::Vector>> centres;
  for (int o=0; o<cphfSolver_->nocc(); ++o) {
       centres.push_back(std::make_shared<psi::Vector>(*(cphfSolver_->lmo_centroid(o))));
  }
  par->set_centres(centres);

  // Return
  return par;
}
//-- UnitaryTransformedMOPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::UnitaryTransformedMOPolarGEFactory::UnitaryTransformedMOPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::AbInitioPolarGEFactory(wfn, opt)
{

}
oepdev::UnitaryTransformedMOPolarGEFactory::~UnitaryTransformedMOPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::UnitaryTransformedMOPolarGEFactory::compute()
{
  // Sizing
  int npoints = 1;
  int nbf = wfn_->basisset()->nbf();
  int nocc = cphfSolver_->nocc();
  size_t n1 = nocc;
  size_t n2 = nocc*nocc;
  size_t n3 = nocc*nocc*nocc;
  size_t n4 = n3*nocc;
  size_t n5 = n3*n2;
  size_t n6 = n3*n3;

  // Allocate
  std::map<int, char> mm;
  mm[0] = 'X'; mm[1] = 'Y'; mm[2] = 'Z';
  std::shared_ptr<MatrixFactory> matrixFactory = std::make_shared<psi::MatrixFactory>();
  matrixFactory->init_with(1, &nbf, &nbf);
  psi::IntegralFactory integrals(wfn_->basisset());
  std::shared_ptr<psi::Matrix> Sao = std::make_shared<psi::Matrix>("Overlap AO Integrals", nbf, nbf);
  std::shared_ptr<psi::Matrix> X   = std::make_shared<psi::Matrix>("Lowdin Symmetric Orthogonalizer: S^{-1/2}", nbf, nbf);
  std::shared_ptr<psi::Matrix> Y   = std::make_shared<psi::Matrix>("Lowdin Symmetric Deorthogonalizer: S^{+1/2}", nbf, nbf);
  std::vector<std::shared_ptr<psi::Matrix>> Mao, Mao_bar, Mmo_ao_bar, L, G;
  for (int z = 0; z<3; ++z) {
       Mao.push_back(std::make_shared<psi::Matrix>("Dipole AO Ints", nbf, nbf));
       G  .push_back(std::make_shared<psi::Matrix>("Left Inverse L Tensor", nocc, nbf));
  }
  std::shared_ptr<psi::Matrix> eigvec   = matrixFactory->create_shared_matrix("Eigenvectors of Sao");
  std::shared_ptr<psi::Matrix> eigval_m = matrixFactory->create_shared_matrix("Eigenvalues of Sao");
  std::shared_ptr<psi::Matrix> temp     = matrixFactory->create_shared_matrix("Temporary");
  std::shared_ptr<psi::Matrix> proj     = matrixFactory->create_shared_matrix("Projector in Orthogonal AO Space");
  std::shared_ptr<psi::Vector> eigval    (matrixFactory->create_vector());

  // Compute Dipole Integrals
  std::shared_ptr<psi::OneBodyAOInt> dipInt(integrals.ao_dipole());
  dipInt->compute(Mao);

  // Compute Overlap Integrals
  std::shared_ptr<psi::OneBodyAOInt> ovlInt(integrals.ao_overlap());
  ovlInt->compute(Sao);

  // Compute Lowdin Symmetric Orthogonalizer
  Sao->diagonalize(eigvec, eigval);
  double min_S = fabs(eigval->get(0,0));
  for (int i=0; i<nbf; ++i) {
       if (min_S > eigval->get(0, i)) 
           min_S = eigval->get(0, i);
       double v = 1.0 / sqrt(eigval->get(0, i));
       eigval->set(0, i, v);
  }
  outfile->Printf("  Info: Minimum eigenvalue in the overlap matrix is %14.10E.\n", min_S);

  eigval_m->set_diagonal(eigval);
  temp->gemm(false, true, 1.0, eigval_m, eigvec, 0.0);
  X->gemm(false, false, 1.0, eigvec, temp, 0.0);

  // Compute Deorthogonalizer
  Y->copy(X);
  Y->invert();

  // LCAO-LMO Coefficients in Non-Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> U;
  if (options_.get_bool("CPHF_LOCALIZE") == true) {
    U = psi::Matrix::doublet(wfn_->Ca_subset("AO","OCC"), cphfSolver_->localizer()->U(), false, false);
  } else {
    U = wfn_->Ca_subset("AO","OCC");
  }
  // LCAO-LMO Coefficients in Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> Ubar = psi::Matrix::doublet(Y, U, false, false);
  Ubar->set_name("LCAO-LMO Coefficients in Orthogonal AO Basis");
  double** Ca = Ubar->pointer();

  // One-Particle Density Matrix in Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> Dbar = psi::Matrix::doublet(Ubar, Ubar, false, true);
  Dbar->set_name("One-Particle Density Matrix in Orthogonal AO Basis");
  double** Da = Dbar->pointer();

  // Transform Dipole Integrals to Orthogonal Basis
  for (int z = 0; z<3; ++z) {
       //Mao[z]->scale(-1.0); // This is removed to be consistent with the paper
       Mao_bar.push_back(psi::Matrix::triplet(X, Mao[z], X, true, false, false));
  }

  // Transform Left Axis of Dipole AO Integrals to LMO Basis
  for (int z = 0; z<3; ++z) {
       Mmo_ao_bar.push_back(psi::Matrix::doublet(Ubar, Mao_bar[z], true, false));
  }

  // Compute L Vector of Matrices
  proj->identity();
  proj->subtract(Dbar);
  for (int z = 0; z<3; ++z) {
       L.push_back(psi::Matrix::doublet(Mmo_ao_bar[z], proj, false, false));
  }

  // Compute Collection of Left L Inverse Matrces (G matrices)
  for (int o = 0; o < nocc; ++o) {
       std::shared_ptr<psi::Matrix> Li = std::make_shared<psi::Matrix>("", nbf, 3);
       for (int z = 0; z < 3; ++z) {
            for (int n = 0; n < nbf; ++n) {
                 double v = L[z]->get(o, n);
                 Li->set(n, z, v);
            }
       }
       std::shared_ptr<psi::Matrix> Q  = psi::Matrix::doublet(Li, Li, true, false);
       Q->invert();
       std::shared_ptr<psi::Matrix> Li_left = psi::Matrix::doublet(Q, Li, false, true);
       for (int z = 0; z < 3; ++z) {
            for (int n = 0; n < nbf; ++n) {
                 double v = Li_left->get(z, n);
                 G[z]->set(o, n, v);
            }
       }
  }

  // ===> Compute the X unitary matrix <=== //

  // --> Set up the trial electric fields and compute perturbed density matrices <-- //
  std::vector<std::shared_ptr<psi::Vector>> fields;
  std::vector<std::shared_ptr<psi::Matrix>> dmats;
  for (int N=0; N<npoints; ++N) {
       fields.push_back(draw_field());
       cout << oepdev::string_sprintf(" Computation for N=%2d F=[%14.4f, %14.4f, %14.4f]\n",N+1,
                                        fields[N]->get(0), fields[N]->get(1), fields[N]->get(2));

       std::shared_ptr<psi::Matrix> dmat = std::make_shared<psi::Matrix>("",nbf,nbf);
       dmat->copy(perturbed_state(fields[N])->Da());
       std::shared_ptr<psi::Matrix> dmat_diff_bar = psi::Matrix::triplet(Y, dmat, Y, false, false, false);
       dmat_diff_bar->subtract(Dbar);
       dmats.push_back(dmat_diff_bar);
  }

  // --> Allocate data <-- //
  double* R = nullptr;
  double* P = nullptr;
  try {
    R  = new double[n6];
    P  = new double[n3];
  } catch (std::bad_alloc &e) {
    psi::outfile->Printf("Error allocating 6-rank and 3-rd rank tensors \n%s\n", e.what());
    exit(EXIT_FAILURE);
  }

  // --> Compute the least-squares R and P tensors <-- //
  psi::timer_on(" Computation of R tensor");
  for (int i=0; i<nocc; ++i) {
  for (int j=0; j<nocc; ++j) {
  for (int k=0; k<nocc; ++k) {
  for (int l=0; l<nocc; ++l) {
  for (int m=0; m<nocc; ++m) {
  for (int n=0; n<nocc; ++n) {
       double v = 0.0;
       for (int N=0; N<npoints; ++N) {
            for (int ia=0; ia<nbf; ++ia) {
             for (int ib=0; ib<nbf; ++ib) {
              double cai = Ca[ia][i];
              double cbi = Ca[ib][i];
              double caj = Ca[ia][j];
              double cbj = Ca[ib][j];
              std::vector<double> vg1;
              std::vector<double> vg2;
              for (int z=0; z<3; ++z) {
                 vg1.push_back(0.0);
                 vg2.push_back(0.0);
              }
              for (int ic=0; ic<nbf; ++ic) {
                 double cd1 = Da[ia][ic] * cbi + Da[ib][ic] * cai;
                 double cd2 = Da[ia][ic] * cbj + Da[ib][ic] * caj;
                 for (int z=0; z<3; ++z) {
                   vg1[z] += cd1 * G[z]->get(i,ic);
                   vg2[z] += cd2 * G[z]->get(j,ic);
                 }
              }
              for (int z1=0; z1<3; ++z1) {
                   double gia_z1 = G[z1]->get(i,ia);
                   double gib_z1 = G[z1]->get(i,ib);
              for (int z2=0; z2<3; ++z2) { 
                   double gja_z2 = G[z2]->get(j,ia);
                   double gjb_z2 = G[z2]->get(j,ib);
              for (int w1=0; w1<3; ++w1) {
              for (int w2=0; w2<3; ++w2) {
                 v += cphfSolver_->polarizability(k,m)->get(z1,w1) * fields[N]->get(w1) *
                      cphfSolver_->polarizability(l,n)->get(z2,w2) * fields[N]->get(w2) * 
                     (cai * gib_z1 + cbi * gia_z1 - vg1[z1]) * 
                     (caj * gjb_z2 + cbj * gja_z2 - vg2[z2]);
              }}}}
             }
            }
       }
       v /= 16.0;
       R[IDX_6(i,j,k,l,m,n)] = v;
       cout << " R= "<< i << j << k << l << m << n << "   " << v << endl;
  }}}}}}
  psi::timer_off(" Computation of R tensor");

  psi::timer_on (" Computation of P tensor");
  for (int i=0; i<nocc; ++i) {
  for (int j=0; j<nocc; ++j) {
  for (int k=0; k<nocc; ++k) {
       double v = 0.0;
       for (int N=0; N<npoints; ++N) {
            for (int ia=0; ia<nbf; ++ia) {
             for (int ib=0; ib<nbf; ++ib) {
              double d_ref = dmats[N]->get(ia,ib);
              double cai = Ca[ia][i];
              double cbi = Ca[ib][i];
              std::vector<double> vg;
              for (int z=0; z<3; ++z) vg.push_back(0.0);
              for (int ic=0; ic<nbf; ++ic) {
                for (int z=0; z<3; ++z) vg[z] += (Da[ia][ic] * cbi + Da[ib][ic] * cai) * G[z]->get(i,ic);
              }
              for (int z=0; z<3; ++z) {
                   double gia_z = G[z]->get(i,ia);
                   double gib_z = G[z]->get(i,ib);
                   double gz    = cai * gib_z + cbi * gia_z - vg[z];
              for (int w=0; w<3; ++w) {
                   v += d_ref * cphfSolver_->polarizability(j,k)->get(z,w) * gz * fields[N]->get(w);
              }
              }
             }
            }
       }
       v /=-2.0;
       P[IDX_3(i,j,k)] = v;
       cout << " P= "<< i << j << k << "   " << v << endl;
  }}}
  psi::timer_off(" Computation of P tensor");

  double Z_0 = 0.0;
  for (int N=0; N<npoints; ++N) {
      for (int ia=0; ia<nbf; ++ia) {
       for (int ib=0; ib<nbf; ++ib) {
            double d_ref = dmats[N]->get(ia,ib);
            Z_0 += d_ref * d_ref;
  }}}

  // --> Perform the optimization <-- //
  std::shared_ptr<psi::Matrix> Xt;
  double Z;
  {
     oepdev::UnitaryOptimizer_4_2 optimizer(R, P, nocc, 1.0e-9, 10, true); 
     optimizer.minimize();
     Xt = optimizer.X();
     Z  = optimizer.Z();
  }
  psi::outfile->Printf("\n @Optimizer: Z_0     = %14.6E\n", Z_0);
  psi::outfile->Printf(  " @Optimizer: Z       = %14.6E\n", Z  );
  psi::outfile->Printf(  " @Optimizer: Z + Z_0 = %14.6E\n\n", Z+Z_0);

  Xt->set_name("Unitary Transformation Xt");
  Xt->print();

  // Clean up
  delete[] R;
  delete[] P;
  
  // --> Compute the ab-initio b tensors (associated with LCAO-MO coefficients)
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> b_susc;
  for (int i=0; i<nocc; ++i) {
       std::vector<std::shared_ptr<psi::Matrix>> bz;
       for (int z=0; z<3; ++z) {
            bz.push_back(std::make_shared<psi::Matrix>(
            oepdev::string_sprintf("LCAO-MO Susceptibility: b[%c](%d) in Orthogonal AO Basis", mm[z], i+1), nbf, 1));
       }
       b_susc.push_back(bz);
  }
  for (int w=0; w<3; ++w) {
       for (int i=0; i<nocc; ++i) {
            for (int ia=0; ia<nbf; ++ia) {
                 double v = 0.0;
                 for (int j=0; j<nocc; ++j) {
                      for (int k=0; k<nocc; ++k) {
                           double xji = Xt->get(j,i);
                           double xki = Xt->get(k,i);
                           for (int z=0; z<3; ++z) {
                                v += xji * xki * cphfSolver_->polarizability(j,k)->get(z,w) * G[z]->get(i,ia);
                           }
                      }
                 }
                 v /= 4.0;
                 b_susc[i][w]->set(ia,0,v);
            }
       }
  }

  // --> Compute the ab-initio B tensors (associated with density matrix elements)
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> B_susc;
  for (int i=0; i<nocc; ++i) {
       std::vector<std::shared_ptr<psi::Matrix>> Bz;
       for (int z=0; z<3; ++z) {
            Bz.push_back(std::make_shared<psi::Matrix>(
            oepdev::string_sprintf("Density Matrix Susceptibility: B[%c](%d) in Orthogonal AO Basis", mm[z], i+1), nbf, nbf));
       }
       B_susc.push_back(Bz);
  }
  //
  for (int ia=0; ia<nbf; ++ia) {
       for (int ib=0; ib<nbf; ++ib) {
            for (int i=0; i<nocc; ++i) {
                 double cai = Ca[ia][i];
                 double cbi = Ca[ib][i];
                 for (int z=0; z<3; ++z) {
                      double vbz  = cai * b_susc[i][z]->get(ib,0) + cbi * b_susc[i][z]->get(ia,0);
                      for (int ic=0; ic<nbf; ++ic) {
                           vbz -=(cbi * Da[ia][ic] + cai * Da[ib][ic]) * b_susc[i][z]->get(ic,0);
                      }
                      B_susc[i][z]->set(ia,ib,vbz);
                 }
            }
       }
  }

  // Transform the susceptibility to non-orthogonal basis
  for (int i=0; i<nocc; ++i) {
       for (int z=0; z<3; ++z) {
            std::shared_ptr<psi::Matrix> Biz_t = psi::Matrix::triplet(X, B_susc[i][z], X, false, false, false);
            Biz_t->set_name(oepdev::string_sprintf("Density Matrix Susceptibility: B[%c](%d) in Non-Orthogonal AO Basis", mm[z], i+1));
            B_susc[i][z] = std::make_shared<psi::Matrix>(*Biz_t);
       }
  }

  // Construct The Effective Fragment Parameters Object
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  par->set_dipole_polarizability(B_susc);

  psi::SharedMatrix geom = std::make_shared<psi::Matrix>(wfn_->molecule()->geometry());
  par->set_matrix("pos", geom);

  return par;
}
