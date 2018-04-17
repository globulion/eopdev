#include <iostream>
#include <random>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"
#include "../libutil/unitary_optimizer.h"
#include "../libutil/scf_perturb.h"

using namespace std;

#define IDX_3(i,j,k) (n2*(i)+n1*(j)+(k))
#define IDX_6(i,j,k,l,m,n) (n5*(i)+n4*(j)+n3*(k)+n2*(l)+n1*(m)+(n))

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
 randomDistribution_(std::uniform_real_distribution<double>(-1.0, 1.0))
{
   double pad = 10.0; // Put to options
   // populate vdwRadius
   vdwRadius_["C"] = 3.000; // Put to options
   vdwRadius_["H"] = 1.000;
   vdwRadius_["N"] = 2.400;
   vdwRadius_["O"] = 2.600;
   vdwRadius_["F"] = 2.300;
   vdwRadius_["Cl"]= 2.900;

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

//-- PolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::PolarGEFactory::PolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::GenEffParFactory(cphf->wfn(), opt),
 cphfSolver_(cphf)
{

}
oepdev::PolarGEFactory::PolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::GenEffParFactory(cphf->wfn(), cphf->options()),
 cphfSolver_(cphf)
{

}
oepdev::PolarGEFactory::~PolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::PolarGEFactory::compute()
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
  par->set_susceptibility(densityMatrixSusceptibility);

  // Return
  return par;
}
std::shared_ptr<psi::Vector> oepdev::PolarGEFactory::draw_field()
{
  std::shared_ptr<psi::Vector> field = std::make_shared<psi::Vector>("", 3);
  const double scale = 0.001; // Put to options
  double fx = random_double() * scale;
  double fy = random_double() * scale;
  double fz = random_double() * scale;
  field->set(0, fx);
  field->set(1, fy);
  field->set(2, fz);
  return field;
}
std::shared_ptr<psi::Matrix> oepdev::PolarGEFactory::perturbed_dmat(const std::shared_ptr<psi::Vector>& field)
{
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, 
                  oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());
  scf->set_perturbation(field);
  scf->compute_energy();
  std::shared_ptr<psi::Matrix> dmat = scf->Da();
  return dmat;
}
std::shared_ptr<psi::Matrix> oepdev::PolarGEFactory::perturbed_dmat(const std::shared_ptr<psi::Vector>& pos, const double& q)
{
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, 
                  oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());
  scf->set_perturbation(pos, q);
  scf->compute_energy();
  std::shared_ptr<psi::Matrix> dmat = scf->Da();
  return dmat;
}

//-- MOScaledPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::MOScaledPolarGEFactory::MOScaledPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::MOScaledPolarGEFactory::MOScaledPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::MOScaledPolarGEFactory::~MOScaledPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::MOScaledPolarGEFactory::compute()
{
  // Sizing
  int npoints = 150;
  int nbf = wfn_->basisset()->nbf();
  int no = cphfSolver_->nocc();
  size_t n1 = no;
  size_t n2 = no*no;
  size_t n3 = no*no*no;
  size_t n4 = n3*no;
  size_t n5 = n3*n2;
  size_t n6 = n3*n3;

  // Compute the ab-initio B tensors (X = 1)
  std::shared_ptr<oepdev::PolarGEFactory> factory_0 = shared_from_this();
  std::shared_ptr<oepdev::GenEffPar> par_0 = factory_0->compute();

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
  double** Ca = Ubar->pointer();

  // One-Particle Density Matrix in Orthogonal AO Basis
  std::shared_ptr<psi::Matrix> Dbar = psi::Matrix::doublet(Ubar, Ubar, false, true);
  Dbar->set_name("One-Particle Density Matrix in Orthogonal AO Basis");
  double** Da = Dbar->pointer();

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

       std::shared_ptr<psi::Matrix> dmat = std::make_shared<psi::Matrix>("",nbf,nbf);
       dmat->copy(perturbed_dmat(fields[N]));
       std::shared_ptr<psi::Matrix> dmat_bar = psi::Matrix::triplet(Y, dmat, Y, false, false, false);
       dmat_diff->subtract(Dbar);
       dmats.push_back(dmat_diff);
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
  for (int i=0; i<no; ++i) {
  for (int j=0; j<no; ++j) {
  for (int k=0; k<no; ++k) {
  for (int l=0; l<no; ++l) {
  for (int m=0; m<no; ++m) {
  for (int n=0; n<no; ++n) {
       double v = 0.0;
       for (int N=0; N<npoints; ++N) {
            for (int ia=0; ia<nbf; ++ia) {
             for (int ib=0; ib<nbf; ++ib) {
              double cai = Ca[ia][i];
              double cbi = Ca[ib][i];
              double cal = Ca[ia][l];
              double cbl = Ca[ib][l];
              std::vector<double> vg1;
              std::vector<double> vg2;
              for (int z=0; z<3; ++z) {
                 vg1.push_back(0.0);
                 vg2.push_back(0.0);
              }
              for (int ic=0; ic<nbf; ++ic) {
                 double cd1 = Da[ia][ic] * cbi + Da[ib][ic] * cai;
                 double cd2 = Da[ia][ic] * cbl + Da[ib][ic] * cal;
                 for (int z=0; z<3; ++z) {
                   vg1[z] += cd1 * G[z]->get(i,ic);
                   vg2[z] += cd2 * G[z]->get(l,ic);
                 }
              }
              for (int z1=0; z1<3; ++z1) {
                   double gia_z1 = G[z1]->get(i,ia);
                   double gib_z1 = G[z1]->get(i,ib);
              for (int z2=0; z2<3; ++z2) { 
                   double gla_z2 = G[z2]->get(l,ia);
                   double glb_z2 = G[z2]->get(l,ib);
              for (int w1=0; w1<3; ++w1) {
              for (int w2=0; w2<3; ++w2) {
                 v += cphfSolver_->polarizability(j,k)->get(z1,w1) * fields[N]->get(w1) *
                      cphfSolver_->polarizability(m,n)->get(z2,w2) * fields[N]->get(w2) * 
                     (cai * gib_z1 + cbi * gia_z1 - vg1[z1]) * 
                     (cal * glb_z2 + cbl * gla_z2 - vg2[z2]);
              }}}}
             }
            }
       }
       v /= 16.0;
       R[IDX_6(i,j,k,l,m,n)] = v;
  }}}}}}
  psi::timer_off(" Computation of R tensor");

  psi::timer_on (" Computation of P tensor");
  for (int i=0; i<no; ++i) {
  for (int j=0; j<no; ++j) {
  for (int k=0; k<no; ++k) {
       double v = 0.0;
       for (int N=0; N<npoints; ++N) {
            for (int ia=0; ia<nbf; ++ia) {
             for (int ib=0; ib<nbf; ++ib) {
              double d_ref = dmats[N]->get(ai,bi);
              double cai = Ca[ia][i];
              double cbi = Ca[ib][i];
              std::vector<double> vg;
              for (int z=0; z<3; ++z) vg.push_back(0.0);
              for (int ic=0; ic<nbf; ++ic) {
                   vg[z] += (Da[ia][ic] * cbi + Da[ib][ic] * cai) * G[z]->get(i,ic);
              }
              for (int z=0; z<3; ++z) {
                   double gia_z = G[z]->get(i,ia);
                   double gib_z = G[z]->get(i,ib);
                   double gz    = cai * gib_z + cbi * gia_z - vg[z]
              for (int w=0; w<3; ++w) {
                   v += d_ref* cphfSolver_->polarizability(j,k)->get(z,w) * gz * fields[N]->get(w);
              }
              }
             }
            }
       }
       v /=-2.0;
       P[IDX_3(i,j,k)] = v;
  }}}
  psi::timer_off(" Computation of P tensor");

  // --> Perform the optimization <-- //
  std::shared_ptr<psi::Matrix> Xt;
  {
     oepdev::UnitaryOptimizer_4_2 optimizer(R, P, no); 
     optimizer.minimize();
     Xt = optimizer.X();
  }

  // Clean up
  delete[] R;
  delete[] P;
  
  // --> Compute the ab-initio b tensors (associated with LCAO-MO coefficients)
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> b_susc;
  for (int i=0; i<no; ++i) {
       std::vector<std::shared_ptr<psi::Matrix>> bz;
       for (int z=0; z<3; ++z) {
            bz.push_back(std::make_shared<psi::Matrix>("LCAO-MO Susceptibility", nbf, 1));
       }
       b_susc.push_back(bz);
  }
  for (int w=0; w<3; ++w) {
       for (int i=0; i<no; ++i) {
            for (int ia=0; ia<nbf; ++ia) {
                 double v = 0.0;
                 for (int j=0; j<no; ++j) {
                      for (int k=0; k<no; ++k) {
                           double xji = Xt->get(j,i);
                           double xki = Xt->get(k,i);
                           for (int z=0; z<3; ++z) {
                                v += xji * xki * cphfSolver_->polarizability(j,k)->get(u,w) * G[u]->get(i,ia);
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
  for (int i=0; i<no; ++i) {
       std::vector<std::shared_ptr<psi::Matrix>> Bz;
       for (int z=0; z<3; ++z) {
            Bz.push_back(std::make_shared<psi::Matrix>("Density Matrix Susceptibility", nbf, nbf));
       }
       B_susc.push_back(Bz);
  }
  //
  for (int ia=0; ia<nbf; ++ia) {
       for (int ib=0; ib<nbf; ++ib) {
            for (int i=0; i<no; ++i) {
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

  // Construct The Effective Fragment Parameters Object
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  par->set_susceptibility(B_susc);

  return par;
}
//-- FieldScaledPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::FieldScaledPolarGEFactory::FieldScaledPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::FieldScaledPolarGEFactory::FieldScaledPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::FieldScaledPolarGEFactory::~FieldScaledPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::FieldScaledPolarGEFactory::compute()
{
  // Sizing
  int npoints = 150;
  int nbf = wfn_->basisset()->nbf();
  int no = cphfSolver_->nocc();

  // Compute the ab-initio (unscaled) B tensors
  std::shared_ptr<oepdev::PolarGEFactory> factory_0 = shared_from_this();
  std::shared_ptr<oepdev::GenEffPar> par_0 = factory_0->compute();

  // Compute the scales:

  // --> Set up the trial electric fields and compute perturbed density matrices <-- //
  std::vector<std::shared_ptr<psi::Vector>> fields;
  std::vector<std::shared_ptr<psi::Matrix>> dmats;
  for (int i=0; i<npoints; ++i) {
       fields.push_back(draw_field());

       std::shared_ptr<psi::Matrix> dmat = std::make_shared<psi::Matrix>("",nbf,nbf);
       dmat->copy(perturbed_dmat(fields[i]));
       dmat->subtract(wfn_->Da());
       dmats.push_back(dmat);
  }

  // --> Allocate data <-- //
  std::shared_ptr<psi::Matrix> A = std::make_shared<psi::Matrix>("", 3*no, 3*no);
  std::shared_ptr<psi::Matrix> a = std::make_shared<psi::Matrix>("", 1, 3*no);

  // --> Compute the least-squares matrices <-- //
  // TODO 

  // --> Perform the fit <-- //
  A->invert();
  std::shared_ptr<psi::Matrix> s = psi::Matrix::doublet(a, A, false, false);

  // --> Scale the ab-initio B tensors by scale values
  // TODO 

  return par_0;
}
