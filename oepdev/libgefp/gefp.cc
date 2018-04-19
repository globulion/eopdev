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
   vdwRadius_["H"] = 4.000;
   vdwRadius_["N"] = 2.400;
   vdwRadius_["O"] = 5.600;
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
  const double scale = options_.get_double("DMATPOL_FIELD_SCALE");
  double fx = random_double() * scale;
  double fy = random_double() * scale;
  double fz = random_double() * scale;
  field->set(0, fx);
  field->set(1, fy);
  field->set(2, fz);
  return field;
}
std::shared_ptr<psi::Vector> oepdev::PolarGEFactory::field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges,
                                                                  const double& x, const double& y, const double& z)
{
  std::shared_ptr<psi::Vector> field = std::make_shared<psi::Vector>("", 3);
  double fx = 0.0;
  double fy = 0.0; 
  double fz = 0.0;
  for (int np=0; np<charges->nrow(); ++np) {
       double rx = charges->get(np, 0);
       double ry = charges->get(np, 1);
       double rz = charges->get(np, 2);
       double q  = charges->get(np, 3);
       double drx= x - rx;
       double dry= y - ry;
       double drz= z - rz;
       double r = sqrt(drx*drx+dry*dry+drz*drz);
       double r3 = q/(r*r*r);
       fx += r3 * drx;
       fy += r3 * dry;
       fz += r3 * drz;
  }
  field->set(0, fx);
  field->set(1, fy);
  field->set(2, fz);
  return field;
}
std::shared_ptr<psi::Vector> oepdev::PolarGEFactory::field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges,
                                                                          const std::shared_ptr<psi::Vector>& pos)
{
  return this->field_due_to_charges(charges, pos->get(0), pos->get(1), pos->get(2));
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
std::shared_ptr<psi::Matrix> oepdev::PolarGEFactory::perturbed_dmat(const std::shared_ptr<psi::Matrix>& charges)
{
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_,
                  oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());
  for (int i=0; i<charges->nrow(); ++i) {
       double x = charges->get(i, 0);
       double y = charges->get(i, 1);
       double z = charges->get(i, 2);
       double q = charges->get(i, 3);
       scf->set_perturbation(x, y, z, q);
  }
  scf->compute_energy();
  std::shared_ptr<psi::Matrix> dmat = scf->Da();
  return dmat;
}

//-- UnitaryTransformedMOPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::UnitaryTransformedMOPolarGEFactory::UnitaryTransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::UnitaryTransformedMOPolarGEFactory::UnitaryTransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
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

  // Compute the ab-initio B tensors (X = 1)
  //std::shared_ptr<oepdev::PolarGEFactory> factory_0 = shared_from_this();
  //std::shared_ptr<oepdev::GenEffPar> par_0 = factory_0->compute();

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
       dmat->copy(perturbed_dmat(fields[N]));
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
  par->set_susceptibility(B_susc);

  return par;
}
//-- ScaledXYZPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::ScaledXYZPolarGEFactory::ScaledXYZPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::ScaledXYZPolarGEFactory::ScaledXYZPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::ScaledXYZPolarGEFactory::~ScaledXYZPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::ScaledXYZPolarGEFactory::compute()
{
  // Sizing
  int npoints = 150;
  int nbf = wfn_->basisset()->nbf();
  int no = cphfSolver_->nocc();

  // Compute the ab-initio (unscaled) B tensors
  std::shared_ptr<oepdev::PolarGEFactory> factory_0 = std::make_shared<oepdev::PolarGEFactory>(cphfSolver_, options_);
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
//-- ScaledAOPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::ScaledAOPolarGEFactory::ScaledAOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::ScaledAOPolarGEFactory::ScaledAOPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::ScaledAOPolarGEFactory::~ScaledAOPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::ScaledAOPolarGEFactory::compute()
{
  // Sizing
  int npoints = options_.get_int("DMATPOL_NSAMPLES");
  int nbf = wfn_->basisset()->nbf();
  int nocc = cphfSolver_->nocc();

  // Compute the ab-initio (unscaled) B tensors
  std::shared_ptr<oepdev::PolarGEFactory> factory_0 = std::make_shared<oepdev::PolarGEFactory>(cphfSolver_, options_);
  std::shared_ptr<oepdev::GenEffPar> par_0 = factory_0->compute();

  // Compute the scales:

  // --> Set up the trial electric fields and compute perturbed density matrices <-- //
  std::vector<std::shared_ptr<psi::Vector>> fields;
  std::vector<std::shared_ptr<psi::Matrix>> dmats;
  for (int N=0; N<npoints; ++N) {
       fields.push_back(draw_field());
       cout << oepdev::string_sprintf(" Computation for N=%2d F=[%14.4f, %14.4f, %14.4f]\n",N+1,
                                        fields[N]->get(0), fields[N]->get(1), fields[N]->get(2));

       std::shared_ptr<psi::Matrix> dmat = std::make_shared<psi::Matrix>("",nbf,nbf);
       dmat->copy(perturbed_dmat(fields[N]));
       dmat->subtract(wfn_->Da());
       dmats.push_back(dmat);
  }

  // --> Allocate data <-- //
  std::shared_ptr<psi::Matrix> A = std::make_shared<psi::Matrix>("A Matrix", nbf, nbf);
  std::shared_ptr<psi::Matrix> C = std::make_shared<psi::Matrix>("C Matrix", nbf, nbf);
  //std::shared_ptr<psi::Matrix> g = std::make_shared<psi::Matrix>("Gradient", nbf* nbf, 1);
  //std::shared_ptr<psi::Matrix> H = std::make_shared<psi::Matrix>("Hessian Inverted", nbf*nbf, nbf*nbf);
  std::shared_ptr<psi::Matrix> X = std::make_shared<psi::Matrix>("AO Scale Matrix", nbf, nbf);
  double** Ap = A->pointer();
  double** Cp = C->pointer();
  //double** gp = g->pointer();
  //double** Hp = H->pointer();
  double** Xp = X->pointer();

  // --> Compute the least-squares matrices <-- //
  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            double va = 0.0;
            double vc = 0.0;
            for (int N=0; N<npoints; ++N) {
                 for (int o1=0; o1<nocc; ++o1) {
                      for (int z1=0; z1<3; ++z1) {
                           vc += dmats[N]->get(i,j) * par_0->susceptibility(o1, z1)->get(i,j) * fields[N]->get(z1);
                           for (int o2=0; o2<nocc; ++o2) {
                                for (int z2=0; z2<3; ++z2) {
                                     va += par_0->susceptibility(o1, z1)->get(i,j) * fields[N]->get(z1) *
                                           par_0->susceptibility(o2, z2)->get(i,j) * fields[N]->get(z2);
                                }
                           }
                      }
                 }
            }
            Ap[i][j] = va;
            Cp[i][j] = vc;
            //if (std::abs(va) > 0.0000000000001) Xp[i][j] = vc/va;
            //if (std::abs(va) > 0.0) Xp[i][j] = vc/va;
            Xp[i][j] = vc/va;
       }
  }
  A->scale( 2.0);
  C->scale(-2.0);
  X->print();

  // --> Compute Hessian <-- //
  //for (int i=0; i<nbf; ++i) {
  //for (int j=0; j<nbf; ++j) {
  //     int ij = nbf*i + j;
  //     Hp[ij][ij] = Ap[i][j];
  //}
  //}

  A->print();
  C->print();
  //H->print();


  // --> Perform the fit <-- //
  // nothing to do here

  // --> Scale the ab-initio B tensors by scale values
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> B_scaled;
  for (int o=0; o<nocc; ++o) {
       std::vector<std::shared_ptr<psi::Matrix>> Bo;
       B_scaled.push_back(Bo);
       for (int z=0; z<3; ++z) {
            B_scaled[o].push_back(par_0->susceptibility(o,z));
            for (int i=0; i<nbf; ++i) { 
            for (int j=0; j<nbf; ++j) {
                 double v = B_scaled[o][z]->get(i,j) * X->get(i,j);
                 B_scaled[o][z]->set(i,j,v);
            }
            }
       }
  }
  par->set_susceptibility(B_scaled);
 
  return par;
}
//-- TransformedMOPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::TransformedMOPolarGEFactory::TransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::TransformedMOPolarGEFactory::TransformedMOPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::TransformedMOPolarGEFactory::~TransformedMOPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::TransformedMOPolarGEFactory::compute()
{
  // Sizing
  int nsamples = options_.get_int("DMATPOL_NSAMPLES");
  int nbf = wfn_->basisset()->nbf();
  int nocc = cphfSolver_->nocc();
  double q = options_.get_double("DMATPOL_TEST_CHARGE");
  int nq   = options_.get_double("DMATPOL_NTEST_CHARGE");

  // Compute the ab-initio (unscaled) B tensors
  std::shared_ptr<oepdev::PolarGEFactory> factory_0 = std::make_shared<oepdev::PolarGEFactory>(cphfSolver_, options_);
  std::shared_ptr<oepdev::GenEffPar> par_0 = factory_0->compute();

  // Compute the scales:

  // --> Set up the trial electric fields and compute perturbed density matrices <-- //
  std::vector<std::shared_ptr<psi::Matrix>> pointCharges;
  std::vector<std::shared_ptr<psi::Matrix>> dmats;
  for (int N=0; N<nsamples; ++N) {
       std::shared_ptr<psi::Matrix> charges = std::make_shared<psi::Matrix>("Q", nq, 4);
       for (int i=0; i<nq; ++i) {
            std::shared_ptr<psi::Vector> point = draw_random_point();
            charges->set(i, 0, point->get(0));
            charges->set(i, 1, point->get(1));
            charges->set(i, 2, point->get(2));
            charges->set(i, 3, q);
       }
       pointCharges.push_back(charges);
       for (int i=0; i<nq; ++i) cout << oepdev::string_sprintf(" Computation for N=%2d X=[%14.4f, %14.4f, %14.4f]\n",N+1,
                                        pointCharges[N]->get(i,0), pointCharges[N]->get(i,1), pointCharges[N]->get(i,2));
       std::shared_ptr<psi::Matrix> dmat = std::make_shared<psi::Matrix>("",nbf,nbf);
       dmat->copy(perturbed_dmat(pointCharges[N]));
       dmat->subtract(wfn_->Da());
       dmats.push_back(dmat);
  }

  // --> Allocate data <-- //
  std::shared_ptr<psi::Matrix> A = std::make_shared<psi::Matrix>("A Matrix", nocc*nocc, nocc*nocc);
  std::shared_ptr<psi::Matrix> C = std::make_shared<psi::Matrix>("C Matrix", nocc*nocc, 1);
  std::shared_ptr<psi::Matrix> g = std::make_shared<psi::Matrix>("Gradient", nocc*nocc, 1);
  std::shared_ptr<psi::Matrix> H = std::make_shared<psi::Matrix>("Hessian Inverted", nocc*nocc, nocc*nocc);
  std::shared_ptr<psi::Matrix> x = std::make_shared<psi::Matrix>("MO Transform Vector", nocc*nocc, 1);
  double** Ap = A->pointer();
  double** Cp = C->pointer();
  double** gp = g->pointer();
  double** Hp = H->pointer();
  double** xp = x->pointer();

  // --> Compute the least-squares matrices <-- //
  for (int N=0; N<nsamples; ++N) {
  for (int a=0; a<nbf; ++a) {
       for (int b=0; b<nbf; ++b) {
            double dabN = dmats[N]->get(a,b);
            for (int i=0; i<nocc; ++i) {
                 for (int j=0; j<nocc; ++j) {
                      int ij = nocc*i + j;
                      for (int z1=0; z1<3; ++z1) {
                           double bf1 = par_0->susceptibility(i,z1)->get(a,b) *
                                            field_due_to_charges(pointCharges[N], cphfSolver_->lmo_centroid(j))->get(z1);
                           Cp[ij][0] += dabN * bf1;
                           for (int k=0; k<nocc; ++k) {
                                for (int l=0; l<nocc; ++l) {
                                     int kl = nocc*k + l;
                                     for (int z2=0; z2<3; ++z2) {
                                          Ap[ij][kl] += bf1 *
                                                        par_0->susceptibility(k,z2)->get(a,b) * 
                                            field_due_to_charges(pointCharges[N], cphfSolver_->lmo_centroid(l))->get(z2);
                                     }
                                }
                           }
                      }
                 }
            }
       }
  }
  }
  C->scale(-2.0);

  //for (int i=0; i<nocc; ++i) {
  //     for (int j=0; j<nocc; ++j) {
  //          int ij = nocc*i + j;
  //          double vc = 0.0;
  //          for (int N=0; N<nsamples; ++N) {
  //               for (int a=0; a<nbf; ++a) {
  //                    for (int b=0; b<nbf; ++b) {
  //                         double dabN = dmats[N]->get(a,b);
  //                         for (int z1=0; z1<3; ++z1) {
  //                              vc += dabN * par_0->susceptibility(i,z1)->get(a,b) *  
  //                                           field_due_to_charges(pointCharges[N], cphfSolver_->lmo_centroid(j))->get(z1);
  //                         }
  //                    }
  //               }
  //          }
  //          Cp[ij][0] = vc;
  //     }
  //}
  //C->scale(-2.0);

  //for (int i=0; i<nocc; ++i) {
  //     for (int j=0; j<nocc; ++j) {
  //          int ij = nocc*i + j;
  //          for (int k=0; k<nocc; ++k) {
  //               for (int l=0; l<nocc; ++l) {
  //                    int kl = k*nocc + l; 
  //                    double va = 0.0;
  //                    for (int N=0; N<nsamples; ++N) {
  //                    for (int a=0; a<nbf; ++a) { 
  //                    for (int b=0; b<nbf; ++b) {
  //                         double dabN = dmats[N]->get(a,b);
  //                         for (int z1=0; z1<3; ++z1) {
  //                         for (int z2=0; z2<3; ++z2) {
  //                              va += par_0->susceptibility(i,z1)->get(a,b) *
  //                                    par_0->susceptibility(k,z2)->get(a,b) *
  //                                    field_due_to_charges(pointCharges[N], cphfSolver_->lmo_centroid(j))->get(z1) *
  //                                    field_due_to_charges(pointCharges[N], cphfSolver_->lmo_centroid(l))->get(z2);
  //                         }
  //                         }
  //                    }
  //                    }
  //                    }
  //                    Ap[ij][kl] = va;
  //               }
  //          }
  //     }
  //}

  // Compute initial Z
  double Z_0 = 0.0;
  for (int N=0; N<nsamples; ++N) {
  for (int a=0; a<nbf; ++a) { 
  for (int b=0; b<nbf; ++b) {
       double dabN = dmats[N]->get(a,b);
       Z_0 += dabN * dabN;
  }}}

  double Z_init = Z_0;
  for (int i=0; i<nocc; ++i) {
       int ii = nocc*i + i;
       Z_init += Cp[ii][0];
       for (int j=0; j<nocc; ++j) {
            int jj = nocc*j + j;
            Z_init += Ap[ii][jj];
       }
  }



  //double sc = 1000.0/(A->absmax());
  //A->scale(sc);
  std::shared_ptr<psi::Matrix> Ai = std::make_shared<psi::Matrix>(*A);
  Ai->set_name("INVERSE!!!!!");
  Ai->invert();
  A->print();
  C->print();
  Ai->print();

  std::shared_ptr<psi::Matrix> S = psi::Matrix::doublet(C, Ai, true, false);
  S->scale(-0.500);
  S->set_name("SCALES!!!!!!!!!!!!!!!");
  S->print();

  // Compute final Z
  double Z = Z_0 + S->vector_dot(C);
  std::shared_ptr<psi::Matrix> u = psi::Matrix::triplet(S, A, S, false, false, true);  
  u->set_name("BUAUAUAUAUAU!!!");
  u->print();
  Z += u->get(0,0);

  psi::outfile->Printf("\n\n Z_0    = %15.9f\n", Z_0);
  psi::outfile->Printf(    " Z_init = %15.9f\n", Z_init);
  psi::outfile->Printf(    " Z_fin  = %15.9f\n\n\n", Z);



  // --> Compute Hessian and Invert <-- //
  //for (int i=0; i<nocc; ++i) {
  //for (int j=0; j<nocc; ++j) {
  //     int ij = nocc*i + j;
  //     for (int k=0; k<nocc; ++k) {
  //     for (int l=0; l<nocc; ++l) {
  //         int kl = nocc*k + l;
  //         Hp[ij][kl] = Ap[ij][kl] + Ap[kl][ij];
  //     }
  //     }
  //}
  //}
  ////H->scale(2.0); 
  //H->invert();
  //H->print();


  // --> Perform the fit <-- //
  // nothing to do here

  // --> Scale the ab-initio B tensors by scale values
  std::shared_ptr<oepdev::GenEffPar> par = std::make_shared<oepdev::GenEffPar>("Polarization");
  std::vector<std::vector<std::shared_ptr<psi::Matrix>>> B_scaled;
  for (int o=0; o<nocc; ++o) {
       std::vector<std::shared_ptr<psi::Matrix>> Bo;
       B_scaled.push_back(Bo);
       for (int z=0; z<3; ++z) {
            B_scaled[o].push_back(std::make_shared<psi::Matrix>("",nbf,nbf));
       }
  }
  for (int i=0; i<nbf; ++i) {
  for (int j=0; j<nbf; ++j) {
  for (int z=0; z<3; ++z) {
       for (int o1=0; o1<nocc; ++o1) { 
            double v = 0.0;
            for (int o2=0; o2<nocc; ++o2) {
               v += par_0->susceptibility(o2,z)->get(i,j) * S->get(0, o2*nocc+o1);
            }
            B_scaled[o1][z]->set(i,j,v);
       }
  }
  }}
  par->set_susceptibility(B_scaled);
 
  return par;
}
void oepdev::TransformedMOPolarGEFactory::gradient(std::shared_ptr<psi::Matrix> A, std::shared_ptr<psi::Matrix> C, std::shared_ptr<psi::Matrix> S)
{
}
