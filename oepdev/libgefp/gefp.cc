#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

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
 options_(opt)
{

}
oepdev::GenEffParFactory::~GenEffParFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::GenEffParFactory::compute()
{

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
            Biz_t->set_name(oepdev::string_sprintf("B(%d)[%d] in Non-Orthogonal AO Basis", o+1, z+1));
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
//-- ScaledPolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::ScaledPolarGEFactory::ScaledPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt) :
 oepdev::PolarGEFactory(cphf, opt)
{

}
oepdev::ScaledPolarGEFactory::ScaledPolarGEFactory(std::shared_ptr<CPHF> cphf) :
 oepdev::PolarGEFactory(cphf)
{

}
oepdev::ScaledPolarGEFactory::~ScaledPolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::ScaledPolarGEFactory::compute()
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
std::shared_ptr<psi::Vector> oepdev::ScaledPolarGEFactory::draw_field()
{
  std::shared_ptr<psi::Vector> field = std::make_shared<psi::Vector>("", 3);
  //TODO Add here random pick up of the field values!
  return field;
}
std::shared_ptr<psi::Matrix> oepdev::ScaledPolarGEFactory::perturbed_dmat(const std::shared_ptr<psi::Vector>& field)
{
  std::shared_ptr<psi::Matrix> dmat;
  // TODO
  return dmat;
}

