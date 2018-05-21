#include "gefp.h"
#include "psi4/libmints/matrix.h"

using namespace std;

oepdev::GeneralizedPolarGEFactory::GeneralizedPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::PolarGEFactory(cphf, opt),
   nBlocks_(0),
   nSites_(0),
   nParameters_(0),
   nParametersBlock_({}),
   Gradient_(std::make_shared<psi::Matrix>()),
   Hessian_(std::make_shared<psi::Matrix>()),
   Parameters_(std::make_shared<psi::Matrix>()),
   PolarizationSusceptibilities_(std::make_shared<oepdev::GenEffPar>("Polarization")),
   hasDipolePolarizability_(false),
   hasDipoleDipoleHyperpolarizability_(false),
   hasQuadrupolePolarizability_(false),
   nSamples_(options_.get_int("DMATPOL_NSAMPLES")),
   mField_(options_.get_double("DMATPOL_SCALE_1")),
   Zinit_(-1.0),
   Z_(-1.0),
   referenceStatisticalSet_({{},{},{},{},{}}),
   modelStatisticalSet_({{},{},{},{}}),
   VMatrixSet_({}),
   electricFieldSet_({}),
   electricFieldGradientSet_({}),
   electricFieldSumSet_({}),
   electricFieldGradientSumSet_({}),
   jk_(nullptr)
{
   // Allocate memory for matrix sets
   for (int n=0; n<nSamples_; ++n) {
        VMatrixSet_                                         .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.DensityMatrixSet           .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.JKMatrixSet                .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.InducedDipoleSet           .push_back(std::make_shared<psi::Vector>("", 3));
        referenceStatisticalSet_.InducedQuadrupoleSet       .push_back(std::make_shared<psi::Matrix>("", 3, 3));
        referenceStatisticalSet_.InducedInteractionEnergySet.push_back(0.0);
        modelStatisticalSet_.DensityMatrixSet               .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        modelStatisticalSet_.JKMatrixSet                    .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        modelStatisticalSet_.InducedDipoleSet               .push_back(std::make_shared<psi::Vector>("", 3));
        modelStatisticalSet_.InducedQuadrupoleSet           .push_back(std::make_shared<psi::Matrix>("", 3, 3));
        modelStatisticalSet_.InducedInteractionEnergySet    .push_back(0.0);
   }
   // Construct the JK object
   jk_ = psi::JK::build_JK(wfn_->basisset(), psi::BasisSet::zero_ao_basis_set(), options_);
   jk_->set_memory((options_.get_double("SCF_MEM_SAFETY_FACTOR")*(psi::Process::environment.get_memory() / 8L)));
   jk_->initialize();
   jk_->print_header();
}
oepdev::GeneralizedPolarGEFactory::GeneralizedPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::GeneralizedPolarGEFactory(cphf, cphf->options())
{
}
oepdev::GeneralizedPolarGEFactory::~GeneralizedPolarGEFactory()
{
   //if (!options_.get_bool("SAVE_JK")) {
   //    jk_.reset();
   //}
}
std::shared_ptr<oepdev::GenEffPar> oepdev::GeneralizedPolarGEFactory::compute(void)
{
   compute_samples();
   compute_parameters();
   compute_statistics();
   return PolarizationSusceptibilities_; 
}
// protected methods
void oepdev::GeneralizedPolarGEFactory::invert_hessian(void)
{
   std::shared_ptr<psi::Matrix> Htemp= std::make_shared<psi::Matrix>(Hessian_);
   Hessian_->invert();
   Hessian_->set_name("\nInverse Hessian");
   Hessian_->print();
   // Perform the identity test
   std::shared_ptr<psi::Matrix> I = psi::Matrix::doublet(Hessian_, Htemp, false, false);
   double s = 0.0;
   for (int a=0; a<nParameters_; ++a) {
        for (int b=0; b<nParameters_; ++b) {
             s += sqrt( I->get(a,b) * I->get(a,b) );
        }
   }
   s /= (double)(nParameters_);
   cout << " Hessian Matrix Identity Test = " << s << endl;
   if (std::abs(s-1.0)>0.000001) cout << " ----> Warning!! Hessian inverse has numerical error! <----\n";
   else {cout << " Hessian Matrix Identity Test: Passed (threshold: 0.000001) ! :)\n";}
}
void oepdev::GeneralizedPolarGEFactory::compute_parameters(void)
{
   if (hasDipoleDipoleHyperpolarizability_) compute_electric_field_sums();
   if (hasQuadrupolePolarizability_)        compute_electric_field_gradient_sums();
   compute_hessian(); 
   Hessian_->print();
   invert_hessian();
   for (int i=0; i<nbf_; ++i) {
        for (int j=0; j<nbf_; ++j) {
             compute_gradient(i, j);
             fit();
             save(i, j);
        }
   }
}
void oepdev::GeneralizedPolarGEFactory::fit(void)
{
   Parameters_ = psi::Matrix::doublet(Hessian_, Gradient_, false, false);
   Parameters_->scale(-1);
}
void oepdev::GeneralizedPolarGEFactory::save(int i, int j)
{
  if (!hasDipolePolarizability_) throw psi::PSIEXCEPTION("Dipole polarizability must be set!");

  // Un-Pack the parameters from Parameters_ vector into dipole and quadrupole (hyper)polarizabilities
  if (hasDipolePolarizability_ and !hasDipoleDipoleHyperpolarizability_ and !hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) {
           for (int z=0; z<3; ++z) {
                double val = Parameters_->get(3*n + z, 0);
                PolarizationSusceptibilities_->dipole_polarizability(n, z)->set(i, j, val);
           }
      }
  } else if (hasDipolePolarizability_ and hasDipoleDipoleHyperpolarizability_ and !hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           for (int z1=0; z1<3; ++z1) {
                double val = Parameters_->get(3*n + z1, 0);
                PolarizationSusceptibilities_->dipole_polarizability(n, z1)->set(i, j, val);

                int iz1 = nSites_ * 3 + 3 * n + z1; // Second block
                for (int z2 = 0; z2<3; ++z2) {
                     int iz2 = nSites_ * 3 + 3 * n + z2; // Second block

                     val = Parameters_->get(iz1, 0) + Parameters_->get(iz2, 0);
                     PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, z1*3+z2)->set(i, j, val*mField_);
                }
           }
      }
  } else if (hasDipolePolarizability_ and !hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           for (int z1=0; z1<3; ++z1) {
                double val = Parameters_->get(3*n + z1, 0);
                PolarizationSusceptibilities_->dipole_polarizability(n, z1)->set(i, j, val);

                int iz1 = nSites_ * 3 + 3 * n + z1; // Second block
                for (int z2 = 0; z2<3; ++z2) {
                     int iz2 = nSites_ * 3 + 3 * n + z2; // Second block

                     val = Parameters_->get(iz1, 0) + Parameters_->get(iz2, 0);
                     PolarizationSusceptibilities_->quadrupole_polarizability(n, z1*3+z2)->set(i, j, val);
                }
           }
      }
  } else if (hasDipolePolarizability_ and hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           for (int z1=0; z1<3; ++z1) {
                double val = Parameters_->get(3*n + z1, 0);
                PolarizationSusceptibilities_->dipole_polarizability(n, z1)->set(i, j, val);

                int i1z1 = nSites_ * 3 + 3 * n + z1; // Second block
                int i2z1 = nSites_ * 6 + 3 * n + z1; // Third block
                for (int z2 = 0; z2<3; ++z2) {
                     int i1z2 = nSites_ * 3 + 3 * n + z2; // Second block
                     int i2z2 = nSites_ * 6 + 3 * n + z2; // Third block

                     val = Parameters_->get(i1z1, 0) + Parameters_->get(i1z2, 0);
                     PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, z1*3+z2)->set(i, j, val*mField_);
                     val = Parameters_->get(i2z1, 0) + Parameters_->get(i2z2, 0);
                     PolarizationSusceptibilities_->quadrupole_polarizability(n, z1*3+z2)->set(i, j, val);
                }
           }
      }
  } else {
    throw psi::PSIEXCEPTION("The model you refer to is not implemented.");
  }
  
}
void oepdev::GeneralizedPolarGEFactory::allocate(void)
{
  // Blocks
  nBlocks_  = (int)hasDipolePolarizability_ + (int)hasDipoleDipoleHyperpolarizability_ + (int)hasQuadrupolePolarizability_;
  nParameters_ = nBlocks_*(nSites_ * 3);
  for (int z=0; z<nBlocks_; ++z) nParametersBlock_.push_back(nSites_ * 3);
  // Parameter spaces
  Gradient_   = std::make_shared<psi::Matrix>("Gradient"  , nParameters_, 1);
  Hessian_    = std::make_shared<psi::Matrix>("Hessian"   , nParameters_, nParameters_);
  Parameters_ = std::make_shared<psi::Matrix>("Parameters", nParameters_, 1);
  // Susceptibilities
  if (hasDipolePolarizability_           ) PolarizationSusceptibilities_->allocate(1, 0, nSites_, nbf_);
  if (hasDipoleDipoleHyperpolarizability_) PolarizationSusceptibilities_->allocate(2, 0, nSites_, nbf_);
  if (hasQuadrupolePolarizability_       ) PolarizationSusceptibilities_->allocate(0, 1, nSites_, nbf_);
}
void oepdev::GeneralizedPolarGEFactory::compute_electric_field_sums(void) {
  for (int n=0; n<nSamples_; ++n) {
     std::vector<double> sum;
     for (int i=0; i<nSites_; ++i) {
          std::shared_ptr<psi::Vector> field = electricFieldSet_[n][i];
          sum.push_back(field->get(0) + field->get(1) + field->get(2));
     }
     electricFieldSumSet_.push_back(sum);
  }
}
void oepdev::GeneralizedPolarGEFactory::compute_statistics(void) {
   cout << " Statistical evaluation ...\n";

   // Grab some unperturbed quantities
   std::shared_ptr<psi::Matrix> Hcore  = wfn_->H();
   std::shared_ptr<psi::Matrix> Fock   = wfn_->Fa();
   std::shared_ptr<psi::Matrix> D      = wfn_->Da();

   // Compute model difference density matrices
   if (hasQuadrupolePolarizability_) {
      for (int n=0; n<nSamples_; ++n) {                                                                          
           modelStatisticalSet_.DensityMatrixSet[n]->copy(PolarizationSusceptibilities_->compute_density_matrix(
                     electricFieldSet_[n], electricFieldGradientSet_[n]));
      }
   } else {
      for (int n=0; n<nSamples_; ++n) {                                                                          
           modelStatisticalSet_.DensityMatrixSet[n]->copy(PolarizationSusceptibilities_->compute_density_matrix(
                     electricFieldSet_[n]));
      }
   }

   // Compute changes in the Coulomb and exchange HF matrices
   std::vector<psi::SharedMatrix>& C_left = jk_->C_left();
   const std::vector<std::shared_ptr<psi::Matrix>>& J = jk_->J();
   const std::vector<std::shared_ptr<psi::Matrix>>& K = jk_->K();

   C_left.clear();
   for (int n=0; n<nSamples_; ++n) {
        std::shared_ptr<psi::Matrix> C = referenceStatisticalSet_.DensityMatrixSet[n]->partial_cholesky_factorize();
        C_left.push_back(C);
   }
   jk_->compute();
   for (int n=0; n<nSamples_; ++n) {
        referenceStatisticalSet_.JKMatrixSet[n]->add(J[n]);
        referenceStatisticalSet_.JKMatrixSet[n]->add(K[n]);
   }
   C_left.clear();
   for (int n=0; n<nSamples_; ++n) {
       std::shared_ptr<psi::Matrix> C = modelStatisticalSet_.DensityMatrixSet[n]->partial_cholesky_factorize();
       C_left.push_back(C);
   }
   jk_->compute();
   for (int n=0; n<nSamples_; ++n) {
        modelStatisticalSet_.JKMatrixSet[n]->add(J[n]);
        modelStatisticalSet_.JKMatrixSet[n]->add(K[n]);
   }

   for (int n=0; n<nSamples_; ++n) {

        // ---> Compute interaction energy <--- //

        double e_nuc = 0.0; // TODO: compute this 
        double eint_refer = 2.0 * D->vector_dot(VMatrixSet_[n]) + e_nuc;
        double eint_model = eint_refer;

        std::shared_ptr<psi::Matrix> A = Hcore->clone(); 
        A->add(VMatrixSet_[n]);
        A->add(VMatrixSet_[n]);
        A->add(Fock);
        eint_refer += referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(A);
        eint_model +=     modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(A);
 
        std::shared_ptr<psi::Matrix> Dr = D->clone(); Dr->add(referenceStatisticalSet_.DensityMatrixSet[n]);
        std::shared_ptr<psi::Matrix> Dm = D->clone(); Dm->add(    modelStatisticalSet_.DensityMatrixSet[n]);
        eint_refer += Dr->vector_dot(referenceStatisticalSet_.JKMatrixSet[n]);
        eint_model += Dm->vector_dot(    modelStatisticalSet_.JKMatrixSet[n]);
        //std::shared_ptr<psi::Matrix> Cr = referenceStatisticalSet_.DensityMatrixSet[n]->partial_cholesky_factorize();
        //std::shared_ptr<psi::Matrix> Cm =     modelStatisticalSet_.DensityMatrixSet[n]->partial_cholesky_factorize();
 
        referenceStatisticalSet_.InducedInteractionEnergySet[n] = eint_refer;
            modelStatisticalSet_.InducedInteractionEnergySet[n] = eint_model;
        cout << eint_refer << " " << eint_model << endl;

        // ---> Compute induced dipole moments <--- //
        // TODO

   }
   // Finish with JK object
   //jk_->finalize();
}
void oepdev::GeneralizedPolarGEFactory::compute_electric_field_gradient_sums(void) {
  // TODO
}
// abstract methods
void oepdev::GeneralizedPolarGEFactory::compute_samples(void)
{
}
void oepdev::GeneralizedPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::GeneralizedPolarGEFactory::compute_hessian(void)
{
}
