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
   referenceDensityMatrixSet_({}),
   guessDensityMatrixSet_({}),
   electricFieldSet_({}),
   electricFieldGradientSet_({}),
   electricFieldSumSet_({}),
   electricFieldGradientSumSet_({})
{
   // Allocate memory for density matrix sets
   for (int n=0; n<nSamples_; ++n) {
        referenceDensityMatrixSet_.push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        guessDensityMatrixSet_    .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
   }
}
oepdev::GeneralizedPolarGEFactory::GeneralizedPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::GeneralizedPolarGEFactory(cphf, cphf->options())
{
}
oepdev::GeneralizedPolarGEFactory::~GeneralizedPolarGEFactory()
{
}
std::shared_ptr<oepdev::GenEffPar> oepdev::GeneralizedPolarGEFactory::compute(void)
{
   compute_samples();
   compute_parameters();
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
