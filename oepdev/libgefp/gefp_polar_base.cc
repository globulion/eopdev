#include "gefp.h"
#include <iostream>
#include <fstream>

using namespace std;

//-- PolarGEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::PolarGEFactory::PolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::GenEffParFactory(wfn, opt)
{

}
oepdev::PolarGEFactory::~PolarGEFactory()
{

}
std::shared_ptr<oepdev::GenEffPar> oepdev::PolarGEFactory::compute()
{

}
double oepdev::PolarGEFactory::draw_charge()
{
  const double scale = options_.get_double("DMATPOL_TEST_CHARGE");
  double q = random_double() * scale;
  return q;
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
std::shared_ptr<psi::Matrix> oepdev::PolarGEFactory::field_gradient_due_to_charges(
                                                                  const std::shared_ptr<psi::Matrix>& charges,
                                                                  const double& x, const double& y, const double& z)
{
  std::shared_ptr<psi::Matrix> gradient = std::make_shared<psi::Matrix>("", 3, 3);
  double gxx = 0.0; double gyx = 0.0; double gzx = 0.0;
  double gxy = 0.0; double gyy = 0.0; double gzy = 0.0;
  double gxz = 0.0; double gyz = 0.0; double gzz = 0.0;
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
       double r5 = 3.0*r3/(r*r);
       gxx += r3 - r5 * drx * drx;
       gxy -=      r5 * drx * dry;
       gxz -=      r5 * drx * drz;
       gyx -=      r5 * dry * drx;
       gyy += r3 - r5 * dry * dry;
       gyz -=      r5 * dry * drz;
       gzx -=      r5 * drz * drx;
       gzy -=      r5 * drz * dry;
       gzz += r3 - r5 * drz * drz;
  }
  gradient->set(0, 0, gxx);
  gradient->set(0, 1, gxy);
  gradient->set(0, 2, gxz);
  gradient->set(1, 0, gyx);
  gradient->set(1, 1, gyy);
  gradient->set(1, 2, gyz);
  gradient->set(2, 0, gzx);
  gradient->set(2, 1, gzy);
  gradient->set(2, 2, gzz);
  return gradient;
}
std::shared_ptr<psi::Vector> oepdev::PolarGEFactory::field_due_to_charges(const std::shared_ptr<psi::Matrix>& charges,
                                                                          const std::shared_ptr<psi::Vector>& pos)
{
  return this->field_due_to_charges(charges, pos->get(0), pos->get(1), pos->get(2));
}
std::shared_ptr<psi::Matrix> oepdev::PolarGEFactory::field_gradient_due_to_charges(
                                                                          const std::shared_ptr<psi::Matrix>& charges,
                                                                          const std::shared_ptr<psi::Vector>& pos)
{
  return this->field_gradient_due_to_charges(charges, pos->get(0), pos->get(1), pos->get(2));
}
std::shared_ptr<oepdev::RHFPerturbed> oepdev::PolarGEFactory::perturbed_state(const std::shared_ptr<psi::Vector>& field)
{
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, 
                  oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());
  scf->set_perturbation(field);
  scf->compute_energy();
  return scf;
}
std::shared_ptr<oepdev::RHFPerturbed> oepdev::PolarGEFactory::perturbed_state(const std::shared_ptr<psi::Vector>& pos, const double& q)
{
  std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn_, 
                  oepdev::create_superfunctional("HF", options_), options_, wfn_->psio());
  scf->set_perturbation(pos, q);
  scf->compute_energy();
  return scf;
}
std::shared_ptr<oepdev::RHFPerturbed> oepdev::PolarGEFactory::perturbed_state(const std::shared_ptr<psi::Matrix>& charges)
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
  return scf;
}
// -- GeneralizedPolarGEFactory ---------------------------------------------------------------------------------------
oepdev::GeneralizedPolarGEFactory::GeneralizedPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::PolarGEFactory(wfn, opt),
   nBlocks_(0),
   nSites_(0),
   nSitesAbInitio_(0),
   nParameters_(0),
   nParametersBlock_({}),
   Gradient_(std::make_shared<psi::Matrix>()),
   Hessian_(std::make_shared<psi::Matrix>()),
   Parameters_(std::make_shared<psi::Matrix>()),
   PolarizationSusceptibilities_(std::make_shared<oepdev::GenEffPar>("Polarization")),
   abInitioPolarizationSusceptibilities_(std::make_shared<oepdev::GenEffPar>("Polarization")),
   hasDipolePolarizability_(false),
   hasDipoleDipoleHyperpolarizability_(false),
   hasQuadrupolePolarizability_(false),
   hasAbInitioDipolePolarizability_(false),
   nSamples_(options_.get_int("DMATPOL_NSAMPLES")),
   mField_(options_.get_double("DMATPOL_SCALE_1")),
   Zinit_(-1.0),
   Z_(-1.0),
   symmetryNumber_{1.0, 2.0, 2.0, 1.0, 2.0, 1.0}, /* XX, XY, XZ, YY, YZ, ZZ */
   referenceStatisticalSet_({{},{},{},{},{}}),
   referenceDpolStatisticalSet_({{},{},{},{},{}}),
   modelStatisticalSet_({{},{},{},{},{}}),
   abInitioModelStatisticalSet_({{},{},{},{},{}}),
   VMatrixSet_({}),
   electricFieldSet_({}),
   electricFieldGradientSet_({}),
   electricFieldSumSet_({}),
   electricFieldGradientSumSet_({}),
   abInitioModelElectricFieldSet_({}),
   jk_(nullptr)
{
   // Allocate memory for matrix sets
   for (int n=0; n<nSamples_; ++n) {
        VMatrixSet_                                         .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.DensityMatrixSet           .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.JKMatrixSet                .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        referenceStatisticalSet_.InducedDipoleSet           .push_back(std::make_shared<psi::Vector>("", 3));
        referenceStatisticalSet_.InducedQuadrupoleSet       .push_back(std::make_shared<psi::Vector>("", 6));
        referenceStatisticalSet_.InducedInteractionEnergySet.push_back(0.0);

        modelStatisticalSet_.DensityMatrixSet               .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        modelStatisticalSet_.JKMatrixSet                    .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
        modelStatisticalSet_.InducedDipoleSet               .push_back(std::make_shared<psi::Vector>("", 3));
        modelStatisticalSet_.InducedQuadrupoleSet           .push_back(std::make_shared<psi::Vector>("", 6));
        modelStatisticalSet_.InducedInteractionEnergySet    .push_back(0.0);

        //if (hasAbInitioDipolePolarizability_) {

           referenceDpolStatisticalSet_.InducedDipoleSet           .push_back(std::make_shared<psi::Vector>("", 3)); 
           referenceDpolStatisticalSet_.InducedQuadrupoleSet       .push_back(std::make_shared<psi::Vector>("", 6));
           referenceDpolStatisticalSet_.InducedInteractionEnergySet.push_back(0.0);

           abInitioModelStatisticalSet_.DensityMatrixSet           .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_)); 
           abInitioModelStatisticalSet_.JKMatrixSet                .push_back(std::make_shared<psi::Matrix>("", nbf_, nbf_));
           abInitioModelStatisticalSet_.InducedDipoleSet           .push_back(std::make_shared<psi::Vector>("", 3));
           abInitioModelStatisticalSet_.InducedQuadrupoleSet       .push_back(std::make_shared<psi::Vector>("", 6));
           abInitioModelStatisticalSet_.InducedInteractionEnergySet.push_back(0.0);
        //}
   }
   // Construct the JK object
   jk_ = psi::JK::build_JK(wfn_->basisset(), psi::BasisSet::zero_ao_basis_set(), options_);
   jk_->set_memory((options_.get_double("SCF_MEM_SAFETY_FACTOR")*(psi::Process::environment.get_memory() / 8L)));
   jk_->initialize();
   jk_->print_header();
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
   set_distributed_centres();
}
void oepdev::GeneralizedPolarGEFactory::set_distributed_centres(void)
{
   // Centres are assumed to be atoms
   std::vector<std::shared_ptr<psi::Vector>> centres;
   for (int s=0; s<nSites_; ++s) {
        centres.push_back(std::make_shared<psi::Vector>(oepdev::string_sprintf("Atomic centre (%d)", s+1), 3));
        centres[s]->set(0, wfn_->molecule()->x(s));
        centres[s]->set(1, wfn_->molecule()->y(s));
        centres[s]->set(2, wfn_->molecule()->z(s));
   }
   PolarizationSusceptibilities_->set_centres(centres);
}
void oepdev::GeneralizedPolarGEFactory::fit(void)
{
   Parameters_ = psi::Matrix::doublet(Hessian_, Gradient_, false, false);
   Parameters_->scale(-1);
}
void oepdev::GeneralizedPolarGEFactory::save(int i, int j)
{
  if (!hasDipolePolarizability_) throw psi::PSIEXCEPTION("Dipole polarizability must be set!");

  // --> Un-Pack the parameters from Parameters_ vector into dipole and quadrupole (hyper)polarizabilities <-- //

  // Dipole Polarizability (always)
  for (int n=0; n<nSites_; ++n) {
       for (int z=0; z<3; ++z) {
            double val = Parameters_->get(3*n + z, 0);
            PolarizationSusceptibilities_->dipole_polarizability(n, z)->set(i, j, val);
       }
  }
  // Only Dipole-Dipole Hyperpolarizability
  if (hasDipoleDipoleHyperpolarizability_ and !hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           double val_xx = Parameters_->get(nSites_*3 + 6*n + 0, 0);
           double val_xy = Parameters_->get(nSites_*3 + 6*n + 1, 0);
           double val_xz = Parameters_->get(nSites_*3 + 6*n + 2, 0);
           double val_yy = Parameters_->get(nSites_*3 + 6*n + 3, 0);
           double val_yz = Parameters_->get(nSites_*3 + 6*n + 4, 0);
           double val_zz = Parameters_->get(nSites_*3 + 6*n + 5, 0);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 0)->set(i, j, val_xx);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 1)->set(i, j, val_xy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 2)->set(i, j, val_xz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 3)->set(i, j, val_xy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 4)->set(i, j, val_yy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 5)->set(i, j, val_yz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 6)->set(i, j, val_xz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 7)->set(i, j, val_yz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 8)->set(i, j, val_zz);
      }
  // Only Quadrupole Polarizability
  } else if (!hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           double val_xx = Parameters_->get(nSites_*3 + 6*n + 0, 0);
           double val_xy = Parameters_->get(nSites_*3 + 6*n + 1, 0);
           double val_xz = Parameters_->get(nSites_*3 + 6*n + 2, 0);
           double val_yy = Parameters_->get(nSites_*3 + 6*n + 3, 0);
           double val_yz = Parameters_->get(nSites_*3 + 6*n + 4, 0);
           double val_zz = Parameters_->get(nSites_*3 + 6*n + 5, 0);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 0)->set(i, j, val_xx);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 1)->set(i, j, val_xy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 2)->set(i, j, val_xz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 3)->set(i, j, val_xy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 4)->set(i, j, val_yy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 5)->set(i, j, val_yz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 6)->set(i, j, val_xz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 7)->set(i, j, val_yz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 8)->set(i, j, val_zz);
     }
  // Both Dipole-Dipole Hyperpolarizability and Quadrupole Polarizability
  } else if (hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) { 
           double val_xx = Parameters_->get(nSites_*3 + 6*n + 0, 0);
           double val_xy = Parameters_->get(nSites_*3 + 6*n + 1, 0);
           double val_xz = Parameters_->get(nSites_*3 + 6*n + 2, 0);
           double val_yy = Parameters_->get(nSites_*3 + 6*n + 3, 0);
           double val_yz = Parameters_->get(nSites_*3 + 6*n + 4, 0);
           double val_zz = Parameters_->get(nSites_*3 + 6*n + 5, 0);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 0)->set(i, j, val_xx);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 1)->set(i, j, val_xy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 2)->set(i, j, val_xz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 3)->set(i, j, val_xy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 4)->set(i, j, val_yy);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 5)->set(i, j, val_yz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 6)->set(i, j, val_xz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 7)->set(i, j, val_yz);
           PolarizationSusceptibilities_->dipole_dipole_hyperpolarizability(n, 8)->set(i, j, val_zz);
           val_xx = Parameters_->get(nSites_*9 + 6*n + 0, 0);
           val_xy = Parameters_->get(nSites_*9 + 6*n + 1, 0);
           val_xz = Parameters_->get(nSites_*9 + 6*n + 2, 0);
           val_yy = Parameters_->get(nSites_*9 + 6*n + 3, 0);
           val_yz = Parameters_->get(nSites_*9 + 6*n + 4, 0);
           val_zz = Parameters_->get(nSites_*9 + 6*n + 5, 0);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 0)->set(i, j, val_xx);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 1)->set(i, j, val_xy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 2)->set(i, j, val_xz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 3)->set(i, j, val_xy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 4)->set(i, j, val_yy);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 5)->set(i, j, val_yz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 6)->set(i, j, val_xz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 7)->set(i, j, val_yz);
           PolarizationSusceptibilities_->quadrupole_polarizability(n, 8)->set(i, j, val_zz);
      }
  } //else {
  //  throw psi::PSIEXCEPTION("The model you refer to is not implemented.");
  //}
}
void oepdev::GeneralizedPolarGEFactory::allocate(void)
{
  // Blocks
  nBlocks_  = (int)hasDipolePolarizability_ + (int)hasDipoleDipoleHyperpolarizability_ + (int)hasQuadrupolePolarizability_;
  nParametersBlock_.push_back(nSites_ * 3);
  if (hasDipoleDipoleHyperpolarizability_) nParametersBlock_.push_back(nSites_ * 6);
  if (hasQuadrupolePolarizability_) nParametersBlock_.push_back(nSites_ * 6);
  for (int b=0; b<nParametersBlock_.size(); ++b) nParameters_ += nParametersBlock_[b];

  // Parameter spaces
  Gradient_   = std::make_shared<psi::Matrix>("Gradient"  , nParameters_, 1);
  Hessian_    = std::make_shared<psi::Matrix>("Hessian"   , nParameters_, nParameters_);
  Parameters_ = std::make_shared<psi::Matrix>("Parameters", nParameters_, 1);

  // Ab initio model
  hasAbInitioDipolePolarizability_ = options_.get_bool("DMATPOL_DO_AB_INITIO");

  // Susceptibilities
  if (hasDipolePolarizability_           ) PolarizationSusceptibilities_->allocate(1, 0, nSites_, nbf_);
  if (hasDipoleDipoleHyperpolarizability_) PolarizationSusceptibilities_->allocate(2, 0, nSites_, nbf_);
  if (hasQuadrupolePolarizability_       ) PolarizationSusceptibilities_->allocate(0, 1, nSites_, nbf_);
  if (hasAbInitioDipolePolarizability_   ) {
     abInitioPolarizationSusceptibilities_->allocate(1, 0, nSitesAbInitio_, nbf_);
     if ((options_.get_int("DMATPOL_FIELD_RANK")==2) && (options_.get_bool("DMATPOL_FF_AB_INITIO")))
     abInitioPolarizationSusceptibilities_->allocate(2, 0, nSitesAbInitio_, nbf_);
  }

  // Ab Initio Model
  if (hasAbInitioDipolePolarizability_   ) compute_ab_initio();
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
  for (int n=0; n<nSamples_; ++n) {
     std::vector<std::shared_ptr<psi::Vector>> sum;
     for (int i=0; i<nSites_; ++i) {
          std::shared_ptr<psi::Matrix> gradient = electricFieldGradientSet_[n][i];
          std::shared_ptr<psi::Vector> sum_i = std::make_shared<psi::Vector>("", 3);
          double sx = gradient->get(0, 0) + gradient->get(1, 0) + gradient->get(2, 0);
          double sy = gradient->get(0, 1) + gradient->get(1, 1) + gradient->get(2, 1);
          double sz = gradient->get(0, 2) + gradient->get(1, 2) + gradient->get(2, 2);
          sum_i->set(0, sx);
          sum_i->set(1, sy);
          sum_i->set(2, sz);
          sum.push_back(sum_i);
     }
     electricFieldGradientSumSet_.push_back(sum);
  }
}
void oepdev::GeneralizedPolarGEFactory::compute_ab_initio(void) {
  std::shared_ptr<oepdev::GenEffParFactory> factory;
  if (options_.get_bool("DMATPOL_FF_AB_INITIO") == true) {
      factory = std::make_shared<oepdev::FFAbInitioPolarGEFactory>(wfn_, options_);
  } else {
      factory = std::make_shared<oepdev::AbInitioPolarGEFactory>(wfn_, options_);
  }
  std::shared_ptr<oepdev::GenEffPar> parameters = factory->compute();
  abInitioPolarizationSusceptibilities_->set_dipole_polarizability(parameters->dipole_polarizability());
  if ((options_.get_int("DMATPOL_FIELD_RANK")==2) && options_.get_bool("DMATPOL_FF_AB_INITIO"))
  abInitioPolarizationSusceptibilities_->set_dipole_dipole_hyperpolarizability(parameters->dipole_dipole_hyperpolarizability());
  abInitioPolarizationSusceptibilities_->set_centres(parameters->centres());
  abInitioPolarizationSusceptibilitiesFactory_ = factory;
}
void oepdev::GeneralizedPolarGEFactory::compute_statistics(void) {
   cout << " Statistical evaluation ...\n";

   // ---> Grab some unperturbed quantities <--- //
   std::shared_ptr<psi::Matrix> H = wfn_->Fa()->clone(); H->add(wfn_->H());
   std::shared_ptr<psi::Matrix> D = wfn_->Da()->clone();

   // ---> Compute model difference density matrices <--- //
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
   if (hasAbInitioDipolePolarizability_) {
      for (int n=0; n<nSamples_; ++n) {
           abInitioModelStatisticalSet_.DensityMatrixSet[n]->copy(abInitioPolarizationSusceptibilities_->compute_density_matrix(abInitioModelElectricFieldSet_[n]));
      }
   }

   // ---> Compute changes in the Coulomb and exchange HF matrices <--- //
   std::vector<psi::SharedMatrix>& C_left = jk_->C_left();
   std::vector<psi::SharedMatrix>& C_right= jk_->C_right();
   const std::vector<std::shared_ptr<psi::Matrix>>& J = jk_->J();
   const std::vector<std::shared_ptr<psi::Matrix>>& K = jk_->K();

   C_left.clear(); C_right.clear();
   for (int n=0; n<nSamples_; ++n) {
        std::shared_ptr<psi::Matrix> Cl = referenceStatisticalSet_.DensityMatrixSet[n]->clone();
        std::shared_ptr<psi::Matrix> Cr= Cl->clone(); Cr->identity();
        C_left.push_back(Cl);
        C_right.push_back(Cr);
   }
   jk_->compute();
   for (int n=0; n<nSamples_; ++n) {
        referenceStatisticalSet_.JKMatrixSet[n]->axpy(2.0,J[n]);
        referenceStatisticalSet_.JKMatrixSet[n]->subtract(K[n]);
   }
   C_left.clear(); C_right.clear();
   for (int n=0; n<nSamples_; ++n) {
        std::shared_ptr<psi::Matrix> Cl = modelStatisticalSet_.DensityMatrixSet[n]->clone();
        std::shared_ptr<psi::Matrix> Cr= Cl->clone(); Cr->identity();
        C_left.push_back(Cl);
        C_right.push_back(Cr);
   }
   jk_->compute();
   for (int n=0; n<nSamples_; ++n) {
        modelStatisticalSet_.JKMatrixSet[n]->axpy(2.0,J[n]);
        modelStatisticalSet_.JKMatrixSet[n]->subtract(K[n]);
   }
   if (hasAbInitioDipolePolarizability_) {
       C_left.clear(); C_right.clear();                                                          
       for (int n=0; n<nSamples_; ++n) {
            std::shared_ptr<psi::Matrix> Cl = abInitioModelStatisticalSet_.DensityMatrixSet[n]->clone();
            std::shared_ptr<psi::Matrix> Cr= Cl->clone(); Cr->identity();
            C_left.push_back(Cl);
            C_right.push_back(Cr);
       }
       jk_->compute();
       for (int n=0; n<nSamples_; ++n) {
            abInitioModelStatisticalSet_.JKMatrixSet[n]->axpy(2.0,J[n]);
            abInitioModelStatisticalSet_.JKMatrixSet[n]->subtract(K[n]);
       }
   }

   // ---> Compute the dipole and the quadrupole integrals <--- //
   std::vector<std::shared_ptr<psi::Matrix>> dipInts, qadInts;
   for (int z=0; z<3; ++z) dipInts.push_back(std::make_shared<psi::Matrix>("",wfn_->basisset()->nbf(), wfn_->basisset()->nbf()));
   for (int z=0; z<6; ++z) qadInts.push_back(std::make_shared<psi::Matrix>("",wfn_->basisset()->nbf(), wfn_->basisset()->nbf()));
   std::shared_ptr<psi::OneBodyAOInt> dipInt(wfn_->integral()->ao_dipole());
   std::shared_ptr<psi::OneBodyAOInt> qadInt(wfn_->integral()->ao_quadrupole());
   dipInt->set_origin(wfn_->molecule()->center_of_mass()); dipInt->compute(dipInts);
   qadInt->set_origin(wfn_->molecule()->center_of_mass()); qadInt->compute(qadInts);


   // ---> Printer(s) <--- //
   //std::filebuf fb, fb2;
   //fb.open(options_.get_str("DMATPOL_OUT_STATS"), std::ios::out);
   //std::shared_ptr<std::ostream> os = std::make_shared<std::ostream>(&fb);
   //std::shared_ptr<psi::PsiOutStream> printer = std::make_shared<psi::PsiOutStream>(os), printer_abini;
   std::shared_ptr<psi::PsiOutStream> printer = std::make_shared<psi::PsiOutStream>(options_.get_str("DMATPOL_OUT_STATS")); 
   std::shared_ptr<psi::PsiOutStream> printer_abini;
   printer->Printf("# Eint(Ref)    Eint(Model)\n");
   if (hasAbInitioDipolePolarizability_) {
       //fb2.open(options_.get_str("DMATPOL_OUT_STATS_AB_INITIO"), std::ios::out);
       //std::shared_ptr<std::ostream> os2 = std::make_shared<std::ostream>(&fb2);               
       //printer_abini = std::make_shared<psi::PsiOutStream>(os2);
       printer_abini = std::make_shared<psi::PsiOutStream>(options_.get_str("DMATPOL_OUT_STATS_AB_INITIO"));
       printer_abini->Printf("# Eint(Ref)    Eint(Model)\n");
   }

   for (int n=0; n<nSamples_; ++n) {

        // ---> Compute interaction energy <--- //

        double tr_rest_r = referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(H) + 
                           referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(referenceStatisticalSet_.JKMatrixSet[n]) +
                           referenceStatisticalSet_.JKMatrixSet[n]->vector_dot(D);
        double tr_rest_m = modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(H) + 
                           modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(modelStatisticalSet_.JKMatrixSet[n]) +
                           modelStatisticalSet_.JKMatrixSet[n]->vector_dot(D);
        double tr_2v_r = 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(VMatrixSet_[n]);
        double tr_2v_m = 2.0 *     modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(VMatrixSet_[n]);

            referenceStatisticalSet_.InducedInteractionEnergySet[n] = tr_2v_r + tr_rest_r;
                modelStatisticalSet_.InducedInteractionEnergySet[n] = tr_2v_m + tr_rest_m;


        // ---> Compute induced dipole and quadrupole moments <--- //
        double dr_av, dm_av, da_av;
        double qr_av, qm_av, qa_av;

        double dr_x = 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[0]);
        double dr_y = 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[1]);
        double dr_z = 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[2]);
        double qr_xx= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[0]);
        double qr_xy= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[1]);
        double qr_xz= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[2]);
        double qr_yy= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[3]);
        double qr_yz= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[4]);
        double qr_zz= 2.0 * referenceStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[5]);

        double dm_x = 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[0]);
        double dm_y = 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[1]);
        double dm_z = 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[2]);
        double qm_xx= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[0]);
        double qm_xy= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[1]);
        double qm_xz= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[2]);
        double qm_yy= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[3]);
        double qm_yz= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[4]);
        double qm_zz= 2.0 * modelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[5]);

        referenceStatisticalSet_.InducedDipoleSet    [n]->set(0, dr_x);
        referenceStatisticalSet_.InducedDipoleSet    [n]->set(1, dr_y);
        referenceStatisticalSet_.InducedDipoleSet    [n]->set(2, dr_z);
        modelStatisticalSet_.InducedDipoleSet    [n]->set(0, dm_x);
        modelStatisticalSet_.InducedDipoleSet    [n]->set(1, dm_y);
        modelStatisticalSet_.InducedDipoleSet    [n]->set(2, dm_z);

        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(0, qr_xx);
        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(1, qr_xy);
        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(2, qr_xz);
        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(3, qr_yy);
        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(4, qr_yz);
        referenceStatisticalSet_.InducedQuadrupoleSet[n]->set(5, qr_zz);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(0, qm_xx);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(1, qm_xy);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(2, qm_xz);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(3, qm_yy);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(4, qm_yz);
        modelStatisticalSet_.InducedQuadrupoleSet[n]->set(5, qm_zz);

        dr_av = oepdev::average_moment(referenceStatisticalSet_.InducedDipoleSet[n]);
        qr_av = oepdev::average_moment(referenceStatisticalSet_.InducedQuadrupoleSet[n]);
        dm_av = oepdev::average_moment(modelStatisticalSet_.InducedDipoleSet[n]);
        qm_av = oepdev::average_moment(modelStatisticalSet_.InducedQuadrupoleSet[n]);


        if (hasAbInitioDipolePolarizability_) {
        double da_x = 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[0]);
        double da_y = 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[1]);
        double da_z = 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(dipInts[2]);
        double qa_xx= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[0]);
        double qa_xy= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[1]);
        double qa_xz= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[2]);
        double qa_yy= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[3]);
        double qa_yz= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[4]);
        double qa_zz= 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(qadInts[5]);

        abInitioModelStatisticalSet_.InducedDipoleSet    [n]->set(0, da_x);
        abInitioModelStatisticalSet_.InducedDipoleSet    [n]->set(1, da_y);
        abInitioModelStatisticalSet_.InducedDipoleSet    [n]->set(2, da_z);

        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(0, qa_xx);
        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(1, qa_xy);
        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(2, qa_xz);
        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(3, qa_yy);
        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(4, qa_yz);
        abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]->set(5, qa_zz);

        da_av = oepdev::average_moment(abInitioModelStatisticalSet_.InducedDipoleSet[n]);
        qa_av = oepdev::average_moment(abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]);
        }


        // ---> Print all to the output file <--- //
        printer->Printf("%20.8E %20.8E %20.8E %20.8E %20.8E %20.8E\n", 
                                           referenceStatisticalSet_.InducedInteractionEnergySet[n],
                                               modelStatisticalSet_.InducedInteractionEnergySet[n],
                                                      dr_av, dm_av, qr_av, qm_av);


        // ---> Optional calculations for the Ab Initio Model <--- //
        if (hasAbInitioDipolePolarizability_) {

        double tr_rest_a = abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(H) + 
                           abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(abInitioModelStatisticalSet_.JKMatrixSet[n]) +
                           abInitioModelStatisticalSet_.JKMatrixSet[n]->vector_dot(D);
        double tr_2v_a = 2.0 * abInitioModelStatisticalSet_.DensityMatrixSet[n]->vector_dot(VMatrixSet_[n]);

        abInitioModelStatisticalSet_.InducedInteractionEnergySet[n]  = tr_2v_a + tr_rest_a;

        std::shared_ptr<psi::Vector> inddip = std::make_shared<psi::Vector>("Induced Dipole (DPOL)", 3);
        std::shared_ptr<psi::Vector> indqad = std::make_shared<psi::Vector>("Induced Quadrupole (DPOL)", 6);
        double mmx = 0.0; double qqxx = 0.0; double qqyz = 0.0;
        double mmy = 0.0; double qqxy = 0.0; double qqyy = 0.0;
        double mmz = 0.0; double qqxz = 0.0; double qqzz = 0.0;
        double edipind = 0.0;

        if (options_.get_bool("DMATPOL_FF_AB_INITIO") == false) {
           for (int o=0; o<nSitesAbInitio_; ++o) {                                                                               
                std::shared_ptr<psi::Matrix> pol = abInitioPolarizationSusceptibilitiesFactory_->cphf_solver()->polarizability(o);
                std::shared_ptr<psi::Vector> lmoc= abInitioPolarizationSusceptibilitiesFactory_->cphf_solver()->lmo_centroid(o);

                double r0x= wfn_->molecule()->center_of_mass().get(0);
                double r0y= wfn_->molecule()->center_of_mass().get(1);
                double r0z= wfn_->molecule()->center_of_mass().get(2);
                double rox= lmoc->get(0);
                double roy= lmoc->get(1);
                double roz= lmoc->get(2);

                double fx = abInitioModelElectricFieldSet_[n][o]->get(0);
                double fy = abInitioModelElectricFieldSet_[n][o]->get(1);
                double fz = abInitioModelElectricFieldSet_[n][o]->get(2);

                double dmmx = pol->get(0,0)*fx + pol->get(0,1)*fy + pol->get(0,2)*fz;
                double dmmy = pol->get(1,0)*fx + pol->get(1,1)*fy + pol->get(1,2)*fz;
                double dmmz = pol->get(2,0)*fx + pol->get(2,1)*fy + pol->get(2,2)*fz;

                mmx += dmmx;
                mmy += dmmy;
                mmz += dmmz;

                qqxx-= 2.0 * dmmx * (r0x - rox);
                qqyy-= 2.0 * dmmy * (r0y - roy);
                qqzz-= 2.0 * dmmz * (r0z - roz);
                qqxy-= dmmx * (r0y - roy) + dmmy * (r0x - rox);
                qqxz-= dmmx * (r0z - roz) + dmmz * (r0x - rox);
                qqyz-= dmmy * (r0z - roz) + dmmz * (r0y - roy);

                edipind -= dmmx * fx + dmmy * fy + dmmz * fz;
           }
        } else {
          // nothing to do
        }
        inddip->set(0, mmx); indqad->set(0, qqxx); indqad->set(3, qqyy);  
        inddip->set(1, mmy); indqad->set(1, qqxy); indqad->set(4, qqyz);
        inddip->set(2, mmz); indqad->set(2, qqxz); indqad->set(5, qqzz);

        referenceDpolStatisticalSet_.InducedDipoleSet[n]->copy(*inddip);
        referenceDpolStatisticalSet_.InducedQuadrupoleSet[n]->copy(*indqad);
        referenceDpolStatisticalSet_.InducedInteractionEnergySet[n] = 0.0;
        referenceDpolStatisticalSet_.InducedInteractionEnergySet[n] += (1.0/2.0) * edipind;

        double drd_av = oepdev::average_moment(inddip);
        double qrd_av = oepdev::average_moment(indqad);

        printer_abini->Printf("%20.8E %20.8E %20.8E %20.8E %20.8E %20.8E %20.8E %20.8E %20.8E\n", 
                                                 referenceStatisticalSet_.InducedInteractionEnergySet[n],
                                             abInitioModelStatisticalSet_.InducedInteractionEnergySet[n],
                                             referenceDpolStatisticalSet_.InducedInteractionEnergySet[n],
                                                     dr_av, da_av, drd_av, qr_av, qa_av, qrd_av);
        }
   }
   // ---> Finish with the JK object (clear but not destroy it) <--- //
   jk_->finalize();

   // ---> Close the statistical output files <--- //
   //fb.close();
   //if (hasAbInitioDipolePolarizability_) fb2.close();

   // ---> Compute RMS indicators <--- //
   double rmse = 0.0; double r2e = 1.0; double sre = 0.0; double ste = 0.0;
   double rmsd = 0.0; double r2d = 1.0; double srd = 0.0; double std = 0.0;
   double rmsq = 0.0; double r2q = 1.0; double srq = 0.0; double stq = 0.0;
   double rmse1= 0.0; double r2e1= 1.0; double sre1= 0.0;
   double rmsd1= 0.0; double r2d1= 1.0; double srd1= 0.0;
   double rmsq1= 0.0; double r2q1= 1.0; double srq1= 0.0;

   double rmse2= 0.0; double r2e2= 1.0; double sre2= 0.0; 
   double rmsd2= 0.0; double r2d2= 1.0; double srd2= 0.0; 
   double rmsq2= 0.0; double r2q2= 1.0; double srq2= 0.0; 
  

   double ave = 0.0, avd = 0.0, avq = 0.0;
   for (int n=0; n<nSamples_; ++n) {
        ave += referenceStatisticalSet_.InducedInteractionEnergySet[n];
        avd += oepdev::average_moment(referenceStatisticalSet_.InducedDipoleSet[n]);
        avq += oepdev::average_moment(referenceStatisticalSet_.InducedQuadrupoleSet[n]);
   }
   ave /= (double)nSamples_;
   avd /= (double)nSamples_;
   avq /= (double)nSamples_;
   for (int n=0; n<nSamples_; ++n) {
        double er = referenceStatisticalSet_.InducedInteractionEnergySet[n];
        double em =     modelStatisticalSet_.InducedInteractionEnergySet[n];
        double dr = oepdev::average_moment(referenceStatisticalSet_.InducedDipoleSet[n]);
        double dm = oepdev::average_moment(    modelStatisticalSet_.InducedDipoleSet[n]);
        double qr = oepdev::average_moment(referenceStatisticalSet_.InducedQuadrupoleSet[n]);
        double qm = oepdev::average_moment(    modelStatisticalSet_.InducedQuadrupoleSet[n]);
        rmse += pow(er-em,2.0); sre += pow(er-em,2.0); ste += pow(er - ave,2.0);
        rmsd += pow(dr-dm,2.0); srd += pow(dr-dm,2.0); std += pow(dr - avd,2.0);
        rmsq += pow(qr-qm,2.0); srq += pow(qr-qm,2.0); stq += pow(qr - avq,2.0);
        if (hasAbInitioDipolePolarizability_) {
        double em1= abInitioModelStatisticalSet_.InducedInteractionEnergySet[n];
        double dm1= oepdev::average_moment(abInitioModelStatisticalSet_.InducedDipoleSet[n]);
        double qm1= oepdev::average_moment(abInitioModelStatisticalSet_.InducedQuadrupoleSet[n]);
        double em2= referenceDpolStatisticalSet_.InducedInteractionEnergySet[n];
        double dm2= oepdev::average_moment(referenceDpolStatisticalSet_.InducedDipoleSet[n]);
        double qm2= oepdev::average_moment(referenceDpolStatisticalSet_.InducedQuadrupoleSet[n]);
        rmse1+= pow(er-em1,2.0); sre1+= pow(er-em1,2.0);
        rmsd1+= pow(dr-dm1,2.0); srd1+= pow(dr-dm1,2.0);
        rmsq1+= pow(qr-qm1,2.0); srq1+= pow(qr-qm1,2.0);
        rmse2+= pow(er-em2,2.0); sre2+= pow(er-em2,2.0);
        rmsd2+= pow(dr-dm2,2.0); srd2+= pow(dr-dm2,2.0);
        rmsq2+= pow(qr-qm2,2.0); srq2+= pow(qr-qm2,2.0);
        }
   }
   rmse /= (double)nSamples_; rmse = sqrt(rmse); r2e -= sre/ste;
   rmsd /= (double)nSamples_; rmsd = sqrt(rmsd); r2d -= srd/std;
   rmsq /= (double)nSamples_; rmsq = sqrt(rmsq); r2q -= srq/stq;
   if (hasAbInitioDipolePolarizability_) {
   rmse1/= (double)nSamples_; rmse1= sqrt(rmse1); r2e1-= sre1/ste;
   rmsd1/= (double)nSamples_; rmsd1= sqrt(rmsd1); r2d1-= srd1/std;
   rmsq1/= (double)nSamples_; rmsq1= sqrt(rmsq1); r2q1-= srq1/stq;
   rmse2/= (double)nSamples_; rmse2= sqrt(rmse2); r2e2-= sre2/ste;
   rmsd2/= (double)nSamples_; rmsd2= sqrt(rmsd2); r2d2-= srd2/std;
   rmsq2/= (double)nSamples_; rmsq2= sqrt(rmsq2); r2q2-= srq2/stq;
   }

   // ---> Compute average electric fields as well as interaction energies <--- //
   double fx_av = 0.0;
   double fy_av = 0.0;
   double fz_av = 0.0;
   double eint_av = 0.0;
   for (int n=0; n<nSamples_; ++n) {
        double fx_av_n = 0.0; 
        double fy_av_n = 0.0;
        double fz_av_n = 0.0;
        for (int s=0; s<nSites_; ++s) {
             fx_av_n += std::abs(electricFieldSet_[n][s]->get(0));
             fy_av_n += std::abs(electricFieldSet_[n][s]->get(1));
             fz_av_n += std::abs(electricFieldSet_[n][s]->get(2));
        }
        fx_av += fx_av_n / (double)nSites_;
        fy_av += fy_av_n / (double)nSites_;
        fz_av += fz_av_n / (double)nSites_;
        //
        eint_av += std::abs(referenceStatisticalSet_.InducedInteractionEnergySet[n]);
   }
   fx_av   /= (double)nSamples_;
   fy_av   /= (double)nSamples_;
   fz_av   /= (double)nSamples_;
   eint_av /= (double)nSamples_;


   // ---> Print to output file <--- //
   const double c1 = 627.509;
   const double c2 = 1.0;
   const double c3 = 1.0;
   psi::outfile->Printf(" \n ===> Statistical Results: Generalized Susceptibility <===\n\n");
   psi::outfile->Printf(" RMSE = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmse, rmse*c1, r2e);
   psi::outfile->Printf(" RMSD = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsd, rmsd*c2, r2d);
   psi::outfile->Printf(" RMSQ = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsq, rmsq*c3, r2q);
   if (hasAbInitioDipolePolarizability_) {
   psi::outfile->Printf(" \n ===> Statistical Results: Ab Initio Susceptibility <===\n\n");
   psi::outfile->Printf(" RMSE = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmse1, rmse1*c1, r2e1);
   psi::outfile->Printf(" RMSD = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsd1, rmsd1*c2, r2d1);
   psi::outfile->Printf(" RMSQ = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsq1, rmsq1*c3, r2q1);
   psi::outfile->Printf(" \n ===> Statistical Results: Reference Distributed Polarizability Model <===\n\n");
   psi::outfile->Printf(" RMSE = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmse2, rmse2*c1, r2e2);
   psi::outfile->Printf(" RMSD = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsd2, rmsd2*c2, r2d2);
   psi::outfile->Printf(" RMSQ = %14.8f [A.U.]  %14.8f [--------]  R^2=%8.6f\n", rmsq2, rmsq2*c3, r2q2);
   }
   psi::outfile->Printf(" \n\n ===> Average Fields and Interaction Energies <===\n\n");
   psi::outfile->Printf(" <|FX|> = %13.4f\n", fx_av);
   psi::outfile->Printf(" <|FY|> = %13.4f\n", fy_av);
   psi::outfile->Printf(" <|FZ|> = %13.4f\n", fz_av);
   psi::outfile->Printf(" <|dE|> = %13.4f [A.U.]  %13.4f [kcal/mole]\n", eint_av, eint_av*c1);

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
