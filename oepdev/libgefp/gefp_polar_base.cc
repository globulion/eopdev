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
   if (hasQuadrupolePolarizability_       ) compute_electric_field_gradient_sums();
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
  if (!hasDipoleDipoleHyperpolarizability_ and !hasQuadrupolePolarizability_) {
      for (int n=0; n<nSites_; ++n) {
           for (int z=0; z<3; ++z) {
                double val = Parameters_->get(3*n + z, 0);
                PolarizationSusceptibilities_->dipole_polarizability(n, z)->set(i, j, val);
           }
      }
  } else if (hasDipoleDipoleHyperpolarizability_ and !hasQuadrupolePolarizability_) {
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
  } else if (!hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
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
  } else if (hasDipoleDipoleHyperpolarizability_ and hasQuadrupolePolarizability_) {
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
void oepdev::GeneralizedPolarGEFactory::compute_statistics(void) {
   cout << " Statistical evaluation ...\n";

   // Grab some unperturbed quantities
   std::shared_ptr<psi::Matrix> H = wfn_->Fa(); H->add(wfn_->H());
   std::shared_ptr<psi::Matrix> D = wfn_->Da();

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

   // ---> Printer <--- //
   std::filebuf fb;
   fb.open(options_.get_str("DMATPOL_OUT_STATS"), std::ios::out);
   std::shared_ptr<std::ostream> os = std::make_shared<std::ostream>(&fb);
   std::shared_ptr<psi::PsiOutStream> printer = std::make_shared<psi::PsiOutStream>(os);
   printer->Printf("# Eint(Ref)    Eint(Model)\n");

   for (int n=0; n<nSamples_; ++n) {

        // ---> Compute interaction energy <--- // nuclear part is already computed in void compute_samples(void)

        double eint_refer = H->vector_dot(referenceStatisticalSet_.DensityMatrixSet[n]);
        double eint_model = H->vector_dot(    modelStatisticalSet_.DensityMatrixSet[n]);

        std::shared_ptr<psi::Matrix> Dn1= D->clone(); Dn1->add(referenceStatisticalSet_.DensityMatrixSet[n]);
        std::shared_ptr<psi::Matrix> Dn2= D->clone(); Dn2->add(    modelStatisticalSet_.DensityMatrixSet[n]);
        std::shared_ptr<psi::Matrix> G = D->clone(); G->zero();
        G->add(VMatrixSet_[n]); G->scale(2.0);
        std::shared_ptr<psi::Matrix> G1 = G->clone(); G1->add(referenceStatisticalSet_.JKMatrixSet[n]);
        std::shared_ptr<psi::Matrix> G2 = G->clone(); G2->add(    modelStatisticalSet_.JKMatrixSet[n]);
        eint_refer+= G1->vector_dot(Dn1);
        eint_model+= G2->vector_dot(Dn2);

        referenceStatisticalSet_.InducedInteractionEnergySet[n] += eint_refer;
            modelStatisticalSet_.InducedInteractionEnergySet[n] += eint_model;


        // ---> Compute induced dipole moments <--- //
        // TODO

        // ---> Print all to the output file <--- //
        printer->Printf("%20.8E %20.8E\n", referenceStatisticalSet_.InducedInteractionEnergySet[n],
                                               modelStatisticalSet_.InducedInteractionEnergySet[n]);

   }
   // Finish with the JK object (clear but not destroy it)
   jk_->finalize();

   // ---> Compute RMS indicators <--- //
   double rmse = 0.0; double r2e = 1.0; double sre = 0.0; double ste = 0.0;
   double rmsd = 0.0; double r2d = 1.0;
   double rmsq = 0.0; double r2q = 1.0;
   double ave = 0.0;
   for (int n=0; n<nSamples_; ++n) {
        ave += referenceStatisticalSet_.InducedInteractionEnergySet[n];
   }
   ave /= (double)nSamples_;
   for (int n=0; n<nSamples_; ++n) {
        double er = referenceStatisticalSet_.InducedInteractionEnergySet[n];
        double em =     modelStatisticalSet_.InducedInteractionEnergySet[n];
        double dr = 0.0; // TODO
        double dm = 0.0; // TODO
        double qr = 0.0; // TODO
        double qm = 0.0; // TODO
        rmse += pow(er-em,2.0); sre += pow(er-em,2.0); ste += pow(er - ave,2.0);
        rmsd += pow(dr-dm,2.0);
        rmsq += pow(qr-qm,2.0);
   }
   rmse /= (double)nSamples_; rmse = sqrt(rmse); r2e -= sre/ste;
   rmsd /= (double)nSamples_; rmsd = sqrt(rmsd); 
   rmsq /= (double)nSamples_; rmsq = sqrt(rmsq);

   // ---> Print to output file <--- //
   const double c1 = 627.509;
   const double c2 = 1.0;
   const double c3 = 1.0;
   psi::outfile->Printf(" RMSE = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmse, rmse*c1, r2e);
   psi::outfile->Printf(" RMSD = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmsd, rmsd*c2, r2d);
   psi::outfile->Printf(" RMSQ = %14.8f [A.U.]  %14.8f [kcal/mol]  R^2=%8.6f\n", rmsq, rmsq*c3, r2q);
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
