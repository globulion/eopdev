#include "util.h"
#include "quambo.h"
#include "wavefunction_union.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace oepdev{

using namespace psi;
using namespace std;


WavefunctionUnion::WavefunctionUnion(SharedWavefunction ref_wfn, Options& options) 
   : Wavefunction(options), hasLocalizedOrbitals_(false), integrals_(nullptr)
{
   // Primary Sanity-check
   if (ref_wfn->molecule()->nfragments()==1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. Molecule has to have minimum two fragments.\n");
   if (ref_wfn->nirrep()>1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. C1 symmetry must be enforced.\n");

   // Secondary Sanity-check (unimplemented features)
   if (ref_wfn->molecule()->multiplicity()!=1 || !ref_wfn->same_a_b_orbs() || !ref_wfn->same_a_b_dens()) throw
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only CLOSE SHELLS are supported!\n");
   if (ref_wfn->molecule()->nfragments()>2) throw
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only DIMERS (nfrag=2) are supported!\n");
   
   // If passed sanity-checks, then initialize the object
   shallow_copy(ref_wfn);

   // Extract the monomers
   SharedSuperFunctional functional = create_superfunctional("HF", options_);
   SharedMolecule molecule_1 = extract_monomer(molecule_, 1);
   SharedMolecule molecule_2 = extract_monomer(molecule_, 2);
   molecule_1->set_name("Monomer 1");
   molecule_2->set_name("Monomer 2");
   molecule_ ->set_name("Aggregate (Dimer)");
   SharedBasisSet primary_1  = basissets_["BASIS_1"];
   SharedBasisSet primary_2  = basissets_["BASIS_2"];
   SharedBasisSet auxiliary_1  = basissets_["BASIS_DF_OEP_1"];
   SharedBasisSet auxiliary_2  = basissets_["BASIS_DF_OEP_2"];
   SharedBasisSet auxiliary_df_1  = basissets_["BASIS_DF_SCF_1"];
   SharedBasisSet auxiliary_df_2  = basissets_["BASIS_DF_SCF_2"];
   SharedBasisSet intermediate_1  = basissets_["BASIS_INT_OEP_1"];
   SharedBasisSet intermediate_2  = basissets_["BASIS_INT_OEP_2"];
   SharedBasisSet guess_1  = basissets_["BASIS_GUESS_1"];
   SharedBasisSet guess_2  = basissets_["BASIS_GUESS_2"];
   SharedWavefunction wfn_1  = solve_scf(molecule_1, primary_1, auxiliary_df_1, guess_1, functional, options_, psi::PSIO::shared_object(), false);
   SharedWavefunction wfn_2  = solve_scf(molecule_2, primary_2, auxiliary_df_2, guess_2, functional, options_, psi::PSIO::shared_object(), false);

   // Finish initialize
   common_init(ref_wfn, molecule_1, molecule_2, 
		        primary_1, primary_2, 
			auxiliary_1, auxiliary_2, 
			auxiliary_df_1, auxiliary_df_2, 
			intermediate_1, intermediate_2, 
                        guess_1, guess_2,
			wfn_1, wfn_2);
}

WavefunctionUnion::WavefunctionUnion(
		SharedMolecule dimer,
		SharedBasisSet primary,
		SharedBasisSet auxiliary_df,
                SharedBasisSet guess,
		SharedBasisSet primary_1,
		SharedBasisSet primary_2,
		SharedBasisSet auxiliary_1,
		SharedBasisSet auxiliary_2,
		SharedBasisSet auxiliary_df_1,
		SharedBasisSet auxiliary_df_2,
		SharedBasisSet intermediate_1,
		SharedBasisSet intermediate_2,
		SharedBasisSet guess_1,
		SharedBasisSet guess_2,
		SharedWavefunction wfn_1,
		SharedWavefunction wfn_2,
		Options& options
		)
 : Wavefunction(options),
   hasLocalizedOrbitals_(false),
   integrals_(nullptr)
{
   // Primary Sanity-check
   if (dimer->nfragments()==1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. Molecule has to have minimum two fragments.\n");
   if (wfn_1->nirrep()>1 || wfn_2->nirrep()>1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. C1 symmetry must be enforced in monomer wavefunctons!.\n");

   // Secondary Sanity-check (unimplemented features)
   if (wfn_1->molecule()->multiplicity()!=1 || !wfn_1->same_a_b_orbs() || !wfn_1->same_a_b_dens() ||
       wfn_2->molecule()->multiplicity()!=1 || !wfn_2->same_a_b_orbs() || !wfn_2->same_a_b_dens()
		   ) throw
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only CLOSE SHELLS are supported!\n");
   if (dimer->nfragments()>2) throw
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only DIMERS (nfrag=2) are supported!\n");

   // Create full dimer wavefunction: compute_mints is set to true to clean-up scratch because this constructor
   // is used only for interfacing with Python (building WavefunctionUnion from Python level)
   SharedWavefunction ref_wfn = oepdev::solve_scf(dimer, primary, auxiliary_df, guess,
   		                create_superfunctional("HF", options_), options_, psi::PSIO::shared_object(), true);
   shallow_copy(ref_wfn);
   set_basisset("BASIS_DF_SCF", auxiliary_df);
   set_basisset("BASIS_GUESS", guess);

   // Extract monomers
   SharedMolecule molecule_1 = extract_monomer(molecule_, 1);
   SharedMolecule molecule_2 = extract_monomer(molecule_, 2);
   molecule_1->set_name("Monomer 1");
   molecule_2->set_name("Monomer 2");
   molecule_ ->set_name("Aggregate (Dimer)");

   // Finish initialize
   common_init(ref_wfn, molecule_1, molecule_2, 
		        primary_1, primary_2, 
			auxiliary_1, auxiliary_2, 
			auxiliary_df_1, auxiliary_df_2, 
			intermediate_1, intermediate_2, 
                        guess_1, guess_2,
			wfn_1, wfn_2);
}

WavefunctionUnion::~WavefunctionUnion() 
{
}
void WavefunctionUnion::common_init(
		SharedWavefunction ref_wfn,
		SharedMolecule molecule_1,
		SharedMolecule molecule_2,
		SharedBasisSet primary_1,
		SharedBasisSet primary_2,
		SharedBasisSet auxiliary_1,
		SharedBasisSet auxiliary_2,
		SharedBasisSet auxiliary_df_1,
		SharedBasisSet auxiliary_df_2,
		SharedBasisSet intermediate_1,
		SharedBasisSet intermediate_2,
		SharedBasisSet guess_1,
		SharedBasisSet guess_2,
		SharedWavefunction wfn_1,
		SharedWavefunction wfn_2
		)
{
   // Compute QUAMBOs if required                                                
   psi::SharedMatrix Ca_occ_1, Ca_occ_2, Ca_vir_1, Ca_vir_2;
   psi::SharedVector eps_a_occ_1, eps_a_occ_2, eps_a_vir_1, eps_a_vir_2;

   if (options_.get_bool("OEPDEV_USE_VVO")) {
       // VVOs
       std::shared_ptr<QUAMBO> solver_1 = std::make_shared<QUAMBO>(wfn_1, options_.get_bool("QUAMBO_ACBS"));
       solver_1->compute();
       std::shared_ptr<QUAMBO> solver_2 = std::make_shared<QUAMBO>(wfn_2, options_.get_bool("QUAMBO_ACBS"));
       solver_2->compute();
                                                                                 
       Ca_occ_1 = solver_1->Ca_subset("AO","OCC");
       Ca_vir_1 = solver_1->Ca_subset("AO","VIR");
       eps_a_occ_1 = solver_1->epsilon_a_subset("MO","OCC");
       eps_a_vir_1 = solver_1->epsilon_a_subset("MO","VIR");
                                                                                 
       Ca_occ_2 = solver_2->Ca_subset("AO","OCC");
       Ca_vir_2 = solver_2->Ca_subset("AO","VIR");
       eps_a_occ_2 = solver_2->epsilon_a_subset("MO","OCC");
       eps_a_vir_2 = solver_2->epsilon_a_subset("MO","VIR");
                                                                                 
       nmo_ = solver_1->nbas() + solver_2->nbas();
     //nso_ = this->basisset()->nbf(); --> no need to change this
       nmopi_[0] = nmo_;
     //nsopi_[0] = nso_; --> no need to change this

   } else {
       // RHF canonical orbitals
       Ca_occ_1 = wfn_1->Ca_subset("AO","OCC");  
       Ca_occ_2 = wfn_2->Ca_subset("AO","OCC");
       Ca_vir_1 = wfn_1->Ca_subset("AO","VIR");
       Ca_vir_2 = wfn_2->Ca_subset("AO","VIR");
       eps_a_occ_1 = wfn_1->epsilon_a_subset("MO","OCC");
       eps_a_occ_2 = wfn_2->epsilon_a_subset("MO","OCC");
       eps_a_vir_1 = wfn_1->epsilon_a_subset("MO","VIR");
       eps_a_vir_2 = wfn_2->epsilon_a_subset("MO","VIR");
   }

   // Sizing of the union
   SharedSuperFunctional functional = create_superfunctional("HF", options_);
   const int nocc_1                = Ca_occ_1->ncol();
   const int nocc_2                = Ca_occ_2->ncol();
   const int nvir_1                = Ca_vir_1->ncol();
   const int nvir_2                = Ca_vir_2->ncol();
   const int nmo_1                 = nocc_1 + nvir_1;
   const int nmo_2                 = nocc_2 + nvir_2;
   const int nbf_union             = this->basisset()->nbf();
   const int nmo_union             = nmo_1 + nmo_2;

   // Properties of the union
   nIsolatedMolecules_ = ref_wfn->molecule()->nfragments();
   energy_             = wfn_1->reference_energy() + wfn_2->reference_energy();
   efzc_               = wfn_1->efzc() + wfn_2->efzc();

   // Properties of isolated monomers organize in vectors
   l_energy_        .push_back(wfn_1->reference_energy()); l_energy_        .push_back(wfn_2->reference_energy()   );  
   l_molecule_      .push_back(molecule_1               ); l_molecule_      .push_back(molecule_2                  );   
   l_wfn_           .push_back(wfn_1                    ); l_wfn_           .push_back(wfn_2                       );
   l_primary_       .push_back(primary_1                ); l_primary_       .push_back(primary_2                   );
   l_auxiliary_     .push_back(auxiliary_1              ); l_auxiliary_     .push_back(auxiliary_2                 );
   l_intermediate_  .push_back(intermediate_1           ); l_intermediate_  .push_back(intermediate_2              );
   l_guess_         .push_back(guess_1                  ); l_guess_         .push_back(guess_2                     );
   l_name_          .push_back(wfn_1->name()            ); l_name_          .push_back(wfn_2->name()               );
   l_nbf_           .push_back(primary_1->nbf()         ); l_nbf_           .push_back(primary_2->nbf()            );
   l_nmo_           .push_back(nmo_1                    ); l_nmo_           .push_back(nmo_2                       );
   l_nso_           .push_back(nmo_1                    ); l_nso_           .push_back(nmo_2                       );
   l_ndocc_         .push_back(nocc_1                   ); l_ndocc_         .push_back(nocc_2                      );
   l_nvir_          .push_back(nvir_1                   ); l_nvir_          .push_back(nvir_2                      );
   l_efzc_          .push_back(wfn_1->efzc()            ); l_efzc_          .push_back(wfn_2->efzc()               );
   l_density_fitted_.push_back(wfn_1->density_fitted()  ); l_density_fitted_.push_back(wfn_2->density_fitted()     );
   l_nalpha_        .push_back(wfn_1->nalpha()          ); l_nalpha_        .push_back(wfn_2->nalpha()             );  
   l_nbeta_         .push_back(wfn_1->nbeta()           ); l_nbeta_         .push_back(wfn_2->nbeta()              );
   l_nfrzc_         .push_back(wfn_1->nfrzc()           ); l_nfrzc_         .push_back(wfn_2->nfrzc()              );
   l_noffs_ao_      .push_back(0                        ); l_noffs_ao_      .push_back(primary_1->nbf()            );
   l_ca_occ_        .push_back(Ca_occ_1                 ); l_ca_occ_        .push_back(Ca_occ_2                    );
   l_ca_vir_        .push_back(Ca_vir_1                 ); l_ca_vir_        .push_back(Ca_vir_2                    );
   l_eps_a_occ_     .push_back(eps_a_occ_1              ); l_eps_a_occ_     .push_back(eps_a_occ_2                 );
   l_eps_a_vir_     .push_back(eps_a_vir_1              ); l_eps_a_vir_     .push_back(eps_a_vir_2                 );


  
   //dimer_wavefunction_     = ...; // store the original wavefunction
   //reference_wavefunction_ = reference_wavefunction();
   //reference_wavefunction_   = ref_wfn->reference_wavefunction();
   //reference_wavefunction_ = shared_from_this();

   /* Now, the matrices of the original wavefunction will be 
    * overriden by new matrices
    */

   // <==== Sanity Checks ====> //
   if (!epsilon_a_ || !wfn_1->epsilon_a()) throw 
                                   PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion Epsilon. Orbital Energies not available!");
   if (!Ca_ || !wfn_1->Ca()) throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion LCAO-MO Coeffs. Wavefunction Coefficients not available!");
   if (!Da_ || !wfn_1->Da()) throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion Density. Density Matrix not available!"); 
   if (!Fa_ || !wfn_1->Fa()) throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion Fock. Fock Matrix not available!");
   if (!H_  || !wfn_1->H ()) throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion CoreH. Core Hamiltonian not available!");
   if (!S_  || !wfn_1->S ()) throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion Overlap. Overlap Matrix not available");


   // <---- Wavefunction Coefficients (LCAO-MO Matrices) ----> //
   // <---- Orbital Energies                             ----> //
   epsilon_a_ = std::make_shared<psi::Vector>("Union Alpha Orbital Energies", nmo_union);
   epsilon_b_ = std::make_shared<psi::Vector>("Union Beta Orbital Energies", nmo_union);
   Ca_ = std::make_shared<psi::Matrix>("Union Alpha Orbitals", nbf_union, nmo_union);
   Cb_ = std::make_shared<psi::Matrix>("Union Beta Orbitals", nbf_union, nmo_union);

   double** pCa = Ca_->pointer(); 
   double** pCb = Cb_->pointer();
   double*  pea = epsilon_a_->pointer();
   double*  peb = epsilon_b_->pointer();

   std::vector<std::vector<int>> orbitals_occ, orbitals_vir;
   std::vector<int> orbitals_occ_all, orbitals_vir_all;
   std::vector<int> dummy;
   /* implementation below is for CLOSED SHELLS only! */
   int nbf, nmo, nmo_occ, nmo_vir;
   int nmo_occ_all = 0; for (int nf=0;nf<nIsolatedMolecules_;++nf) nmo_occ_all += l_ndocc_[nf];
   int nOffsetAO = 0, nOffsetMOOcc = 0, nOffsetMOVir = nmo_occ_all;
   //
   for (int nf=0; nf<nIsolatedMolecules_; nf++) {
        nbf      = l_wfn_[nf]->basisset()->nbf();
        nmo      = l_nmo_[nf];
        nmo_occ  = l_ndocc_[nf];
        nmo_vir  = nmo - nmo_occ;
        for (int i=0; i<nbf; i++) {      
             // <--- Occupied Orbitals ---> //
             for (int jo=0; jo<nmo_occ; jo++) {
                  pCa[i+nOffsetAO][jo+nOffsetMOOcc] = l_ca_occ_[nf]->get(0, i, jo);
             }
             // <--- Virtual Orbitals ---> //
             for (int jv=0; jv<nmo_vir; jv++) {
                  pCa[i+nOffsetAO][jv+nOffsetMOVir] = l_ca_vir_[nf]->get(0, i, jv);
             }
        }
        // <--- Occupied orbital Energies and Indices ---> //
        for (int jo=0; jo<nmo_occ; jo++) {
             pea[jo+nOffsetMOOcc] = l_eps_a_occ_[nf]->get(0, jo);
             dummy.push_back(jo+nOffsetMOOcc);
             orbitals_occ_all.push_back(jo+nOffsetMOOcc);
        }
        orbitals_occ.push_back(dummy); dummy.clear();

        // <--- Virtual orbital Energies and Indices ---> //
        for (int jv=0; jv<nmo_vir; jv++) {
             pea[jv+nOffsetMOVir] = l_eps_a_vir_[nf]->get(0, jv);
             dummy.push_back(jv+nOffsetMOVir);
             orbitals_vir_all.push_back(jv+nOffsetMOVir);
        }
        orbitals_vir.push_back(dummy); dummy.clear();

        // 
        nOffsetAO    += nbf;
        nOffsetMOOcc += nmo_occ;
        nOffsetMOVir += nmo_vir;
   }
   Cb_       ->copy( Ca_       ); // Because implemented for closed-shell only!
   epsilon_b_->copy(*epsilon_a_); // 

   // <---- AO One-Particle Density Matrices (OPDM's) ----> //
   // <---- AO Fock Matrices                          ----> //
   // <---- AO Core One-Electron Hamiltonian          ----> //
   // <---- AO Overlap Matrix                         ----> //

   Da_->zero(); Db_->zero();                                                      
   Fa_->zero(); Fb_->zero();                                                      
   H_ ->zero(); S_ ->zero();

   double** pDa = Da_->pointer(); 
   double** pDb = Db_->pointer();
   double** pFa = Fa_->pointer(); 
   double** pFb = Fb_->pointer();
   double** pH  = H_ ->pointer();
   double** pS  = S_ ->pointer(); 

   nOffsetAO = 0;
   for (int nf = 0; nf<nIsolatedMolecules_; nf++) {
        nbf  = l_wfn_[nf]->basisset()->nbf();
        for (int i=0; i<nbf; i++) {
             for (int j=0; j<nbf; j++) {
                  pDa[i+nOffsetAO][j+nOffsetAO] = l_wfn_[nf]->Da()->get(0, i, j);
                  pFa[i+nOffsetAO][j+nOffsetAO] = l_wfn_[nf]->Fa()->get(0, i, j);
                  pH [i+nOffsetAO][j+nOffsetAO] = l_wfn_[nf]->H ()->get(0, i, j);
                  pS [i+nOffsetAO][j+nOffsetAO] = l_wfn_[nf]->S ()->get(0, i, j);
             }
        }
        nOffsetAO += nbf;
   }
   Db_->copy(Da_);
   Fb_->copy(Fa_);

   // <---- MO spaces of the union ----> //
   /* as for now for two fragments only */
   if (nIsolatedMolecules_>2) throw 
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only DIMERS (nfrag=2) are supported!\n");

   SharedMOSpace space_1_occ = std::make_shared<MOSpace>('I', orbitals_occ[0], dummy);
   SharedMOSpace space_2_occ = std::make_shared<MOSpace>('J', orbitals_occ[1], dummy);
   SharedMOSpace space_1_vir = std::make_shared<MOSpace>('X', orbitals_vir[0], dummy);
   SharedMOSpace space_2_vir = std::make_shared<MOSpace>('Y', orbitals_vir[1], dummy);
   std::map<const std::string, SharedMOSpace> spaces_1, spaces_2;
   spaces_1["OCC"] = space_1_occ;   spaces_1["VIR"] = space_1_vir;
   spaces_2["OCC"] = space_2_occ;   spaces_2["VIR"] = space_2_vir;
   l_mospace_.push_back(spaces_1);
   l_mospace_.push_back(spaces_2);
   if (false) { //debugging
       for (int i=0; i<orbitals_occ[0].size(); ++i) { std::cout << orbitals_occ[0][i] << "  " ;} std::cout << std::endl; 
       for (int i=0; i<orbitals_occ[1].size(); ++i) { std::cout << orbitals_occ[1][i] << "  " ;} std::cout << std::endl;
       for (int i=0; i<orbitals_vir[0].size(); ++i) { std::cout << orbitals_vir[0][i] << "  " ;} std::cout << std::endl;
       for (int i=0; i<orbitals_vir[1].size(); ++i) { std::cout << orbitals_vir[1][i] << "  " ;} std::cout << std::endl;
   }
   mospacesUnion_["OCC"] = std::make_shared<MOSpace>('F', orbitals_occ_all, dummy); // (F)illed
   mospacesUnion_["VIR"] = std::make_shared<MOSpace>('E', orbitals_vir_all, dummy); // (E)mpty


   // <==== Compute One-Electron Property object ====> //
   oeprop_ = std::make_shared<OEProp>(ref_wfn);
   oeprop_->set_title(" One-Electron Properties of Wavefunction Union");
   oeprop_->set_Da_ao(Da_);
   oeprop_->add("DIPOLE");
   oeprop_->add("QUADRUPOLE");
   oeprop_->add("MULLIKEN CHARGES");
   oeprop_->add("ESP AT NUCLEI");
   oeprop_->compute();
}

SharedIntegralTransform WavefunctionUnion::integrals() const { 
   if (integrals_) {return integrals_;}
   else 
   {
       throw PSIEXCEPTION(" OEPDEV: Error. WavefunctionUnion IntegralTransform. Integral Transform was not created!");
   }
}

SharedLocalizer WavefunctionUnion::l_localizer(int n) const { 
   if (hasLocalizedOrbitals_) {return l_localizer_[n];}
   else 
   {
       throw PSIEXCEPTION(" OEPDEV: Error. WavefunctionUnion Localizer. Union orbitals were not localized yet!");
   }
}



double WavefunctionUnion::compute_energy() {}

void WavefunctionUnion::localize_orbitals() {
  // ===> Replace canonical occupied orbitals with localized ones <=== //
  /* 
     Note: Updated orbitals are:
           Matrix::doublet(wfn->Ca_subset("AO","OCC"), localizer->U(), false, false)  
           which is exactly equal to localizer->L() 
           Orbitals of the monomers are also changed to localized ones.
  */
  double** pCa = Ca_->pointer();
  int nbf, nmo_occ;
  int nOffsetAO = 0, nOffsetMOOcc = 0;
  for (int nf = 0; nf < nIsolatedMolecules_; ++nf) {
       l_localizer_.push_back(Localizer::build(options_.get_str("WFN_UNION_LOCALIZER"), 
                           l_primary_[nf], l_ca_occ_[nf], options_));
       l_localizer_[nf]->set_maxiter(options_.get_int("OEPDEV_LOCALIZER_MAXITER"));
       l_localizer_[nf]->localize();
       //
       nbf      = l_wfn_[nf]->basisset()->nbf();
       nmo_occ  = l_ndocc_[nf];
       //
       double** pca = l_ca_occ_[nf]->pointer();
     //double** pcb = l_ca->pointer();
       // 
       for (int i=0; i<nbf; ++i) {
            for (int jo=0; jo<nmo_occ; ++jo) {
                 pCa[i+nOffsetAO][jo+nOffsetMOOcc] = l_localizer_[nf]->L()->get(0, i, jo);
                 pca[i][jo] = l_localizer_[nf]->L()->get(0, i, jo);
               //pcb[i][jo] = pca[i][jo];
            }
       }
       nOffsetAO    += nbf;
       nOffsetMOOcc += nmo_occ;

  }
  Cb_->copy(Ca_);
  hasLocalizedOrbitals_ = true;
}

void WavefunctionUnion::transform_integrals() 
{
  clock_t t_time = -clock(); // Clock BEGIN

    SharedMOSpaceVector spaces;
    SharedMOSpace space_1o = l_mospace(0,"OCC");
    SharedMOSpace space_2o = l_mospace(1,"OCC");
    SharedMOSpace space_1v = l_mospace(0,"VIR");
    SharedMOSpace space_2v = l_mospace(1,"VIR");
    //SharedMOSpace space_12o=   mospace(  "OCC");//MOSpace::occ; //--> equivalent
    SharedMOSpace space_12o= MOSpace::occ;
    spaces.push_back(space_1o);
    spaces.push_back(space_2o);
    spaces.push_back(space_1v);
    spaces.push_back(space_2v);
    spaces.push_back(space_12o);
    integrals_ = std::make_shared<IntegralTransform>(shared_from_this(), 
    //int nbf = this->basisset()->nbf();
    //psi::SharedMatrix c = std::make_shared<psi::Matrix>("", nbf, 0); // -> frozen core
    //psi::SharedMatrix v = std::make_shared<psi::Matrix>("", nbf, 0); // -> frozen virtuals
    //psi::SharedMatrix i = this->Ca_subset("AO","OCC");               // -> active occupied
    //psi::SharedMatrix a = this->Ca_subset("AO","VIR");               // -> active virtual
    //cout << "WW " << c->ncol() << " " << i->ncol() << " " << a->ncol() << " " << v->ncol() << " " << Ca_->ncol() << endl;
    //cout << "WW " << c->nrow() << " " << i->nrow() << " " << a->nrow() << " " << v->nrow() << " " << Ca_->nrow() << endl;
  //integrals_ = std::make_shared<IntegralTransform>(this->H_, c, i, a, v,
                                                     spaces,
                                                     IntegralTransform::TransformationType::Restricted,
                                                     IntegralTransform::OutputType::DPDOnly,
                                                     IntegralTransform::MOOrdering::QTOrder,
                                                     IntegralTransform::FrozenOrbitals::None);

    integrals_->set_keep_dpd_so_ints(true);
    integrals_->set_print(0);

    // Trans (II|II)
    timer_on("Trans (II|II)");
    integrals_->transform_tei(space_1o, space_1o, space_1o, space_1o, IntegralTransform::HalfTrans::MakeAndKeep);
    timer_off("Trans (II|II)");

    // Trans (II|IJ)
    timer_on("Trans (II|IJ)");
    integrals_->transform_tei(space_1o, space_1o, space_1o, space_2o, IntegralTransform::HalfTrans::ReadAndKeep);
    timer_off("Trans (II|IJ)");

    // Trans (II|JJ)
    timer_on("Trans (II|JJ)");
    integrals_->transform_tei(space_1o, space_1o, space_2o, space_2o, IntegralTransform::HalfTrans::ReadAndNuke);
    timer_off("Trans (II|JJ)");

    //// Trans (JJ|II) BBB-T
    //timer_on("Trans (JJ|II)");
    //integrals_->transform_tei(space_2o, space_2o, space_1o, space_1o, IntegralTransform::HalfTrans::MakeAndNuke);
    //timer_off("Trans (JJ|II)");

    // Trans (IJ|IJ)
    timer_on("Trans (IJ|IJ)");
    integrals_->transform_tei(space_1o, space_2o, space_1o, space_2o, IntegralTransform::HalfTrans::MakeAndKeep);
    timer_off("Trans (IJ|IJ)");

    // Trans (IJ|JJ)
    timer_on("Trans (IJ|JJ)");
    integrals_->transform_tei(space_1o, space_2o, space_2o, space_2o, IntegralTransform::HalfTrans::ReadAndNuke);
    timer_off("Trans (IJ|JJ)");

    // Trans (JJ|JJ)
    timer_on("Trans (JJ|JJ)");
    integrals_->transform_tei(space_2o, space_2o, space_2o, space_2o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (JJ|JJ)");

    // Trans (OO|OO)
    timer_on("Trans (OO|OO)");
    integrals_->transform_tei(space_12o, space_12o, space_12o, space_12o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (OO|OO)");

    // Trans (XJ|II)
    timer_on("Trans (XJ|II)");
    integrals_->transform_tei(space_1v, space_2o, space_1o, space_1o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (XJ|II)");

    // Trans (XI|JJ)
    timer_on("Trans (XI|JJ)");
    integrals_->transform_tei(space_1v, space_1o, space_2o, space_2o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (XI|JJ)");

    // Trans (XI|JI)
    timer_on("Trans (XI|JI)");
    integrals_->transform_tei(space_1v, space_1o, space_2o, space_1o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (XI|JI)");

    // Trans (YI|JJ)
    timer_on("Trans (YI|JJ)");
    integrals_->transform_tei(space_2v, space_1o, space_2o, space_2o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (YI|JJ)");

    // Trans (YJ|II)
    timer_on("Trans (YJ|II)");
    integrals_->transform_tei(space_2v, space_2o, space_1o, space_1o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (YJ|II)");

    // Trans (YJ|IJ)
    timer_on("Trans (YJ|IJ)");
    integrals_->transform_tei(space_2v, space_2o, space_1o, space_2o, IntegralTransform::HalfTrans::MakeAndNuke);
    timer_off("Trans (YJ|IJ)");


  t_time += clock(); // Clock END
  cout << " o TIME Integral transformation (WFN UNION)  : " << ((double)t_time/CLOCKS_PER_SEC) << endl;

}

void WavefunctionUnion::print_mo_integrals(void) { 

    if (!integrals_) 
        throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion MO IntegralsPrint. No MO integrals were created yet!\n");
    if (nIsolatedMolecules_ != 2)
        throw PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion MO IntegralsPrint. Only 2 fragments are now supported.\n");

    std::shared_ptr<PSIO> psio = PSIO::shared_object();

    dpd_set_default(integrals_->get_dpd_id());
    dpdbuf4 buf_IIJJ, buf_IJJJ, buf_IIIJ, buf_IJIJ, buf_YIJJ, buf_XIJJ, buf_YJII, buf_XJII;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    if (options_.get_int("PRINT") > 2) psio->tocprint(PSIF_LIBTRANS_DPD);

    global_dpd_->buf4_init(&buf_IIIJ, PSIF_LIBTRANS_DPD, 0, 
                           integrals_->DPD_ID("[I,I]"  ), integrals_->DPD_ID("[I,J]"  ),
                           integrals_->DPD_ID("[I>=I]+"), integrals_->DPD_ID("[I,J]"  ), 0, "MO Ints (II|IJ)");
    global_dpd_->buf4_init(&buf_IIJJ, PSIF_LIBTRANS_DPD, 0, 
                           integrals_->DPD_ID("[I,I]"  ), integrals_->DPD_ID("[J,J]"  ),
                           integrals_->DPD_ID("[I>=I]+"), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (II|JJ)");
    global_dpd_->buf4_init(&buf_IJIJ, PSIF_LIBTRANS_DPD, 0, 
                           integrals_->DPD_ID("[I,J]"  ), integrals_->DPD_ID("[I,J]"  ),
                           integrals_->DPD_ID("[I,J]"  ), integrals_->DPD_ID("[I,J]"  ), 0, "MO Ints (IJ|IJ)");
    global_dpd_->buf4_init(&buf_IJJJ, PSIF_LIBTRANS_DPD, 0,
                           integrals_->DPD_ID("[I,J]"  ), integrals_->DPD_ID("[J,J]"  ),
                           integrals_->DPD_ID("[I,J]"  ), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (IJ|JJ)");

    // YIJJ = (V_2 O_1 | O_2 O_2) = (ni|jj)
    global_dpd_->buf4_init(&buf_YIJJ, PSIF_LIBTRANS_DPD, 0,
                           integrals_->DPD_ID("[Y,I]"  ), integrals_->DPD_ID("[J,J]"  ),
                           integrals_->DPD_ID("[Y,I]"  ), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (YI|JJ)");
    // XIJJ = (V_1 O_1 | O_2 O_2) = (mi|jj)
    global_dpd_->buf4_init(&buf_XIJJ, PSIF_LIBTRANS_DPD, 0,
          	         integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J,J]"  ),
          		 integrals_->DPD_ID("[X,I]"  ), integrals_->DPD_ID("[J>=J]+"), 0, "MO Ints (XI|JJ)");
    // YJII = (V_2 O_2 | O_1 O_1) = (nj|kk)
    global_dpd_->buf4_init(&buf_YJII, PSIF_LIBTRANS_DPD, 0,
                           integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I,I]"  ),
                           integrals_->DPD_ID("[Y,J]"  ), integrals_->DPD_ID("[I>=I]+"), 0, "MO Ints (YJ|II)");
    // XJII = (V_1 O_2 | O_1 O_1) = (mj|kk)
    global_dpd_->buf4_init(&buf_XJII, PSIF_LIBTRANS_DPD, 0,
          	         integrals_->DPD_ID("[X,J]"  ), integrals_->DPD_ID("[I,I]"  ),
          		 integrals_->DPD_ID("[X,J]"  ), integrals_->DPD_ID("[I>=I]+"), 0, "MO Ints (XJ|II)");


    psi::outfile->Printf("\n <=== buf_IIIJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_IIIJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_IIIJ, h);
         for (int pq = 0; pq < buf_IIIJ.params->rowtot[h]; ++pq) {
              int p = buf_IIIJ.params->roworb[h][pq][0];
              int q = buf_IIIJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_IIIJ.params->coltot[h]; ++rs) {
                   int r = buf_IIIJ.params->colorb[h][rs][0];
                   int s = buf_IIIJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_IIIJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_IIIJ, h);
    }

    psi::outfile->Printf("\n <=== buf_IIJJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_IIJJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_IIJJ, h);
         for (int pq = 0; pq < buf_IIJJ.params->rowtot[h]; ++pq) {
              int p = buf_IIJJ.params->roworb[h][pq][0];
              int q = buf_IIJJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_IIJJ.params->coltot[h]; ++rs) {
                   int r = buf_IIJJ.params->colorb[h][rs][0];
                   int s = buf_IIJJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_IIJJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_IIJJ, h);
    }

    psi::outfile->Printf("\n <=== buf_IJJJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_IJJJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_IJJJ, h);
         for (int pq = 0; pq < buf_IJJJ.params->rowtot[h]; ++pq) {
              int p = buf_IJJJ.params->roworb[h][pq][0];
              int q = buf_IJJJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_IJJJ.params->coltot[h]; ++rs) {
                   int r = buf_IJJJ.params->colorb[h][rs][0];
                   int s = buf_IJJJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_IJJJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_IJJJ, h);
    }

    psi::outfile->Printf("\n <=== buf_IJIJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_IJIJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_IJIJ, h);
         for (int pq = 0; pq < buf_IJIJ.params->rowtot[h]; ++pq) {
              int p = buf_IJIJ.params->roworb[h][pq][0];
              int q = buf_IJIJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_IJIJ.params->coltot[h]; ++rs) {
                   int r = buf_IJIJ.params->colorb[h][rs][0];
                   int s = buf_IJIJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_IJIJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_IJIJ, h);
    }

    psi::outfile->Printf("\n <=== buf_XIJJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_XIJJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_XIJJ, h);
         for (int pq = 0; pq < buf_XIJJ.params->rowtot[h]; ++pq) {
              int p = buf_XIJJ.params->roworb[h][pq][0];
              int q = buf_XIJJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_XIJJ.params->coltot[h]; ++rs) {
                   int r = buf_XIJJ.params->colorb[h][rs][0];
                   int s = buf_XIJJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_XIJJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_XIJJ, h);
    }
    psi::outfile->Printf("\n <=== buf_YIJJ MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_YIJJ, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_YIJJ, h);
         for (int pq = 0; pq < buf_YIJJ.params->rowtot[h]; ++pq) {
              int p = buf_YIJJ.params->roworb[h][pq][0];
              int q = buf_YIJJ.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_YIJJ.params->coltot[h]; ++rs) {
                   int r = buf_YIJJ.params->colorb[h][rs][0];
                   int s = buf_YIJJ.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_YIJJ.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_YIJJ, h);
    }
    psi::outfile->Printf("\n <=== buf_XJII MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_XJII, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_XJII, h);
         for (int pq = 0; pq < buf_XJII.params->rowtot[h]; ++pq) {
              int p = buf_XJII.params->roworb[h][pq][0];
              int q = buf_XJII.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_XJII.params->coltot[h]; ++rs) {
                   int r = buf_XJII.params->colorb[h][rs][0];
                   int s = buf_XJII.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_XJII.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_XJII, h);
    }
    psi::outfile->Printf("\n <=== buf_YJII MO Integrals ===>\n\n");
    for (int h = 0; h < shared_from_this()->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_YJII, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_YJII, h);
         for (int pq = 0; pq < buf_YJII.params->rowtot[h]; ++pq) {
              int p = buf_YJII.params->roworb[h][pq][0];
              int q = buf_YJII.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_YJII.params->coltot[h]; ++rs) {
                   int r = buf_YJII.params->colorb[h][rs][0];
                   int s = buf_YJII.params->colorb[h][rs][1];
                   psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_YJII.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_YJII, h);
    }

    //
    global_dpd_->buf4_close(&buf_IIJJ);
    global_dpd_->buf4_close(&buf_IJJJ);
    global_dpd_->buf4_close(&buf_IIIJ);
    global_dpd_->buf4_close(&buf_IJIJ);
    global_dpd_->buf4_close(&buf_XIJJ);
    global_dpd_->buf4_close(&buf_YIJJ);
    global_dpd_->buf4_close(&buf_XJII);
    global_dpd_->buf4_close(&buf_YJII);

    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
}

void WavefunctionUnion::print_header() {
    psi::outfile->Printf("\n  ===> Wavefunction Union <===\n\n");
    psi::outfile->Printf("        nfrag = %4d\n", nIsolatedMolecules_);
    for (int nf=0; nf<nIsolatedMolecules_; ++nf) {
    psi::outfile->Printf("        mol %d: nbf = %4d;  nmo = %4d;  docc = %4d;  nvir = %4d\n", 
                                   nf+1, l_nbf_[nf], l_nmo_[nf], l_ndocc_[nf], l_nvir_[nf] );
    }
    psi::outfile->Printf("\n");
}

double WavefunctionUnion::nuclear_repulsion_interaction_energy()
{
  double e_nuc_nuc = 0.0;
  for (int x=0; x < l_molecule_[0]->natom(); ++x) {
       for (int y=0; y < l_molecule_[1]->natom(); ++y) {
            double rxy = sqrt(pow(l_molecule_[0]->x(x) - l_molecule_[1]->x(y), 2.0) +
                              pow(l_molecule_[0]->y(x) - l_molecule_[1]->y(y), 2.0) +
                              pow(l_molecule_[0]->z(x) - l_molecule_[1]->z(y), 2.0) );
            e_nuc_nuc += (double)l_molecule_[0]->Z(x) *
                         (double)l_molecule_[1]->Z(y) / rxy;
       }
  }
  return e_nuc_nuc;
}

// Copied from psi::Wavefunction
SharedMatrix WavefunctionUnion::Ca_subset(const std::string &basis, const std::string &subset)
{
    return C_subset_helper(Ca_, nalphapi_, epsilon_a_, basis, subset);
}
// Copied from psi::Wavefunction
SharedMatrix WavefunctionUnion::Cb_subset(const std::string &basis, const std::string &subset)
{
    return C_subset_helper(Cb_, nalphapi_, epsilon_b_, basis, subset);
}

// Copied from psi::Wavefunction. One line (sorting of orbitals according to energy) was removed.
SharedMatrix WavefunctionUnion::C_subset_helper(SharedMatrix C, const Dimension &noccpi, SharedVector epsilon, const std::string &basis, const std::string &subset)
{
    std::vector <std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int) positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }
    SharedMatrix C2(new Matrix("C " + basis + " " + subset, nsopi_, nmopi));
    for (int h = 0; h < (int) positions.size(); h++) {
        for (int i = 0; i < (int) positions[h].size(); i++) {
            C_DCOPY(nsopi_[h], &C->pointer(h)[0][positions[h][i]], nmopi_[h], &C2->pointer(h)[0][i], nmopi[h]);
        }
    }

    if (basis == "AO") {

        SharedMatrix C3(new Matrix("C " + basis + " " + subset, nso_, nmopi.sum()));
        std::swap(C2, C3);

        std::vector <std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        // The following line of code has to be erased - sorting cannot be used in WavefunctionUnion!
        // std::sort(order.begin(), order.end(), std::less < std::tuple < double, int, int > > ());

        for (int index = 0; index < (int) order.size(); index++) {
            int i = std::get<1>(order[index]);
            int h = std::get<2>(order[index]);

            int nao = nso_;
            int nso = nsopi_[h];

            if (!nso) continue;

            C_DGEMV('N', nao, nso, 1.0, AO2SO_->pointer(h)[0], nso, &C3->pointer(h)[0][i], nmopi[h], 0.0, &C2->pointer()[0][index], nmopi.sum());
        }

    } else if (basis == "SO" || basis == "MO") {
        // Already done
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}

// Copied from psi::Wavefunction. One line (sorting of orbitals according to energy) was removed.
SharedVector WavefunctionUnion::epsilon_subset_helper(SharedVector epsilon, const Dimension &noccpi, const std::string &basis, const std::string &subset)
{
    std::vector <std::vector<int>> positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < (int) positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }

    SharedVector C2;

    if (basis == "AO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi.sum()));

        std::vector <std::tuple<double, int, int>> order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                order.push_back(std::tuple<double, int, int>(epsilon->get(h, positions[h][i]), i, h));
            }
        }

        // This line was removed because WavefunctionUnion cannot have sorted orbital energies 
        // std::sort(order.begin(), order.end(), std::less < std::tuple < double, int, int > > ());

        for (int index = 0; index < (int) order.size(); index++) {
            C2->set(0, index, std::get<0>(order[index]));
        }

    } else if (basis == "SO" || basis == "MO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi));
        for (int h = 0; h < (int) positions.size(); h++) {
            for (int i = 0; i < (int) positions[h].size(); i++) {
                C2->set(h, i, epsilon->get(h, positions[h][i]));
            }
        }

    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}
void WavefunctionUnion::clear_dpd() {
  // Destruct the IntegralTransform object
  this->integrals_.reset();
}



} // EndNameSpace oepdev
