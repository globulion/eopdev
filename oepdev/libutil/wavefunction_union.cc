#include "util.h"

namespace oepdev_libutil{

using namespace psi;
using namespace std;


WavefunctionUnion::WavefunctionUnion(SharedWavefunction ref_wfn, Options& options) 
   : Wavefunction(options)
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
   common_init(ref_wfn);
}

WavefunctionUnion::~WavefunctionUnion() 
{
}

void WavefunctionUnion::common_init(SharedWavefunction ref_wfn) {
   SharedSuperFunctional functional = create_superfunctional("HF", options_);
   SharedMolecule molecule_1 = extract_monomer(molecule_, 1);
   SharedMolecule molecule_2 = extract_monomer(molecule_, 2);
   molecule_1->set_name("Monomer 1");
   molecule_2->set_name("Monomer 2");
   molecule_ ->set_name("Aggregate (Dimer)");
   SharedBasisSet primary_1  = basissets_["BASIS_1"];
   SharedBasisSet primary_2  = basissets_["BASIS_2"];
   SharedWavefunction wfn_1  = solve_scf(molecule_1, primary_1, functional, options_, psio_);
   SharedWavefunction wfn_2  = solve_scf(molecule_2, primary_2, functional, options_, psio_);

   // Properties of the union
   nIsolatedMolecules_ = ref_wfn->molecule()->nfragments();
   energy_             = wfn_1->reference_energy() + wfn_2->reference_energy();
   efzc_               = wfn_1->efzc() + wfn_2->efzc();

   // Properties of isolated monomers organize in vectors
   l_energy_        .push_back(wfn_1->reference_energy()); l_energy_        .push_back(wfn_2->reference_energy()   );  
   l_molecule_      .push_back(molecule_1               ); l_molecule_      .push_back(molecule_2                  );   
   l_wfn_           .push_back(wfn_1                    ); l_wfn_           .push_back(wfn_2                       );
   l_primary_       .push_back(primary_1                ); l_primary_       .push_back(primary_2                   );
 //l_auxiliary_     .push_back(auxiliary_1              ); l_auxiliary_     .push_back(auxiliary_2                 );
   l_name_          .push_back(wfn_1->name()            ); l_name_          .push_back(wfn_2->name()               );
   l_nmo_           .push_back(wfn_1->nmo()             ); l_nmo_           .push_back(wfn_2->nmo()                );
   l_nso_           .push_back(wfn_1->nso()             ); l_nso_           .push_back(wfn_2->nso()                );
   l_efzc_          .push_back(wfn_1->efzc()            ); l_efzc_          .push_back(wfn_2->efzc()               );
   l_density_fitted_.push_back(wfn_1->density_fitted()  ); l_density_fitted_.push_back(wfn_2->density_fitted()     );
   l_nalpha_        .push_back(wfn_1->nalpha()          ); l_nalpha_        .push_back(wfn_2->nalpha()             );  
   l_nbeta_         .push_back(wfn_1->nbeta()           ); l_nbeta_         .push_back(wfn_2->nbeta()              );
   l_nfrzc_         .push_back(wfn_1->nfrzc()           ); l_nfrzc_         .push_back(wfn_2->nfrzc()              );

  
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
   epsilon_a_->zero() ; epsilon_b_->zero();
   Ca_->zero() ; Cb_->zero();
   double** pCa = Ca_->pointer(); 
   double** pCb = Cb_->pointer();
   double*  pea = epsilon_a_->pointer();
   double*  peb = epsilon_b_->pointer();

   std::vector<std::vector<int>> orbitals_occ, orbitals_vir;
   std::vector<int> dummy;
   /* implementation below is for CLOSED SHELLS only! */
   int nbf, nmo, nmo_occ, nmo_vir;
   int nmo_occ_all = ref_wfn->doccpi()[0];
   int nOffsetAO = 0, nOffsetMOOcc = 0, nOffsetMOVir = nmo_occ_all;
   //
   for (int nf=0; nf<nIsolatedMolecules_; nf++) {
        nbf      = l_wfn_[nf]->basisset()->nbf();
        nmo      = l_wfn_[nf]->nmo();
        nmo_occ  = l_wfn_[nf]->doccpi()[0];
        nmo_vir  = nmo - nmo_occ;
        for (int i=0; i<nbf; i++) {      
             // <--- Occupied Orbitals ---> //
             for (int jo=0; jo<nmo_occ; jo++) {
                  pCa[i+nOffsetAO][jo+nOffsetMOOcc] = l_wfn_[nf]->Ca_subset("AO", "OCC")->get(0, i, jo);
             }
             // <--- Virtual Orbitals ---> //
             for (int jv=0; jv<nmo_vir; jv++) {
                  pCa[i+nOffsetAO][jv+nOffsetMOVir] = l_wfn_[nf]->Ca_subset("AO", "VIR")->get(0, i, jv);
             }
        }
        // <--- Occupied orbital Energies and Indices ---> //
        for (int jo=0; jo<nmo_occ; jo++) {
             pea[jo+nOffsetMOOcc] = l_wfn_[nf]->epsilon_a_subset("AO", "OCC")->get(0, jo);
             dummy.push_back(jo+nOffsetMOOcc);
        }
        orbitals_occ.push_back(dummy); dummy.clear();

        // <--- Virtual orbital Energies and Indices ---> //
        for (int jv=0; jv<nmo_vir; jv++) {
             pea[jv+nOffsetMOVir] = l_wfn_[nf]->epsilon_a_subset("AO", "VIR")->get(0, jv);
             dummy.push_back(jv+nOffsetMOVir);
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
   SharedMOSpace space_1_occ = std::make_shared<MOSpace>('1', orbitals_occ[0], dummy);
   SharedMOSpace space_2_occ = std::make_shared<MOSpace>('2', orbitals_occ[1], dummy);
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

   // <==== Compute One-Electron Property object ====> //
   //oeprop_ = std::make_shared<OEProp>(static_cast<SharedWavefunction>(shared_from_this()));
   oeprop_ = std::make_shared<OEProp>(ref_wfn);
   oeprop_->set_title(" One-Electron Properties of Wavefunction Union");
   oeprop_->set_Da_ao(Da_);
   oeprop_->add("DIPOLE");
   oeprop_->add("QUADRUPOLE");
   oeprop_->add("MULLIKEN CHARGES");
   oeprop_->compute();
}

double WavefunctionUnion::compute_energy() {}

void WavefunctionUnion::localize_orbitals() {
  throw PSIEXCEPTION(" OEPDEV: NotImplementedError. Localization is not yet implemented.\n");
}

void WavefunctionUnion::transform_integrals() 
{
    SharedMOSpaceVector spaces;
    SharedMOSpace space_1o = l_mospace(0,"OCC");
    SharedMOSpace space_2o = l_mospace(1,"OCC");
    spaces.push_back(space_1o);
    spaces.push_back(space_2o);
    integrals_ = std::make_shared<IntegralTransform>(shared_from_this(), spaces,
                                                     IntegralTransform::Restricted,
                                                     IntegralTransform::DPDOnly,
                                                     IntegralTransform::QTOrder,
                                                     IntegralTransform::None);

    integrals_->set_keep_dpd_so_ints(true);

    // Trans (11|11)
    timer_on("Trans (11|11)");
    integrals_->transform_tei(space_1o, space_1o, space_1o, space_1o, IntegralTransform::MakeAndKeep);
    timer_off("Trans (11|11)");

    // Trans (11|12)
    timer_on("Trans (11|12)");
    integrals_->transform_tei(space_1o, space_1o, space_1o, space_2o, IntegralTransform::ReadAndKeep);
    timer_off("Trans (11|12)");

    // Trans (11|22)
    timer_on("Trans (11|22)");
    integrals_->transform_tei(space_1o, space_1o, space_2o, space_2o, IntegralTransform::ReadAndNuke);
    timer_off("Trans (11|22)");

    // Trans (12|12)
    timer_on("Trans (12|12)");
    integrals_->transform_tei(space_1o, space_2o, space_1o, space_2o, IntegralTransform::MakeAndKeep);
    timer_off("Trans (12|12)");

    // Trans (12|22)
    timer_on("Trans (12|22)");
    integrals_->transform_tei(space_1o, space_2o, space_2o, space_2o, IntegralTransform::ReadAndNuke);
    timer_off("Trans (12|22)");

    // Trans (22|22)
    timer_on("Trans (22|22)");
    integrals_->transform_tei(space_2o, space_2o, space_2o, space_2o, IntegralTransform::MakeAndNuke);
    timer_off("Trans (22|22)");
}

} // EndNameSpace oepdev_libutil
