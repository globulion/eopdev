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
   if (ref_wfn->molecule()->multiplicity()!=1) throw
       PSIEXCEPTION(" OEPDEV: NotImplementedError Wavefunction init. So far only CLOSE SHELLS are supported.\n");
   // Add sanity check for multiplicity (closed-shells only!)

   // Initialize the object
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

   /* below implementation is for CLOSED SHELLS only! */
   int nOffsetAO = 0, nOffsetMO = 0;
   int nbf, nmo;
   for (int nf=0; nf<nIsolatedMolecules_; nf++) {
        nbf  = l_wfn_[nf]->basisset()->nbf();
        nmo  = l_wfn_[nf]->nmo();
        for (int i=0; i<nbf; i++) {      
             for (int j=0; j<nmo; j++) {
                  pCa[i+nOffsetAO][j+nOffsetMO] = l_wfn_[nf]->Ca()->get(0, i, j);
             }
        }
        nOffsetAO += nbf;
        nOffsetMO += nmo;
   }
   Cb_       ->copy( Ca_       ); // Because implemented for closed-shell only!
   epsilon_b_->copy(*epsilon_a_); // 
   //outfile->Printf(" Number of MOs: %d\n", l_wfn_[0]->nmo());
   //outfile->Printf(" Number of BFs: %d\n", l_wfn_[0]->basisset()->nbf());
   //l_wfn_[0]->Ca()->print();
   //Ca_->print();

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
   //l_wfn_[1]->Db()->print();
   //Db_->print();
}

double WavefunctionUnion::compute_energy() {}

} // EndNameSpace oepdev_libutil
