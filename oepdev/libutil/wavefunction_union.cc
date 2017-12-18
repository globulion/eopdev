#include "util.h"

namespace oepdev_libutil{

using namespace psi;
using namespace std;


WavefunctionUnion::WavefunctionUnion(SharedWavefunction ref_wfn, Options& options) 
   : Wavefunction(options)
{
   // Sanity-check
   if (ref_wfn->molecule()->nfragments()==1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. Molecule has to have minimum two fragments.\n");
   if (ref_wfn->nirrep()>1) throw 
       PSIEXCEPTION(" OEPDEV: Error WavefunctionUnion init. C1 symmetry must be enforced.\n");
   // Add sanity check for multiplicity (closed-shells only!)

   // Initialize the object
   shallow_copy(ref_wfn);
   common_init();
}

WavefunctionUnion::~WavefunctionUnion() 
{
}

void WavefunctionUnion::common_init() {
   SharedSuperFunctional functional = create_superfunctional("HF", options_);
   SharedMolecule molecule_1 = extract_monomer(molecule_, 1);
   SharedMolecule molecule_2 = extract_monomer(molecule_, 2);
   molecule_1->set_name("Monomer 1");
   molecule_2->set_name("Monomer 2");
   molecule  ->set_name("Aggregate (Dimer)");
   SharedBasisSet primary_1  = basissets_["BASIS_1"];
   SharedBasisSet primary_2  = basissets_["BASIS_2"];
   SharedWavefunction wfn_1  = solve_scf(molecule_1, primary_1, functional, options_, psio_);
   SharedWavefunction wfn_2  = solve_scf(molecule_2, primary_2, functional, options_, psio_);

   // Properties of the union
   nIsolatedMolecules_ = ref_wfn->molecule()->nfragments();
   energy_             = wfn_1->energy() + wfn_2->energy();
   nfzc_               = wfn_1->efzc() + wfn_2->efzc();

   // Properties of isolated monomers organize in vectors
   l_energy_        .push_back(wfn_1->energy()         ); l_energy_        .push_back(wfn_2->energy()             );  
   l_molecule_      .push_back(molecule_1              ); l_molecule_      .push_back(molecule_2                  );
   l_wfn_           .push_back(wfn_1                   ); l_wfn_           .push_back(wfn_2                       );
   l_primary_       .push_back(primary_1               ); l_primary_       .push_back(primary_2                   );
 //l_auxiliary_     .push_back(auxiliary_1             ); l_auxiliary_     .push_back(auxiliary_2                 );
   l_name_          .push_back(wfn_1->name()           ); l_name_          .push_back(wfn_2->name()               );
   l_nmo_           .push_back(wfn_1->nmo()            ); l_nmo_           .push_back(wfn_2->nmo()                );
   l_nso_           .push_back(wfn_1->nso()            ); l_nso_           .push_back(wfn_2->nso()                );
   l_efzc_          .push_back(wfn_1->efzc()           ); l_efzc_          .push_back(wfn_2->efzc()               );
   l_density_fitted_.push_back(wfn_1->density_fitted() ); l_density_fitted_.push_back(wfn_2->density_fitted()     );
   l_nalpha_        .push_back(wfn_1->nalpha()         ); l_nalpha_        .push_back(wfn_2->nalpha()             );  
   l_nbeta_         .push_back(wfn_1->nbeta()          ); l_nbeta_         .push_back(wfn_2->nbeta()              );
   l_nfrzc_         .push_back(wfn_1->nfrzc()          ); l_nfrzc_         .push_back(wfn_2->nfrzc()              );
   
  
   //dimer_wavefunction_     = ...; // store the original wavefunction

   /* Set the reference wavefunction to original wavefunction
    * Next, the matrices of the original wavefunction will be 
    * overriden by new matrices
    */
   reference_wavefunction_ = reference_wavefunction();

   // <==== Replace the vectors and martrices ====> //
   
   // <---- Wavefunction Coefficients ----> //
   if (Ca_ && wfn_1->Ca()) {
       Ca_ = Ca()  ; Cb_ = Cb() ;
       Ca_->zero() ; Cb_->zero();
       double** pCa = Ca_->pointer(); 
       double** pCb = Cb_->pointer();

       int nOffsetAO, nOffsetMO = 0, 0;
       //int nbf, ndocc, nvir, nmo;
       int nbf, nmo;
       for (int nf=0; nf<nIsolatedMolecules_; nf++) {
            // below implementation is for CLOSED SHELLS only!
            nbf  = l_wfn_[nf]->basisset()->nbf();
            nmo  = l_wfn_[nf]->nmo()/2;
            //ndocc= l_wfn_[nf]->doccpi()[0];
            //nvir = l_wfn_[nf]->nmo() - ndocc;
            for (int i=0; i<nbf; i++) {      
                 for (int j=0; j<nmo; j++) {
                      pCa[i][j] = l_wfn_[nf]->Ca()->get(0, i+nOffsetAO, j+nOffsetMO);
                 }
            }
            nOffsetAO += nbf;
            nOffsetMO += nmo;
       }
       Cb_->copy(Ca_); // Because implemented for closed-shell only!
   }
}

void Wavefunction::insert_block_(double** p) 
{
}


double WavefunctionUnion::compute_energy() {}

} // EndNameSpace oepdev_libutil
