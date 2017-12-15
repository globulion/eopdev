#include "util.h"

namespace oepdev_libutil{

using namespace psi;
using namespace std;


WavefunctionUnion::WavefunctionUnion(SharedWavefunction ref_wfn, Options& options) 
   : Wavefunction(options)
{
   shallow_copy(ref_wfn);
   common_init();
}

WavefunctionUnion::~WavefunctionUnion() 
{
}

void WavefunctionUnion::common_init() {
   if (nirrep()>1) throw PSIEXCEPTION(" ERROR: WavefunctionUnion: ONLY C1 SYMMETRY IS SUPPORTED!\n");
   SharedSuperFunctional functional = create_superfunctional("HF", options_);
   molecule_1_ = extract_monomer(molecule_, 1);
   molecule_2_ = extract_monomer(molecule_, 2);
   molecule_1_->set_name("Monomer 1");
   molecule_2_->set_name("Monomer 2");
   molecule_  ->set_name("Aggregate (Dimer)");
   primary_1_  = basissets_["BASIS_1"];
   primary_2_  = basissets_["BASIS_2"];
   wfn_1_  = solve_scf(molecule_1_, primary_1_, functional, options_, psio_);
   wfn_2_  = solve_scf(molecule_2_, primary_2_, functional, options_, psio_);

   // Monomeric properties (implement all later!)
   nmo_1_       = wfn_1_->nmo()               ; nmo_2_       = wfn_2_->nmo()                ;
   nso_1_       = wfn_1_->nso()               ; nso_2_       = wfn_2_->nso()                ;
   nalpha_1_    = wfn_1_->nalpha()            ; nalpha_2_    = wfn_2_->nalpha()             ;
   nbeta_1_     = wfn_1_->nbeta()             ; nbeta_2_     = wfn_2_->nbeta()              ;
   energy_1_    = wfn_1_->reference_energy()  ; energy_2_    = wfn_2_->reference_energy()   ;

   // Properties of the union
   energy_      = energy_1_ + energy_2_;

   //dimer_wavefunction_     = ...; // store the original wavefunction

   // Set the reference wavefunction to original wavefunction
   // Next, the matrices of the original wavefunction will be 
   // overriden by new matrices
   reference_wavefunction_ = reference_wavefunction();
   // replace the vectors and martrices

}

double WavefunctionUnion::compute_energy() {}

} // EndNameSpace oepdev_libutil
