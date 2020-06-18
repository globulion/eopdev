#include <iostream>
#include "gefp.h"

using namespace std;

oepdev::FragmentedSystem::FragmentedSystem(
             std::vector<SharedGenEffFrag> bsm, std::vector<int> ind) 
 : aggregate_({}), 
   basis_prim_({}),
   basis_aux_({}),
   nfrag_(ind.size()), bsm_(bsm), ind_(ind), fragments_({}) 
{
  // Initialize all the fragments
  for (int i=0; i<this->nfrag_; ++i) {
       SharedGenEffFrag fragment = this->bsm_[this->ind_[i]]->clone();
       fragments_.push_back(fragment);
  }
}
oepdev::FragmentedSystem::~FragmentedSystem(){}

oepdev::SharedFragmentedSystem oepdev::FragmentedSystem::build(std::vector<std::shared_ptr<oepdev::GenEffFrag>> bsm, std::vector<int> ind) {
    oepdev::SharedFragmentedSystem system = std::make_shared<oepdev::FragmentedSystem>(bsm, ind);
    return system;
}

void oepdev::FragmentedSystem::superimpose() {
   if (this->aggregate_.empty()) {
       throw psi::PSIEXCEPTION("You want to compute total energy in FragmentedSystem but have provided no target molecules!\n");
   }
   /* TODO The below part is a temporary handling of molecule and basis sets.
    * Once rotation and translation of those objects is available in Psi4 or OEPDev,
    * the three first steps in the loop will be no longer required leaving only fragment superimposition,
    * which stores all these information.
    */
   for (int i=0; i<this->nfrag_; ++i) {
        fragments_[i]->set_molecule(aggregate_[i]);
        if (!(basis_prim_.empty())) fragments_[i]->set_basisset("primary", basis_prim_[i]);
        if (!(basis_aux_.empty())) fragments_[i]->set_basisset("auxiliary", basis_aux_[i]);
        fragments_[i]->superimpose();
   }
}

double oepdev::FragmentedSystem::compute_energy_term(std::string theory, bool manybody) {
   this->superimpose();
   double energy = oepdev::GenEffFrag::compute_energy_term(theory, fragments_, manybody);
   return energy;
}
double oepdev::FragmentedSystem::compute_energy(std::string theory) {
   this->superimpose();
   double energy = oepdev::GenEffFrag::compute_energy(theory, fragments_);
   return energy;
}

