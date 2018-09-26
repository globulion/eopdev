#include "dmtp.h"

namespace oepdev{

using namespace std;

CAMM::CAMM(psi::SharedWavefunction wfn, int n) 
 : DMTPole(wfn, n)
{
  // Distribution centers and origins
  nCentres_ = mol_->natom();
  nOrigins_ = nCentres_;
  // Available multipoles
  hasCharges_ = true;
  hasDipoles_ = true;
  hasQuadrupoles_ = true;
  hasOctupoles_ = true;
  hasHexadecapoles_ = true;
  // Allocate memory
  this->allocate();
  this->set_sites();
  // Compute necessary integrals
  this->compute_integrals();
}
CAMM::~CAMM()
{

}
void CAMM::set_sites(void) 
{
  for (int n=0; n<mol_->natom(); ++n) {
       double x = mol_->x(n);
       double y = mol_->y(n);
       double z = mol_->z(n);
       centres_->set(n, 0, x);  origins_->set(n, 0, x);
       centres_->set(n, 1, y);  origins_->set(n, 1, y);
       centres_->set(n, 2, z);  origins_->set(n, 2, z);
  }
}
void CAMM::compute(psi::SharedMatrix D, bool transition, int i) {
 // TODO
}

} // EndNameSpace oepdev
