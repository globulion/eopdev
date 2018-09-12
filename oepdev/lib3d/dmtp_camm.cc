#include "dmtp.h"

namespace oepdev{

using namespace std;

CAMM::CAMM(psi::SharedWavefunction wfn, int n) 
 : DMTPole(wfn, n)
{
  this->allocate();
}
CAMM::~CAMM()
{

}

void CAMM::compute(psi::SharedMatrix D, int n) {
 // TODO
}

} // EndNameSpace oepdev
