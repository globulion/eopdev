#include "dmtp.h"


namespace oepdev{

using namespace std;

DMTPole::DMTPole(psi::SharedMolecule mol, int n) 
 : mol_(mol), wfn_(nullptr), nDMTPs_(n),
   name_("none"),
   nCentres_(0),
   nOrigins_(0),
   hasHexadecapoles_(false),
   centres_(nullptr),
   origins_(nullptr),
   charges_({}),
   dipoles_({}),
   quadrupoles_({}),
   octupoles_({}),
   hexadecapoles_({})
{

}
DMTPole::DMTPole(psi::SharedWavefunction wfn, int n) 
 : DMTPole(wfn->molecule(), n)
{
  wfn_ = wfn;
}
DMTPole::~DMTPole()
{

}
std::shared_ptr<DMTPole> DMTPole::build(std::shared_ptr<psi::Wavefunction> wfn, 
                                        const std::string& type, 
                                        int n)
{
  std::shared_ptr<DMTPole> dmtp;
  if (type == "CAMM") dmtp = std::make_shared<CAMM>(wfn, n);
  else throw psi::PSIEXCEPTION("Invalid DMTP type requested.");
  return dmtp;
}

void DMTPole::allocate()
{
  for (int i=0; i<nDMTPs_; ++i) {
       charges_      .push_back( std::make_shared<psi::Matrix>("DMTP 0-th order tensor", nCentres_  , 1 ) );
       dipoles_      .push_back( std::make_shared<psi::Matrix>("DMTP 1-th order tensor", nCentres_  , 3 ) );
       quadrupoles_  .push_back( std::make_shared<psi::Matrix>("DMTP 2-th order tensor", nCentres_  , 6 ) );
       octupoles_    .push_back( std::make_shared<psi::Matrix>("DMTP 3-th order tensor", nCentres_  , 10) );
       if (hasHexadecapoles_) 
       hexadecapoles_.push_back( std::make_shared<psi::Matrix>("DMTP 4-th order tensor", nCentres_  , 15) );
  }
}

void DMTPole::compute(std::vector<psi::SharedMatrix> D) {
 if (D.size() != nDMTPs_) throw psi::PSIEXCEPTION("The number of OED's does not match the allocated size of DMTP object!");
 for (int i=0; i<nDMTPs_; ++i) compute(D.at(i), i);
}
// abstract methods
void DMTPole::compute(psi::SharedMatrix D, int n) {
}

} // EndNameSpace oepdev
