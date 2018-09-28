#include "dmtp.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/multipolesymmetry.h"


namespace oepdev{

using namespace std;

DMTPole::DMTPole(psi::SharedWavefunction wfn, int n) 
 : mol_(wfn->molecule()), 
   wfn_(wfn), 
   primary_(wfn->basisset()),
   nDMTPs_(n),
   name_("none"),
   order_(0),
   nCentres_(0),
   nOrigins_(0),
   hasCharges_(false),
   hasDipoles_(false),
   hasQuadrupoles_(false),
   hasOctupoles_(false),
   hasHexadecapoles_(false),
   centres_(nullptr),
   origins_(nullptr),
   charges_({}),
   dipoles_({}),
   quadrupoles_({}),
   octupoles_({}),
   hexadecapoles_({}),
   mpInts_({})
{

}
DMTPole::~DMTPole()
{

}
std::shared_ptr<DMTPole> DMTPole::build(const std::string& type, 
                                        std::shared_ptr<psi::Wavefunction> wfn, 
                                        int n)
{
  std::shared_ptr<DMTPole> dmtp;
  if (type == "CAMM") dmtp = std::make_shared<oepdev::CAMM>(wfn, n);
  else throw psi::PSIEXCEPTION("Invalid DMTP type requested.");
  return dmtp;
}
void DMTPole::allocate()
{
  for (int i=0; i<nDMTPs_; ++i) {
       if (hasCharges_      )  charges_      .push_back( std::make_shared<psi::Matrix>("DMTP 0-th order tensor", nCentres_  , 1 ) );
       if (hasDipoles_      )  dipoles_      .push_back( std::make_shared<psi::Matrix>("DMTP 1-th order tensor", nCentres_  , 3 ) );
       if (hasQuadrupoles_  )  quadrupoles_  .push_back( std::make_shared<psi::Matrix>("DMTP 2-th order tensor", nCentres_  , 6 ) );
       if (hasOctupoles_    )  octupoles_    .push_back( std::make_shared<psi::Matrix>("DMTP 3-th order tensor", nCentres_  , 10) );
       if (hasHexadecapoles_)  hexadecapoles_.push_back( std::make_shared<psi::Matrix>("DMTP 4-th order tensor", nCentres_  , 15) );
  }
  centres_ = std::make_shared<psi::Matrix>("DMTP Centres", nCentres_, 3);
  origins_ = std::make_shared<psi::Matrix>("DMTP Origins", nOrigins_, 3);
}
void DMTPole::compute(std::vector<psi::SharedMatrix> D, std::vector<bool> transition) {
 if (D.size() != nDMTPs_) throw psi::PSIEXCEPTION("The number of OED's does not match the allocated size of DMTP object!");
 for (int i=0; i<nDMTPs_; ++i) this->compute(D.at(i), transition.at(i), i);
}
void DMTPole::compute(void) {
  psi::SharedMatrix D = wfn_->Da(); D->add(wfn_->Db());
  this->compute(D, false, 0);
}
void DMTPole::compute_order(void) {
 order_ = (int)hasCharges_ + (int)hasDipoles_ + (int)hasQuadrupoles_ + (int)hasOctupoles_ + (int)hasHexadecapoles_ - 1;
}
void DMTPole::compute_integrals(void) {
 this->compute_order(); // computes order_
 std::shared_ptr<psi::IntegralFactory> integral = std::make_shared<psi::IntegralFactory>(wfn_->basisset());
 psi::MultipoleSymmetry mpsymm(order_, mol_, integral, wfn_->matrix_factory());
 mpInts_ = mpsymm.create_matrices("Multipole Integrals", true);
 std::shared_ptr<psi::OneBodyAOInt> aompOBI(integral->ao_multipoles(order_));
 aompOBI->compute(mpInts_);
}
void DMTPole::recenter(psi::SharedMatrix new_origins, int i)
{
 //TODO
}
void DMTPole::recenter(psi::SharedMatrix new_origins)
{
 for (int i=0; i<nDMTPs_; ++i) this->recenter(new_origins, i);
}
std::vector<double> DMTPole::energy(std::shared_ptr<DMTPole> other, const std::string& type)
{
  std::vector<double> energies;
  for (int i=0; i<nDMTPs_; ++i) energies.push_back(0.0);
  //TODO
  // ... if (type == "R-5") ...
  // Return
  return energies;
}


// abstract methods
void DMTPole::compute(psi::SharedMatrix D, bool transition, int i) {/* nothing to implement here */}

} // EndNameSpace oepdev
