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
 aompOBI->set_origin(psi::Vector3());
 aompOBI->compute(mpInts_);
}
void DMTPole::recenter(psi::SharedMatrix new_origins, int i)
{
 psi::SharedMatrix dipoles_new = std::make_shared<psi::Matrix>(dipoles_[i]);
 psi::SharedMatrix qdpoles_new = std::make_shared<psi::Matrix>(quadrupoles_[i]);
 psi::SharedMatrix ocpoles_new = std::make_shared<psi::Matrix>(octupoles_[i]);
 psi::SharedMatrix hdpoles_new = std::make_shared<psi::Matrix>(hexadecapoles_[i]);

 double** cp = charges_[i]->pointer();
 double** mp = dipoles_[i]->pointer();
 double** mp_= dipoles_new->pointer();
 double** qp = quadrupoles_[i]->pointer();
 double** qp_= qdpoles_new->pointer();
 double** op = octupoles_[i]->pointer();
 double** op_= ocpoles_new->pointer();
 double** hd = hexadecapoles_[i]->pointer();
 double** hd_= hdpoles_new->pointer();

 // Recenter
 for (int ic=0; ic<nOrigins_; ++ic) {
      double rx_o = origins_->get(ic, 0);
      double ry_o = origins_->get(ic, 1);
      double rz_o = origins_->get(ic, 2);
      double rx_n = new_origins->get(ic, 0);
      double ry_n = new_origins->get(ic, 1);
      double rz_n = new_origins->get(ic, 2);

      double d1x = rx_n - rx_o;
      double d1y = ry_n - ry_o;
      double d1z = rz_n - rz_o;
      double d2xx = rx_n * rx_n - rx_o * rx_o;
      double d2xy = rx_n * ry_n - rx_o * ry_o;
      double d2xz = rx_n * rz_n - rx_o * rz_o;
      double d2yy = ry_n * ry_n - ry_o * ry_o;
      double d2yz = ry_n * rz_n - ry_o * rz_o;
      double d2zz = rz_n * rz_n - rz_o * rz_o;

      double c = cp[ic][0];

      double mx = mp[ic][0] - c * d1x;
      double my = mp[ic][1] - c * d1y;
      double mz = mp[ic][2] - c * d1z;

      double qxx = qp[ic][0] + c * d2xx - 2.0 * mp[ic][0] * d1x;
      double qxy = qp[ic][1] + c * d2xy -       mp[ic][0] * d1y - mp[ic][1] * d1x;
      double qxz = qp[ic][2] + c * d2xz -       mp[ic][0] * d1z - mp[ic][2] * d1x;
      double qyy = qp[ic][3] + c * d2yy - 2.0 * mp[ic][1] * d1y;
      double qyz = qp[ic][4] + c * d2yz -       mp[ic][1] * d1z - mp[ic][2] * d1y; 
      double qzz = qp[ic][5] + c * d2zz - 2.0 * mp[ic][2] * d1z;

      // Collect
      mp_[ic][0] = mx;
      mp_[ic][1] = my;
      mp_[ic][2] = mz;

      qp_[ic][0] = qxx;
      qp_[ic][1] = qxy;
      qp_[ic][2] = qxz;
      qp_[ic][3] = qyy;
      qp_[ic][4] = qyz;
      qp_[ic][5] = qzz;
 }

 // Save
 dipoles_[i]->copy(dipoles_new);
 quadrupoles_[i]->copy(qdpoles_new);
}
void DMTPole::recenter(psi::SharedMatrix new_origins)
{
 for (int i=0; i<nDMTPs_; ++i) this->recenter(new_origins, i);
 origins_->copy(new_origins);
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
