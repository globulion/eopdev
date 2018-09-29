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
   nSites_(0),
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
       if (hasCharges_      )  {charges_      .push_back( std::make_shared<psi::Matrix>("DMTP 0-th order tensor", nSites_  , 1 ) );}
       else {charges_      .push_back( std::make_shared<psi::Matrix>() );}
       if (hasDipoles_      )  {dipoles_      .push_back( std::make_shared<psi::Matrix>("DMTP 1-th order tensor", nSites_  , 3 ) );}
       else {dipoles_      .push_back( std::make_shared<psi::Matrix>() );}
       if (hasQuadrupoles_  )  {quadrupoles_  .push_back( std::make_shared<psi::Matrix>("DMTP 2-th order tensor", nSites_  , 6 ) );}
       else {quadrupoles_  .push_back( std::make_shared<psi::Matrix>() );}
       if (hasOctupoles_    )  {octupoles_    .push_back( std::make_shared<psi::Matrix>("DMTP 3-th order tensor", nSites_  , 10) );}
       else {octupoles_    .push_back( std::make_shared<psi::Matrix>() );}
       if (hasHexadecapoles_)  {hexadecapoles_.push_back( std::make_shared<psi::Matrix>("DMTP 4-th order tensor", nSites_  , 15) );}
       else {hexadecapoles_.push_back( std::make_shared<psi::Matrix>() );}
  }
  centres_ = std::make_shared<psi::Matrix>("DMTP Centres", nSites_, 3);
  origins_ = std::make_shared<psi::Matrix>("DMTP Origins", nSites_, 3);
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
 double** hp = hexadecapoles_[i]->pointer();
 double** hp_= hdpoles_new->pointer();

 // Recenter
 for (int ic=0; ic<nSites_; ++ic) {
      if (hasDipoles_) {
      double rx_o = origins_->get(ic, 0);
      double ry_o = origins_->get(ic, 1);
      double rz_o = origins_->get(ic, 2);
      double rx_n = new_origins->get(ic, 0);
      double ry_n = new_origins->get(ic, 1);
      double rz_n = new_origins->get(ic, 2);

      double d1x = rx_n - rx_o;
      double d1y = ry_n - ry_o;
      double d1z = rz_n - rz_o;

      double c = cp[ic][0];

      double mx_o = mp[ic][0];
      double my_o = mp[ic][1];
      double mz_o = mp[ic][2];

      double mx_n = mx_o - c * d1x;
      double my_n = my_o - c * d1y;
      double mz_n = mz_o - c * d1z;

      // Collect
      mp_[ic][0] = mx_n;
      mp_[ic][1] = my_n;
      mp_[ic][2] = mz_n;

      if (hasQuadrupoles_) {
      double rxx_n = rx_n * rx_n;  double rxx_o = rx_o * rx_o;
      double rxy_n = rx_n * ry_n;  double rxy_o = rx_o * ry_o;
      double rxz_n = rx_n * rz_n;  double rxz_o = rx_o * rz_o;
      double ryy_n = ry_n * ry_n;  double ryy_o = ry_o * ry_o;
      double ryz_n = ry_n * rz_n;  double ryz_o = ry_o * rz_o;
      double rzz_n = rz_n * rz_n;  double rzz_o = rz_o * rz_o;
      
      double d2xx = rxx_n - rxx_o;
      double d2xy = rxy_n - rxy_o;
      double d2xz = rxz_n - rxz_o;
      double d2yy = ryy_n - ryy_o;
      double d2yz = ryz_n - ryz_o;
      double d2zz = rzz_n - rzz_o;

      double qxx_o = qp[ic][0];
      double qxy_o = qp[ic][1];
      double qxz_o = qp[ic][2];
      double qyy_o = qp[ic][3];
      double qyz_o = qp[ic][4];
      double qzz_o = qp[ic][5];

      double qxx_n = qxx_o + c * d2xx - 2.0 * mx_o * d1x;
      double qxy_n = qxy_o + c * d2xy -       mx_o * d1y - my_o * d1x;
      double qxz_n = qxz_o + c * d2xz -       mx_o * d1z - mz_o * d1x;
      double qyy_n = qyy_o + c * d2yy - 2.0 * my_o * d1y;
      double qyz_n = qyz_o + c * d2yz -       my_o * d1z - mz_o * d1y; 
      double qzz_n = qzz_o + c * d2zz - 2.0 * mz_o * d1z;

      // Collect
      qp_[ic][0] = qxx_n;
      qp_[ic][1] = qxy_n;
      qp_[ic][2] = qxz_n;
      qp_[ic][3] = qyy_n;
      qp_[ic][4] = qyz_n;
      qp_[ic][5] = qzz_n;

      if (hasOctupoles_) {
      double rxxx_n = rxx_n * rx_n;  double rxxx_o = rxx_o * rx_o;  
      double rxxy_n = rxx_n * ry_n;  double rxxy_o = rxx_o * ry_o;
      double rxxz_n = rxx_n * rz_n;  double rxxz_o = rxx_o * rz_o;
      double rxyy_n = rxy_n * ry_n;  double rxyy_o = rxy_o * ry_o;
      double rxyz_n = rxy_n * rz_n;  double rxyz_o = rxy_o * rz_o;
      double rxzz_n = rxz_n * rz_n;  double rxzz_o = rxz_o * rz_o;
      double ryyy_n = ryy_n * ry_n;  double ryyy_o = ryy_o * ry_o;
      double ryyz_n = ryy_n * rz_n;  double ryyz_o = ryy_o * rz_o;
      double ryzz_n = ryz_n * rz_n;  double ryzz_o = ryz_o * rz_o;
      double rzzz_n = rzz_n * rz_n;  double rzzz_o = rzz_o * rz_o;

      double d3xxx = rxxx_n - rxxx_o; 
      double d3xxy = rxxy_n - rxxy_o;
      double d3xxz = rxxz_n - rxxz_o;
      double d3xyy = rxyy_n - rxyy_o;
      double d3xyz = rxyz_n - rxyz_o;
      double d3xzz = rxzz_n - rxzz_o;
      double d3yyy = ryyy_n - ryyy_o;
      double d3yyz = ryyz_n - ryyz_o;
      double d3yzz = ryzz_n - ryzz_o;
      double d3zzz = rzzz_n - rzzz_o;

      double oxxx_o = op[ic][0];
      double oxxy_o = op[ic][1];
      double oxxz_o = op[ic][2];
      double oxyy_o = op[ic][3];
      double oxyz_o = op[ic][4];
      double oxzz_o = op[ic][5];
      double oyyy_o = op[ic][6];
      double oyyz_o = op[ic][7];
      double oyzz_o = op[ic][8];
      double ozzz_o = op[ic][9];

      double oxxx_n = oxxx_o  - c * d3xxx  + mx_o * d2xx + mx_o * d2xx + mx_o * d2xx  - qxx_o * d1x - qxx_o * d1x - qxx_o * d1x;   
      double oxxy_n = oxxy_o  - c * d3xxy  + mx_o * d2xy + my_o * d2xx + mx_o * d2xy  - qxy_o * d1x - qxx_o * d1y - qxy_o * d1x;
      double oxxz_n = oxxz_o  - c * d3xxz  + mx_o * d2xz + mz_o * d2xx + mx_o * d2xz  - qxz_o * d1x - qxx_o * d1z - qxz_o * d1x;
      double oxyy_n = oxyy_o  - c * d3xyy  + mx_o * d2yy + my_o * d2xy + my_o * d2xy  - qxy_o * d1y - qxy_o * d1y - qyy_o * d1x;
      double oxyz_n = oxyz_o  - c * d3xyz  + mx_o * d2yz + mz_o * d2xy + my_o * d2xz  - qxz_o * d1y - qxy_o * d1z - qyz_o * d1x;
      double oxzz_n = oxzz_o  - c * d3xzz  + mx_o * d2zz + mz_o * d2xz + mz_o * d2xz  - qxz_o * d1z - qxz_o * d1z - qzz_o * d1x;
      double oyyy_n = oyyy_o  - c * d3yyy  + my_o * d2yy + my_o * d2yy + my_o * d2yy  - qyy_o * d1y - qyy_o * d1y - qyy_o * d1y;
      double oyyz_n = oyyz_o  - c * d3yyz  + my_o * d2yz + mz_o * d2yy + my_o * d2yz  - qyz_o * d1y - qyy_o * d1z - qyz_o * d1y;
      double oyzz_n = oyzz_o  - c * d3yzz  + my_o * d2zz + mz_o * d2yz + mz_o * d2yz  - qyz_o * d1z - qyz_o * d1z - qzz_o * d1y;
      double ozzz_n = ozzz_o  - c * d3zzz  + mz_o * d2zz + mz_o * d2zz + mz_o * d2zz  - qzz_o * d1z - qzz_o * d1z - qzz_o * d1z;

      // Collect
      op_[ic][0] = oxxx_n;
      op_[ic][1] = oxxy_n;
      op_[ic][2] = oxxz_n;
      op_[ic][3] = oxyy_n;
      op_[ic][4] = oxyz_n;
      op_[ic][5] = oxzz_n;
      op_[ic][6] = oyyy_n;
      op_[ic][7] = oyyz_n;
      op_[ic][8] = oyzz_n;
      op_[ic][9] = ozzz_n;

      if (hasHexadecapoles_) {


      }}}} // EndIFHasMultipoles
 }

 // Save
 if (hasDipoles_      ) dipoles_      [i]->copy(dipoles_new);
 if (hasQuadrupoles_  ) quadrupoles_  [i]->copy(qdpoles_new);
 if (hasOctupoles_    ) octupoles_    [i]->copy(ocpoles_new);
 if (hasHexadecapoles_) hexadecapoles_[i]->copy(hdpoles_new);
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
