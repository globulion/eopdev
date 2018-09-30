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
      double rxxxx_n = rxxx_n * rx_n;  double rxxxx_o = rxxx_o * rx_o;
      double rxxxy_n = rxxx_n * ry_n;  double rxxxy_o = rxxx_o * ry_o;
      double rxxxz_n = rxxx_n * rz_n;  double rxxxz_o = rxxx_o * rz_o;
      double rxxyy_n = rxxy_n * ry_n;  double rxxyy_o = rxxy_o * ry_o;
      double rxxyz_n = rxxy_n * rz_n;  double rxxyz_o = rxxy_o * rz_o;
      double rxxzz_n = rxxz_n * rz_n;  double rxxzz_o = rxxz_o * rz_o;
      double rxyyy_n = rxyy_n * ry_n;  double rxyyy_o = rxyy_o * ry_o;
      double rxyyz_n = rxyy_n * rz_n;  double rxyyz_o = rxyy_o * rz_o;
      double rxyzz_n = rxyz_n * rz_n;  double rxyzz_o = rxyz_o * rz_o;
      double rxzzz_n = rxzz_n * rz_n;  double rxzzz_o = rxzz_o * rz_o;
      double ryyyy_n = ryyy_n * ry_n;  double ryyyy_o = ryyy_o * ry_o;
      double ryyyz_n = ryyy_n * rz_n;  double ryyyz_o = ryyy_o * rz_o;
      double ryyzz_n = ryyz_n * rz_n;  double ryyzz_o = ryyz_o * rz_o;
      double ryzzz_n = ryzz_n * rz_n;  double ryzzz_o = ryzz_o * rz_o;
      double rzzzz_n = rzzz_n * rz_n;  double rzzzz_o = rzzz_o * rz_o;

      double d4xxxx = rxxxx_n - rxxxx_o;
      double d4xxxy = rxxxy_n - rxxxy_o;
      double d4xxxz = rxxxz_n - rxxxz_o;
      double d4xxyy = rxxyy_n - rxxyy_o;
      double d4xxyz = rxxyz_n - rxxyz_o;
      double d4xxzz = rxxzz_n - rxxzz_o;
      double d4xyyy = rxyyy_n - rxyyy_o;
      double d4xyyz = rxyyz_n - rxyyz_o;
      double d4xyzz = rxyzz_n - rxyzz_o;
      double d4xzzz = rxzzz_n - rxzzz_o;
      double d4yyyy = ryyyy_n - ryyyy_o;
      double d4yyyz = ryyyz_n - ryyyz_o;
      double d4yyzz = ryyzz_n - ryyzz_o;
      double d4yzzz = ryzzz_n - ryzzz_o;
      double d4zzzz = rzzzz_n - rzzzz_o;

      double hxxxx_o = hp[ic][ 0];
      double hxxxy_o = hp[ic][ 1];
      double hxxxz_o = hp[ic][ 2];
      double hxxyy_o = hp[ic][ 3];
      double hxxyz_o = hp[ic][ 4];
      double hxxzz_o = hp[ic][ 5];
      double hxyyy_o = hp[ic][ 6];
      double hxyyz_o = hp[ic][ 7];
      double hxyzz_o = hp[ic][ 8];
      double hxzzz_o = hp[ic][ 9];
      double hyyyy_o = hp[ic][10];
      double hyyyz_o = hp[ic][11];
      double hyyzz_o = hp[ic][12];
      double hyzzz_o = hp[ic][13];
      double hzzzz_o = hp[ic][14];

      double hxxxx_n = hxxxx_o + c * d4xxxx - mx_o * d3xxx - mx_o * d3xxx - mx_o * d3xxx - mx_o * d3xxx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx + qxx_o * d2xx - oxxx_o * d1x - oxxx_o * d1x - oxxx_o * d1x - oxxx_o * d1x; 
      double hxxxy_n = hxxxy_o + c * d4xxxy - mx_o * d3xxy - mx_o * d3xxy - mx_o * d3xxy - my_o * d3xxx + qxx_o * d2xy + qxx_o * d2xy + qxy_o * d2xx + qxy_o * d2xx + qxy_o * d2xx + qxx_o * d2xy - oxxx_o * d1y - oxxy_o * d1x - oxxy_o * d1x - oxxy_o * d1x;
      double hxxxz_n = hxxxz_o + c * d4xxxz - mx_o * d3xxz - mx_o * d3xxz - mx_o * d3xxz - mz_o * d3xxx + qxx_o * d2xz + qxx_o * d2xz + qxz_o * d2xx + qxz_o * d2xx + qxz_o * d2xx + qxx_o * d2xz - oxxx_o * d1z - oxxz_o * d1x - oxxz_o * d1x - oxxz_o * d1x;
      double hxxyy_n = hxxyy_o + c * d4xxyy - mx_o * d3xyy - mx_o * d3xyy - my_o * d3xxy - my_o * d3xxy + qxx_o * d2yy + qxy_o * d2xy + qxy_o * d2xy + qyy_o * d2xx + qxy_o * d2xy + qxy_o * d2xy - oxxy_o * d1y - oxxy_o * d1y - oxyy_o * d1x - oxyy_o * d1x;
      double hxxyz_n = hxxyz_o + c * d4xxyz - mx_o * d3xyz - mx_o * d3xyz - my_o * d3xxz - mz_o * d3xxy + qxx_o * d2yz + qxy_o * d2xz + qxz_o * d2xy + qyz_o * d2xx + qxz_o * d2xy + qxy_o * d2xz - oxxy_o * d1z - oxxz_o * d1y - oxyz_o * d1x - oxyz_o * d1x;
      double hxxzz_n = hxxzz_o + c * d4xxzz - mx_o * d3xzz - mx_o * d3xzz - mz_o * d3xxz - mz_o * d3xxz + qxx_o * d2zz + qxz_o * d2xz + qxz_o * d2xz + qzz_o * d2xx + qxz_o * d2xz + qxz_o * d2xz - oxxz_o * d1z - oxxz_o * d1z - oxzz_o * d1x - oxzz_o * d1x;
      double hxyyy_n = hxyyy_o + c * d4xyyy - mx_o * d3yyy - my_o * d3xyy - my_o * d3xyy - my_o * d3xyy + qxy_o * d2yy + qxy_o * d2yy + qyy_o * d2xy + qyy_o * d2xy + qxy_o * d2yy + qyy_o * d2xy - oxyy_o * d1y - oxyy_o * d1y - oxyy_o * d1y - oyyy_o * d1x;
      double hxyyz_n = hxyyz_o + c * d4xyyz - mx_o * d3yyz - my_o * d3xyz - my_o * d3xyz - mz_o * d3xyy + qxy_o * d2yz + qxy_o * d2yz + qyz_o * d2xy + qyz_o * d2xy + qxz_o * d2yy + qyy_o * d2xz - oxyy_o * d1z - oxyz_o * d1y - oxyz_o * d1y - oyyz_o * d1x;
      double hxyzz_n = hxyzz_o + c * d4xyzz - mx_o * d3yzz - my_o * d3xzz - mz_o * d3xyz - mz_o * d3xyz + qxy_o * d2zz + qxz_o * d2yz + qyz_o * d2xz + qzz_o * d2xy + qxz_o * d2yz + qyz_o * d2xz - oxyz_o * d1z - oxyz_o * d1z - oxzz_o * d1y - oyzz_o * d1x;
      double hxzzz_n = hxzzz_o + c * d4xzzz - mx_o * d3zzz - mz_o * d3xzz - mz_o * d3xzz - mz_o * d3xzz + qxz_o * d2zz + qxz_o * d2zz + qzz_o * d2xz + qzz_o * d2xz + qxz_o * d2zz + qzz_o * d2xz - oxzz_o * d1z - oxzz_o * d1z - oxzz_o * d1z - ozzz_o * d1x;
      double hyyyy_n = hyyyy_o + c * d4yyyy - my_o * d3yyy - my_o * d3yyy - my_o * d3yyy - my_o * d3yyy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy + qyy_o * d2yy - oyyy_o * d1y - oyyy_o * d1y - oyyy_o * d1y - oyyy_o * d1y;
      double hyyyz_n = hyyyz_o + c * d4yyyz - my_o * d3yyz - my_o * d3yyz - my_o * d3yyz - mz_o * d3yyy + qyy_o * d2yz + qyy_o * d2yz + qyz_o * d2yy + qyz_o * d2yy + qyz_o * d2yy + qyy_o * d2yz - oyyy_o * d1z - oyyz_o * d1y - oyyz_o * d1y - oyyz_o * d1y;
      double hyyzz_n = hyyzz_o + c * d4yyzz - my_o * d3yzz - my_o * d3yzz - mz_o * d3yyz - mz_o * d3yyz + qyy_o * d2zz + qyz_o * d2yz + qyz_o * d2yz + qzz_o * d2yy + qyz_o * d2yz + qyz_o * d2yz - oyyz_o * d1z - oyyz_o * d1z - oyzz_o * d1y - oyzz_o * d1y;
      double hyzzz_n = hyzzz_o + c * d4yzzz - my_o * d3zzz - mz_o * d3yzz - mz_o * d3yzz - mz_o * d3yzz + qyz_o * d2zz + qyz_o * d2zz + qzz_o * d2yz + qzz_o * d2yz + qyz_o * d2zz + qzz_o * d2yz - oyzz_o * d1z - oyzz_o * d1z - oyzz_o * d1z - ozzz_o * d1y;
      double hzzzz_n = hzzzz_o + c * d4zzzz - mz_o * d3zzz - mz_o * d3zzz - mz_o * d3zzz - mz_o * d3zzz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz + qzz_o * d2zz - ozzz_o * d1z - ozzz_o * d1z - ozzz_o * d1z - ozzz_o * d1z;

      // Collect
      hp_[ic][ 0] = hxxxx_n; 
      hp_[ic][ 1] = hxxxy_n;
      hp_[ic][ 2] = hxxxz_n;
      hp_[ic][ 3] = hxxyy_n;
      hp_[ic][ 4] = hxxyz_n;
      hp_[ic][ 5] = hxxzz_n;
      hp_[ic][ 6] = hxyyy_n;
      hp_[ic][ 7] = hxyyz_n;
      hp_[ic][ 8] = hxyzz_n;
      hp_[ic][ 9] = hxzzz_n;
      hp_[ic][10] = hyyyy_n;
      hp_[ic][11] = hyyyz_n;
      hp_[ic][12] = hyyzz_n;
      hp_[ic][13] = hyzzz_n;
      hp_[ic][14] = hzzzz_n;

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
std::vector<double> DMTPole::potential(std::shared_ptr<DMTPole> other, const std::string& type)
{
  std::vector<double> potentials;
  for (int i=0; i<nDMTPs_; ++i) potentials.push_back(0.0);
  //TODO
  // ... if (type == "R-5") ...
  // Return
  return potentials;
}

// abstract methods
void DMTPole::compute(psi::SharedMatrix D, bool transition, int i) {/* nothing to implement here */}

} // EndNameSpace oepdev
