#include "dmtp.h"

namespace oepdev{

using namespace std;

CAMM::CAMM(psi::SharedWavefunction wfn, int n) 
 : DMTPole(wfn, n)
{
  name_ = "Cumulative Atomic Multipole Moments";
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
       // Centres are set to atomic positions. Origins are initially set to zero.
       centres_->set(n, 0, x);  origins_->set(n, 0, 0.0);
       centres_->set(n, 1, y);  origins_->set(n, 1, 0.0);
       centres_->set(n, 2, z);  origins_->set(n, 2, 0.0);
  }
}
void CAMM::compute(psi::SharedMatrix D, bool transition, int i) {

 double** Dp = D->pointer();
 double** Sp = wfn_->S()->pointer();

 for (int ia = 0; ia < mol_->natom(); ++ia) {

      // Nuclear contribution
      double rx = mol_->x(ia);
      double ry = mol_->y(ia);
      double rz = mol_->z(ia);

      double q = (double)mol_->Z(ia); if (transition) q = 0.0;

      double mx = q * rx;
      double my = q * ry;
      double mz = q * rz;

      double qxx = mx * rx;
      double qxy = mx * ry;
      double qxz = mx * rz;
      double qyy = my * ry;
      double qyz = my * rz;
      double qzz = mz * rz;

      double oxxx = qxx * rx;
      double oxxy = qxx * ry;
      double oxxz = qxx * rz;
      double oxyy = qxy * ry;
      double oxyz = qxy * rz;
      double oxzz = qxz * rz;
      double oyyy = qyy * ry;
      double oyyz = qyy * rz;
      double oyzz = qyz * rz;
      double ozzz = qzz * rz;

      double hxxxx = oxxx * rx;
      double hxxxy = oxxx * ry;
      double hxxxz = oxxx * rz;
      double hxxyy = oxxy * ry;
      double hxxyz = oxxy * rz;
      double hxxzz = oxxz * rz;
      double hxyyy = oxyy * ry;
      double hxyyz = oxyy * rz;
      double hxyzz = oxyz * rz;
      double hxzzz = oxzz * rz;
      double hyyyy = oyyy * ry;
      double hyyyz = oyyy * rz;
      double hyyzz = oyyz * rz;
      double hyzzz = oyzz * rz;
      double hzzzz = ozzz * rz;

      // Electronic contribution
      for (int I = 0; I < primary_->nbf(); ++I) {
           if (ia == primary_->function_to_center(I)) { 
               for (int J = 0; J < primary_->nbf(); ++J) {

                    double  d = Dp[I][J];
                    q    -= d * Sp[I][J];              
                                                       
                    mx   += d * mpInts_[ 0]->get(I,J); 
                    my   += d * mpInts_[ 1]->get(I,J);
                    mz   += d * mpInts_[ 2]->get(I,J);
                                                        
                    qxx  += d * mpInts_[ 3]->get(I,J);
                    qxy  += d * mpInts_[ 4]->get(I,J);
                    qxz  += d * mpInts_[ 5]->get(I,J);
                    qyy  += d * mpInts_[ 6]->get(I,J);
                    qyz  += d * mpInts_[ 7]->get(I,J);
                    qzz  += d * mpInts_[ 8]->get(I,J);
                                                         
                    oxxx += d * mpInts_[ 9]->get(I,J);
                    oxxy += d * mpInts_[10]->get(I,J);
                    oxxz += d * mpInts_[11]->get(I,J);
                    oxyy += d * mpInts_[12]->get(I,J);
                    oxyz += d * mpInts_[13]->get(I,J);
                    oxzz += d * mpInts_[14]->get(I,J);
                    oyyy += d * mpInts_[15]->get(I,J);
                    oyyz += d * mpInts_[16]->get(I,J);
                    oyzz += d * mpInts_[17]->get(I,J);
                    ozzz += d * mpInts_[18]->get(I,J);

                    hxxxx+= d * mpInts_[19]->get(I,J);
                    hxxxy+= d * mpInts_[20]->get(I,J);
                    hxxxz+= d * mpInts_[21]->get(I,J);
                    hxxyy+= d * mpInts_[22]->get(I,J);
                    hxxyz+= d * mpInts_[23]->get(I,J);
                    hxxzz+= d * mpInts_[24]->get(I,J);
                    hxyyy+= d * mpInts_[25]->get(I,J);
                    hxyyz+= d * mpInts_[26]->get(I,J);
                    hxyzz+= d * mpInts_[27]->get(I,J);
                    hxzzz+= d * mpInts_[28]->get(I,J);
                    hyyyy+= d * mpInts_[29]->get(I,J);
                    hyyyz+= d * mpInts_[30]->get(I,J);
                    hyyzz+= d * mpInts_[31]->get(I,J);
                    hyzzz+= d * mpInts_[32]->get(I,J);
                    hzzzz+= d * mpInts_[33]->get(I,J);
               }
           }
      }

      // Save the components 
      charges_      [i]->set(ia, 0, q    );
      dipoles_      [i]->set(ia, 0, mx   );
      dipoles_      [i]->set(ia, 1, my   );
      dipoles_      [i]->set(ia, 2, mz   );
      quadrupoles_  [i]->set(ia, 0, qxx  );
      quadrupoles_  [i]->set(ia, 1, qxy  );
      quadrupoles_  [i]->set(ia, 2, qxz  );
      quadrupoles_  [i]->set(ia, 3, qyy  );
      quadrupoles_  [i]->set(ia, 4, qyz  );
      quadrupoles_  [i]->set(ia, 5, qzz  );
      octupoles_    [i]->set(ia, 0, oxxx );
      octupoles_    [i]->set(ia, 1, oxxy );
      octupoles_    [i]->set(ia, 2, oxxz );
      octupoles_    [i]->set(ia, 3, oxyy );
      octupoles_    [i]->set(ia, 4, oxyz );
      octupoles_    [i]->set(ia, 5, oxzz );
      octupoles_    [i]->set(ia, 6, oyyy );
      octupoles_    [i]->set(ia, 7, oyyz );
      octupoles_    [i]->set(ia, 8, oyzz );
      octupoles_    [i]->set(ia, 9, ozzz );
      hexadecapoles_[i]->set(ia, 0, hxxxx);
      hexadecapoles_[i]->set(ia, 1, hxxxy);
      hexadecapoles_[i]->set(ia, 2, hxxxz);
      hexadecapoles_[i]->set(ia, 3, hxxyy);
      hexadecapoles_[i]->set(ia, 4, hxxyz);
      hexadecapoles_[i]->set(ia, 5, hxxzz);
      hexadecapoles_[i]->set(ia, 6, hxyyy);
      hexadecapoles_[i]->set(ia, 7, hxyyz);
      hexadecapoles_[i]->set(ia, 8, hxyzz);
      hexadecapoles_[i]->set(ia, 9, hxzzz);
      hexadecapoles_[i]->set(ia,10, hyyyy);
      hexadecapoles_[i]->set(ia,11, hyyyz);
      hexadecapoles_[i]->set(ia,12, hyyzz);
      hexadecapoles_[i]->set(ia,13, hyzzz);
      hexadecapoles_[i]->set(ia,14, hzzzz);
 } // EndOFAtomLoop

 // Change origins from (0, 0, 0) to atomic positions
 this->recenter(centres_, i);
 origins_->copy(centres_);
}

} // EndNameSpace oepdev
