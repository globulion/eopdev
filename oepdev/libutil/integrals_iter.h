#ifndef _oepdev_libutil_integrals_iter_h
#define _oepdev_libutil_integrals_iter_h

#include<cstdio>

#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

namespace oepdev_libutil{

using namespace psi;
using namespace std;

using SharedBasisSet        = std::shared_ptr<BasisSet>;
using SharedIntegralFactory = std::shared_ptr<IntegralFactory>;

class AllAOShellCombinationsIterator {
 private:
   struct Shell {
      int P, Q, R, S;
   };
   Shell current;
   bool done;
   const int nshell_1, nshell_2, nshell_3, nshell_4;
   int pp, qq, rr, ss;

   SharedBasisSet bs_1_;
   SharedBasisSet bs_2_;
   SharedBasisSet bs_3_;
   SharedBasisSet bs_4_;

 public:
   AllAOShellCombinationsIterator(SharedBasisSet bs_1, SharedBasisSet bs_2, 
                                  SharedBasisSet bs_3, SharedBasisSet bs_4);
   AllAOShellCombinationsIterator(SharedIntegralFactory integrals);
    
   void first();
   void next();
   bool is_done() {return done;}

   int P() const { return current.P; }
   int Q() const { return current.Q; }
   int R() const { return current.R; }
   int S() const { return current.S; }
  
   SharedBasisSet bs_1() const { return bs_1_;}
   SharedBasisSet bs_2() const { return bs_2_;}
   SharedBasisSet bs_3() const { return bs_3_;}
   SharedBasisSet bs_4() const { return bs_4_;}

};

class AllAOIntegralsIterator {
 private:
   struct Bazi {
      int i;
      int j;
      int k;
      int l;
      unsigned int index;
   };
   Bazi current;

   bool done;

   const int nishell_1, nishell_2, nishell_3, nishell_4;
   const int  ishell_1,  ishell_2,  ishell_3,  ishell_4;

   int ii__, jj__, kk__, ll__;
   int index__;
   

 public:
   AllAOIntegralsIterator(const AllAOShellCombinationsIterator& shellIter);
   AllAOIntegralsIterator(std::shared_ptr<AllAOShellCombinationsIterator> shellIter);
   
   void first();
   void next();
   bool is_done() {return done;}

   int i() const { return current.i; }
   int j() const { return current.j; }
   int k() const { return current.k; }
   int l() const { return current.l; }

   int index() const { return current.index;}
};


}      // EndNameSpace oepdev_libutil
#endif //_oepdev_libutil_integrals_iter_h
