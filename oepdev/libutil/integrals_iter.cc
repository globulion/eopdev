#include "integrals_iter.h"

namespace oepdev{

using namespace psi;
using namespace std;


AllAOShellCombinationsIterator::AllAOShellCombinationsIterator(SharedBasisSet bs_1, SharedBasisSet bs_2,
                                                               SharedBasisSet bs_3, SharedBasisSet bs_4) 
 : done(false), nshell_1(bs_1->nshell()-1), nshell_2(bs_2->nshell()-1), 
                nshell_3(bs_3->nshell()-1), nshell_4(bs_4->nshell()-1), 
                bs_1_(bs_1), bs_2_(bs_2), bs_3_(bs_3), bs_4_(bs_4),
                pp(0), qq(0), rr(0), ss(0)
{
 
}

AllAOShellCombinationsIterator::AllAOShellCombinationsIterator(SharedIntegralFactory ints) 
 : AllAOShellCombinationsIterator(ints->basis1(), ints->basis2(), 
                                  ints->basis3(), ints->basis4())
{

}

AllAOShellCombinationsIterator::AllAOShellCombinationsIterator(psi::IntegralFactory ints) 
 : AllAOShellCombinationsIterator(ints.basis1(), ints.basis2(), 
                                  ints.basis3(), ints.basis4())
{

}


AllAOIntegralsIterator::AllAOIntegralsIterator(std::shared_ptr<AllAOShellCombinationsIterator> shellIter) 
 : done(false), nishell_1(shellIter->bs_1()->shell(shellIter->P()).nfunction()-1), 
                nishell_2(shellIter->bs_2()->shell(shellIter->Q()).nfunction()-1), 
                nishell_3(shellIter->bs_3()->shell(shellIter->R()).nfunction()-1), 
                nishell_4(shellIter->bs_4()->shell(shellIter->S()).nfunction()-1), 
                 ishell_1(shellIter->bs_1()->shell(shellIter->P()).function_index()), 
                 ishell_2(shellIter->bs_2()->shell(shellIter->Q()).function_index()), 
                 ishell_3(shellIter->bs_3()->shell(shellIter->R()).function_index()), 
                 ishell_4(shellIter->bs_4()->shell(shellIter->S()).function_index()), 
                ii__(0), jj__(0), kk__(0), ll__(0), index__(0)
{
 current.i = 0;
 current.j = 0;
 current.k = 0;
 current.l = 0;
 current.index = 0;
}

AllAOIntegralsIterator::AllAOIntegralsIterator(const AllAOShellCombinationsIterator& shellIter) 
 : done(false), nishell_1(shellIter.bs_1()->shell(shellIter.P()).nfunction()-1), 
                nishell_2(shellIter.bs_2()->shell(shellIter.Q()).nfunction()-1), 
                nishell_3(shellIter.bs_3()->shell(shellIter.R()).nfunction()-1), 
                nishell_4(shellIter.bs_4()->shell(shellIter.S()).nfunction()-1), 
                 ishell_1(shellIter.bs_1()->shell(shellIter.P()).function_index()), 
                 ishell_2(shellIter.bs_2()->shell(shellIter.Q()).function_index()), 
                 ishell_3(shellIter.bs_3()->shell(shellIter.R()).function_index()), 
                 ishell_4(shellIter.bs_4()->shell(shellIter.S()).function_index()), 
                ii__(0), jj__(0), kk__(0), ll__(0), index__(0)
{
 current.i = 0;
 current.j = 0;
 current.k = 0;
 current.l = 0;
 current.index = 0;
}


void AllAOShellCombinationsIterator::first() {
  current.P = 0;
  current.Q = 0;
  current.R = 0;
  current.S = 0;
}

void AllAOIntegralsIterator::first() {
  current.i = ishell_1;
  current.j = ishell_2;
  current.k = ishell_3;
  current.l = ishell_4;
}


void AllAOShellCombinationsIterator::next() {
  ++ss;
  if(ss > nshell_4) {
     ss = 0;
     ++rr;
     if(rr > nshell_3) {
        rr = 0;
        ++qq;
        if(qq > nshell_2) {
           qq = 0;
           ++pp;
           if(pp > nshell_1) {
              done = true;
           }
        }
     }
  }
  //
  current.P = pp;
  current.Q = qq;
  current.R = rr;
  current.S = ss;
}

void AllAOIntegralsIterator::next() {
   ++ll__; ++index__;
   if(ll__ > nishell_4) {
      ll__ = 0;
      ++kk__;
      if(kk__ > nishell_3) {
         kk__ = 0;
         ++jj__;
         if(jj__ > nishell_2) {
            jj__ = 0;
            ++ii__;
            if(ii__ > nishell_1) {
               done = true;
            }
         }
      }
   }
   //
   //std::cout << ii__ << jj__ << kk__ << ll__ << " " << ishell_1 << ishell_2 << ishell_3 << ishell_4 << "\n";
   //std::cout << " Debug: " << nishell_1 << nishell_2 << nishell_3 << nishell_4 << " "  << ishell_1 << ishell_2 << ishell_3 << ishell_4 << "\n";
   current.i = ii__ + ishell_1;
   current.j = jj__ + ishell_2;
   current.k = kk__ + ishell_3;
   current.l = ll__ + ishell_4;
   // std::cout << current.i << current.j << current.k << current.l << "\n";
   current.index = index__;
}

void AllAOShellCombinationsIterator::compute_shell(SharedTwoBodyAOInt tei) const 
{
   tei->compute_shell(current.P, current.Q, current.R, current.S);
}


} // EndNameSpace oepdev
