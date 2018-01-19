#ifndef _oepdev_libutil_integrals_iter_h
#define _oepdev_libutil_integrals_iter_h

#include<cstdio>

#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

namespace oepdev{

using namespace psi;
using namespace std;

using SharedBasisSet        = std::shared_ptr<BasisSet>;
using SharedIntegralFactory = std::shared_ptr<IntegralFactory>;
using SharedTwoBodyAOInt    = std::shared_ptr<TwoBodyAOInt>;

/** \brief Loop over all possible ERI shells.
 *
 * Constructed by providing shared pointer to IntegralFactory object or
 * shared pointers to four basis set spaces.
 *
 * Suggested usage:
 *
 * \code{.cpp}
 *  SharedIntegralFactory ints(bs1, bs2, bs3, bs4);
 *  SharedTwoBodyAOInt tei = std::make_shared<TwoBodyAOInt>(ints->eri());
 *  AllAOShellCombinationsIterator shellIter(ints);
 *  const doube * buffer = tei->buffer();
 *  for (shellIter.first(); shellIter.is_done()==false; shellIter.next())
 *  {
 *       shellIter.compute_shell(tei);
 *       AllAOIntegralsIterator intsIter(shellIter);
 *       for (intsIter.first(); intsIter.is_done()==false; intsIter.next())
 *       {
 *            // Grab (ij|kl) integrals and indices here
 *            int i = intsIter.i();
 *            int j = intsIter.j();
 *            int k = intsIter.k();
 *            int l = intsIter.l();
 *            double integral = buffer[intsIter.index()];
 *       }
 *  }
 * \endcode
 */
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

   //@{
   /** Grab the current shell indices */
   int P() const { return current.P; }
   int Q() const { return current.Q; }
   int R() const { return current.R; }
   int S() const { return current.S; }
   //@}
 
   SharedBasisSet bs_1() const { return bs_1_;}
   SharedBasisSet bs_2() const { return bs_2_;}
   SharedBasisSet bs_3() const { return bs_3_;}
   SharedBasisSet bs_4() const { return bs_4_;}

   void compute_shell(SharedTwoBodyAOInt tei) const;

};

/** \brief Loop over all possible ERI within a particular shell.
 *
 * Constructed by providing a const reference or shared pointer 
 * to an AllAOShellCombinationsIterator object.
 *
 * Suggested usage:
 *
 * \code{.cpp}
 *  SharedIntegralFactory ints(bs1, bs2, bs3, bs4);
 *  SharedTwoBodyAOInt tei = std::make_shared<TwoBodyAOInt>(ints->eri());
 *  AllAOShellCombinationsIterator shellIter(ints);
 *  const doube * buffer = tei->buffer();
 *  for (shellIter.first(); shellIter.is_done()==false; shellIter.next())
 *  {
 *       shellIter.compute_shell(tei);
 *       AllAOIntegralsIterator intsIter(shellIter);
 *       for (intsIter.first(); intsIter.is_done()==false; intsIter.next())
 *       {
 *            // Grab (ij|kl) integrals and indices here
 *            int i = intsIter.i();
 *            int j = intsIter.j();
 *            int k = intsIter.k();
 *            int l = intsIter.l();
 *            double integral = buffer[intsIter.index()];
 *       }
 *  }
 * \endcode
 */
class AllAOIntegralsIterator {
 private:
   struct Integral {
      int i;
      int j;
      int k;
      int l;
      unsigned int index;
   };
   Integral current;

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

   //@{
   /** Grab the current integral indices */
   int i() const { return current.i; }
   int j() const { return current.j; }
   int k() const { return current.k; }
   int l() const { return current.l; }
   //@}

   /** Grab the current index of integral value stored in the buffer */
   int index() const { return current.index;}
};


}      // EndNameSpace oepdev
#endif //_oepdev_libutil_integrals_iter_h
