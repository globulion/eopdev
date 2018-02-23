#ifndef _oepdev_libutil_integrals_iter_h
#define _oepdev_libutil_integrals_iter_h
/** @file integrals_iter.h */

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
/** \addtogroup OEPDEV_INTEGRAL_HELPERS
 * @{
 */

/** \brief Loop over all possible ERI shells.
 *
 * Constructed by providing shared pointer to IntegralFactory object or
 * shared pointers to four basis set spaces.
 *
 * Suggested usage:
 *
 * \code{.cpp}
 *  SharedIntegralFactory ints = std::make_shared<IntegralFactory>(bs1, bs2, bs3, bs4);
 *  SharedTwoBodyAOInt tei(ints->eri());
 *  AllAOShellCombinationsIterator shellIter(ints);
 *  const double * buffer = tei->buffer();
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
 public:
 
   /**\brief Construct by providing basis sets for each axis.
    *  The basis sets must be defined for the same molecule.
    *  @param bs_1 - basis set of axis 1
    *  @param bs_2 - basis set of axis 2
    *  @param bs_3 - basis set of axis 3
    *  @param bs_4 - basis set of axis 4
    */  
   AllAOShellCombinationsIterator(SharedBasisSet bs_1, SharedBasisSet bs_2, 
                                  SharedBasisSet bs_3, SharedBasisSet bs_4);

   /**\brief Construct by providing integral factory.
    *  
    *  @param integrals - integral factory object
    */  
   AllAOShellCombinationsIterator(SharedIntegralFactory integrals);
   AllAOShellCombinationsIterator(psi::IntegralFactory integrals);
   
   /// First iteration 
   void first();
   /// Next iteration
   void next();
   /// Check status of iterations
   bool is_done() {return done;}

   /// Grab the current shell *P* index
   int P() const { return current.P; }
   /// Grab the current shell *Q* index
   int Q() const { return current.Q; }
   /// Grab the current shell *R* index
   int R() const { return current.R; }
   /// Grab the current shell *S* index
   int S() const { return current.S; }

   /// Grab the basis set of axis 1
   SharedBasisSet bs_1() const { return bs_1_;}
   /// Grab the basis set of axis 2
   SharedBasisSet bs_2() const { return bs_2_;}
   /// Grab the basis set of axis 3
   SharedBasisSet bs_3() const { return bs_3_;}
   /// Grab the basis set of axis 4
   SharedBasisSet bs_4() const { return bs_4_;}

   /**\brief Compute ERI's for the current shell.
    * The eris are stored in the buffer of the argument object.
    * @param tei - two electron AO integral
    */
   void compute_shell(SharedTwoBodyAOInt tei) const;

 private:
   /// Integral Shell: stores indices of a shell (PQ|RS)
   struct Shell {
      int P, Q, R, S;
   };
   /// Current shell
   Shell current;
   /// Status of iterator
   bool done;
   const int nshell_1, nshell_2, nshell_3, nshell_4;
   int pp, qq, rr, ss;

   /// Basis set of axis 1
   SharedBasisSet bs_1_;
   /// Basis set of axis 2
   SharedBasisSet bs_2_;
   /// Basis set of axis 3
   SharedBasisSet bs_3_;
   /// Basis set of axis 4
   SharedBasisSet bs_4_;

};

/** \brief Loop over all possible ERI within a particular shell.
 *
 * Constructed by providing a const reference or shared pointer 
 * to an AllAOShellCombinationsIterator object.
 *
 * Suggested usage:
 *
 * \see AllAOShellCombinationsIterator
 */
class AllAOIntegralsIterator {
 public:

   /**\brief Construct by shell iterator (const object)
    *  @param shellIter - shell iterator object
    */  
   AllAOIntegralsIterator(const AllAOShellCombinationsIterator& shellIter);

   /**\brief Construct by shell iterator (pointed by shared pointer)
    *  @param shellIter - shell iterator object
    */  
   AllAOIntegralsIterator(std::shared_ptr<AllAOShellCombinationsIterator> shellIter);
  
   /// First iteration 
   void first();
   /// Next iteration
   void next();
   /// Check status of iterations
   bool is_done() {return done;}

   /// Grab the current integral *i* index
   int i() const { return current.i; }
   /// Grab the current integral *j* index
   int j() const { return current.j; }
   /// Grab the current integral *k* index
   int k() const { return current.k; }
   /// Grab the current integral *l* index
   int l() const { return current.l; }

   /** Grab the current index of integral value stored in the buffer */
   int index() const { return current.index;}

 private:
   /// Integral (stores current indices and buffer index)
   struct Integral {
      int i;
      int j;
      int k;
      int l;
      unsigned int index;
   };
   /// Current integral
   Integral current;

   /// Status of iteration procedure
   bool done;

   const int nishell_1, nishell_2, nishell_3, nishell_4;
   const int  ishell_1,  ishell_2,  ishell_3,  ishell_4;

   int ii__, jj__, kk__, ll__;
   int index__;
};

/** \example example_integrals_iter.cc
 * Iterations over electron repulsion integrals in AO basis.
 * This is an example of how to use
 *   - the AllAOShellCombinationsIterator class
 *   - the AllAOIntegralsIterator class.
 */

/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_integrals_iter_h
