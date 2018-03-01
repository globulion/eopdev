#ifndef _oepdev_libutil_integrals_iter_h
#define _oepdev_libutil_integrals_iter_h
/** @file integrals_iter.h */

#include<cstdio>

#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

#include "../libpsi/integral.h"

namespace oepdev{

using namespace std;

using SharedBasisSet        = std::shared_ptr<psi::BasisSet>;
using SharedIntegralFactory = std::shared_ptr<IntegralFactory>;
using SharedTwoBodyAOInt    = std::shared_ptr<TwoBodyAOInt>;

class AOIntegralsIterator;

/** \addtogroup OEPDEV_INTEGRAL_HELPERS
 * @{
 */

/** \brief Iterator for Shell Combinations. Abstract Base.
 *
 *  \date 2018/03/01 17:22:00
 */
class ShellCombinationsIterator
{
  public:
   /** \brief Constructor
    *  @param nshell - number of shells this iterator is for
    */
   ShellCombinationsIterator(int nshell);

   /// Destructor
   virtual ~ShellCombinationsIterator();

   /** \brief Build shell iterator from oepdev::IntegralFactory.
    *  @param ints    - integral factory
    *  @param mode    - mode of iteration (either `ALL` or `UNIQUE`)
    *  @param nshell  - number of shells to iterate through
    *  @return shell iterator
    */
   static std::shared_ptr<ShellCombinationsIterator> build(const IntegralFactory& ints, 
                                                           std::string mode = "ALL", int nshell = 4);
   /*!
    * \overload
    */
   static std::shared_ptr<ShellCombinationsIterator> build(std::shared_ptr<IntegralFactory> ints, 
                                                           std::string mode = "ALL", int nshell = 4);
   /** \brief Build shell iterator from psi::IntegralFactory.
    *  @param ints    - integral factory
    *  @param mode    - mode of iteration (either `ALL` or `UNIQUE`)
    *  @param nshell  - number of shells to iterate through
    *  @return shell iterator
    */
   static std::shared_ptr<ShellCombinationsIterator> build(const psi::IntegralFactory& ints, 
                                                           std::string mode = "ALL", int nshell = 4);
   /*!
    * \overload
    */
   static std::shared_ptr<ShellCombinationsIterator> build(std::shared_ptr<psi::IntegralFactory> ints, 
                                                           std::string mode = "ALL", int nshell = 4);

   /// First iteration
   virtual void first(void) = 0;

   /// Next iteration
   virtual void next(void) = 0;

   //!@{ 
   /** Compute integrals in a current shell.
    *  Works both for oepdev::TwoBodyAOInt and psi::TwoBodyAOInt
    *  @param tei - two body integral object
    */
   virtual void compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const = 0;
   virtual void compute_shell(std::shared_ptr<psi   ::TwoBodyAOInt> tei) const = 0;
   //!@}

   /// Grab the basis set of axis 1
   virtual std::shared_ptr<psi::BasisSet> bs_1(void) const {return bs_1_;}
   /// Grab the basis set of axis 2
   virtual std::shared_ptr<psi::BasisSet> bs_2(void) const {return bs_2_;}
   /// Grab the basis set of axis 3
   virtual std::shared_ptr<psi::BasisSet> bs_3(void) const {return bs_3_;}
   /// Grab the basis set of axis 4
   virtual std::shared_ptr<psi::BasisSet> bs_4(void) const {return bs_4_;}

   /// Grab the current shell *P* index
   virtual int P(void) const;
   /// Grab the current shell *Q* index
   virtual int Q(void) const;
   /// Grab the current shell *R* index
   virtual int R(void) const;
   /// Grab the current shell *S* index
   virtual int S(void) const;

   /// Return status of an iterator
   virtual bool is_done(void) {return done;}

   /// Return number of shells this iterator is for
   virtual const int nshell(void) const {return nshell_;}

   /** Make an AO integral iterator based on current shell
    * @param mode - either "ALL" or "UNIQUE" (iterate over all or unique integrals)
    * @return iterator over AO integrals */
   virtual std::shared_ptr<AOIntegralsIterator> ao_iterator(std::string mode = "ALL") const;

  protected:
   /// Basis set of axis 1
   SharedBasisSet bs_1_;
   /// Basis set of axis 2
   SharedBasisSet bs_2_;
   /// Basis set of axis 3
   SharedBasisSet bs_3_;
   /// Basis set of axis 4
   SharedBasisSet bs_4_;

   /// Number of shells this iterator is for
   const int nshell_;
   /// Status of an iterator
   bool done;
};

/** \brief Iterator for AO Integrals. Abstract Base.
 *
 */
class AOIntegralsIterator
{
  public:
   /// Base Constructor
   AOIntegralsIterator();
   /// Base Destructor
   virtual ~AOIntegralsIterator();
   /** Build AO integrals iterator from current state of iterator over shells
    *  @param shellIter - iterator over shells
    *  @mode - either "ALL" or "UNIQUE" (iterate over all or unique integrals)
    *  @return iterator over AO integrals
    */
   static std::shared_ptr<AOIntegralsIterator> build(const ShellCombinationsIterator* shellIter, std::string mode = "ALL");
   /**
    \overload
    */
   static std::shared_ptr<AOIntegralsIterator> build(std::shared_ptr<ShellCombinationsIterator> shellIter,  
                                                     std::string mode = "ALL");

   /// Do the first iteration
   virtual void first(void) = 0;
   /// Do the next iteration
   virtual void next(void) = 0;

   /// Grab *i*-th index
   virtual int i(void) const;
   /// Grab *j*-th index
   virtual int j(void) const;
   /// Grab *k*-th index
   virtual int k(void) const;
   /// Grab *l*-th index
   virtual int l(void) const;
   /// Grab index in the integral buffer
   virtual int index(void) const = 0;

   /// Returns the status of an iterator
   virtual bool is_done(void) {return done;}

  protected:
   /// The status of an iterator
   bool done;
};


/** \brief Loop over all possible ERI shells in a shell quartet.
 *
 * Constructed by providing IntegralFactory object or
 * shared pointers to four basis set spaces.
 */
class AllAOShellCombinationsIterator_4 : public ShellCombinationsIterator
{
 public:
 
   /**\brief Iterate over shell quartets. 
    *  Construct by providing basis sets for each axis.
    *  The basis sets must be defined for the same molecule.
    *  @param bs_1 - basis set of axis 1
    *  @param bs_2 - basis set of axis 2
    *  @param bs_3 - basis set of axis 3
    *  @param bs_4 - basis set of axis 4
    */  
   AllAOShellCombinationsIterator_4(SharedBasisSet bs_1, SharedBasisSet bs_2, 
                                    SharedBasisSet bs_3, SharedBasisSet bs_4);

   /**\brief Construct by providing integral factory.
    *  
    *  @param integrals - OepDev integral factory object
    */  
   AllAOShellCombinationsIterator_4(std::shared_ptr<IntegralFactory> integrals);
   /**
    * \overload
    */
   AllAOShellCombinationsIterator_4(const IntegralFactory& integrals);
   /**\brief Construct by providing integral factory.
    *  
    *  @param integrals - OepDev integral factory object
    */  
   AllAOShellCombinationsIterator_4(std::shared_ptr<psi::IntegralFactory> integrals);
   /**
    * \overload
    */
   AllAOShellCombinationsIterator_4(const psi::IntegralFactory& integrals);

   /// Do the first iteration   
   void first();
   /// Do the next iteration   
   void next();
   void compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const;
   void compute_shell(std::shared_ptr<psi   ::TwoBodyAOInt> tei) const;

   /// Grab the current shell *P* index
   int P() const { return current.P; }
   /// Grab the current shell *Q* index
   int Q() const { return current.Q; }
   /// Grab the current shell *R* index
   int R() const { return current.R; }
   /// Grab the current shell *S* index
   int S() const { return current.S; }

 private:
   // Integral Shell: stores indices of a shell (PQ|RS)
   struct Shell {
      int P, Q, R, S;
   };
   // Current shell
   Shell current;

   const int nshell_1, nshell_2, nshell_3, nshell_4;
   int pp, qq, rr, ss;

};

/** \brief Loop over all possible ERI shells in a shell doublet.
 *
 * Constructed by providing IntegralFactory object or
 * shared pointers to two basis set spaces.
 *
 */
class AllAOShellCombinationsIterator_2 : public ShellCombinationsIterator
{
 public:
 
   /**\brief Iterate over shell doublets. 
    *  Construct by providing basis sets for each axis.
    *  The basis sets must be defined for the same molecule.
    *  @param bs_1 - basis set of axis 1
    *  @param bs_2 - basis set of axis 2
    */  
   AllAOShellCombinationsIterator_2(SharedBasisSet bs_1, SharedBasisSet bs_2);

   /**\brief Construct by providing integral factory.
    *  
    *  @param integrals - integral factory object 
    *  @param ib1       - BasisSet axis in integral factory for *P* shell
    *  @param ib2       - BasisSet axis in integral factory for *Q* shell
    */  
   AllAOShellCombinationsIterator_2(std::shared_ptr<IntegralFactory> integrals);
   AllAOShellCombinationsIterator_2(const IntegralFactory& integrals);
   AllAOShellCombinationsIterator_2(std::shared_ptr<psi::IntegralFactory> integrals);
   AllAOShellCombinationsIterator_2(const psi::IntegralFactory& integrals);
  
   /// First iteration 
   void first();
   /// Next iteration
   void next();
   /**\brief Compute ERI's for the current shell.
    * The eris are stored in the buffer of the argument object.
    * @param tei - two electron AO integral
    */
   void compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const;
   void compute_shell(std::shared_ptr<   psi::TwoBodyAOInt> tei) const;

   /// Grab the current shell *P* index
   int P() const { return current.P; }
   /// Grab the current shell *Q* index
   int Q() const { return current.Q; }

 private:
   // Integral Shell: stores indices of a shell (P|Q)
   struct Shell {
      int P, Q;
   };
   // Current shell
   Shell current;

   const int nshell_1, nshell_2;
   int pp, qq;

};


/** \brief Loop over all possible ERI within a particular shell quartet.
 *
 * Constructed by providing a const reference or shared pointer 
 * to an AllAOShellCombinationsIterator object.
 *
 * \see AllAOShellCombinationsIterator_4
 */
class AllAOIntegralsIterator_4 : public AOIntegralsIterator
{
 public:

   /**\brief Construct by shell iterator (const object)
    *  @param shellIter - shell iterator object
    */  
   //AllAOIntegralsIterator(const AllAOShellCombinationsIterator& shellIter);
   AllAOIntegralsIterator_4(const ShellCombinationsIterator* shellIter);

   /**\brief Construct by shell iterator (pointed by shared pointer)
    *  @param shellIter - shell iterator object
    */  
   //AllAOIntegralsIterator(std::shared_ptr<AllAOShellCombinationsIterator> shellIter);
   AllAOIntegralsIterator_4(std::shared_ptr<ShellCombinationsIterator> shellIter);
  
   /// First iteration 
   void first();
   /// Next iteration
   void next();

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
   // Integral (stores current indices and buffer index)
   struct Integral {
      int i;
      int j;
      int k;
      int l;
      unsigned int index;
   };
   // Current integral
   Integral current;

   const int nishell_1, nishell_2, nishell_3, nishell_4;
   const int  ishell_1,  ishell_2,  ishell_3,  ishell_4;

   int ii__, jj__, kk__, ll__;
   int index__;
};

/** \brief Loop over all possible ERI within a particular shell doublet.
 *
 * Constructed by providing a const reference or shared pointer 
 * to an AllAOShellCombinationsIterator object.
 *
 * \see AllAOShellCombinationsIterator_2
 */
class AllAOIntegralsIterator_2 : public AOIntegralsIterator
{
 public:

   /**\brief Construct by shell iterator (const object)
    *  @param shellIter - shell iterator object
    */  
   AllAOIntegralsIterator_2(const ShellCombinationsIterator* shellIter);

   /**\brief Construct by shell iterator (pointed by shared pointer)
    *  @param shellIter - shell iterator object
    */  
   AllAOIntegralsIterator_2(std::shared_ptr<ShellCombinationsIterator> shellIter);
  
   /// First iteration 
   void first();
   /// Next iteration
   void next();

   /// Grab the current integral *i* index
   int i() const { return current.i; }
   /// Grab the current integral *j* index
   int j() const { return current.j; }

   /** Grab the current index of integral value stored in the buffer */
   int index() const { return current.index;}

 private:
   // Integral (stores current indices and buffer index)
   struct Integral {
      int i;
      int j;
      unsigned int index;
   };
   // Current integral
   Integral current;

   const int nishell_1, nishell_2;
   const int  ishell_1,  ishell_2;

   int ii__, jj__;
   int index__;
};

// ---> Quick use typedefs <--- //

/// Iterator over shells as shared pointer
using SharedShellsIterator = std::shared_ptr<ShellCombinationsIterator>;

/// Iterator over AO integrals as shared pointer
using SharedAOIntsIterator = std::shared_ptr<AOIntegralsIterator>;

/** \example example_integrals_iter.cc
 * Iterations over electron repulsion integrals in AO basis.
 * This is an example of how to use
 *   - the oepdev::ShellCombinationsIterator class
 *   - the oepdev::AOIntegralsIterator class.
 */

/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_integrals_iter_h
