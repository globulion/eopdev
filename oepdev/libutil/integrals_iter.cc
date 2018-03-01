#include "integrals_iter.h"

using namespace oepdev;
using namespace std;

ShellCombinationsIterator::ShellCombinationsIterator(int nshell) :
 done(false), 
 nshell_(nshell), 
 bs_1_(psi::BasisSet::zero_ao_basis_set()),
 bs_2_(psi::BasisSet::zero_ao_basis_set()),
 bs_3_(psi::BasisSet::zero_ao_basis_set()),
 bs_4_(psi::BasisSet::zero_ao_basis_set())
{
}
ShellCombinationsIterator::~ShellCombinationsIterator() {}
void ShellCombinationsIterator::first() {}
void ShellCombinationsIterator::next() {}
void ShellCombinationsIterator::compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const {}
void ShellCombinationsIterator::compute_shell(std::shared_ptr<psi   ::TwoBodyAOInt> tei) const {}
int ShellCombinationsIterator::P() const { throw psi::PSIEXCEPTION("Wrong usage of ShellCombinationsIterator!");}
int ShellCombinationsIterator::Q() const { throw psi::PSIEXCEPTION("Wrong usage of ShellCombinationsIterator!");}
int ShellCombinationsIterator::R() const { throw psi::PSIEXCEPTION("Wrong usage of ShellCombinationsIterator!");}
int ShellCombinationsIterator::S() const { throw psi::PSIEXCEPTION("Wrong usage of ShellCombinationsIterator!");}
std::shared_ptr<AOIntegralsIterator> ShellCombinationsIterator::ao_iterator(std::string mode) const
{
  return AOIntegralsIterator::build(this, mode);
}
std::shared_ptr<ShellCombinationsIterator> ShellCombinationsIterator::build(const oepdev::IntegralFactory& ints,
                                                                            std::string mode, int nshell)
{
    std::shared_ptr<ShellCombinationsIterator> iterator;
    if (mode == "ALL") {
        if      (nshell == 4) iterator = std::make_shared<AllAOShellCombinationsIterator_4>(ints);
        else if (nshell == 2) iterator = std::make_shared<AllAOShellCombinationsIterator_2>(ints);
        else {
           throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator: Invalid number of shells specified");
        }
    }
    else if (mode == "UNIQUE") {
       throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator: Unique mode is not implemented yet!");
    }
    else {
       throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator: Invalid mode used.");
    }
    return iterator;
}
std::shared_ptr<ShellCombinationsIterator> ShellCombinationsIterator::build(const psi::IntegralFactory& ints,
                                                                            std::string mode, int nshell)
{
    std::shared_ptr<ShellCombinationsIterator> iterator;
    if (mode == "ALL") {
        if      (nshell == 4) iterator = std::make_shared<AllAOShellCombinationsIterator_4>(ints);
        else {
          throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator does not work for psi::IntegralFactory and shell doublet!");
        }
    }
    else if (mode == "UNIQUE") {
       throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator: Unique mode is not implemented yet!");
    }
    else {
       throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator: Invalid mode used.");
    }
    return iterator;
}
std::shared_ptr<ShellCombinationsIterator> ShellCombinationsIterator::build(std::shared_ptr<oepdev::IntegralFactory> ints,
                                                                            std::string mode, int nshell)
{
   return ShellCombinationsIterator::build(*ints, mode, nshell);
}
std::shared_ptr<ShellCombinationsIterator> ShellCombinationsIterator::build(std::shared_ptr<psi::IntegralFactory> ints,
                                                                            std::string mode, int nshell)
{
   return ShellCombinationsIterator::build(*ints, mode, nshell);
}
AllAOShellCombinationsIterator_4::AllAOShellCombinationsIterator_4(SharedBasisSet bs_1, SharedBasisSet bs_2,
                                                                   SharedBasisSet bs_3, SharedBasisSet bs_4) 
 : ShellCombinationsIterator(4),
   nshell_1(bs_1->nshell()-1), nshell_2(bs_2->nshell()-1), 
   nshell_3(bs_3->nshell()-1), nshell_4(bs_4->nshell()-1), 
   pp(0), qq(0), rr(0), ss(0)
{
  bs_1_ = std::make_shared<psi::BasisSet>(*bs_1);
  bs_2_ = std::make_shared<psi::BasisSet>(*bs_2);
  bs_3_ = std::make_shared<psi::BasisSet>(*bs_3);
  bs_4_ = std::make_shared<psi::BasisSet>(*bs_4);
}
AllAOShellCombinationsIterator_4::AllAOShellCombinationsIterator_4(std::shared_ptr<oepdev::IntegralFactory> ints) 
 : AllAOShellCombinationsIterator_4(ints->basis1(), ints->basis2(), 
                                    ints->basis3(), ints->basis4())
{

}
AllAOShellCombinationsIterator_4::AllAOShellCombinationsIterator_4(const oepdev::IntegralFactory& ints) 
 : AllAOShellCombinationsIterator_4(ints.basis1(), ints.basis2(), 
                                    ints.basis3(), ints.basis4())
{

}
AllAOShellCombinationsIterator_4::AllAOShellCombinationsIterator_4(std::shared_ptr<psi::IntegralFactory> ints) 
 : AllAOShellCombinationsIterator_4(ints->basis1(), ints->basis2(), 
                                    ints->basis3(), ints->basis4())
{

}
AllAOShellCombinationsIterator_4::AllAOShellCombinationsIterator_4(const psi::IntegralFactory& ints) 
 : AllAOShellCombinationsIterator_4(ints.basis1(), ints.basis2(), 
                                    ints.basis3(), ints.basis4())
{

}
void AllAOShellCombinationsIterator_4::first() {
  current.P = 0;
  current.Q = 0;
  current.R = 0;
  current.S = 0;
}
void AllAOShellCombinationsIterator_4::next() {
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
void AllAOShellCombinationsIterator_4::compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const 
{
   tei->compute_shell(current.P, current.Q, current.R, current.S);
}
void AllAOShellCombinationsIterator_4::compute_shell(std::shared_ptr<psi::TwoBodyAOInt> tei) const 
{
   tei->compute_shell(current.P, current.Q, current.R, current.S);
}
//-------------------------------------------------------------------------------------------------------
AllAOShellCombinationsIterator_2::AllAOShellCombinationsIterator_2(SharedBasisSet bs_1, SharedBasisSet bs_2)
 : ShellCombinationsIterator(2),
   nshell_1(bs_1->nshell()-1), nshell_2(bs_2->nshell()-1), 
   pp(0), qq(0)
{
  bs_1_ = std::make_shared<psi::BasisSet>(*bs_1);
  bs_2_ = std::make_shared<psi::BasisSet>(*bs_2);
}
AllAOShellCombinationsIterator_2::AllAOShellCombinationsIterator_2(std::shared_ptr<oepdev::IntegralFactory> ints) 
 : AllAOShellCombinationsIterator_2(ints->basis1(), ints->basis2())
{

}
AllAOShellCombinationsIterator_2::AllAOShellCombinationsIterator_2(const oepdev::IntegralFactory& ints) 
 : AllAOShellCombinationsIterator_2(ints.basis1(), ints.basis2())
{

}
AllAOShellCombinationsIterator_2::AllAOShellCombinationsIterator_2(std::shared_ptr<psi::IntegralFactory> ints) 
 : AllAOShellCombinationsIterator_2(ints->basis1(), ints->basis2())
{
  throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator does not work for psi::IntegralFactory and shell doublets!");
}
AllAOShellCombinationsIterator_2::AllAOShellCombinationsIterator_2(const psi::IntegralFactory& ints) 
 : AllAOShellCombinationsIterator_2(ints.basis1(), ints.basis2())
{
  throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsIterator does not work for psi::IntegralFactory and shell doublets!");
}

void AllAOShellCombinationsIterator_2::first() {
  current.P = 0;
  current.Q = 0;
}
void AllAOShellCombinationsIterator_2::next() {
  ++qq;
  if(qq > nshell_2) {
     qq = 0;
     ++pp;
     if(pp > nshell_1) {
        done = true;
     }
  }
  //
  current.P = pp;
  current.Q = qq;
}
void AllAOShellCombinationsIterator_2::compute_shell(std::shared_ptr<oepdev::TwoBodyAOInt> tei) const 
{
   tei->compute_shell(current.P, current.Q);
}
void AllAOShellCombinationsIterator_2::compute_shell(std::shared_ptr<psi::TwoBodyAOInt> tei) const 
{
   throw psi::PSIEXCEPTION(" oepdev::ShellCombinationsFactory_2 does not operate on psi::TwoBodyAOInt!");
}
//-------------------------------------------------------------------------------------------------------
AOIntegralsIterator::AOIntegralsIterator() :
 done(false)
{
}
AOIntegralsIterator::~AOIntegralsIterator() 
{
}
std::shared_ptr<AOIntegralsIterator> AOIntegralsIterator::build(const ShellCombinationsIterator* shellIter, std::string mode)
{
   std::shared_ptr<AOIntegralsIterator> iterator;
   if (mode == "ALL") {
     if      (shellIter->nshell() == 4) iterator = std::make_shared<AllAOIntegralsIterator_4>(shellIter);
     else if (shellIter->nshell() == 2) iterator = std::make_shared<AllAOIntegralsIterator_2>(shellIter);
   } else if (mode == "UNIQUE") {
          throw psi::PSIEXCEPTION(" oepdev::AOIntegralsIterator: Unique iteration is not implemented yet!");
   } else throw psi::PSIEXCEPTION(" oepdev::AOIntegealsIterator: Invalid mode chosen!");
   return iterator;
}
std::shared_ptr<AOIntegralsIterator> AOIntegralsIterator::build(std::shared_ptr<ShellCombinationsIterator> shellIter,
                                                                std::string mode)
{
   return AOIntegralsIterator::build(shellIter.get(), mode);
}
void AOIntegralsIterator::first() {}
void AOIntegralsIterator::next() {}
int AOIntegralsIterator::i() const {}
int AOIntegralsIterator::j() const {}
int AOIntegralsIterator::k() const {}
int AOIntegralsIterator::l() const {}
int AOIntegralsIterator::index() const {}
//-------------------------------------------------------------------------------------------------------
AllAOIntegralsIterator_4::AllAOIntegralsIterator_4(std::shared_ptr<ShellCombinationsIterator> shellIter) 
 : AllAOIntegralsIterator_4(shellIter.get())
{
}
AllAOIntegralsIterator_4::AllAOIntegralsIterator_4(const ShellCombinationsIterator* shellIter) 
 : AOIntegralsIterator(),
                nishell_1(shellIter->bs_1()->shell(shellIter->P()).nfunction()-1), 
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
void AllAOIntegralsIterator_4::first() {
  current.i = ishell_1;
  current.j = ishell_2;
  current.k = ishell_3;
  current.l = ishell_4;
}
void AllAOIntegralsIterator_4::next() {
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
   current.i = ii__ + ishell_1;
   current.j = jj__ + ishell_2;
   current.k = kk__ + ishell_3;
   current.l = ll__ + ishell_4;
   current.index = index__;
}
//-------------------------------------------------------------------------------------------------------
AllAOIntegralsIterator_2::AllAOIntegralsIterator_2(std::shared_ptr<ShellCombinationsIterator> shellIter) 
 : AllAOIntegralsIterator_2(shellIter.get())
{
}
AllAOIntegralsIterator_2::AllAOIntegralsIterator_2(const ShellCombinationsIterator* shellIter) 
 : AOIntegralsIterator(),
   nishell_1(shellIter->bs_1()->shell(shellIter->P()).nfunction()-1), 
   nishell_2(shellIter->bs_2()->shell(shellIter->Q()).nfunction()-1), 
    ishell_1(shellIter->bs_1()->shell(shellIter->P()).function_index()), 
    ishell_2(shellIter->bs_2()->shell(shellIter->Q()).function_index()), 
   ii__(0), jj__(0), index__(0)
{
 current.i = 0;
 current.j = 0;
 current.index = 0;
}
void AllAOIntegralsIterator_2::first() {
  current.i = ishell_1;
  current.j = ishell_2;
}
void AllAOIntegralsIterator_2::next() {
   ++jj__; ++index__;
   if(jj__ > nishell_2) {
      jj__ = 0;
      ++ii__;
      if(ii__ > nishell_1) {
         done = true;
      }
   }
   current.i = ii__ + ishell_1;
   current.j = jj__ + ishell_2;
   current.index = index__;
}
