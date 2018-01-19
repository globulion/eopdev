#include "solver.h"

namespace oepdev{

using namespace psi;
using namespace std;

OEPDevSolver::OEPDevSolver(SharedWavefunctionUnion wfn_union)
 : wfn_union_(wfn_union), methods_oepBased_(), methods_benchmark_()
{
  
}

OEPDevSolver::~OEPDevSolver()
{

}

double OEPDevSolver::compute_oep_based(const std::string& method) {}
double OEPDevSolver::compute_benchmark(const std::string& method) {}


ElectrostaticEnergySolver::ElectrostaticEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  methods_oepBased_.push_back("DEFAULT");
  methods_benchmark_.push_back("DEFAULT");
}

ElectrostaticEnergySolver::~ElectrostaticEnergySolver() 
{

}

double ElectrostaticEnergySolver::compute_oep_based(const std::string& method) {}
double ElectrostaticEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  double e_nuc_nuc = 0.0;
  double e_nuc_el  = 0.0; 
  double e_el_el   = 0.0;

  SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
  SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);

  // ===> Nuc(A) --- Nuc(B) <=== //
  e_nuc_nuc = wfn_union_->nuclear_repulsion_interaction_energy();

  // ===> Nuc(A) --- Nuc(B) <=== //
  // ===> Nuc(A) --- Nuc(B) <=== //
  // ===> Nuc(A) --- Nuc(B) <=== //

  e = e_nuc_nuc + e_nuc_el + e_el_el;
  return e;
}



std::shared_ptr<OEPDevSolver> OEPDevSolver::build(const std::string& target, SharedWavefunctionUnion wfn_union)
{
   std::shared_ptr<OEPDevSolver> solver;
   if     (target == "ELECTROSTATIC ENERGY") solver = std::make_shared<ElectrostaticEnergySolver>(wfn_union);
   else throw psi::PSIEXCEPTION("OEPDEV: Error. OEPDevSolver: Incorrect targer property chosen!\n");
   return solver;
}

} // EndNameSpace oepdev
