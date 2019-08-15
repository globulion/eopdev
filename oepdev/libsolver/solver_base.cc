#include "solver.h"

using namespace std;
using namespace psi;
using namespace oepdev;

OEPDevSolver::OEPDevSolver(SharedWavefunctionUnion wfn_union)
 : wfn_union_(wfn_union), methods_oepBased_(), methods_benchmark_(), options_(wfn_union->options())
{
  
}
OEPDevSolver::~OEPDevSolver() {}
double OEPDevSolver::compute_oep_based(const std::string& method) {}
double OEPDevSolver::compute_benchmark(const std::string& method) {}

// Build: factory static method 
std::shared_ptr<OEPDevSolver> OEPDevSolver::build(const std::string& target, SharedWavefunctionUnion wfn_union)
{
   std::shared_ptr<OEPDevSolver> solver;
   if      (target == "ELECTROSTATIC ENERGY"  ) solver = std::make_shared< ElectrostaticEnergySolver>(wfn_union);
   else if (target == "REPULSION ENERGY"      ) solver = std::make_shared<     RepulsionEnergySolver>(wfn_union);
   else if (target == "CHARGE TRANSFER ENERGY") solver = std::make_shared<ChargeTransferEnergySolver>(wfn_union);
   else if (target == "EET COUPLING CONSTANT" ) solver = std::make_shared<         EETCouplingSolver>(wfn_union);   
   else throw psi::PSIEXCEPTION("OEPDEV: Error. OEPDevSolver: Incorrect targer property chosen!\n");
   return solver;
}
