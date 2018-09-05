//#include "psi4/libtrans/integraltransform.h"
//#include "psi4/libdpd/dpd.h"

#include "solver.h"

using namespace std;
using namespace psi;
using namespace oepdev;

// CT Solver//
ChargeTransferEnergySolver::ChargeTransferEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Benchmarks
  methods_benchmark_.push_back("OTTO_LADIK"      );
  methods_benchmark_.push_back("EFP2"            );
  // OEP-based
  methods_oepBased_ .push_back("OTTO_LADIK"      );
}
ChargeTransferEnergySolver::~ChargeTransferEnergySolver() {}
double ChargeTransferEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "OTTO_LADIK"      )  e = compute_oep_based_otto_ladik();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for repulsion energy calculations!\n");
  }
  return e;
}
double ChargeTransferEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" ||
           method == "OTTO_LADIK"   ) e = compute_benchmark_otto_ladik();
  else if (method == "EFP2"         ) e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
double ChargeTransferEnergySolver::compute_benchmark_otto_ladik(){}
double ChargeTransferEnergySolver::compute_benchmark_efp2(){}
double ChargeTransferEnergySolver::compute_oep_based_otto_ladik(){}
