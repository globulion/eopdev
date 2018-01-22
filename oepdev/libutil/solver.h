#ifndef _oepdev_libutil_solver_h
#define _oepdev_libutil_solver_h

#include<cstdio>
#include<string>
#include<map>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"

#include "wavefunction_union.h"
#include "integrals_iter.h"
#include "../liboep/oep.h"

namespace oepdev{

using namespace std;
using namespace psi;

using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedWavefunctionUnion  = std::shared_ptr<WavefunctionUnion>;
using SharedOEPotential        = std::shared_ptr<OEPotential>;


/**\brief Solver of properties of molecular aggregates. Abstract base.
 *
 *
 */
class OEPDevSolver : public std::enable_shared_from_this<OEPDevSolver>
{
 public:

   // <--- Constructor and Destructor ---> //
   
   /**\brief Take wavefunction union and initialize the Solver.
    * 
    * @param wfn_union   - wavefunction union of isolated molecular wavefunctions
    */
   OEPDevSolver(SharedWavefunctionUnion wfn_union);

   /// Destroctor
   virtual ~OEPDevSolver();


   // <--- Factory Methods ---> //

   /**\brief Build a solver of a particular property for given molecular cluster.
    *
    * @param target         - target property
    * @param wfn_union      - wavefunction union of isolated molecular wavefunctions
    */
   static std::shared_ptr<OEPDevSolver> build(const std::string& target, SharedWavefunctionUnion wfn_union);

  
   // <--- Computers ---> //

   /**\brief Compute property by using OEP's
    * 
    *  @param method - flavour of OEP model
    */
   virtual double compute_oep_based(const std::string& method = "DEFAULT") = 0;

   /**\brief Compute property by using benchmark method
    *
    *  @param method - benchmark method
    */
   virtual double compute_benchmark(const std::string& method = "DEFAULT") = 0;


 protected:

   /// Wavefunction union
   SharedWavefunctionUnion wfn_union_;

   /// Names of all OEP-based methods implemented for a solver
   std::vector<std::string> methods_oepBased_;

   /// Names of all benchmark methods implemented for a solver
   std::vector<std::string> methods_benchmark_;

};


class ElectrostaticEnergySolver : public OEPDevSolver
{
  public:
    ElectrostaticEnergySolver(SharedWavefunctionUnion wfn_union);
    virtual ~ElectrostaticEnergySolver();

    virtual double compute_oep_based(const std::string& method = "DEFAULT");
    virtual double compute_benchmark(const std::string& method = "DEFAULT");
   
};


}      // EndNameSpace oepdev
#endif //_oepdev_libutil_wavefunction_union_h

