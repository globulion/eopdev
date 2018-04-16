#ifndef _oepdev_libutil_scf_perturb_h
#define _oepdev_libutil_scp_perturb_h
/** @file scf_perturb.h */
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libscf_solver/rhf.h"

namespace oepdev {

using namespace std;

/** \addtogroup OEPDEV_UTILITIES
 * @{
 */

/**\brief Structure to hold perturbing charges
 *
 */
struct PerturbCharges
{
  /// Vector of charge values
  std::vector<double> charges;
  /// Vector of charge position vectors
  std::vector<std::shared_ptr<psi::Vector>> positions;
};

/**\brief RHF theory under electrostatic perturbation.
 *
 * Compute RHF wavefunction under the following conditions:
 *  - external uniform electric field
 *  - set of point charges
 * The mixed conditions can also be used.
 */
class RHFPerturbed : public psi::scf::RHF
{
  public:
    /// Build from wavefunction and superfunctional
    RHFPerturbed(std::shared_ptr<psi::Wavefunction> ref_wfn, std::shared_ptr<psi::SuperFunctional> functional);
    /// Build from wavefunction and superfunctional + options and psio
    RHFPerturbed(std::shared_ptr<psi::Wavefunction> ref_wfn, std::shared_ptr<psi::SuperFunctional> functional,
        psi::Options& options, std::shared_ptr<psi::PSIO> psio);

    /// Clear memory
    virtual ~RHFPerturbed();

    /// Compute total energy
    virtual double compute_energy();

    /// Perturb the system with external electric field
    virtual void set_perturbation(std::shared_ptr<psi::Vector> field);
    /*!
     * \overload
     */
    virtual void set_perturbation(const double& fx, const double& fy, const double& fz);
    /// Perturb the system with a point charge
    virtual void set_perturbation(std::shared_ptr<psi::Vector> position, const double& charge);
    /*!
     * \overload
     */
    virtual void set_perturbation(const double& rx, const double& ry, const double& rz, const double& charge);

  protected:

    /// Perturbing electric field
    std::shared_ptr<psi::Vector> perturbField_;
    /// Perturbing charges
    std::shared_ptr<PerturbCharges> perturbCharges_;
    
    /// Add the electrostatic perturbation to the Hcore matrix
    virtual void perturb_Hcore();
};

/** \example example_scf_perturb.cc
 * Perturb HF Hamiltonian with external electrostatic potential.
 * This is an example of how to use the oepdev::RHFPerturb class.
 */


/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_scf_perturbed_h
