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
 *
 * ##Theory
 * 
 * The electrostatic perturbation is here understood as a distribution of external
 * (generally non-uniform) electric field. It is assumed that this perturbation 
 * is one-electron in nature.
 * Therefore, the one-electron Hamiltonian is changed according to the following
 * \f[
 *   {\bf H}^{\rm core} \rightarrow {\bf H}^{\rm core} + \sum_n q_n {\bf V}^{(n)} - \mathbb{M} \cdot {\bf F}
 * \f]
 * where \f$ q_n \f$ is the external classical point charge, \f$ {\bf V}^{(n)} \f$ is the associated matrix 
 * of potential integrals, \f$ \mathbb{M} \f$ is the vector of dipole integrals and \f$ {\bf F} \f$ is
 * an external uniform electric field.
 * The total energy is then computed by performing an SCF procedure on the above one-electron Hamiltionian.
 * The contribution due to nuclei is included, i.e.,
 * \f[
 *    E_{\rm Nuc} \rightarrow E_{\rm Nuc-Nuc} + \sum_{In} \frac{q_n Z_I}{r_{In}} - {\bf \mu_{\rm Nuc}} \cdot {\bf F} 
 * \f]
 * where \f$ {\bf \mu_{\rm Nuc}} \f$ is the nuclear dipole moment and \f$ Z_I \f$ is the atomic number 
 * of the \f$ I \f$th nucleus.
 * It is added in the nuclear repulsion energy \f$ E_{\rm Nuc-Nuc} \f$ (note that the resulting energy can be negative as well
 * depending on the electric field direction and configuration of point charges.
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

    /// Get a copy of the perturbation potential one-electron matrix
    std::shared_ptr<psi::Matrix> Vpert() const {return std::make_shared<psi::Matrix>(Vpert_);}

  protected:

    /// Perturbing electric field
    std::shared_ptr<psi::Vector> perturbField_;
    /// Perturbing charges
    std::shared_ptr<PerturbCharges> perturbCharges_;

    /// Perturbation potential one-electron matrix
    std::shared_ptr<psi::Matrix> Vpert_;
    
    /// Add the electrostatic perturbation to the Hcore matrix
    virtual void perturb_Hcore();

  private:

    /// Initialize the Vpert_ as zero matrix
    void common_init();
};

/** \example example_scf_perturb.cc
 * Perturb HF Hamiltonian with external electrostatic potential.
 * This is an example of how to use the oepdev::RHFPerturbed class.
 */


/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_scf_perturbed_h
