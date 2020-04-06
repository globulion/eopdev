#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(wfn, opt)
{
   // Atoms are assumed to be distributed centres
   nSites_ = wfn_->molecule()->natom();

   // Set the number of sites for Ab Initio model
   nSitesAbInitio_= wfn_->doccpi()[0];
}
oepdev::NonUniformEFieldPolarGEFactory::~NonUniformEFieldPolarGEFactory()
{

}
// implementations of abstract methods from base
void oepdev::NonUniformEFieldPolarGEFactory::compute_samples(void)
{
   int nq   = options_.get_double("DMATPOL_NTEST_CHARGE");

   std::shared_ptr<psi::Matrix> D = wfn_->Da()->clone();

   for (int n=0; n<nSamples_; ++n) {
        std::shared_ptr<psi::Matrix> charges = std::make_shared<psi::Matrix>("Q", nq, 4);
        for (int i=0; i<nq; ++i) {
             std::shared_ptr<psi::Vector> point = draw_random_point();
             charges->set(i, 0, point->get(0));
             charges->set(i, 1, point->get(1));
             charges->set(i, 2, point->get(2));
             charges->set(i, 3, draw_charge());
        }
        std::vector<std::shared_ptr<psi::Vector>> fields_n;
        std::vector<std::shared_ptr<psi::Matrix>> grads_n;
        for (int o=0; o<nSites_; ++o) {
             double x = wfn_->molecule()->x(o);
             double y = wfn_->molecule()->y(o);
             double z = wfn_->molecule()->z(o);
             std::shared_ptr<psi::Vector> field = field_due_to_charges(charges, x, y, z);
             std::shared_ptr<psi::Matrix> grad  = field_gradient_due_to_charges(charges, x, y, z);
             fields_n.push_back(field);
             grads_n .push_back(grad);
             cout << oepdev::string_sprintf(" Computation for N=%2d S=%2d F=[%14.4f, %14.4f, %14.4f]\n",n+1, o+1,
                                              field->get(0), field->get(1), field->get(2));
             cout << oepdev::string_sprintf("                      grad F=[[%14.8f, %14.8f, %14.8f]\n",
                                              grad->get(0,0), grad->get(0,1), grad->get(0,2));
             cout << oepdev::string_sprintf("                              [%14.8f, %14.8f, %14.8f]\n",
                                              grad->get(1,0), grad->get(1,1), grad->get(1,2));
             cout << oepdev::string_sprintf("                              [%14.8f, %14.8f, %14.8f]]\n",
                                              grad->get(2,0), grad->get(2,1), grad->get(2,2));
        }

        std::shared_ptr<oepdev::RHFPerturbed> pert = perturbed_state(charges);
        referenceStatisticalSet_.DensityMatrixSet[n]->copy(pert->Da());
        referenceStatisticalSet_.DensityMatrixSet[n]->subtract(wfn_->Da());
        VMatrixSet_[n]->copy(pert->Vpert());


        electricFieldSet_        .push_back(fields_n);
        electricFieldGradientSet_.push_back(grads_n);

        referenceStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();
            modelStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();

        cout << oepdev::string_sprintf(" Interaction Energy = %15.6f\n", pert->reference_energy() - wfn_->reference_energy());

        if (hasAbInitioDipolePolarizability_) {
            std::vector<std::shared_ptr<psi::Vector>> fields_abini_n;
            for (int o=0; o<nSitesAbInitio_; ++o) {
                 std::shared_ptr<psi::Vector> field = field_due_to_charges(charges, abInitioPolarizationSusceptibilities_->centre(o));
                 fields_abini_n.push_back(field);
            }
            abInitioModelElectricFieldSet_.push_back(fields_abini_n);
            abInitioModelStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();
            referenceDpolStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy() 
                         + 2.0*D->vector_dot(VMatrixSet_[n]);
        }

   }
 
}
// abstract methods
void oepdev::NonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::NonUniformEFieldPolarGEFactory::compute_hessian(void)
{
}
