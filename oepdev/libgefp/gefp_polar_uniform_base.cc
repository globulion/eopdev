#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::UniformEFieldPolarGEFactory::UniformEFieldPolarGEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(wfn, opt)
{
   // Set the number of sites to one since the electric field is uniform in space
   nSites_ = 1;
}
oepdev::UniformEFieldPolarGEFactory::~UniformEFieldPolarGEFactory()
{

}
// implementations of abstract methods from base
void oepdev::UniformEFieldPolarGEFactory::compute_samples(void)
{
   std::shared_ptr<psi::Matrix> D = wfn_->Da();

   for (int n=0; n<nSamples_; ++n) {
        std::vector<std::shared_ptr<psi::Vector>> fields;
        std::shared_ptr<psi::Vector> field = draw_field();

        fields.push_back(field); // only one element since field is uniform in space (nSites_ = 1)

        cout << oepdev::string_sprintf(" Computation for N=%2d F=[%14.4f, %14.4f, %14.4f]\n",n+1, 
                                         field->get(0), field->get(1), field->get(2));

        std::shared_ptr<oepdev::RHFPerturbed> pert = perturbed_state(field);
        referenceStatisticalSet_.DensityMatrixSet[n]->copy(pert->Da());
        referenceStatisticalSet_.DensityMatrixSet[n]->subtract(wfn_->Da());
        VMatrixSet_[n]->copy(pert->Vpert());

        electricFieldSet_.push_back(fields);

            referenceStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();
                modelStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();


        cout << oepdev::string_sprintf(" Interaction Energy = %15.6f\n", pert->reference_energy() - wfn_->reference_energy());

        if (hasAbInitioDipolePolarizability_) {
            for (int o=1; o<wfn_->doccpi()[0]; ++o) fields.push_back(field);
            abInitioModelElectricFieldSet_.push_back(fields);
            abInitioModelStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy();
            referenceDpolStatisticalSet_.InducedInteractionEnergySet[n] = pert->nuclear_interaction_energy() 
                         + 2.0*D->vector_dot(VMatrixSet_[n]);
        }
   }
}
// abstract methods
void oepdev::UniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::UniformEFieldPolarGEFactory::compute_hessian(void)
{
}
