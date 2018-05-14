#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::UniformEFieldPolarGEFactory::UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(cphf, opt)
{
   // Set the number of sites to one since the electric field is uniform in space
   nSites_ = 1;
}
oepdev::UniformEFieldPolarGEFactory::UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::UniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::UniformEFieldPolarGEFactory::~UniformEFieldPolarGEFactory()
{
}
// implementations of abstract methods from base
void oepdev::UniformEFieldPolarGEFactory::compute_samples(void)
{
   for (int n=0; n<nSamples_; ++n) {
        std::vector<std::shared_ptr<psi::Vector>> fields;
        std::shared_ptr<psi::Vector> field = draw_field();

        fields.push_back(field); // only one element since field is uniform in space (nSites_ = 1)

        cout << oepdev::string_sprintf(" Computation for N=%2d F=[%14.4f, %14.4f, %14.4f]\n",n+1, 
                                         field->get(0), field->get(1), field->get(2));

        referenceDensityMatrixSet_[n]->copy(perturbed_dmat(field));
        referenceDensityMatrixSet_[n]->subtract(wfn_->Da());

        electricFieldSet_         .push_back(fields);
   }
}
// abstract methods
void oepdev::UniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::UniformEFieldPolarGEFactory::compute_hessian(void)
{
}
