#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(cphf, opt)
{
   // Atoms are assumed to be distributed centres
   nSites_ = wfn_->molecule()->natom();
}
oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, cphf->options())
{
}
oepdev::NonUniformEFieldPolarGEFactory::~NonUniformEFieldPolarGEFactory()
{
}
// implementations of abstract methods from base
void oepdev::NonUniformEFieldPolarGEFactory::compute_samples(void)
{
   int nq   = options_.get_double("DMATPOL_NTEST_CHARGE");
   double q = options_.get_double("DMATPOL_TEST_CHARGE");

   for (int n=0; n<nSamples_; ++n) {
        std::shared_ptr<psi::Matrix> charges = std::make_shared<psi::Matrix>("Q", nq, 4);
        for (int i=0; i<nq; ++i) {
             std::shared_ptr<psi::Vector> point = draw_random_point();
             charges->set(i, 0, point->get(0));
             charges->set(i, 1, point->get(1));
             charges->set(i, 2, point->get(2));
             charges->set(i, 3, q);
        }
        std::vector<std::shared_ptr<psi::Vector>> fields_n;
        for (int o=0; o<nSites_; ++o) {
             double x = wfn_->molecule()->x(o);
             double y = wfn_->molecule()->y(o);
             double z = wfn_->molecule()->z(o);
             std::shared_ptr<psi::Vector> field = field_due_to_charges(charges, x, y, z);
             fields_n.push_back(field);
             cout << oepdev::string_sprintf(" Computation for N=%2d S=%2d F=[%14.4f, %14.4f, %14.4f]\n",n+1, o+1,
                                              field->get(0), field->get(1), field->get(2));
        }

        referenceDensityMatrixSet_[n]->copy(perturbed_dmat(charges));
        referenceDensityMatrixSet_[n]->subtract(wfn_->Da());

        electricFieldSet_         .push_back(fields_n);
   }
 
}
// abstract methods
void oepdev::NonUniformEFieldPolarGEFactory::compute_gradient(int i, int j)
{
}
void oepdev::NonUniformEFieldPolarGEFactory::compute_hessian(void)
{
}
