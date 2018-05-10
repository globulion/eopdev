#include <iostream>
#include <random>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"
#include "../libutil/unitary_optimizer.h"
#include "../libutil/scf_perturb.h"

using namespace std;

oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(cphf, opt)
{
}
oepdev::NonUniformEFieldPolarGEFactory::NonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::GeneralizedPolarGEFactory(cphf)
{
}
oepdev::NonUniformEFieldPolarGEFactory::~NonUniformEFieldPolarGEFactory()
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::NonUniformEFieldPolarGEFactory::compute(void)
//{
//}
