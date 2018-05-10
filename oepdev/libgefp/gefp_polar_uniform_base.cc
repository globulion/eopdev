#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::UniformEFieldPolarGEFactory::UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::GeneralizedPolarGEFactory(cphf, opt)
{
}
oepdev::UniformEFieldPolarGEFactory::UniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::GeneralizedPolarGEFactory(cphf)
{
}
oepdev::UniformEFieldPolarGEFactory::~UniformEFieldPolarGEFactory()
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::UniformEFieldPolarGEFactory::compute(void)
//{
//}
