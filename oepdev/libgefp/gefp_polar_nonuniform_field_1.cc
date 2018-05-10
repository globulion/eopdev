#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearNonUniformEFieldPolarGEFactory::LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
}
oepdev::LinearNonUniformEFieldPolarGEFactory::LinearNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf)
{
}
oepdev::LinearNonUniformEFieldPolarGEFactory::~LinearNonUniformEFieldPolarGEFactory()
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::LinearNonUniformEFieldPolarGEFactory::compute(void)
//{
//}
