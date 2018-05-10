#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
}
oepdev::LinearGradientNonUniformEFieldPolarGEFactory::LinearGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf)
{
}
oepdev::LinearGradientNonUniformEFieldPolarGEFactory::~LinearGradientNonUniformEFieldPolarGEFactory()
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::LinearGradientNonUniformEFieldPolarGEFactory::compute(void)
//{
//}
