#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::LinearUniformEFieldPolarGEFactory::LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::UniformEFieldPolarGEFactory(cphf, opt)
{
}
oepdev::LinearUniformEFieldPolarGEFactory::LinearUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::UniformEFieldPolarGEFactory(cphf)
{
}
oepdev::LinearUniformEFieldPolarGEFactory::~LinearUniformEFieldPolarGEFactory()
{
}
std::shared_ptr<oepdev::GenEffPar> oepdev::LinearUniformEFieldPolarGEFactory::compute(void)
{
}
