#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticNonUniformEFieldPolarGEFactory::QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
}
oepdev::QuadraticNonUniformEFieldPolarGEFactory::QuadraticNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf)
{
}
oepdev::QuadraticNonUniformEFieldPolarGEFactory::~QuadraticNonUniformEFieldPolarGEFactory()
{
}
//std::shared_ptr<oepdev::GenEffPar> oepdev::QuadraticNonUniformEFieldPolarGEFactory::compute(void)
//{
//}
