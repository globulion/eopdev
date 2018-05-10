#include <iostream>
#include "psi4/libmints/matrix.h"
#include "gefp.h"
#include "../libutil/util.h"

using namespace std;

oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf, psi::Options& opt)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf, opt)
{
}
oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::QuadraticGradientNonUniformEFieldPolarGEFactory(std::shared_ptr<CPHF> cphf)
 : oepdev::NonUniformEFieldPolarGEFactory(cphf)
{
}
oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::~QuadraticGradientNonUniformEFieldPolarGEFactory()
{
}
std::shared_ptr<oepdev::GenEffPar> oepdev::QuadraticGradientNonUniformEFieldPolarGEFactory::compute(void)
{
}
