#include "cis.h"
#include <iostream>

namespace oepdev{

R_CISComputer::R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt, psi::IntegralTransform::TransformationType::Restricted)
{
}

R_CISComputer::~R_CISComputer() {}

void R_CISComputer::print_excited_state_character_(int I) {
//TODO
}



} // EndNameSpace oepdev
