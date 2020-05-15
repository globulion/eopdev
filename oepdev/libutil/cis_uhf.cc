#include "cis.h"
#include <iostream>

namespace oepdev{

U_CISComputer::U_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt, psi::IntegralTransform::TransformationType::Unrestricted)
{
}

U_CISComputer::~U_CISComputer() {}

void U_CISComputer::print_excited_state_character_(int I) {
//TODO
}



} // EndNameSpace oepdev
