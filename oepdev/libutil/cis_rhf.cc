#include "cis.h"

namespace oepdev{

R_CISComputer::R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt)
{
}

R_CISComputer::~R_CISComputer() {}

void R_CISComputer::set_beta_(void) {//TODO
 Fb_oo_ = Fa_oo_;
 Fb_vv_ = Fa_vv_;
}


} // EndNameSpace oepdev
