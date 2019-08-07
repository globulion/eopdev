#include "cis.h"

namespace oepdev{

U_CISComputer::U_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt)
{
}

U_CISComputer::~U_CISComputer() {}

void U_CISComputer::set_beta_(void) {//TODO
 Fb_oo_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","OCC"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","OCC"));
 Fb_vv_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","VIR"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","VIR"));
}


} // EndNameSpace oepdev
