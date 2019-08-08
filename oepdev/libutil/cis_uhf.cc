#include "cis.h"

namespace oepdev{

U_CISComputer::U_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt, psi::IntegralTransform::TransformationType::Unrestricted)
{
}

U_CISComputer::~U_CISComputer() {}

void U_CISComputer::set_beta_(void) {//TODO
// Fb_oo_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","OCC"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","OCC"), true, false, false);
// Fb_vv_ = psi::Matrix::triplet(ref_wfn_->Cb_subset("AO","VIR"), ref_wfn_->Fb(), ref_wfn_->Cb_subset("AO","VIR"), true, false, false);
 eps_b_o_ = ref_wfn_->epsilon_b_subset("MO", "OCC");
 eps_b_v_ = ref_wfn_->epsilon_b_subset("MO", "VIR");
}

void U_CISComputer::build_hamiltonian_(void) {//TODO
}

} // EndNameSpace oepdev
