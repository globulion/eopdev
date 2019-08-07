#include "cis.h"

namespace oepdev{

const std::vector<std::string> CISComputer::reference_types = {"RHF", "UHF"};

CISComputer::CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
          ref_wfn_(wfn),
          options_(opt),
          Fa_oo_(nullptr),
          Fa_vv_(nullptr),
          Fb_oo_(nullptr),
          Fb_vv_(nullptr),
          E_(nullptr),
          U_(nullptr),
          H_(nullptr),
          nmo_(ref_wfn_->nmo()),
          naocc_(ref_wfn_->nalpha()),
          nbocc_(ref_wfn_->nbeta()),
          navir_(nmo_ - naocc_),
          nbvir_(nmo_ - nbocc_)
{
  this->common_init(); 
}

CISComputer::~CISComputer() {}

void CISComputer::compute(void) {
 this->prepare_for_cis_();
 this->build_hamiltonian_();
 this->diagonalize_hamiltonian_(); 
}

void CISComputer::prepare_for_cis_(void) {//TODO
 Fa_oo_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","OCC"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","OCC"));
 Fa_vv_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","VIR"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","VIR"));
 this->set_beta_();
}

void CISComputer::build_hamiltonian_(void) {//TODO
}

void CISComputer::diagonalize_hamiltonian_(void) {//TODO
}

std::shared_ptr<CISComputer> CISComputer::build(const std::string& type, 
                                                std::shared_ptr<psi::Wavefunction> ref_wfn, 
                                                psi::Options& opt, const std::string& reference) {

  // Determine reference if not specified
  std::string ref = reference;
  if (ref.empty()) {
      ref += "RHF";
      if (!ref_wfn->same_a_b_orbs() && !ref_wfn->same_a_b_dens()) ref += "UHF";
  }

  // Sanity checks
  bool b = false;
  for (auto &refc : CISComputer::reference_types) {
       if (ref == refc) {b = true; break;}
  }
  if (!b) {throw psi::PSIEXCEPTION("Incorrect reference wavefunction type chosen. Only RHF and UHF are available");}

  if (ref =="RHF" and ref_wfn->molecule()->multiplicity() != 1)
   throw psi::PSIEXCEPTION("RHF reference cannot be set for open-shell system!");

  // Create
  std::shared_ptr<CISComputer> cis;

  if ((ref_wfn->molecule()->multiplicity() != 1) || (ref == "UHF")) 
     { cis = std::make_shared<U_CISComputer>(ref_wfn, opt); }
  else cis = std::make_shared<R_CISComputer>(ref_wfn, opt);
  
  // Return 
  return cis;
}

void CISComputer::set_beta_(void) {}

void CISComputer::common_init(void) {//TODO
 ndets_ = naocc_ * navir_ + nbocc_ * nbvir_;
 H_ = std::make_shared<psi::Matrix>("CIS Excited State Hamiltonian", ndets_, ndets_);
}


} // EndNameSpace oepdev
