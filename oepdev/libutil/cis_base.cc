#include "cis.h"
#include<iostream>

//#include "psi4/libpsi4util/PsiOutStream.h"

namespace oepdev{

const std::vector<std::string> CISComputer::reference_types = {"RHF", "UHF"};

CISComputer::CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt, 
                         psi::IntegralTransform::TransformationType trans_type):
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
          nbvir_(nmo_ - nbocc_),
          transformation_type_(trans_type),
          inttrans_(nullptr)
{
  this->common_init(); 
}

CISComputer::~CISComputer() {}

void CISComputer::compute(void) {
 this->prepare_for_cis_();
 this->build_hamiltonian_();
 this->diagonalize_hamiltonian_(); 
 E_->print_out();
}

void CISComputer::prepare_for_cis_(void) {//TODO
 Fa_oo_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","OCC"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","OCC"), true, false, false);
 Fa_vv_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","VIR"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","VIR"), true, false, false);
 this->set_beta_();
 this->transform_integrals_();
}

void CISComputer::transform_integrals_(void) {
 SharedMOSpaceVector spaces;
 spaces.push_back(psi::MOSpace::occ);
 spaces.push_back(psi::MOSpace::vir);
 inttrans_ = std::make_shared<psi::IntegralTransform>(ref_wfn_, spaces, this->transformation_type_, 
                                                   psi::IntegralTransform::OutputType::DPDOnly,
                                                   psi::IntegralTransform::MOOrdering::QTOrder,
                                                   psi::IntegralTransform::FrozenOrbitals::None); 
 inttrans_->set_keep_dpd_so_ints(true);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::occ);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::vir, psi::MOSpace::vir);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::vir, psi::MOSpace::occ, psi::MOSpace::vir);

 //inttrans_ = std::shared_ptr<psi::IntegralTransform>(&inttrans);
}

void CISComputer::diagonalize_hamiltonian_(void) {
 H_->diagonalize(U_, E_);
}

std::shared_ptr<CISComputer> CISComputer::build(const std::string& type, 
                                                std::shared_ptr<psi::Wavefunction> ref_wfn, 
                                                psi::Options& opt, const std::string& reference) {

  // Determine reference if not specified
  std::string ref = reference;
  if (ref.empty()) {
      if (!ref_wfn->same_a_b_orbs() && !ref_wfn->same_a_b_dens()) {ref += "UHF";}
      else { ref += "RHF";}
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
 U_ = std::make_shared<psi::Matrix>("CIS Eigenvectors", ndets_, ndets_);
 E_ = std::make_shared<psi::Vector>("CIS Eigenvalues", ndets_);
}


} // EndNameSpace oepdev
