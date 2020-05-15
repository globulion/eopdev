#include "cis.h"
#include <iostream>
#include "psi4/libpsi4util/process.h"

namespace oepdev{

R_CISComputer::R_CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt):
 CISComputer(wfn, opt, psi::IntegralTransform::TransformationType::Restricted)
{
}

R_CISComputer::~R_CISComputer() {}

void R_CISComputer::print_excited_state_character_(int I) {
 const double amplitude_print_threshold = this->options_.get_double("OEPDEV_AMPLITUDE_PRINT_THRESHOLD");

 for (int i=0; i<this->naocc_; ++i) {
 for (int a=0; a<this->navir_; ++a) {
      int ia = this->navir_*i + a;
      double t_ia = this->U_->get(ia,I);

      if (std::abs(t_ia) >= amplitude_print_threshold) {

      int h = this->naocc_ - i - 1;
      int l = a;

      psi::outfile->Printf("       T(H-%3d ---> L+%3d)= %9.6f\n", h, l, t_ia);

      }

 }
 }
 psi::outfile->Printf("\n");
}



} // EndNameSpace oepdev
