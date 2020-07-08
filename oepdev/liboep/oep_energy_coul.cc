#include "oep.h"
#include "../lib3d/esp.h"
#include "../lib3d/dmtp.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;
using SharedDMTPole = std::shared_ptr<oepdev::DMTPole>;

// <============== Electrostatic Energy (demo class) ==============> //
SharedOEPotential ElectrostaticEnergyOEPotential::clone(void) const {
    SharedOEPotential temp = std::make_shared<ElectrostaticEnergyOEPotential>(this);
    return temp;
}
ElectrostaticEnergyOEPotential::ElectrostaticEnergyOEPotential(const ElectrostaticEnergyOEPotential* f) : OEPotential(f) {}
ElectrostaticEnergyOEPotential::ElectrostaticEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

ElectrostaticEnergyOEPotential::~ElectrostaticEnergyOEPotential() {}
void ElectrostaticEnergyOEPotential::common_init() 
{
   // In this demo OEPotential, there is only one OEP type. In addition,
   // only one OEP exists within this type: the electrostatic potential of entire molecule
   SharedMatrix mat = std::make_shared<psi::Matrix>("OEP V", wfn_->molecule()->natom(), 1);
   OEPType type = OEPType("V", false, 1, mat, nullptr, nullptr);
   oepTypes_["V"] = type;
}

void ElectrostaticEnergyOEPotential::compute(const std::string& oepType) 
{
  if (oepType == "V" || oepType == "TOTAL") {

      //psi::timer_on("OEP    E(Coul)                  ");

      // ESP charges
      SharedField3D potential = oepdev::Field3D::build("ELECTROSTATIC POTENTIAL", 
                                                   options_.get_int   ("ESP_NPOINTS_PER_ATOM") * wfn_->molecule()->natom(), 
                                                   options_.get_double("ESP_PAD_SPHERE"      ), 
                                                   wfn_, options_);
      potential->compute();
      oepdev::ESPSolver esp(potential);
      esp.compute();
     
      for (int i=0; i<esp.charges()->nrow(); ++i) {
           oepTypes_["V"].matrix->set(i, 0, esp.charges()->get(i, 0));
      }

      // CAMM
      SharedDMTPole camm = oepdev::DMTPole::build("CAMM", wfn_);
      camm->compute();
      oepTypes_["V"].dmtp = camm;

      //psi::timer_off("OEP    E(Coul)                  ");

  } else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
  }

}
double ElectrostaticEnergyOEPotential::compute_3D_V(const double& x, const double& y, const double& z){
      double val = 0.0;
      // ===> Nuclear contribution <=== //
      for (int i=0; i<wfn_->molecule()->natom(); ++i) {
           val+= (double)wfn_->molecule()->Z(i) /
                     sqrt( pow(wfn_->molecule()->x(i)-x, 2.0) 
                         + pow(wfn_->molecule()->y(i)-y, 2.0) 
                         + pow(wfn_->molecule()->z(i)-z, 2.0) );
      }

      // ===> Electronic contribution <=== //
      // double v;
      potInt_->set_charge_field(x, y, z);
      OEInt_ = potInt_;
      OEInt_->compute(potMat_);
      val += 2.0 * potMat_->vector_dot(wfn_->Da()); 
      //for (int i=0; i<primary_->nbf(); ++i) {
      //     for (int j=0; j<=i; ++j) {
      //          v= (wfn_->Da()->get(i,j) + wfn_->Db()->get(i,j)) * potMat_->get(i,j);
      //          val += 2.0*v;
      //          if (i==j) val -= v;
      //     }
      //}
      potMat_->zero();
 return val;
}
void ElectrostaticEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) 
{
   double val;
   if (oepType == "V" || oepType == "TOTAL") val = compute_3D_V(x, y, z);
   else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
   }

   // Assign final value
   v->set(0, val);
}
void ElectrostaticEnergyOEPotential::print_header(void) const 
{
   psi::outfile->Printf("  ==> OEPotential: %s <==\n\n", name_.c_str());
   oepTypes_.at("V").matrix->print();
}
void ElectrostaticEnergyOEPotential::rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) {
  // Potential "V"
  oepTypes_.at("V").dmtp->rotate(r);

  this->rotate_basic(r, R_prim, R_aux);
}
void ElectrostaticEnergyOEPotential::translate_oep(psi::SharedVector t) {
  // Potential "V"
  oepTypes_.at("V").dmtp->translate(t);

  this->translate_basic(t);
}
