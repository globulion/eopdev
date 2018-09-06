#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

// <============== CT Energy ===============> //

ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   throw psi::PSIEXCEPTION("OEPDEV: Construction of Charge Transfer Energy OEP requires auxiliary basis set!\n");
   common_init();
}

ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, 
                                                                 SharedBasisSet auxiliary, 
                                                                 SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

ChargeTransferEnergyOEPotential::~ChargeTransferEnergyOEPotential() {}
void ChargeTransferEnergyOEPotential::common_init() 
{
    int n1 = primary_->nbf();
    int n2 = auxiliary_->nbf();
    int n3 = wfn_->molecule()->natom();

    SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G ", n2, n1);
    SharedMatrix mat_2 = std::make_shared<psi::Matrix>("Q1", n3, n1);
    SharedMatrix mat_3 = std::make_shared<psi::Matrix>("Q2", n3, 1);

    OEPType type_1 = {"Otto-Ladik.V1", true , n1, mat_1};
    OEPType type_2 = {"Otto-Ladik.V2", false, n1, mat_2};
    OEPType type_3 = {"Otto-Ladik.V3", false,  1, mat_3};

    oepTypes_[type_1.name] = type_1;
    oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
}

void ChargeTransferEnergyOEPotential::compute(const std::string& oepType) {}
void ChargeTransferEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void ChargeTransferEnergyOEPotential::print_header(void) const {}

void ChargeTransferEnergyOEPotential::compute_otto_ladik_v1() {}
void ChargeTransferEnergyOEPotential::compute_otto_ladik_v2() {}
void ChargeTransferEnergyOEPotential::compute_otto_ladik_v3() {}
