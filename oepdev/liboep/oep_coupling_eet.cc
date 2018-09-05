#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

// <============== EET Coupling ==============> //

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, 
                                               SharedBasisSet auxiliary, 
                                               SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

EETCouplingOEPotential::~EETCouplingOEPotential() {}
void EETCouplingOEPotential::common_init() 
{
    int n1 = primary_->nbf();
    int n2 = auxiliary_->nbf();
    int n3 = wfn_->molecule()->natom();

    // Implement these matrices (here provide correct dimensions)
    SharedMatrix mat_1 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_2 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  
    SharedMatrix mat_3 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_4 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  
    SharedMatrix mat_5 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_6 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  

    // Provide correct classes (DF- or ESP-based) and number of OEP's involved in each type
    OEPType type_1 = {"Fujimoto.ET1", false,  0, mat_1};
    OEPType type_2 = {"Fujimoto.ET2", false,  0, mat_2};
    OEPType type_3 = {"Fujimoto.HT1", false,  0, mat_3};
    OEPType type_4 = {"Fujimoto.HT2", false,  0, mat_4};
    OEPType type_5 = {"Fujimoto.CT1", false,  0, mat_5};
    OEPType type_6 = {"Fujimoto.CT2", false,  0, mat_6};

    oepTypes_[type_1.name] = type_1;
    oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
    oepTypes_[type_4.name] = type_4;
    oepTypes_[type_5.name] = type_5;
    oepTypes_[type_6.name] = type_6;
}

void EETCouplingOEPotential::compute(const std::string& oepType) {}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void EETCouplingOEPotential::print_header(void) const {}
