#include "gefp.h"
#include <iostream>
#include <fstream>

using namespace std;

//-- OEP_EFP2_GEFactory --////////////////////////////////////////////////////////////////////////////////
oepdev::OEP_EFP2_GEFactory::OEP_EFP2_GEFactory(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt) :
 oepdev::EFP2_GEFactory(wfn, opt) {}
oepdev::OEP_EFP2_GEFactory::~OEP_EFP2_GEFactory() {}
std::shared_ptr<oepdev::GenEffPar> oepdev::OEP_EFP2_GEFactory::compute()
{

}
