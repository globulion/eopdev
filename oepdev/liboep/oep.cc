#include "oep.h"
namespace oepdev_liboep{

OEPotential::OEPotential(const SharedWavefunction& wfn, Options& options) 
   : _wfn(new Wavefunction(*wfn)),
     _options(options),
     is_density_fitted(false),
     is_esp_based(true)
{
}

OEPotential::OEPotential(const SharedWavefunction& wfn, const SharedBasisSet& auxiliary, Options& options) 
   : _wfn(new Wavefunction(*wfn)),
     _auxiliary(new BasisSet(*auxiliary)),
     _options(options),
     is_density_fitted(true),
     is_esp_based(false)
{
}



//OEPotential::compute_3D(std::string oepType, std::string fileName) {
//  std::shared_ptr<CubicScalarGrid> grid(new CubicScalarGrid(_auxiliary, _options));
//  grid->write_cube_file(v_oep, oepType);
//}

} // EndNameSpace oepdev_liboep
