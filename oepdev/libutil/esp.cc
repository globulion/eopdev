#include "esp.h"

namespace oepdev{

using namespace psi;
using namespace std;


Potential3D::Potential3D(SharedMatrix pot) : pot_(std::make_shared<Matrix>(pot)) 
{
  common_init();
}

Potential3D::Potential3D(SharedWavefunction wfn) 
{
  common_init();
}

Potential3D::Potential3D(SharedOEPotential oep, const std::string& oepType)
{
  common_init();
}

void Potential3D::common_init()
{

}

Potential3D::~Potential3D() 
{

}

std::shared_ptr<Potential3D> Potential3D::build(SharedMatrix pot)
{
   std::shared_ptr<Potential3D> potential = std::make_shared<BarePotential3D>(pot);
   return potential;
}
std::shared_ptr<Potential3D> Potential3D::build(const std::string& type, SharedWavefunction wfn)
{
   std::shared_ptr<Potential3D> potential;
   if     (type == "ELECTROSTATIC") potential = std::make_shared<ElectrostaticPotential3D>(wfn);
   else throw PSIEXCEPTION("OEPDEV: Potential3D init. Unknown type of potential specified.");
   return potential;
}
std::shared_ptr<Potential3D> Potential3D::build(SharedOEPotential oep, const std::string& oepType)
{
   std::shared_ptr<Potential3D> potential = std::make_shared<OEPotential3D>(oep, oepType);
   return potential;
}

void Potential3D::compute() 
{

}

/// BarePotential3D
BarePotential3D::BarePotential3D(SharedMatrix pot)
 : Potential3D(pot)
{
  common_init();
}

void BarePotential3D::common_init()
{
  
}

void BarePotential3D::compute() 
{

}


/// ElectrostaticPotential3D
ElectrostaticPotential3D::ElectrostaticPotential3D(SharedWavefunction wfn)
 : Potential3D(wfn)
{
  common_init();
}

void ElectrostaticPotential3D::common_init()
{

}

void ElectrostaticPotential3D::compute() 
{

}


/// OEPotential3D
OEPotential3D::OEPotential3D(SharedOEPotential oep, const std::string& oepType) 
 : Potential3D(oep, oepType)
{
  common_init();
}

void OEPotential3D::common_init()
{

}

void OEPotential3D::compute() 
{

}


} // EndNameSpace oepdev
