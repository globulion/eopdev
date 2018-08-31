#include "oep_gdf.h"

using namespace oepdev;

// GeneralizedDensityFit - Abstract Base
GeneralizedDensityFit::GeneralizedDensityFit()
 : G_(nullptr), 
   H_(nullptr),
   bs_a_(nullptr),
   bs_i_(nullptr),
   ints_aa_(nullptr),
   ints_ai_(nullptr)
{
}
GeneralizedDensityFit::~GeneralizedDensityFit() {}
GeneralizedDensityFit::compute() {}
// SingleGeneralizedDensityFit
SingleGeneralizedDensityFit::(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                              std::shared_ptr<psi::Matrix> v_vector)
 : GeneralizedDensityFit()
{
// TODO
}
SingleGeneralizedDensityFit::~SingleGeneralizedDensityFit() {}
SingleGeneralizedDensityFit::compute()
{
// TODO
}
// DoubleGeneralizedDensityFit
DoubleGeneralizedDensityFit::(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                              std::shared_ptr<psi::BasisSet> bs_intermediate,
                              std::shared_ptr<psi::Matrix> v_vector)
 : GeneralizedDensityFit()
{
// TODO
}
DoubleGeneralizedDensityFit::~DoubleGeneralizedDensityFit() {}
DoubleGeneralizedDensityFit::compute()
{
// TODO
}
// Static Factory Methods
std::shared_ptr<GeneralizedDensityFit> 
GeneralizedDensityFit::build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                             std::shared_ptr<psi::Matrix> v_vector)
{
  return std::make_shared<SingleGeneralizedDensityFit>(bs_auxiliary, v_vector);
}
std::shared_ptr<GeneralizedDensityFit> 
GeneralizedDensityFit::build(std::shared_ptr<psi::BasisSet> bs_auxiliary,
                             std::shared_ptr<psi::BasisSet> bs_intermediate,
                             std::shared_ptr<psi::Matrix> v_vector)
{
  return std::make_shared<DoubleGeneralizedDensityFit>(bs_auxiliary, bs_intermediate, v_vector);
}

