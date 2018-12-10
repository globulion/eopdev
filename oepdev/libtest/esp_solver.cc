#include <iostream>
#include "test.h"
#include "../lib3d/esp.h"
#include "../lib3d/space3d.h"

using namespace std;

double oepdev::test::Test::test_esp_solver(void)
{
  // Reference data for MeNH2 at RHF/6-311++G** 
  const double charge_ref[7] = { -1.043052,  0.458551,   0.372558,
                                  0.374973, -0.030429,  -0.103094,
                                 -0.029507};


  // ESP settings
  const double pad = options_.get_double("ESP_PAD_SPHERE");
  const int npoints = options_.get_double("ESP_NPOINTS_PER_ATOM") * wfn_->molecule()->natom();

  // Perform ESP
  psi::timer_on("Test: ESPSolver   Calculation              ");
  std::shared_ptr<ElectrostaticPotential3D> field = std::make_shared<oepdev::ElectrostaticPotential3D>(npoints, pad, wfn_, options_);
  field->compute();
  std::shared_ptr<ESPSolver> esp = std::make_shared<oepdev::ESPSolver>(field);
  esp->compute();
  psi::timer_off("Test: ESPSolver   Calculation              ");

  esp->charges()->print();

  // Accumulate errors
  double r_sum = 0.0;
  double** c = esp->charges()->pointer();
  const double*  r = charge_ref;

  
  for (int i=0; i<7; ++i) {
       r_sum += pow(c[0][i] - *r, 2.0);
       r++;
  }
  r_sum = sqrt(r_sum);

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;


}
