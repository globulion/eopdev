#include <iostream>
#include "test.h"
#include "../libutil/cphf.h"


using namespace std;

double oepdev::test::Test::test_cphf(void)
{
  // Reference data for MeNH2 at RHF/6-311++G** (6D)
  const double pol_ref[9] = { 22.375420,  -0.545125,   0.630563,
                              -0.545355,  20.939757,  -0.263902,  
                               0.630892,  -0.263922,  20.901605};

  std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn_, options_);
  solver->compute();
  for (int i=0; i<solver->nocc(); i++) {
       solver->lmo_centroid(i)->print();
       solver->polarizability(i)->print();
  }
  solver->polarizability()->print();

  // Accumulate errors
  double r_sum = 0.0;
  double** p = solver->polarizability()->pointer();
  const double*  r = pol_ref;
  for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
            r_sum += pow(p[i][j] - *r, 2.0);
            r++;
       }
  }

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}

