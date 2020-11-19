#include <iostream>
#include "test.h"
#include "../libutil/unitary_optimizer.h"

using namespace std;


double oepdev::test::Test::test_unitaryOptimizer_2()
{
  double result = 0.0;

  double P[27] = {  -5.7237E-02,
    -4.8214E-01,  -8.0553E-01,  -1.5206E-01,  -7.4836E-01,  -2.9927E-01,  -4.5974E-01,  -5.1164E-02,
    -3.7566E-01,  -1.6202E-01,  -9.9207E-01,  -1.0660E-02,  -9.2229E-01,  -6.7787E-01,  -5.3848E-02,
    -9.9106E-01,  -1.7727E-01,  -1.3879E-01,  -5.6017E-01,  -7.4425E-01,  -1.9731E-01,  -5.2214E-01,
    -8.6566E-01,  -7.2151E-02,  -1.0403E-01,  -5.0845E-01,  -1.4330E-01 };

  // the reference Xmin and Xmax values were set by checking the output of the UnitaryOptimizer_2
  // for the given P tensor. To fully ensure the correct result, revise this test.
  

  const double Xmin_ref[9] = {  -0.44555350,   0.58229464,   0.680011039,
                                 0.42752510,   0.80575731,  -0.409850497,
                                -0.78657762,   0.10811146,  -0.607953579};
  
  const double Xmax_ref[9] = {   0.73412257,   0.55674902,  -0.388708844,
                                 0.65062607,  -0.74056654,   0.168068170,
                                -0.19429297,  -0.37628674,  -0.905902049};

  oepdev::UnitaryOptimizer_2 optimizer(P, 3, 1.0e-8, 100, true);

  psi::outfile->Printf(" ==> Unitary Minimization 2 in 3 Dimensions <==\n");
  bool success_min = optimizer.minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer.X();

  psi::outfile->Printf(" ==> Unitary Maximization 2 in 3 Dimensions <==\n");
  bool success_max = optimizer.maximize();
  std::shared_ptr<psi::Matrix> X_max = optimizer.X();

  // Accumulate errors
  double** xmin = X_min->pointer();
  double** xmax = X_max->pointer();
  const double*  xmin_ref = Xmin_ref;
  const double*  xmax_ref = Xmax_ref;
  for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
            result += pow(std::abs(xmin[i][j]) - std::abs(*xmin_ref), 2.0);
            result += pow(std::abs(xmax[i][j]) - std::abs(*xmax_ref), 2.0);
            xmin_ref++; 
            xmax_ref++;
       }
  }
  X_min->print();
  X_max->print();

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
