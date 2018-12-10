#include <iostream>
#include "test.h"
#include "../libutil/unitary_optimizer.h"

using namespace std;


double oepdev::test::Test::test_unitaryOptimizer()
{
  double result = 0.0;
  const double P[5] = {-0.2563, -0.95325, 1.535253, 0.0553, 2.04};
  const double R[25]= {0.242412, -0.3324216, -0.61421416, -0.9023525, 0.072185, 
                       -1.4562215367, -1.52156, -0.58757, 0.305301,
                        0.45, 0.9, 1.1, -0.5, 0.01, 0.4, 0.723, -0.7512, 0.04, -0.1121, 
                        -1.533, 0.553, -0.04566, 0.05654, 0.1432, -1.04};
  const double Xmin_ref[25] = { -8.74267861e-01,  -3.17344968e-01,  -3.67352484e-01,  -6.75418301e-05,
                                 1.14274386e-05,   2.71298173e-01,  -7.21888732e-01,  -2.19297103e-02,    
                                -6.27836730e-01,   1.03024980e-01,  -3.39007266e-01,  -1.42396112e-01,
                                 9.29823205e-01,  -1.48347999e-02,   2.47344382e-03,   2.13741091e-01,
                                -5.89011051e-01,   2.07712921e-06,   7.77760943e-01,   4.96864340e-02,
                                -3.79713378e-02,   1.04681613e-01,  -3.67038294e-05,   2.62477372e-02,
                                 9.93433921e-01};
  const double Xmax_ref[25] = {-0.96081674, -0.17487789,  0.08508092, -0.03368997,  0.19461540, -0.07491281,
                               -0.40602029, -0.82986219, -0.01744834, -0.37491303, -0.18717579,  0.87781524,  
                               -0.43742598, -0.00431968,  0.05518537, -0.05791162,  0.01419595,  0.02906414,
                                0.99123309, -0.11426684, -0.18119255,  0.18385879,  0.33452004, -0.12648636,
                               -0.89747569};

  std::shared_ptr<psi::Matrix> Rm = std::make_shared<psi::Matrix>("R matrix",5,5);
  std::shared_ptr<psi::Vector> Pm = std::make_shared<psi::Vector>("P vector",5);
  for (int i=0; i< 5; ++i) {
       Pm->set(i, P[i]);
       for (int j=0; j<5; ++j) {
            Rm->set(i,j,R[5*i+j]);
       }
  }
  oepdev::UnitaryOptimizer optimizer(Rm, Pm, 1.0e-8, 100, true);

  psi::outfile->Printf(" ==> Unitary Minimization in 5 Dimensions <==\n");
  bool success_min = optimizer.minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer.X();

  psi::outfile->Printf(" ==> Unitary Maximization in 5 Dimensions <==\n");
  bool success_max = optimizer.maximize();
  std::shared_ptr<psi::Matrix> X_max = optimizer.X();

  // Accumulate errors
  double** xmin = X_min->pointer();
  double** xmax = X_max->pointer();
  const double*  xmin_ref = Xmin_ref;
  const double*  xmax_ref = Xmax_ref;
  for (int i=0; i<5; ++i) {
       for (int j=0; j<5; ++j) {
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
