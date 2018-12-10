#include <iostream>
#include "test.h"


using namespace std;

double oepdev::test::Test::test_basic(void)
{
  // Read reference data
  const double ref_energy = -74.965002974633;
  const double ref_D[49]  = {     2.10706290,
      -0.45100674,     0.08310605,     0.06647334,    -0.02146005,    -0.02536451,
      -0.02539307,    -0.45100674,     1.98400112,    -0.46570102,    -0.37261261,
       0.12028659,    -0.04161382,    -0.04174240,     0.08310605,    -0.46570102,
       1.01664358,     0.23894684,    -0.01948050,     0.71844075,     0.02334713,
       0.06647334,    -0.37261261,     0.23894684,     0.99697621,     0.25666998,
      -0.10823398,     0.70131230,    -0.02146005,     0.12028659,    -0.01948050,
       0.25666998,     1.93244525,    -0.00345464,    -0.18802698,    -0.02536451,
      -0.04161382,     0.71844075,    -0.10823398,    -0.00345464,     0.61777967,
      -0.18193560,    -0.02539307,    -0.04174240,     0.02334713,     0.70131230,
      -0.18802698,    -0.18193560,     0.61752184};
  
  const double ref_S[49]  = {     1.00000002,
       0.23670394,     0.00000000,    -0.00000000,    -0.00000000,     0.04988317,
       0.04995228,     0.23670394,     1.00000003,     0.00000000,     0.00000000,
      -0.00000000,     0.45329197,     0.45366441,     0.00000000,     0.00000000,
       1.00000002,     0.00000000,    -0.00000000,     0.37260686,    -0.01487029,
      -0.00000000,     0.00000000,     0.00000000,     1.00000002,     0.00000000,
      -0.08261482,     0.36889805,    -0.00000000,    -0.00000000,    -0.00000000,
       0.00000000,     1.00000002,     0.00526845,    -0.09768319,     0.04988317,
       0.45329197,     0.37260686,    -0.08261482,     0.00526845,     0.99999999,
       0.23348072,     0.04995228,     0.45366441,    -0.01487029,     0.36889805,
      -0.09768319,     0.23348072,     0.99999999};
  
  const double ref_F[49]  = {   -20.23721180,
      -5.16293310,    -0.02178888,    -0.01744041,     0.00562972,    -1.10550392,
      -1.10702115,    -5.16293310,    -2.43580923,    -0.08638505,    -0.06927099,
       0.02235340,    -0.96807077,    -0.96899727,    -0.02178888,    -0.08638505,
      -0.31709770,    -0.01787751,     0.00151548,    -0.52223226,    -0.02449398,
      -0.01744041,    -0.06927099,    -0.01787751,    -0.31581677,    -0.01880895,
       0.07126681,    -0.50871736,     0.00562972,     0.02235340,     0.00151548,
      -0.01880895,    -0.38443448,     0.00448625,     0.13673079,    -1.10550392,
      -0.96807077,    -0.52223226,     0.07126681,     0.00448625,    -0.53571599,
      -0.35502764,    -1.10702115,    -0.96899727,    -0.02449398,    -0.50871736,
       0.13673079,    -0.35502764,    -0.53644448};
  
  const double ref_H[49]  = {   -32.68421101,
      -7.60411655,    -0.01359058,    -0.01088508,     0.00351329,    -1.61526323,
      -1.61749754,    -7.60411655,    -9.30126043,    -0.16169184,    -0.12947179,
       0.04179037,    -3.53236969,    -3.53577701,    -0.01359058,    -0.16169184,
      -7.54061793,     0.02880600,    -0.00271303,    -2.45623567,     0.04537431,
      -0.01088508,    -0.12947179,     0.02880600,    -7.54370603,     0.02897248,
       0.49288256,    -2.42256523,     0.00351329,     0.04179037,    -0.00271303,
       0.02897248,    -7.43777932,    -0.02094154,     0.64383980,    -1.61526323,
      -3.53236969,    -2.45623567,     0.49288256,    -0.02094154,    -4.94164328,
      -1.46364788,    -1.61749754,    -3.53577701,     0.04537431,    -2.42256523,
       0.64383980,    -1.46364788,    -4.94376112};

  // Error accumulator
  double r_sum = 0.0;

  int nbf = wfn_->basisset()->nbf();
  // Matrices
  std::shared_ptr<psi::Matrix> C = wfn_->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> D = wfn_->Da();
  std::shared_ptr<psi::Matrix> F = wfn_->Fa();
  D->scale(2.0); // Total density matrix

  // 1-Electron Integrals
  psi::IntegralFactory fact(wfn_->basisset());

  std::shared_ptr<psi::Matrix> S = std::make_shared<psi::Matrix>("S", nbf, nbf);
  std::shared_ptr<psi::Matrix> T = std::make_shared<psi::Matrix>("T", nbf, nbf);
  std::shared_ptr<psi::Matrix> V = std::make_shared<psi::Matrix>("V", nbf, nbf);

  std::shared_ptr<psi::OneBodyAOInt> sI(fact.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> tI(fact.ao_kinetic());
  std::shared_ptr<psi::OneBodyAOInt> vI(fact.ao_potential());

  sI->compute(S);
  tI->compute(T);
  vI->compute(V);
  T->add(V); // Construct the core Hamiltonian

  // Accumulate the error
  r_sum = sqrt(pow(wfn_->reference_energy() - ref_energy, 2.0));
  double** pD = D->pointer();
  double** pF = F->pointer();
  double** pS = S->pointer();
  double** pH = T->pointer();
  const double* pref_D = ref_D;
  const double* pref_F = ref_F;
  const double* pref_S = ref_S;
  const double* pref_H = ref_H;

  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            r_sum += pow(pD[i][j] - *pref_D, 2.0);
            r_sum += pow(pF[i][j] - *pref_F, 2.0);
            r_sum += pow(pS[i][j] - *pref_S, 2.0);
            r_sum += pow(pH[i][j] - *pref_H, 2.0);
            pref_D++;
            pref_F++;
            pref_S++;
            pref_H++;
       }
  }

  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test SCF energy= " << wfn_->reference_energy() << std::endl;
  std::cout << " Ref  SCF energy= " <<             ref_energy   << std::endl;
 
  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;
}
