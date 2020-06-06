#include <iostream>
#include "test.h"
#include "../lib3d/dmtp.h"
#include "../libutil/wavefunction_union.h"

using namespace std;

double oepdev::test::Test::test_dmtp_superimposition(void) {
  // Sanity check for multimer test
  if (options_.get_str("OEPDEV_TEST_MODE") != "DIMER") 
     throw psi::PSIEXCEPTION("Monomer test mode cannot be used for this test. Set the OEPDEV_TEST_MODE to DIMER");

  // This test is for H2O dimer at HF/STO-3G (PyQuante-mod basis)
  // Rotation: mol_1 --> mol_2

  // Reference multipoles
  const double d_ref[9]  = {-3.850832E-02, 1.363762E-01, 1.714598E-01,
                            -2.937028E-02, 1.733716E-02, 1.012935E-02,
                             1.450756E-02, 1.517663E-02, 2.804010E-02};
  const double Q_ref[27] = {-3.420835E+00,-2.853787E-02, 4.404615E-02,-3.559506E+00, 3.795026E-02,-3.522726E+00,
                            -2.932453E-01, 3.363217E-02, 1.318864E-01,-4.449207E-01, 5.330635E-02,-3.346298E-01,
                            -2.231463E-01,-9.898882E-02,-1.694343E-02,-4.245383E-01, 2.529205E-02,-4.494886E-01};
  const double O_ref[81] = { 1.950784E-01,-2.185210E-01,-2.056819E-01, 6.251934E-02,-2.109664E-02,-1.484525E-01, 
                             2.734087E-03,-3.602960E-02,-7.984671E-02,-1.177391E-01,-4.956953E-01,-5.473696E-02,
                            -2.576221E-01,-1.013417E-01,-8.830077E-02,-2.919527E-01,-1.433590E-02,-7.563961E-02,
                            -1.116672E-01,-3.421141E-01, 7.485620E-01,-2.511674E-01,-4.635358E-02, 1.838533E-01, 
                             6.164565E-02, 1.124052E-01,-9.429049E-02,-2.160236E-02,-3.059421E-02, 4.263874E-02};
 const double H_ref[243] = {-6.570445E+00, 1.435883E-01, 1.437984E-02,-2.190978E+00,-1.128851E-01,-2.216735E+00, 
                             4.428762E-02, 8.239104E-03,-4.983846E-02,-1.236806E-01,-6.276044E+00,-4.164814E-02,
                            -2.121586E+00,-9.346266E-02,-6.418189E+00, 2.053049E-01, 1.306208E-01, 6.425287E-01,
                            -1.780228E-01, 1.454769E-01, 1.582921E-01, 9.699597E-02, 1.664590E-01, 1.648321E-01, 
                             6.498213E-01,-1.065813E+00, 1.476136E-01,-2.004955E-01, 2.897397E-01,-1.797658E-01, 
                             1.081619E+00,-6.714522E-01,-1.300105E-01, 6.813061E-02, 1.273056E-01,-7.915023E-02,
                            -3.516791E-01,-7.436144E-02,-1.149948E-01,-3.130480E-02,-9.105143E-01, 9.076238E-02,
                            -3.239872E-01, 6.272127E-02,-1.076322E+00};

  // Create WFN Union
  std::shared_ptr<oepdev::WavefunctionUnion> wfn_union = std::make_shared<oepdev::WavefunctionUnion>(wfn_, options_);

  // Compute CAMM's for the first monomer only
  std::shared_ptr<DMTPole> dmtp = oepdev::DMTPole::build("CAMM", wfn_union->l_wfn(0));
  dmtp->compute();

  // Superimpose the first DMTP onto the second molecule
  psi::SharedMatrix final_xyz = std::make_shared<psi::Matrix>(wfn_union->l_molecule(1)->geometry());
  dmtp->superimpose(final_xyz, {});

  dmtp->print();

  // Accumulate error stored in result
  double result = 0.0;
  double** pd = dmtp->dipoles(0)->pointer();
  double** pQ = dmtp->quadrupoles(0)->pointer();
  double** pO = dmtp->octupoles(0)->pointer();
  double** pH = dmtp->hexadecapoles(0)->pointer();
  const double* pref_d = d_ref;
  const double* pref_Q = Q_ref;
  const double* pref_O = O_ref;
  const double* pref_H = H_ref;

  for (int i=0; i<dmtp->n_sites(); ++i) {
       for (int z=0; z<3; ++z) {
            result += pow(pd[i][z] - *pref_d, 2.0);
            pref_d++;
       }
       for (int z=0; z<6; ++z) {
            result += pow(pQ[i][z] - *pref_Q, 2.0);
            pref_Q++;
       }
       for (int z=0; z<10; ++z) {
            result += pow(pO[i][z] - *pref_O, 2.0);
            pref_O++;
       }
       for (int z=0; z<15; ++z) {
          //result += pow(pH[i][z] - *pref_H, 2.0);
            pref_H++;
       }
  }
  result = sqrt(result);

  // Print result
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << result << std::endl;

  return result;
}
