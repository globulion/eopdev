#include <iostream>
#include "test.h"
#include "../libutil/unitary_optimizer.h"

using namespace std;


double oepdev::test::Test::test_unitaryOptimizer_4_2()
{
  double result = 0.0;

  double R[729] = {  -5.8298E-01,
    -2.7968E-01,  -9.9989E-01,  -6.9767E-01,  -8.5324E-01,  -9.0766E-01,  -8.1374E-01,  -6.5444E-01,
    -6.0323E-01,  -4.6118E-01,  -5.8081E-01,  -3.1478E-01,  -7.9555E-01,  -1.2188E-01,  -9.7261E-01,
    -3.2953E-01,  -5.8270E-01,  -4.4131E-01,  -8.5961E-01,  -8.0190E-01,  -1.9926E-01,  -3.1738E-02,
    -6.8658E-01,  -3.0768E-01,  -1.2361E-01,  -1.0539E-01,  -9.1496E-01,  -9.6095E-01,  -8.3017E-01,
    -1.2186E-01,  -9.0165E-01,  -5.7889E-01,  -4.2110E-02,  -4.6683E-01,  -3.0812E-01,  -6.8448E-01,
    -3.1350E-01,  -1.6537E-01,  -9.8171E-01,  -2.4986E-01,  -1.1139E-02,  -2.5183E-01,  -7.1956E-01,
    -2.1072E-01,  -8.9677E-01,  -5.5211E-01,  -9.1404E-02,  -7.0639E-01,  -7.1222E-01,  -8.6997E-01,
    -9.8063E-01,  -3.2116E-01,  -7.8837E-01,  -7.3445E-01,  -5.0843E-01,  -9.4664E-01,  -4.2588E-01,
    -8.5327E-01,  -4.1069E-01,  -3.0024E-01,  -8.9767E-01,  -5.8594E-01,  -3.0560E-01,  -5.8582E-01,
    -9.5005E-01,  -4.6410E-01,  -3.3621E-01,  -4.8511E-01,  -5.5405E-02,  -4.1344E-01,  -9.6598E-02,
    -8.6253E-01,  -8.6072E-01,  -1.9261E-01,  -6.0232E-01,  -8.3465E-01,  -7.2491E-02,  -6.5223E-01,
    -2.4919E-01,  -2.7400E-01,  -1.1669E-01,  -3.7633E-01,  -2.4906E-01,  -6.5110E-01,  -7.3007E-01,
    -1.0411E-01,  -5.7191E-01,  -3.5160E-02,  -3.3656E-01,  -3.7830E-01,  -8.8525E-01,  -5.0511E-02,
    -5.5009E-01,  -4.2161E-01,  -5.9186E-01,  -7.6297E-01,  -9.6620E-02,  -4.2632E-01,  -9.9713E-01,
    -3.8286E-01,  -6.7336E-01,  -4.7294E-01,  -1.1406E-01,  -6.4273E-01,  -9.1465E-02,  -3.7664E-01,
    -9.8418E-01,  -7.0563E-02,  -3.0910E-01,  -2.6771E-03,  -8.2766E-01,  -8.6286E-01,  -6.7405E-02,
    -3.0318E-01,  -9.3400E-01,  -2.4454E-01,  -2.4612E-01,  -7.6975E-02,  -2.8848E-01,  -8.7573E-01,
    -9.8012E-01,  -9.7379E-01,  -9.7169E-01,  -7.5379E-01,  -1.3997E-01,  -4.6117E-01,  -4.4718E-01,
    -1.5797E-01,  -8.7583E-01,  -7.2082E-01,  -4.1424E-01,  -3.0404E-02,  -4.3897E-01,  -9.8135E-01,
    -1.9937E-01,  -7.6703E-01,  -1.9289E-01,  -6.1214E-01,  -1.3646E-01,  -2.5288E-01,  -4.4376E-01,
    -8.6354E-01,  -9.4008E-01,  -8.7866E-01,  -9.5545E-01,  -8.9251E-01,  -7.7429E-01,  -2.8701E-01,
    -4.4028E-01,  -9.8744E-01,  -9.2803E-01,  -3.2724E-02,  -4.3190E-01,  -7.9671E-01,  -7.4767E-01,
    -2.5617E-01,  -8.0457E-01,  -4.1864E-01,  -2.9980E-02,  -1.5317E-01,  -7.6015E-01,  -5.0623E-01,
    -3.8004E-01,  -1.7102E-01,  -8.4321E-01,  -9.8142E-01,  -9.2998E-01,  -5.1365E-01,  -3.9367E-01,
    -4.3115E-01,  -6.8264E-01,  -1.1384E-02,  -4.2025E-01,  -6.1986E-01,  -4.4905E-01,  -2.5467E-01,
    -3.3077E-01,  -7.3508E-01,  -9.3367E-01,  -6.2992E-01,  -3.7028E-01,  -7.8983E-01,  -2.4724E-01,
    -9.3346E-01,  -7.3968E-01,  -1.9525E-01,  -8.0657E-01,  -3.6054E-01,  -4.7533E-01,  -7.5192E-02,
    -7.3670E-01,  -9.3404E-01,  -2.6493E-01,  -2.2782E-01,  -9.2184E-02,  -6.8028E-02,  -9.8605E-01,
    -7.6564E-01,  -3.8322E-01,  -5.0984E-02,  -4.9824E-02,  -4.4335E-01,  -8.4394E-02,  -3.5843E-01,
    -6.0999E-01,  -5.1401E-01,  -3.9569E-01,  -4.5045E-01,  -7.3819E-02,  -8.1267E-02,  -6.0512E-01,
    -3.6737E-02,  -8.2604E-01,  -8.7367E-01,  -8.6492E-01,  -4.9434E-01,  -9.7848E-01,  -5.2030E-02,
    -1.7288E-01,  -9.8498E-01,  -8.2380E-01,  -6.6794E-01,  -8.6900E-01,  -1.9051E-01,  -6.5526E-01,
    -5.9893E-02,  -4.1799E-01,  -1.2117E-01,  -1.5527E-01,  -9.4608E-02,  -5.4012E-01,  -4.5365E-01,
    -2.0140E-01,  -7.1428E-01,  -5.0975E-01,  -4.0089E-01,  -9.8447E-01,  -4.0652E-01,  -5.6632E-01,
    -1.9264E-01,  -6.8476E-01,  -1.0711E-01,  -4.2214E-01,  -8.1599E-01,  -2.1207E-01,  -3.8797E-01,
    -9.4609E-01,  -5.7981E-01,  -3.2093E-01,  -8.1398E-02,  -9.9960E-01,  -2.3241E-02,  -6.2342E-01,
    -2.6216E-02,  -3.9528E-01,  -1.7115E-01,  -4.2529E-01,  -3.7192E-01,  -7.1442E-01,  -4.1317E-01,
    -2.4998E-01,  -1.4169E-01,  -2.4492E-01,  -3.0194E-01,  -1.3552E-01,  -6.7732E-01,  -3.2921E-01,
    -5.4913E-01,  -6.1790E-01,  -5.8919E-01,  -5.9852E-01,  -6.8262E-01,  -3.7808E-01,  -5.6975E-01,
    -2.6198E-02,  -3.2220E-01,  -8.0143E-01,  -5.7330E-01,  -6.5665E-01,  -2.0236E-01,  -1.2000E-01,
    -9.6158E-02,  -3.3728E-01,  -7.2979E-01,  -7.4763E-01,  -1.4510E-01,  -4.7229E-01,  -1.9784E-01,
    -4.2751E-01,  -2.6686E-01,  -4.8099E-01,  -2.2912E-01,  -4.3114E-01,  -5.3429E-01,  -6.5731E-01,
    -9.3179E-01,  -6.2208E-01,  -9.2037E-01,  -1.7183E-02,  -8.1839E-01,  -1.8814E-01,  -1.2504E-01,
    -3.1159E-01,  -4.3051E-01,  -8.3903E-01,  -5.3312E-01,  -6.5483E-01,  -7.7496E-01,  -4.0749E-01,
    -6.8773E-01,  -8.3694E-02,  -9.0364E-02,  -7.4288E-01,  -8.8911E-01,  -8.0704E-01,  -5.0042E-01,
    -2.7141E-01,  -7.9181E-01,  -7.5197E-01,  -1.4833E-01,  -5.8415E-01,  -3.8331E-01,  -7.6633E-01,
    -8.9803E-01,  -4.8414E-01,  -5.2286E-01,  -8.4733E-01,  -3.7819E-01,  -4.5599E-01,  -3.4586E-01,
    -8.5545E-01,  -2.4847E-01,  -7.7795E-01,  -4.8065E-01,  -2.1470E-01,  -9.7767E-01,  -6.7564E-01,
    -1.2708E-01,  -1.5529E-01,  -4.6156E-01,  -1.3339E-01,  -5.0194E-02,  -1.7359E-01,  -1.4588E-01,
    -9.0126E-01,  -3.4870E-01,  -2.9648E-01,  -3.8976E-01,  -2.0038E-01,  -9.6543E-01,  -2.2976E-01,
    -2.6827E-01,  -7.4030E-01,  -7.4293E-01,  -3.6770E-01,  -6.5470E-01,  -2.0341E-01,  -5.5385E-01,
    -2.1725E-01,  -9.5282E-03,  -6.9975E-01,  -8.5699E-01,  -9.8692E-02,  -4.5844E-01,  -2.5260E-02,
    -3.6340E-01,  -6.0870E-03,  -4.5393E-01,  -4.7357E-01,  -8.6457E-01,  -6.4429E-01,  -9.7378E-01,
    -8.3960E-01,  -2.5436E-01,  -9.6960E-01,  -6.3346E-01,  -1.3765E-01,  -3.0732E-01,  -3.0906E-01,
    -8.1136E-01,  -5.5810E-01,  -4.1842E-01,  -1.0248E-02,  -7.9609E-01,  -7.5227E-01,  -7.3783E-01,
    -2.4983E-01,  -5.4302E-01,  -9.4307E-01,  -4.9148E-01,  -7.8804E-01,  -2.0140E-01,  -7.0267E-01,
    -9.7239E-01,  -4.0657E-01,  -1.5616E-01,  -6.1898E-01,  -2.5014E-01,  -4.8886E-01,  -4.5905E-01,
    -4.0566E-02,  -1.9604E-01,  -9.6768E-01,  -2.9061E-01,  -5.3500E-01,  -5.2451E-02,  -7.7857E-01,
    -7.3293E-01,  -9.1853E-01,  -5.7138E-01,  -8.9098E-01,  -3.6621E-01,  -1.9704E-01,  -3.0320E-01,
    -2.3379E-01,  -6.5755E-01,  -1.5415E-01,  -5.7123E-01,  -1.7599E-01,  -3.7350E-01,  -8.5658E-01,
    -9.2161E-01,  -9.8167E-01,  -9.3328E-01,  -5.4142E-01,  -8.8666E-01,  -9.7222E-01,  -2.4514E-01,
    -6.0515E-01,  -2.5306E-01,  -5.4760E-01,  -5.4991E-01,  -5.2193E-01,  -5.2600E-01,  -1.9684E-01,
    -5.9761E-01,  -9.5314E-02,  -9.6294E-01,  -2.2613E-01,  -8.7436E-01,  -3.8149E-01,  -9.8964E-01,
    -4.6137E-01,  -9.9698E-01,  -4.8806E-02,  -9.4598E-02,  -2.0403E-01,  -8.4726E-02,  -8.5444E-01,
    -8.4227E-01,  -8.1237E-01,  -3.7750E-01,  -9.4191E-02,  -1.0045E-02,  -2.8888E-01,  -2.6820E-01,
    -9.0707E-02,  -5.9913E-01,  -7.5015E-01,  -8.2657E-01,  -8.8054E-01,  -1.8739E-01,  -8.5321E-01,
    -7.3570E-01,  -1.8091E-01,  -6.8941E-01,  -1.7583E-02,  -7.3336E-01,  -4.6635E-01,  -6.8553E-01,
    -8.9227E-02,  -6.3344E-01,  -5.6641E-01,  -4.8771E-01,  -6.1114E-02,  -9.6905E-01,  -2.8312E-01,
    -1.0898E-01,  -9.7271E-01,  -4.7795E-01,  -6.7401E-01,  -1.4051E-01,  -4.4148E-01,  -3.0977E-01,
    -5.4715E-01,  -3.7169E-01,  -7.0990E-01,  -9.9065E-01,  -4.2324E-01,  -6.8856E-01,  -4.8273E-01,
    -8.3594E-02,  -5.7353E-01,  -7.5260E-01,  -6.2871E-01,  -6.8139E-02,  -6.3132E-02,  -1.5567E-01,
    -7.9793E-02,  -7.7210E-01,  -9.1252E-01,  -7.7269E-01,  -6.8562E-01,  -8.2523E-01,  -3.9291E-01,
    -5.8641E-01,  -1.8365E-01,  -8.1487E-01,  -2.9812E-01,  -7.5964E-01,  -4.2578E-01,  -6.5101E-01,
    -9.4304E-01,  -7.7119E-01,  -3.3590E-01,  -5.0275E-01,  -4.8098E-01,  -8.2528E-01,  -4.2928E-01,
    -3.2466E-03,  -1.8316E-01,  -4.0563E-01,  -2.4011E-02,  -9.8437E-02,  -4.0439E-01,  -9.6757E-01,
    -9.0642E-01,  -9.3463E-01,  -5.4827E-01,  -6.2457E-01,  -2.4650E-02,  -8.3202E-01,  -2.7212E-02,
    -2.3253E-01,  -1.7576E-01,  -3.6738E-01,  -3.3127E-01,  -5.2312E-01,  -9.8686E-01,  -6.4699E-01,
    -5.0793E-01,  -2.6991E-01,  -5.3137E-01,  -5.4260E-01,  -8.6234E-01,  -9.8911E-01,  -2.4172E-01,
    -6.8005E-01,  -1.5617E-02,  -7.7977E-01,  -6.6129E-01,  -4.7610E-01,  -2.4511E-01,  -5.3614E-01,
    -8.7518E-01,  -6.8750E-01,  -4.9548E-01,  -3.2615E-01,  -2.2985E-01,  -8.6966E-01,  -9.7708E-01,
    -4.8092E-01,  -1.9001E-01,  -9.8740E-01,  -3.2753E-01,  -3.1319E-01,  -5.5075E-01,  -8.5211E-02,
    -3.5564E-01,  -9.9476E-01,  -5.1557E-01,  -1.4068E-01,  -1.6960E-01,  -3.5085E-01,  -3.2630E-01,
    -4.2150E-01,  -7.2588E-01,  -4.3947E-01,  -3.2827E-01,  -6.4757E-01,  -1.4417E-01,  -8.0496E-01,
    -2.5268E-01,  -7.1040E-01,  -2.2620E-01,  -5.7226E-01,  -1.9230E-01,  -6.4647E-01,  -7.8631E-01,
    -2.3272E-01,  -6.9136E-01,  -2.6675E-01,  -2.5553E-01,  -7.7860E-01,  -7.8589E-01,  -8.0105E-01,
    -8.5748E-01,  -6.2292E-01,  -9.7337E-01,  -8.8908E-01,  -3.2544E-01,  -2.0022E-01,  -9.1947E-01,
    -7.6830E-01,  -7.9237E-01,  -8.2666E-02,  -2.8869E-01,  -4.4612E-01,  -6.9548E-01,  -1.6515E-01,
    -5.6469E-01,  -7.6544E-02,  -2.9395E-01,  -5.2197E-01,  -8.7379E-01,  -2.3956E-02,  -8.4017E-01,
    -7.9740E-01,  -5.6882E-01,  -5.9580E-01,  -8.5325E-01,  -2.7068E-01,  -8.1125E-01,  -3.5610E-01,
    -2.4569E-01,  -7.8927E-01,  -3.9905E-01,  -2.5107E-01,  -3.6178E-01,  -4.0287E-01,  -7.0452E-01,
    -2.6839E-01,  -5.4692E-02,  -5.7444E-01,  -2.1782E-01,  -9.4386E-01,  -1.6473E-01,  -8.0775E-01,
    -6.0490E-01,  -6.9992E-01,  -9.1990E-01,  -9.5369E-02,  -6.2985E-01,  -4.6930E-01,  -5.0588E-01,
    -8.6784E-01,  -7.9355E-01,  -9.2381E-01,  -4.9208E-01,  -7.3845E-01,  -6.4294E-01,  -8.9193E-01,
    -2.1245E-01,  -8.9342E-01,  -1.4291E-02,  -8.2284E-01,  -4.2759E-01,  -9.5515E-01,  -2.1288E-01,
    -8.1039E-01,  -4.7210E-01,  -2.5992E-01,  -8.5007E-01,  -4.4891E-01,  -7.8338E-01,  -2.4080E-01,
    -2.7708E-01,  -8.2345E-01,  -1.3803E-01,  -9.8022E-01,  -1.3976E-01,  -4.4110E-01,  -5.9678E-01,
    -2.4125E-01,  -2.8307E-01,  -1.2674E-02,  -7.2191E-01,  -9.9621E-01,  -6.6097E-02,  -1.4210E-01,
    -2.7115E-01,  -4.8331E-01,  -2.9304E-01,  -2.1947E-01,  -6.2512E-01,  -2.2968E-01,  -2.4938E-01,
    -3.8679E-01,  -5.9813E-01,  -3.0269E-01,  -9.9689E-01,  -2.2510E-01,  -1.0358E-01,  -7.6068E-01,
    -8.7923E-01,  -7.7972E-01,  -6.9790E-01,  -1.1697E-01,  -4.5683E-01,  -7.1329E-01,  -8.6165E-01,
    -7.0986E-01,  -3.8613E-01,  -6.7586E-01,  -5.4264E-01,  -5.5588E-01,  -1.7186E-01,  -5.7365E-01,
    -6.5430E-01,  -3.2503E-01,  -7.7852E-01,  -5.3275E-01,  -6.8523E-01,  -3.7314E-01,  -1.2264E-01,
    -5.5231E-01,  -2.1554E-01,  -5.4303E-01,  -3.4377E-01,  -8.6816E-01,  -5.6702E-01,  -9.0688E-02,
    -3.9452E-01,  -2.3323E-01,  -4.9530E-01,  -5.0194E-01,  -1.5710E-01,  -9.3219E-01,  -4.2673E-01 };
  
  double P[27] = {  -5.7237E-02,
    -4.8214E-01,  -8.0553E-01,  -1.5206E-01,  -7.4836E-01,  -2.9927E-01,  -4.5974E-01,  -5.1164E-02,
    -3.7566E-01,  -1.6202E-01,  -9.9207E-01,  -1.0660E-02,  -9.2229E-01,  -6.7787E-01,  -5.3848E-02,
    -9.9106E-01,  -1.7727E-01,  -1.3879E-01,  -5.6017E-01,  -7.4425E-01,  -1.9731E-01,  -5.2214E-01,
    -8.6566E-01,  -7.2151E-02,  -1.0403E-01,  -5.0845E-01,  -1.4330E-01 };
  
  const double Xmin_ref[9] = {   4.4346E-01,
    -6.1498E-01,   6.5203E-01,  -7.4667E-01,  -6.5590E-01,  -1.1081E-01,   4.9581E-01,  -4.3771E-01,
    -7.5006E-01 };
  
  const double Xmax_ref[9] = {   6.9716E-01,
     2.5450E-01,  -6.7022E-01,   3.4562E-01,  -9.3837E-01,   3.1831E-03,  -6.2811E-01,  -2.3386E-01,
    -7.4215E-01 };

  oepdev::UnitaryOptimizer_4_2 optimizer(R, P, 3, 1.0e-8, 100, true);

  psi::outfile->Printf(" ==> Unitary Minimization 4_2 in 3 Dimensions <==\n");
  bool success_min = optimizer.minimize();
  std::shared_ptr<psi::Matrix> X_min = optimizer.X();

  psi::outfile->Printf(" ==> Unitary Maximization 4_2 in 3 Dimensions <==\n");
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
