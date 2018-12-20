#include <iostream>
#include "test.h"
#include "../lib3d/space3d.h"

using namespace std;

double oepdev::test::Test::test_points_collection3d(void)
{
  // Reference data for H2O at RHF/STO-3G 
  const double points_ref[27] = {
  -6.64425E-03, -1.11079E-02, -1.21072E-02,
  -7.40143E-03, -1.65748E-02, -1.96386E-02,
  -1.65566E-03, -3.62771E-03, -6.75760E-03,
  -6.93719E-03, -1.50064E-02, -1.65116E-02,
  -9.07012E-03, -7.94549E-02, -1.00238E-01,
   4.82831E-03,  5.62557E-02,  1.34968E-02,
  -6.56789E-04, -1.66398E-04, -4.82657E-04,
   4.16409E-03,  3.55697E-02,  5.69803E-02,
   1.05254E-02,  3.28728E-02,  2.63736E-02};



  //Perform PointsCollection3D
  psi::timer_on("Test: ESPSolver   Calculation              ");
  shared_ptr<Field3D> field = oepdev::Field3D::build("ELECTROSTATIC POTENTIAL", 3, 3, 3, 4.897906, 4.565654, 4.492313, wfn_, options_, 1);
  field->compute();
  psi::timer_off("Test: ESPSolver   Calculation              ");

 
  field->write_cube_file("cube_file");

 // Iterator
  //shared_ptr<Points3DIterator> iter = oepdev::Points3DIterator::build(3, 3, 3, 3.898543, 3.781046, 3.162960, -10.799993, -2.848328, -4.963252);
  // its not good to construct a new iterator again manually: use already existing iterator object:
  shared_ptr<Points3DIterator> iter = field->points_collection()->points_iterator();

  //Accumulate errors
  double r_sum = 0.0;
  const double*  r = points_ref;

  for (iter->first(); iter->is_done()==false; iter->next()) {
       std::shared_ptr<psi::Vector> v = field->compute_xyz(iter->x(), iter->y(), iter->z());
       r_sum += pow(v->get(0) - *r, 2.0);
       r++;
  }


  

  // Print
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << " Test result= " << r_sum << std::endl;

  // Return
  return r_sum;


}
