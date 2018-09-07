#!/usr/bin/python
"""
 -------------------------------------------------------------------
 Create the command for psi::timer code.
 -------------------------------------------------------------------

 It keeps to maintain the 'timer.dat' report file of Psi4 in order.
 
 Usage: [class] [message]

 Notes:
  - for 'class' do not exceed 6 characters. Use: 
    * 'Solver' - when timing inside oepdev::OEPDevSolver class
    * 'OEP'    - when timing inside oepdev::OEPotential class
  - for 'message' do not exceed 25 characters. 

 Examples:
  ./timer.py Solver "E(Paul)" OEP-Based
  creates
  psi::timer_on("Solver E(Paul) OEP-Basedw       ");
 -------------------------------------------------------------------
                             Last Revision: Gundelfingen, 7 Sep 2018
"""
from sys import argv, exit
if len(argv)==1:
   print __doc__; exit()

# read
t_1 = argv[1]
t_2 = " ".join(argv[2:])

# print
cmd = 'psi::timer_on("%7s%25s");' % (t_1.ljust(7), t_2.ljust(25))
print cmd
