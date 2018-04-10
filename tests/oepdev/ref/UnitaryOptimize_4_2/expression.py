#!/usr/bin/python
from sympy import *
A,B,C,D,E,F,G,H,I,J,K,x = symbols("A B C D E F G H I J K x")
w=(A*cos(x)+B*sin(x))*(C*sin(x)+D+E*cos(x))*(F*sin(x)+G+H*cos(x))*(I*sin(x)+J+K*cos(x))
we=expand(w)
#wt=trigsimp(we)
wt =\
A*C*F*I*sin(x)**3*cos(x) + A*C*F*J*sin(x)**2*cos(x) + A*C*F*K*(-cos(4*x) + 1)/8 + A*C*G*I*sin(x)**2*cos(x) + A*C*G*J*sin(2*x)/2 + A*C*G*K*sin(x)*cos(x)**2 + A*C*H*I*(-cos(4*x) + 1)/8 + A*C*H*J*sin(x)*cos(x)**2 + A*C*H*K*sin(x)*cos(x)**3 + A*D*F*I*sin(x)**2*cos(x) + A*D*F*J*sin(2*x)/2 + A*D*F*K*sin(x)*cos(x)**2 + A*D*G*I*sin(2*x)/2 + A*D*G*J*cos(x) + A*D*G*K*cos(x)**2 + A*D*H*I*sin(x)*cos(x)**2 + A*D*H*J*cos(x)**2 + A*D*H*K*cos(x)**3 + A*E*F*I*(-cos(4*x) + 1)/8 + A*E*F*J*sin(x)*cos(x)**2 + A*E*F*K*sin(x)*cos(x)**3 + A*E*G*I*sin(x)*cos(x)**2 + A*E*G*J*cos(x)**2 + A*E*G*K*cos(x)**3 + A*E*H*I*sin(x)*cos(x)**3 + A*E*H*J*cos(x)**3 + A*E*H*K*cos(x)**4 + B*C*F*I*sin(x)**4 + B*C*F*J*sin(x)**3 + B*C*F*K*sin(x)**3*cos(x) + B*C*G*I*sin(x)**3 + B*C*G*J*sin(x)**2 + B*C*G*K*sin(x)**2*cos(x) + B*C*H*I*sin(x)**3*cos(x) + B*C*H*J*sin(x)**2*cos(x) + B*C*H*K*(-cos(4*x) + 1)/8 + B*D*F*I*sin(x)**3 + B*D*F*J*sin(x)**2 + B*D*F*K*sin(x)**2*cos(x) + B*D*G*I*sin(x)**2 + B*D*G*J*sin(x) + B*D*G*K*sin(2*x)/2 + B*D*H*I*sin(x)**2*cos(x) + B*D*H*J*sin(2*x)/2 + B*D*H*K*sin(x)*cos(x)**2 + B*E*F*I*sin(x)**3*cos(x) + B*E*F*J*sin(x)**2*cos(x) + B*E*F*K*(-cos(4*x) + 1)/8 + B*E*G*I*sin(x)**2*cos(x) + B*E*G*J*sin(2*x)/2 + B*E*G*K*sin(x)*cos(x)**2 + B*E*H*I*(-cos(4*x) + 1)/8 + B*E*H*J*sin(x)*cos(x)**2 + B*E*H*K*sin(x)*cos(x)**3

print wt
