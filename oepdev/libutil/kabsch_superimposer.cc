#include "kabsch_superimposer.h"

namespace oepdev{

void KabschSuperimposer::compute(psi::SharedMatrix initial_xyz, psi::SharedMatrix final_xyz) {
  const int n = initial_xyz->nrow();
  if (n!=final_xyz->nrow()) {throw  psi::PSIEXCEPTION(" Superimposition with incompatible initial and final coordinates!");}

  // Clear all previous calculations
  this->clear();

  // Copy coordinates
  psi::SharedMatrix c = std::make_shared<psi::Matrix>("",n,3);
  psi::SharedMatrix r = std::make_shared<psi::Matrix>("",n,3);

  for (int z=0; z<3; ++z) {
       c->set_column(0, z, initial_xyz->get_column(0, z));
       r->set_column(0, z,   final_xyz->get_column(0, z));
  }
  this->initial_xyz = std::make_shared<psi::Matrix>(c);
  this->final_xyz = std::make_shared<psi::Matrix>(r);

  // Subtract average
  double w = 1.0/(double)initial_xyz->nrow();
  psi::SharedMatrix c_av = c->collapse(0); c_av->scale(w);
  psi::SharedMatrix r_av = r->collapse(0); r_av->scale(w);

  double** pc = c->pointer();
  double** pr = r->pointer();
  for (int z=0; z<3; ++z) {
       double dc = c_av->get(z,0);
       double dr = r_av->get(z,0);
       for (int i=0; i<n; ++i) {
            pc[i][z] -= dc;
            pr[i][z] -= dr;
       }
  }

  // Covariance matrix
  psi::SharedMatrix A = psi::Matrix::doublet(c, r, true, false);

  // SVD
  psi::SharedMatrix U = std::make_shared<psi::Matrix>("",3,3);
  psi::SharedMatrix Vt= std::make_shared<psi::Matrix>("",3,3);
  psi::SharedVector S = std::make_shared<psi::Vector>("",3);
  A->svd(U, S, Vt);

  // Compute rotation matrix
  this->rotation = psi::Matrix::doublet(U, Vt, false, false);
  double** rr = this->rotation->pointer();
  double det = rr[0][0]*rr[1][1]*rr[2][2] + rr[0][1]*rr[1][2]*rr[2][0] + rr[0][2]*rr[1][0]*rr[2][1]
              -rr[0][2]*rr[1][1]*rr[2][0] - rr[0][1]*rr[1][0]*rr[2][2] - rr[0][0]*rr[1][2]*rr[2][1]; // Sarrus

  if (det < 0.0) {
      Vt->scale_row(0, 2, -1.0);
      this->rotation = psi::Matrix::doublet(U, Vt, false, false);
  }

  // Compute translation vector
  psi::SharedMatrix T = psi::Matrix::doublet(c_av, this->rotation, true, false);
  T->subtract(r_av);
  T->scale(-1.0);
  this->translation->set(0, T->get(0,0));
  this->translation->set(1, T->get(0,1));
  this->translation->set(2, T->get(0,2));

  this->translation->set_name("Kabsch Translation Vector");
  this->rotation->set_name("Kabsch Rotation Matrix");
}

void KabschSuperimposer::clear(void) {
  this->initial_xyz.reset();
  this->final_xyz.reset();
  this->rotation->identity();
  this->translation->zero();
}

psi::SharedMatrix KabschSuperimposer::get_transformed(void) {
  psi::SharedMatrix xyz = psi::Matrix::doublet(this->initial_xyz, this->rotation, false, false);
  double** p = xyz->pointer();
  for (int z=0; z<3; ++z) {
       double t = this->translation->get(z);
       for (int i=0; i<xyz->nrow(); ++i) {
            p[i][z] += t;
       }
  }
  return xyz;
}

double KabschSuperimposer::rms(void) {
  psi::SharedMatrix xyz = this->get_transformed();
  xyz->subtract(this->final_xyz);
  double rms = xyz->rms();
  return rms;
}

} // EndNameSpace oepdev
