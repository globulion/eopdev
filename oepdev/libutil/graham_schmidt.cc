#include "graham_schmidt.h"

oepdev::GrahamSchmidt::GrahamSchmidt(std::vector<psi::SharedVector> V) 
 : L_(V.size()), 
   V_({})
{
  for (int i=0; i<V.size(); ++i) {
    this->V_.push_back(V[i]);
  }
}

oepdev::GrahamSchmidt::~GrahamSchmidt() {}

void oepdev::GrahamSchmidt::normalize() {
  for (int i=0; i<this->L_; ++i) {
       double d = this->V_[i]->norm();
       this->V_[i]->scale(1.0/d);
  }
}

void oepdev::GrahamSchmidt::append(psi::SharedVector d) {
  this->V_.push_back(d);
  this->L_ += 1;
}

void oepdev::GrahamSchmidt::reset(std::vector<psi::SharedVector> V) {
  this->V_.clear();
  this->L_ = V.size();
  for (int i=0; i<this->L_; ++i) {
    this->V_.push_back(V[i]);
  }
}

void oepdev::GrahamSchmidt::orthonormalize() {
  this->orthogonalize();
  this->normalize();
}

void oepdev::GrahamSchmidt::orthogonalize() {
  for (int k=1; k<this->L_; ++k) {
       for (int i=0; i<k; ++i) {
            this->V_[k]->subtract( this->projection(this->V_[i], this->V_[k]) );
       }
  }
}

psi::SharedVector oepdev::GrahamSchmidt::orthogonalize_vector(psi::SharedVector d, bool normalize) {

 for (int i=0; i<this->L_; ++i) {
      d->subtract( this->projection(this->V_[i], d) );
 }
 if (normalize) {
     double dn = d->norm();
     d->scale(1.0/dn);
 }
 return d;
}

psi::SharedVector oepdev::GrahamSchmidt::projection(psi::SharedVector u, psi::SharedVector v) {
  psi::SharedVector p = std::make_shared<psi::Vector>(*u);
  double un = u->norm();
  double d  = v->vector_dot(u) / (un * un);
  p->scale(d);
  return p;
}
