#include "gram_schmidt.h"

oepdev::GramSchmidt::GramSchmidt(std::vector<psi::SharedVector> V) 
 : L_(V.size()), 
   V_({})
{
  this->reset(V);
}

oepdev::GramSchmidt::GramSchmidt()  
 : L_(0), 
   V_({})
{
}

oepdev::GramSchmidt::~GramSchmidt() {}

void oepdev::GramSchmidt::normalize() {
  for (int i=0; i<this->L_; ++i) {
       double d = this->V_[i]->norm();
       this->V_[i]->scale(1.0/d);
  }
}

void oepdev::GramSchmidt::append(psi::SharedVector d) {
  this->V_.push_back(d);
  this->L_ += 1;
}

void oepdev::GramSchmidt::reset(std::vector<psi::SharedVector> V) {
  this->V_.clear();
  this->L_ = V.size();
  for (int i=0; i<this->L_; ++i) {
    this->V_.push_back(V[i]);
  }
}

void oepdev::GramSchmidt::reset() {
  this->V_.clear();
  this->L_ = 0;
}

void oepdev::GramSchmidt::orthonormalize() {
  this->orthogonalize();
  this->normalize();
}

void oepdev::GramSchmidt::orthogonalize() {
  for (int k=1; k<this->L_; ++k) {
       for (int i=0; i<k; ++i) {
            this->V_[k]->subtract( this->projection(this->V_[i], this->V_[k]) );
       }
  }
}

void oepdev::GramSchmidt::orthogonalize_vector(psi::SharedVector& d, bool normalize) const {

 for (int i=0; i<this->L_; ++i) {
      d->subtract( this->projection(this->V_[i], d) );
 }
 if (normalize) {
     double dn = d->norm();
     d->scale(1.0/dn);
 }
}

psi::SharedVector oepdev::GramSchmidt::projection(psi::SharedVector u, psi::SharedVector v) const {
  psi::SharedVector p = std::make_shared<psi::Vector>(*u);
  double un = u->norm();
  double d  = v->vector_dot(u) / (un * un);
  p->scale(d);
  return p;
}
