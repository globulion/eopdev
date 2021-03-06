#include "diis.h"

namespace std{

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


}

namespace oepdev{


DIISManager::DIISManager (int dim, int na, int nb) :
          _dim(dim),
          _dn(dim+1), 
          _na(na),
          _nb(nb), 
          _ndiis(0), 
          _u(std::make_shared<Matrix>("DIIS Updated Vector", na, nb))
{
  outfile->Printf("\n WARNING@oepdev_libutil::diis::DIISManager: Using Intrinsic DIIS Manager\n\n");
}

void DIISManager::update(std::shared_ptr<Matrix>& other) {
  if (_ndiis == _dim) other->copy(this->_u);
}

DIISManager::~DIISManager() 
 { /* nothing to explicitly do because all memorials are std containers or smart pointers */ }

void DIISManager::put(const std::shared_ptr<const Matrix>& error, 
                  const std::shared_ptr<const Matrix>& vector) {
   _vectors.push_back(std::shared_ptr<Matrix>(new Matrix(*vector)));
   _errors .push_back(std::shared_ptr<Matrix>(new Matrix(*error )));
   if (_ndiis == _dim) {
       _vectors.erase(_vectors.begin());
       _errors .erase(_errors .begin());
   }
   if (_ndiis < _dim) _ndiis++;
}

void DIISManager::compute(void) {
  if (_ndiis == _dim) {
   std::shared_ptr<Matrix> B(new Matrix("DIIS Covariance Matrix"                 , _dn, _dn));
   std::shared_ptr<Vector> c(new Vector("DIIS Coefficients + Lagrange Multiplier", _dn     ));
   double** Bp = B->pointer();
   double*  cp = c->pointer();
   int*     ipiv = init_int_array(_dn);

   // construct B matrix and f vector
   for (int i=0; i<_dim; i++) {
        Bp[i][_dim] = Bp[_dim][i] = 1.0;
        cp[i]        = 0.0;
        for (int j=0; j<_dim; j++) {
             Bp[i][j] = _errors[i]->vector_dot(_errors[j]);
        }
   }
   Bp[_dim][_dim] = 0.0;
   cp[_dim]       = 1.0;

   // solve for c vector
   C_DGESV(_dn, 1, Bp[0], _dn, ipiv, cp, _dn);

   // update vector
   _u->zero();
   double** up = _u->pointer();
   for (int ina=0; ina<_na; ina++) {
        for (int inb=0; inb<_nb; inb++) { 
             for (int k=0; k<_dim; k++) {
                  up[ina][inb] += cp[k] * _vectors[k]->get(ina,inb);
             }
        }
   }

   // free memory
   free(ipiv);
   }
}

} // EndNameSpace oepdev
