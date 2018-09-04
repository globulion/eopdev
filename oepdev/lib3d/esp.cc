#include "esp.h"


namespace oepdev{

using namespace std;

ESPSolver::ESPSolver(SharedScalarField3D field)
 : field_(field), nCentres_(field->wfn()->molecule()->natom()),
   centres_(std::make_shared<psi::Matrix>(field->wfn()->molecule()->geometry()))
{
   common_init();
}

ESPSolver::ESPSolver(SharedScalarField3D field, psi::SharedMatrix centres)
 : field_(field), nCentres_(centres->nrow()),
   centres_(std::make_shared<psi::Matrix>(*centres))
{
   common_init();
}

ESPSolver::~ESPSolver()
{

}

void ESPSolver::common_init()
{
  charges_ = std::make_shared<psi::Vector>("ESP Fitted Charges", nCentres_               );
  bvec__   = std::make_shared<psi::Vector>("ESP b vector"      , nCentres_+1             );
  amat__   = std::make_shared<psi::Matrix>("ESP A matrix"      , nCentres_+1, nCentres_+1);
}

void ESPSolver::compute(void) 
{
   compute_matrices();
   fit();
}

void ESPSolver::compute_matrices(void)
{
    if (!field_->is_computed())
        throw psi::PSIEXCEPTION("OEPDEV Error. Scalar field does not have computed potential data yet!\n");

    int npot = field_->npoints();
    double** fdp = field_->data()->pointer();
    double** ap  = amat__->pointer();
    double*  bp  = bvec__->pointer();
    double** cp  = centres_->pointer();

    double r_mx, r_nx;
    double fdpx, fdpy, fdpz, fdpv; 
    double cpmx, cpmy, cpmz;

    // ===> Compute ESP A matrix and b vector <=== //
    for (int i=0; i < npot; ++i) {
         fdpx = fdp[i][0];
         fdpy = fdp[i][1];
         fdpz = fdp[i][2];
         fdpv = fdp[i][3];
         for (int m=0; m < nCentres_; ++m) {
              cpmx = cp[m][0];
              cpmy = cp[m][1];
              cpmz = cp[m][2];
              r_mx = sqrt( pow(fdpx - cpmx, 2.0) 
                         + pow(fdpy - cpmy, 2.0) 
                         + pow(fdpz - cpmz, 2.0));
              bp[m] = bp[m] + fdpv / r_mx;
              for (int n=0; n <= m; ++n) {
                   r_nx = sqrt( pow(fdpx - cp[n][0], 2.0) 
                              + pow(fdpy - cp[n][1], 2.0) 
                              + pow(fdpz - cp[n][2], 2.0));
                   ap[m][n] = ap[m][n] + 1.0/(r_mx * r_nx);
              }
         }
    }

    for (int m=0; m < nCentres_; ++m) {
         ap[nCentres_][m] = ap[m][nCentres_] = 1.0;
         for (int n=0; n <= m; ++n) {
             ap[n][m] = ap[m][n];
         }
    }
    bp[nCentres_]            = 0.0;
    ap[nCentres_][nCentres_] = 0.0;
}

void ESPSolver::fit(void)
{
    amat__->invert();

    double** ap = amat__->pointer();
    double*  bp = bvec__->pointer();

    double v;
    for (int m=0; m < nCentres_+1; ++m) {
         v = 0.0;
         for (int n=0; n < nCentres_+1; ++n) {
              v += ap[m][n] * bp[n];
         }
         if (m < nCentres_) charges_->set(m, v);
    }
}


} // EndNameSpace oepdev
