#include <iostream>
#include <random>
#include "gefp.h"

using namespace std;

oepdev::GenEffFrag::GenEffFrag(std::string name) : 
  name_(name), frag_(nullptr),
  densityMatrixSusceptibilityGEF_(std::make_shared<oepdev::GenEffPar>("POLARIZATION")) // deprecate this!!!/TODO
//electrostaticEnergyGEF_(nullptr),
//repulsionEnergyGEF_(nullptr),
//chargeTransferEnergyGEF_(nullptr),
//EETCouplingConstantGEF_(nullptr)
{
// Deprecate all below:
//parameters["POLARIZATION"   ] = densityMatrixSusceptibilityGEF_;
//parameters["COULOMBIC   "   ] = electrostaticEnergyGEF_;
//parameters["REPULSION"      ] = repulsionEnergyGEF_;
//parameters["CHARGE_TRANSFER"] = chargeTransferEnergyGEF_;
//parameters["EET_COUPLING"   ] = EETCouplingConstantGEF_;
}
oepdev::GenEffFrag::GenEffFrag(const GenEffFrag* f) {
 name_ = f->name_;
 this->set_molecule(f->frag_);
 for (auto const& x : f->parameters) {
       std::string key = x.first;
       std::shared_ptr<oepdev::GenEffPar> par = x.second->clone();
       this->parameters[key] = par;
  }
  // do not copy basis sets! They must be set externally by basissets["BASIS_KEY"] = some_basis;

 densityMatrixSusceptibilityGEF_ = f->densityMatrixSusceptibilityGEF_;//TODO ->copy it, not pass the pointer
}
oepdev::GenEffFrag::~GenEffFrag() {}
void oepdev::GenEffFrag::set_basisset(std::string key, psi::SharedBasisSet basis) {
 this->basissets[key] = basis;
   for (auto const& x : this->parameters) {
       x.second->set_basisset(key, basis);
  }
}
void oepdev::GenEffFrag::rotate(std::shared_ptr<psi::Matrix> R)
{
  for (auto const& x : this->parameters) {
       x.second->rotate(R);
  }
}
void oepdev::GenEffFrag::translate(std::shared_ptr<psi::Vector> T)
{
  for (auto const& x : this->parameters) {
       x.second->translate(T);
  }
}
void oepdev::GenEffFrag::superimpose(std::shared_ptr<psi::Matrix> targetXYZ, std::vector<int> supList)
{
  for (auto const& x : this->parameters) {
       outfile->Printf(" Superimposing now %s in Fragment %s\n", x.first.c_str(), this->name_.c_str());
       x.second->superimpose(targetXYZ, supList);
  }
}
void oepdev::GenEffFrag::superimpose(psi::SharedMolecule mol, std::vector<int> supList)
{
 psi::SharedMatrix xyz = std::make_shared<psi::Matrix>(mol->geometry());
 this->superimpose(xyz, supList);
}
void oepdev::GenEffFrag::superimpose(void)
{
 psi::SharedMatrix xyz = std::make_shared<psi::Matrix>(this->frag_->geometry());
 this->superimpose(xyz, {});
}
// 
double oepdev::GenEffFrag::energy(std::string theory, std::shared_ptr<GenEffFrag> other) {
 double e_tot = this->compute_pairwise_energy(theory, other);
 return e_tot;
}

double oepdev::GenEffFrag::compute_energy(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments, bool manybody)
{
  double e_tot = 0.0;
  if (manybody) {e_tot += oepdev::GenEffFrag::compute_many_body_energy(theory, fragments);}
  else {

     const int n = fragments.size();                                            
     for (int i=0; i<n; ++i) {
     for (int j=0; j<i; ++j) {
          e_tot += fragments[i]->compute_pairwise_energy(theory, fragments[j]);
     }}
  }
  return e_tot;
}
double oepdev::GenEffFrag::compute_many_body_energy(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments)
{
  double energy;

  if (theory == "EFP2:IND" || theory == "OEPa-EFP2:IND" || theory == "OEPb-EFP2:IND") {
       // Initialize
       oepdev::MultipoleConvergence::ConvergenceLevel clevel = oepdev::DMTPole::determine_dmtp_convergence_level("DMTP_CONVER");
           
       const int n_frag = fragments.size();
      
       std::vector<std::vector<psi::SharedMatrix>> dpols; 
       std::vector<psi::SharedMatrix> rpols; 
       std::vector<std::shared_ptr<oepdev::DMTPole>> dmtps; 
       int nocc_sum = 0;
       for (int n=0; n<n_frag; ++n) {
            dpols.push_back(fragments[n]->parameters["efp2"]->dpol("0"));
            rpols.push_back(fragments[n]->parameters["efp2"]->matrix("lmoc"));
            dmtps.push_back(fragments[n]->parameters["efp2"]->dmtp("camm"));
            nocc_sum += dpols[n].size();
       }

       const int DIM = 3*nocc_sum;
      
       // Calculate D and F matrices
       psi::SharedMatrix D = std::make_shared<psi::Matrix>("D-matrix", DIM, DIM);
       psi::SharedMatrix F = std::make_shared<psi::Matrix>("F-matrix", DIM, 1  );
       double** Dp = D->pointer();
       double** Fp = F->pointer();
      
       int i_pol = 0;
       for (int n=0; n<n_frag; ++n) {
            const int N = dpols[n].size();
            
            i_pol += N;
            for (int a=0; a<N; ++a) {
                 int ix3 = 3*(i_pol - N) + 3*a;
      
                 psi::SharedMatrix A = dpols[n][a]->clone();
                 A->invert();
                 double** Ap = A->pointer();
      
                 Dp[ix3+0][ix3+0] = Ap[0][0];
                 Dp[ix3+0][ix3+1] = Ap[0][1];
                 Dp[ix3+0][ix3+2] = Ap[0][2];
                 Dp[ix3+1][ix3+0] = Ap[1][0];
                 Dp[ix3+1][ix3+1] = Ap[1][1];
                 Dp[ix3+1][ix3+2] = Ap[1][2];
                 Dp[ix3+2][ix3+0] = Ap[2][0];
                 Dp[ix3+2][ix3+1] = Ap[2][1];
                 Dp[ix3+2][ix3+2] = Ap[2][2];
      
                 psi::SharedVector ra = rpols[n]->get_row(0, a);
                 double fx = 0.0; double x = ra->get(0);
                 double fy = 0.0; double y = ra->get(1);
                 double fz = 0.0; double z = ra->get(2);
      
                 int j_pol = 0;
                 for (int m=0; m<n_frag; ++m) {
                      const int M = dpols[m].size();
                      j_pol += M;
      
                      // Consider only other fragments than the n-th one
                      if (n!=m) {
      
                          // Accumulate electric field
                          std::shared_ptr<oepdev::MultipoleConvergence> efield = dmtps[m]->field(x, y, z, clevel);
                          fx += efield->level(clevel)->get(0, 0);
                          fy += efield->level(clevel)->get(0, 1);
                          fz += efield->level(clevel)->get(0, 2);
                       
                          for (int b=0; b<M; ++b) {
                               int jx3 = 3*(j_pol - M) + 3*b;
      
                               psi::SharedVector rb = rpols[m]->get_row(0, b);
      
                               // compute Tab
                               double rab_x = ra->get(0) - rb->get(0);
                               double rab_y = ra->get(1) - rb->get(1);
                               double rab_z = ra->get(2) - rb->get(2);
                               double r = sqrt(rab_x*rab_x+rab_y*rab_y+rab_z*rab_z);
                               double r3 = 1.0/(r*r*r);
                               double r5 = r3/(r*r);
      
                               psi::SharedMatrix Tab = std::make_shared<psi::Matrix>("",3,3); Tab->identity();
                               Tab->scale(-1.0*r3);
                               double** t = Tab->pointer();
                               t[0][0] += 3.0 * r5 * rab_x * rab_x;
                               t[0][1] += 3.0 * r5 * rab_x * rab_y;
                               t[0][2] += 3.0 * r5 * rab_x * rab_z;
                               t[1][1] += 3.0 * r5 * rab_y * rab_y;
                               t[1][2] += 3.0 * r5 * rab_y * rab_z;
                               t[2][2] += 3.0 * r5 * rab_z * rab_z;
                               t[1][0]  = t[0][1];
                               t[2][0]  = t[0][2];
                               t[2][1]  = t[1][2];
      
                               // Set off-diagonals of D with -Tab
                               Dp[ix3+0][jx3+0] = -t[0][0];
                               Dp[ix3+0][jx3+1] = -t[0][1];
                               Dp[ix3+0][jx3+2] = -t[0][2];
                               Dp[ix3+1][jx3+0] = -t[1][0];
                               Dp[ix3+1][jx3+1] = -t[1][1];
                               Dp[ix3+1][jx3+2] = -t[1][2];
                               Dp[ix3+2][jx3+0] = -t[2][0];
                               Dp[ix3+2][jx3+1] = -t[2][1];
                               Dp[ix3+2][jx3+2] = -t[2][2];
                          }
                      }
                 }
      
                 Fp[ix3+0][0] = fx;
                 Fp[ix3+1][0] = fy;
                 Fp[ix3+2][0] = fz;
            }
       }
      
       if (psi::Process::environment.options.get_int("PRINT") > 2) F->print();
       // Compute interaction energy
       D->invert();
       if (psi::Process::environment.options.get_int("PRINT") > 2) D->print();
       psi::SharedMatrix P = psi::Matrix::doublet(D, F, false, false);

       energy = -P->vector_dot(F) / 2.0;

  } else {
     throw psi::PSIEXCEPTION("Wrong manybody theory chosen.");
  }
  return energy;
}
//double oepdev::GenEffFrag::compute_pairwise_total_energy(std::vector<std::shared_ptr<GenEffFrag>> fragments,
//          double (oepdev::GenEffFrag::*pairwise_energy_computer)(void )) {
//  const int n = fragments.size();
//  double e_tot = 0.0;
//  for (int i=0; i<n; ++i) {
//  for (int j=0; j<i; ++j) {
//       e_tot += fragments[i]->
//       e_tot += pairwise_energy_computer(fragments[i], fragments[j]);
//  }
//  }
//  return e_tot;
//}
double oepdev::GenEffFrag::compute_pairwise_energy(std::string theory, std::shared_ptr<GenEffFrag> other) {
  double e;
 if (theory == "EFP2:COUL" || theory == "OEPa-EFP2:COUL" || theory == "OEPb-EFP2:COUL") 
     {e = this->compute_pairwise_energy_efp2_coul(other);}
 else if (theory == "EFP2:IND" || theory == "OEPa-EFP2:IND" || theory == "OEPb-EFP2:IND") 
     {e = this->compute_pairwise_energy_efp2_ind(other);}
 else if (theory == "EFP2:EXREP" || theory == "OEPa-EFP2:EXREP") 
     {e = this->compute_pairwise_energy_efp2_exrep(other);}
 else if (theory == "EFP2:CT")
     {e = this->compute_pairwise_energy_efp2_ct(other);}
 else if (theory == "EFP2:DISP" || theory == "OEPa-EFP2:DISP" || theory == "OEPb-EFP2:DISP")
     {e = this->compute_pairwise_energy_efp2_disp(other);}
 else if (theory == "OEPb-EFP2:EXREP") 
     {e = this->compute_pairwise_energy_oep_efp2_exrep(other);}
 else if (theory == "OEPa-EFP2:CT" || theory == "OEPb-EFP2:CT") 
     {e = this->compute_pairwise_energy_oep_efp2_ct(other);}
 else {
     throw psi::PSIEXCEPTION("Wrong pairwise theory chosen.");
 }
 return e;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_coul(std::shared_ptr<GenEffFrag> other) {
 oepdev::MultipoleConvergence::ConvergenceLevel clevel = oepdev::DMTPole::determine_dmtp_convergence_level("DMTP_CONVER");
 oepdev::SharedDMTPole camm_1 = this->parameters["efp2"]->dmtp("camm");
 oepdev::SharedDMTPole camm_2 =other->parameters["efp2"]->dmtp("camm");
 double e = camm_1->energy(camm_2)->level(clevel)->get(0,0);
 return e;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ind(std::shared_ptr<GenEffFrag> other) {
 std::vector<std::shared_ptr<oepdev::GenEffFrag>> fragments;
 fragments.push_back(shared_from_this());
 fragments.push_back(other);
 double e_ind = oepdev::GenEffFrag::compute_many_body_energy("EFP2:IND", fragments);
 return e_ind;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_disp(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_exrep(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_ct(std::shared_ptr<GenEffFrag> other) {
 //TODO
 return 0.0;
}
