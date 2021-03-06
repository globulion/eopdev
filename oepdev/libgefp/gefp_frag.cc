#include <iostream>
#include <random>
#include "gefp.h"

namespace oepdev {
}

using namespace std;

oepdev::GenEffFrag::GenEffFrag(std::string name) : 
  name_(name), frag_(nullptr), nbf_(0), natom_(0), ndocc_(0),
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
 nbf_  = f->nbf_;
 natom_= f->natom_;
 ndocc_= f->ndocc_;
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
  outfile->Printf(" Superimposing finished\n");
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
double oepdev::GenEffFrag::energy_term(std::string theory, std::shared_ptr<GenEffFrag> other) const {
 double e_tot = this->compute_pairwise_energy(theory, other);
 return e_tot;
}

double oepdev::GenEffFrag::compute_energy(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments)
{
  double e_tot = 0.0;
  /*
   * EFP2 Model.
   */
  if (theory == "EFP2") {
      double e_cou = 0.0;
      double e_exr = 0.0;
      double e_ind = oepdev::GenEffFrag::compute_many_body_energy_term("EFP2:IND", fragments);
      double e_dis = 0.0;
      double e_ct  = 0.0;

      const int n = fragments.size();                                            

      for (int i=0; i<n; ++i) {
      for (int j=0; j<i; ++j) {
           e_cou += fragments[i]->compute_pairwise_energy("EFP2:COUL" , fragments[j]);
           e_exr += fragments[i]->compute_pairwise_energy("EFP2:EXREP", fragments[j]);
           e_dis += fragments[i]->compute_pairwise_energy("EFP2:DISP" , fragments[j]);
           e_ct  += fragments[i]->compute_pairwise_energy("EFP2:CT"   , fragments[j]);
      }}
    //cout << e_cou << " " << e_exr << " " << e_ind << " " << e_ct << endl;
      e_tot = e_cou + e_exr + e_ind + e_dis + e_ct;
  } 
  /*
   * OEPa-EFP2 Model. 
   *
   * EFP2 with CT term replaced with OEP-CT term 
   */
  else if (theory == "OEPa-EFP2") {
      double e_cou = 0.0;
      double e_exr = 0.0;
      double e_ind = oepdev::GenEffFrag::compute_many_body_energy_term("EFP2:IND", fragments);
      double e_dis = 0.0;
      double e_ct  = 0.0;

      const int n = fragments.size();                                            

      for (int i=0; i<n; ++i) {
      for (int j=0; j<i; ++j) {
           e_cou += fragments[i]->compute_pairwise_energy(     "EFP2:COUL" , fragments[j]);
           e_exr += fragments[i]->compute_pairwise_energy(     "EFP2:EXREP", fragments[j]);
           e_dis += fragments[i]->compute_pairwise_energy(     "EFP2:DISP" , fragments[j]);
           e_ct  += fragments[i]->compute_pairwise_energy("OEPb-EFP2:CT"   , fragments[j]);
      }}
      e_tot = e_cou + e_exr + e_ind + e_dis + e_ct;
  }
  /*
   * OEPb-EFP2 Model. 
   *
   * EFP2 with CT term replaced with OEP-CT(Otto-Ladik) term 
   * as well as REP term replaced with OEP-REP(Murrell etal) term
   */
  else if (theory == "OEPb-EFP2") {
      double e_cou = 0.0;
      double e_exr = 0.0;
      double e_ind = oepdev::GenEffFrag::compute_many_body_energy_term("EFP2:IND", fragments);
      double e_dis = 0.0;
      double e_ct  = 0.0;

      const int n = fragments.size();                                            

      for (int i=0; i<n; ++i) {
      for (int j=0; j<i; ++j) {
           e_cou += fragments[i]->compute_pairwise_energy(     "EFP2:COUL" , fragments[j]);
           e_exr += fragments[i]->compute_pairwise_energy("OEPb-EFP2:EXREP", fragments[j]);
           e_dis += fragments[i]->compute_pairwise_energy(     "EFP2:DISP" , fragments[j]);
           e_ct  += fragments[i]->compute_pairwise_energy("OEPb-EFP2:CT"   , fragments[j]);
      }}
    //cout << e_cou << " " << e_exr << " " << e_ind << " " << e_ct << endl;
      e_tot = e_cou + e_exr + e_ind + e_dis + e_ct;
  } else {
    throw psi::PSIEXCEPTION("Unknown theory chosen for GenEffFrag::compute_energy!\n");
  }
  return e_tot;
}

double oepdev::GenEffFrag::compute_energy_term(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments, bool manybody)
{
  double e_term = 0.0;
  if (manybody) {e_term += oepdev::GenEffFrag::compute_many_body_energy_term(theory, fragments);}
  else {

     const int n = fragments.size();                                            
     for (int i=0; i<n; ++i) {
     for (int j=0; j<i; ++j) {
          e_term += fragments[i]->compute_pairwise_energy(theory, fragments[j]);
     }}
  }
  return e_term;
}
double oepdev::GenEffFrag::compute_many_body_energy_term(std::string theory, std::vector<std::shared_ptr<GenEffFrag>> fragments)
{
  double energy;

  // ---> Timer-on <--- //
  clock_t t_time = -clock();

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

      // ---> Timer-off <--- //
      t_time += clock(); // Clock END
      cout << " o TIME EFP2:IND   " << ((double)t_time/CLOCKS_PER_SEC) << endl;

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
double oepdev::GenEffFrag::compute_pairwise_energy(std::string theory, std::shared_ptr<GenEffFrag> other) const {
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
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_coul(std::shared_ptr<GenEffFrag> other) const {
 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 oepdev::MultipoleConvergence::ConvergenceLevel clevel = oepdev::DMTPole::determine_dmtp_convergence_level("DMTP_CONVER");
 oepdev::SharedDMTPole camm_1 = this->parameters.at("efp2")->dmtp("camm");
 oepdev::SharedDMTPole camm_2 =other->parameters.at("efp2")->dmtp("camm");
 double e = camm_1->energy(camm_2)->level(clevel)->get(0,0);

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME EFP2:COUL  " << ((double)t_time/CLOCKS_PER_SEC) << endl;

 /* 
   Coulomb penetration is not included in timing because overlap integrals
   and similar loops can be performed in the EFP2-EXREP routine, which is
   necessary in EFP2 calculations. Adding the below calculation does not
   change the timing noticeably and is neglected here.
  */
 if (psi::Process::environment.options.get_bool("COULOMB_PENETRATION")) {

    const int nbf_1 = this->nbf();
    const int nbf_2 =other->nbf();
    const int nocc_1 = this->ndocc();
    const int nocc_2 =other->ndocc();

    psi::SharedBasisSet primary_1 = this->parameters.at("efp2")->basisset("primary");
    psi::SharedBasisSet primary_2 =other->parameters.at("efp2")->basisset("primary");
    psi::SharedMatrix cmoo_1 =  this->parameters.at("efp2")->matrix("lmoo");
    psi::SharedMatrix cmoo_2 = other->parameters.at("efp2")->matrix("lmoo");
    psi::SharedMatrix lmoc_1 =  this->parameters.at("efp2")->matrix("lmoc");
    psi::SharedMatrix lmoc_2 = other->parameters.at("efp2")->matrix("lmoc");

    psi::SharedMatrix Sao12 = std::make_shared<psi::Matrix>("Sao(1,2)", nbf_1, nbf_2); 
    psi::IntegralFactory fact_12(primary_1, primary_2, primary_1, primary_2);
    std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_12.ao_overlap());
    ovlInt->compute(Sao12);
    psi::SharedMatrix Smo12 = psi::Matrix::triplet(cmoo_1, Sao12 , cmoo_2, true, false, false);

    double** S = Smo12->pointer();
    for (int i=0; i<nocc_1; ++i) {
    for (int j=0; j<nocc_2; ++j) {
         double s = S[i][j];
         double t = 2.0 * log(abs(s));
         double r = sqrt(pow(lmoc_1->get(i,0)-lmoc_2->get(j,0), 2.0) + 
                         pow(lmoc_1->get(i,1)-lmoc_2->get(j,1), 2.0) + 
                         pow(lmoc_1->get(i,2)-lmoc_2->get(j,2), 2.0) );
         e -= 2.0 * sqrt(-1.0 / t) * s*s / r ;
    }
    }
       
 }

 return e;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_exrep(std::shared_ptr<GenEffFrag> other) const {

 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 // Initialize results
 double e_s1 = 0.0, e_s2 = 0.0, e_ex = 0.0, e_tot = 0.0;

 // Extract EFP2 parameters from fragments
 const int nbf_1 = this->nbf();
 const int nbf_2 =other->nbf();
 const int nat_1 = this->natom();
 const int nat_2 =other->natom();
 const int ndocc_1= this->ndocc();
 const int ndocc_2=other->ndocc();
 psi::SharedMolecule mol_1 = this->molecule();
 psi::SharedMolecule mol_2 =other->molecule();
 psi::SharedBasisSet primary_1 = this->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet primary_2 =other->parameters.at("efp2")->basisset("primary");
 psi::SharedMatrix cmoo_1 =  this->parameters.at("efp2")->matrix("lmoo");
 psi::SharedMatrix cmoo_2 = other->parameters.at("efp2")->matrix("lmoo");
 psi::SharedMatrix fock_1 =  this->parameters.at("efp2")->matrix("fock_lmo");
 psi::SharedMatrix fock_2 = other->parameters.at("efp2")->matrix("fock_lmo");
 psi::SharedMatrix lmoc_1 =  this->parameters.at("efp2")->matrix("lmoc");
 psi::SharedMatrix lmoc_2 = other->parameters.at("efp2")->matrix("lmoc");

 // Compute one-electron integrals
 // Cl+
 std::shared_ptr<psi::Matrix> Sao12 = std::make_shared<psi::Matrix>("Sao(1,2)", nbf_1, nbf_2);
 std::shared_ptr<psi::Matrix> Tao12 = std::make_shared<psi::Matrix>("Tao(1,2)", nbf_1, nbf_2);

 t_time += clock();
 psi::IntegralFactory fact_12(primary_1, primary_2, primary_1, primary_2);
 t_time -= clock();

 std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_12.ao_overlap());
 std::shared_ptr<psi::OneBodyAOInt> kinInt(fact_12.ao_kinetic());
 ovlInt->compute(Sao12);
 kinInt->compute(Tao12);
 psi::SharedMatrix Smo12 = psi::Matrix::triplet(cmoo_1, Sao12 , cmoo_2, true, false, false);
 psi::SharedMatrix Tmo12 = psi::Matrix::triplet(cmoo_1, Tao12 , cmoo_2, true, false, false);
 psi::SharedMatrix SF1S  = psi::Matrix::triplet(Smo12, fock_1, Smo12, true, false, false);
 psi::SharedMatrix SF2S  = psi::Matrix::triplet(Smo12, fock_2, Smo12, false, false, true);

 e_s1 = SF1S->trace() + SF2S->trace() - 2.0 * Smo12->vector_dot(Tmo12);
 SF1S.reset(); SF2S.reset();

 // S2
 double** S = Smo12->pointer();
 psi::SharedVector ss1 = std::make_shared<psi::Vector>("s^2 (1)", ndocc_1);
 psi::SharedVector ss2 = std::make_shared<psi::Vector>("s^2 (2)", ndocc_2);
 double* s1 = ss1->pointer();
 double* s2 = ss2->pointer();
 double val, sab;
 for (int a = 0; a<ndocc_1; ++a) {
      val = 0.0;
      for (int b = 0; b<ndocc_2; ++b) {
           sab = S[a][b];
           val += sab * sab;
      }
      s1[a] = val;
 }
 for (int b = 0; b<ndocc_2; ++b) {
      val = 0.0;
      for (int a = 0; a<ndocc_1; ++a) {
           sab = S[a][b];
           val += sab * sab;
      }
      s2[b] = val;
 }

 psi::SharedVector ZR1 = std::make_shared<psi::Vector>("ZR1", ndocc_1);
 psi::SharedVector ZR2 = std::make_shared<psi::Vector>("ZR2", ndocc_2);
 psi::SharedVector IR1 = std::make_shared<psi::Vector>("IR1", ndocc_1);
 psi::SharedVector IR2 = std::make_shared<psi::Vector>("IR2", ndocc_2);

 double* pZR1 = ZR1->pointer();
 double* pZR2 = ZR2->pointer();
 double* pIR1 = IR1->pointer();
 double* pIR2 = IR2->pointer();
 
 double rbx, ray, rad, rbc, rab;

 for (int b=0; b<ndocc_2; ++b){
      val = 0.0;
      for (int x=0; x<nat_1; ++x) { 
           rbx  = sqrt(pow(mol_1->x(x)-lmoc_2->get(b,0), 2.0) + 
                       pow(mol_1->y(x)-lmoc_2->get(b,1), 2.0) + 
                       pow(mol_1->z(x)-lmoc_2->get(b,2), 2.0) );
           val -= (double)mol_1->Z(x) / rbx;
      }
      pZR2[b] = val;
 }
 for (int a=0; a<ndocc_1; ++a){
      val = 0.0;
      for (int y=0; y<nat_2; ++y) { 
           ray  = sqrt(pow(mol_2->x(y)-lmoc_1->get(a,0), 2.0) + 
                       pow(mol_2->y(y)-lmoc_1->get(a,1), 2.0) + 
                       pow(mol_2->z(y)-lmoc_1->get(a,2), 2.0) );
           val -= (double)mol_2->Z(y) / ray;
      }
      pZR1[a] = val;
 }
 for (int a=0; a<ndocc_1; ++a) {
      val = 0.0;
      for (int b=0; b<ndocc_2; ++b) {
           rad  = sqrt(pow(lmoc_1->get(a,0)-lmoc_2->get(b,0), 2.0) + 
                       pow(lmoc_1->get(a,1)-lmoc_2->get(b,1), 2.0) + 
                       pow(lmoc_1->get(a,2)-lmoc_2->get(b,2), 2.0) );
           val += 2.0 / rad;
      }
      pIR1[a] = val;
 }
 for (int b=0; b<ndocc_2; ++b) {
      val = 0.0;
      for (int a=0; a<ndocc_1; ++a) {
           rbc  = sqrt(pow(lmoc_1->get(a,0)-lmoc_2->get(b,0), 2.0) + 
                       pow(lmoc_1->get(a,1)-lmoc_2->get(b,1), 2.0) + 
                       pow(lmoc_1->get(a,2)-lmoc_2->get(b,2), 2.0) );
           val += 2.0 / rbc;
      }
      pIR2[b] = val;
 }
 ZR1->add(IR1);
 ZR2->add(IR2);

 e_s2 = ss1->vector_dot(ZR1) + ss2->vector_dot(ZR2);

 for (int a=0; a<ndocc_1; ++a) {
      for (int b=0; b<ndocc_2; ++b) { 
           sab = S[a][b];
           rab = sqrt(pow(lmoc_1->get(a,0) - lmoc_2->get(b,0), 2.0) +
                      pow(lmoc_1->get(a,1) - lmoc_2->get(b,1), 2.0) +
                      pow(lmoc_1->get(a,2) - lmoc_2->get(b,2), 2.0) );
           e_s2 -= sab * sab / rab;
      }
 }

 // ---> Finalize with repulsion <--- //
 e_s1 *=-2.0;
 e_s2 *= 2.0;

 // EXCH
 for (int a=0; a<ndocc_1; ++a) {
      for (int b=0; b<ndocc_2; ++b) {
           sab = S[a][b];
           rab = sqrt(pow(lmoc_1->get(a,0) - lmoc_2->get(b,0), 2.0) + 
                      pow(lmoc_1->get(a,1) - lmoc_2->get(b,1), 2.0) + 
                      pow(lmoc_1->get(a,2) - lmoc_2->get(b,2), 2.0) );
           e_ex += sqrt(-2.0*log(abs(sab)) / M_PI) * sab * sab / rab;
      }
 }
 e_ex *= -4.0;
 // Cl-

 // Save
 double v_e_rep = psi::Process::environment.globals["EINT REP EFP2 KCAL"];   
 double v_e_exc = psi::Process::environment.globals["EINT EXC EFP2 KCAL"];   
 double v_e_exr = psi::Process::environment.globals["EINT EXR EFP2 KCAL"];   
 double v_e_rs1 = psi::Process::environment.globals["EINT REP EFP2:S1 KCAL"];
 double v_e_rs2 = psi::Process::environment.globals["EINT REP EFP2:S2 KCAL"];
 v_e_rep +=(e_s1 + e_s2        )*OEPDEV_AU_KcalPerMole;
 v_e_exc +=(e_ex               )*OEPDEV_AU_KcalPerMole;
 v_e_exr +=(e_s1 + e_s2 + e_ex )*OEPDEV_AU_KcalPerMole;
 v_e_rs1 +=(e_s1               )*OEPDEV_AU_KcalPerMole;
 v_e_rs2 +=(e_s2               )*OEPDEV_AU_KcalPerMole;

 psi::Process::environment.globals["EINT REP EFP2 KCAL"]    = v_e_rep;
 psi::Process::environment.globals["EINT EXC EFP2 KCAL"]    = v_e_exc;
 psi::Process::environment.globals["EINT EXR EFP2 KCAL"]    = v_e_exr;
 psi::Process::environment.globals["EINT REP EFP2:S1 KCAL"] = v_e_rs1;
 psi::Process::environment.globals["EINT REP EFP2:S2 KCAL"] = v_e_rs2;
 
 e_tot = e_s1 + e_s2 + e_ex;
//cout << e_s1 << " " << e_s2 << " " << e_ex << endl;

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME EFP2:EXREP " << ((double)t_time/CLOCKS_PER_SEC) << endl;

 return e_tot;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ind(std::shared_ptr<GenEffFrag> other) const {

 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 std::vector<std::shared_ptr<oepdev::GenEffFrag>> fragments;
 std::shared_ptr<oepdev::GenEffFrag> t = shared_from_this()->clone();
 fragments.push_back(t);
 fragments.push_back(other);
 double e_ind = oepdev::GenEffFrag::compute_many_body_energy_term("EFP2:IND", fragments);

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME EFP2:IND(P)" << ((double)t_time/CLOCKS_PER_SEC) << endl;

 return e_ind;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_ct(std::shared_ptr<GenEffFrag> other) const {

 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 // Initialize result
 double e_ct  = 0.0;

 // Extract EFP2 parameters from fragments
 const int nbf_1 = this->nbf();
 const int nbf_2 =other->nbf();
 const int ndocc_1= this->ndocc();
 const int ndocc_2=other->ndocc();
 psi::SharedMolecule mol_1 = this->molecule();
 psi::SharedMolecule mol_2 =other->molecule();
 psi::SharedBasisSet primary_1 = this->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet primary_2 =other->parameters.at("efp2")->basisset("primary");

 // V matrices //
 psi::SharedMatrix VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)", nbf_1, nbf_2);
 psi::SharedMatrix VaoB11    = std::make_shared<psi::Matrix>("VaoB(1,1)", nbf_1, nbf_1);
 psi::SharedMatrix VaoA21    = std::make_shared<psi::Matrix>("VaoA(2,1)", nbf_2, nbf_1);
 psi::SharedMatrix VaoA22    = std::make_shared<psi::Matrix>("VaoA(2,2)", nbf_2, nbf_2);

 // S matrices //
 psi::SharedMatrix Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)", nbf_1, nbf_2);

 // T matrices //
 psi::SharedMatrix Tao12     = std::make_shared<psi::Matrix>("Tao(1,2)", nbf_1, nbf_2);
 psi::SharedMatrix Tao11     = std::make_shared<psi::Matrix>("Tao(1,1)", nbf_1, nbf_1);
 psi::SharedMatrix Tao22     = std::make_shared<psi::Matrix>("Tao(2,2)", nbf_2, nbf_2);

 // F matrices //
 psi::SharedMatrix Fao11     = this->parameters.at("efp2")->matrix("fock_ao");
 psi::SharedMatrix Fao22     =other->parameters.at("efp2")->matrix("fock_ao");

 // Ca matrices //
 psi::SharedMatrix cmoo_1 =  this->parameters.at("efp2")->matrix("cmoo"); 
 psi::SharedMatrix cmoo_2 = other->parameters.at("efp2")->matrix("cmoo");
 psi::SharedMatrix cmov_1 =  this->parameters.at("efp2")->matrix("cmov");
 psi::SharedMatrix cmov_2 = other->parameters.at("efp2")->matrix("cmov");
 const int nvir_1 = cmov_1->ncol();
 const int nvir_2 = cmov_2->ncol();

 // IntegralFactory //
 t_time += clock();
 oepdev::IntegralFactory fact_12(primary_1, primary_2, primary_1, primary_2);
 oepdev::IntegralFactory fact_21(primary_2, primary_1, primary_2, primary_1);
 oepdev::IntegralFactory fact_11(primary_1, primary_1, primary_1, primary_1);
 oepdev::IntegralFactory fact_22(primary_2, primary_2, primary_2, primary_2);
 t_time -= clock();

 // Overlap integrals //
 std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_12.ao_overlap());
 ovlInt->compute(Sao12);

 // Kinetic energy integrals //
 std::shared_ptr<psi::OneBodyAOInt> kinInt12(fact_12.ao_kinetic());
 kinInt12->compute(Tao12); 

 //t_time -= clock(); // Clock BEGIN
 std::shared_ptr<psi::OneBodyAOInt> kinInt11(fact_11.ao_kinetic());
 kinInt11->compute(Tao11);

 std::shared_ptr<psi::OneBodyAOInt> kinInt22(fact_22.ao_kinetic());
 kinInt22->compute(Tao22);
 //t_time += clock(); // Clock END

 // Compute potential energy integrals from CAMM and multipole integrals
 std::shared_ptr<oepdev::DMTPole> camm_1 = this->parameters.at("efp2")->dmtp("camm"); 
 std::shared_ptr<oepdev::DMTPole> camm_2 =other->parameters.at("efp2")->dmtp("camm"); 
 //
 const double p1 = 1.0 / 3.0;
 const double p2 = 2.0 / 3.0;
 const double p3 = 1.0 /15.0;
 const double p4 = 3.0 /15.0;
 const double p5 = 6.0 /15.0;
 const double prefacs[20] = {
 /* 0    X    Y    Z    XX  YY  ZZ  XY  XZ  YZ */
    1.0, 1.0, 1.0, 1.0, p1, p1, p1, p2, p2, p2, 
 /*    XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ */
       p3, p3, p3, p4, p4, p4, p4, p4, p4, p5};
 //
 //t_time -= clock(); // Clock BEGIN

 /* __12  ->   V from molecule 2 
    __11  ->   V from molecule 2
    __21  ->   V from molecule 1
    __22  ->   V from molecule 1 */

 // 
 size_t n_multipole_1 = mol_1->natom();
 size_t n_multipole_2 = mol_2->natom();
 assert (n_multipole_1 == camm_1->n_sites());
 assert (n_multipole_2 == camm_2->n_sites());
// t_time += clock();
 auto xyz_1 = this->extract_xyz(mol_1);
 auto xyz_2 = this->extract_xyz(mol_2);
 auto mult_1= this->extract_dmtp(camm_1);
 auto mult_2= this->extract_dmtp(camm_2);
// t_time -= clock();
 //
 int max_k = 3; int kdim = 20;
 if (psi::Process::environment.options.get_bool("EFP2_CT_NO_OCTUPOLES")) { max_k = 2; kdim = 10; }
 //std::shared_ptr<psi::OneBodyAOInt> efp_ints_12(fact_12.ao_efp_multipole_potential());
 //std::shared_ptr<psi::OneBodyAOInt> efp_ints_11(fact_11.ao_efp_multipole_potential());
 //std::shared_ptr<psi::OneBodyAOInt> efp_ints_21(fact_21.ao_efp_multipole_potential());
 //std::shared_ptr<psi::OneBodyAOInt> efp_ints_22(fact_22.ao_efp_multipole_potential());
 std::shared_ptr<psi::OneBodyAOInt> efp_ints_12(fact_12.ao_efp_multipole_potential_new(max_k,0));
 std::shared_ptr<psi::OneBodyAOInt> efp_ints_11(fact_11.ao_efp_multipole_potential_new(max_k,0));
 std::shared_ptr<psi::OneBodyAOInt> efp_ints_21(fact_21.ao_efp_multipole_potential_new(max_k,0));
 std::shared_ptr<psi::OneBodyAOInt> efp_ints_22(fact_22.ao_efp_multipole_potential_new(max_k,0));
 //
 std::vector<psi::SharedMatrix> mats_12, mats_11, mats_21, mats_22;
 for (int i = 0; i < kdim; ++i) {
      psi::SharedMatrix mm12 = std::make_shared<psi::Matrix>("", nbf_1, nbf_2);
      psi::SharedMatrix mm11 = std::make_shared<psi::Matrix>("", nbf_1, nbf_1);
      psi::SharedMatrix mm21 = std::make_shared<psi::Matrix>("", nbf_2, nbf_1);
      psi::SharedMatrix mm22 = std::make_shared<psi::Matrix>("", nbf_2, nbf_2);
      mats_12.push_back(mm12);
      mats_11.push_back(mm11);
      mats_21.push_back(mm21);
      mats_22.push_back(mm22);
      //mats_12.push_back(std::make_shared<psi::Matrix>("", nbf_1, nbf_2));
      //mats_11.push_back(std::make_shared<psi::Matrix>("", nbf_1, nbf_1));
      //mats_21.push_back(std::make_shared<psi::Matrix>("", nbf_2, nbf_1));
      //mats_22.push_back(std::make_shared<psi::Matrix>("", nbf_2, nbf_2));
 }
 //
 double *xyz_1_p = xyz_1->pointer();
 double *xyz_2_p = xyz_2->pointer();
 double *mult_1_p = mult_1->pointer();
 double *mult_2_p = mult_2->pointer();

 // Molecule B CAMM
 for (size_t n2 = 0; n2 < n_multipole_2; n2++) {

      for (int i = 0; i < kdim; ++i) {
           mats_12[i]->zero();
           mats_11[i]->zero();
      }
      psi::Vector3 coords_2(xyz_2_p[n2 * 3], xyz_2_p[n2 * 3 + 1], xyz_2_p[n2 * 3 + 2]);
      efp_ints_12->set_origin(coords_2);
      efp_ints_11->set_origin(coords_2);
      efp_ints_12->compute(mats_12);
      efp_ints_11->compute(mats_11);

      for (int i = 0; i < kdim; ++i) {                            
           mats_12[i]->scale(-prefacs[i] * mult_2_p[20 * n2 + i]);
           mats_11[i]->scale(-prefacs[i] * mult_2_p[20 * n2 + i]);
           VaoB12->add(mats_12[i]);
           VaoB11->add(mats_11[i]);
      }
 }

 // Molecule A CAMM
 for (size_t n1 = 0; n1 < n_multipole_1; n1++) {

      for (int i = 0; i < kdim; ++i) {
           mats_21[i]->zero();
           mats_22[i]->zero();
      }
      psi::Vector3 coords_1(xyz_1_p[n1 * 3], xyz_1_p[n1 * 3 + 1], xyz_1_p[n1 * 3 + 2]);
      efp_ints_21->set_origin(coords_1);
      efp_ints_22->set_origin(coords_1);
      efp_ints_21->compute(mats_21);
      efp_ints_22->compute(mats_22);

      for (int i = 0; i < kdim; ++i) {                            
           mats_21[i]->scale(-prefacs[i] * mult_1_p[20 * n1 + i]);
           mats_22[i]->scale(-prefacs[i] * mult_1_p[20 * n1 + i]);
           VaoA21->add(mats_21[i]);
           VaoA22->add(mats_22[i]);
      }
 }
 //
 //t_time += clock(); // Clock END


 // ---> Transform one electron contributions to MO basis <--- //
  

 // Transform S matrices //
 //t_time -= clock(); // Clock BEGIN
 psi::SharedMatrix Smoij    = psi::Matrix::triplet(cmoo_1, Sao12,  cmoo_2, true, false, false);
 psi::SharedMatrix Smoxn    = psi::Matrix::triplet(cmoo_1, Sao12,  cmov_2, true, false, false); // x \in OCC_A
 psi::SharedMatrix Smoyn    = psi::Matrix::triplet(cmov_1, Sao12,  cmov_2, true, false, false); // y \in VIR_A

 psi::SharedMatrix Smoxm    = psi::Matrix::triplet(cmoo_2, Sao12,  cmov_1, true, true , false); // x \in OCC_B
 psi::SharedMatrix Smoym    = psi::Matrix::triplet(cmov_2, Sao12,  cmov_1, true, true , false); // y \in VIR_B
                                                                                                                                

 // Transform T matrices //
 psi::SharedMatrix Tmonn    = psi::Matrix::triplet(cmov_2, Tao22,  cmov_2, true, false, false);
 psi::SharedMatrix Tmonj    = psi::Matrix::triplet(cmov_2, Tao22,  cmoo_2, true, false, false);                            
 psi::SharedMatrix Tmoxj    = psi::Matrix::triplet(cmoo_1, Tao12,  cmoo_2, true, false, false);  // x \in OCC_A (T_mj)
 psi::SharedMatrix Tmoyj    = psi::Matrix::triplet(cmov_1, Tao12,  cmoo_2, true, false, false);  // y \in VIR_A (T_mj)          
 psi::SharedMatrix Tmomm    = psi::Matrix::triplet(cmov_1, Tao11,  cmov_1, true, false, false);
 psi::SharedMatrix Tmomi    = psi::Matrix::triplet(cmov_1, Tao11,  cmoo_1, true, false, false);  
 psi::SharedMatrix Tmoxi    = psi::Matrix::triplet(cmoo_2, Tao12,  cmoo_1, true, true , false);  // x \in OCC_B (T_ni)
 psi::SharedMatrix Tmoyi    = psi::Matrix::triplet(cmov_2, Tao12,  cmoo_1, true, true , false);  // y \in VIR_B (T_ni)

                                                                                                                                
 // Transform F matrices //
 psi::SharedMatrix Fmoii    = psi::Matrix::triplet(cmoo_1, Fao11,  cmoo_1, true, false, false);
 psi::SharedMatrix Fmojj    = psi::Matrix::triplet(cmoo_2, Fao22,  cmoo_2, true, false, false);


 // Transform V matrices //                                                                                                   
 psi::SharedMatrix VmoBin   = psi::Matrix::triplet(cmoo_1, VaoB12, cmov_2, true, false, false); 
 psi::SharedMatrix VmoBix   = psi::Matrix::triplet(cmoo_1, VaoB11, cmoo_1, true, false, false); // x \in OCC_A (V_im)
 psi::SharedMatrix VmoBiy   = psi::Matrix::triplet(cmoo_1, VaoB11, cmov_1, true, false, false); // y \in VIR_A (V_im)

 psi::SharedMatrix VmoAjm   = psi::Matrix::triplet(cmoo_2, VaoA21, cmov_1, true, false, false); 
 psi::SharedMatrix VmoAjx   = psi::Matrix::triplet(cmoo_2, VaoA22, cmoo_2, true, false, false); // x \in OCC_B (V_jn)
 psi::SharedMatrix VmoAjy   = psi::Matrix::triplet(cmoo_2, VaoA22, cmov_2, true, false, false); // y \in VIR_B (V_jn)


 // ==> Vin_I and Wjm_I <== // 

 // Vin_I = VBin - VBix*Sxn - VBiy*Syn // 
 psi::SharedMatrix Vin_I = VmoBin->clone();
 Vin_I->gemm(false, false, -1.0, VmoBix, Smoxn, 1.0);
 Vin_I->gemm(false, false, -1.0, VmoBiy, Smoyn, 1.0);

 // Wjm_I = VAjm - Vjx*Sxm - Vjy*Sym 
 psi::SharedMatrix Wjm_I = VmoAjm->clone();
 Wjm_I->gemm(false, false, -1.0, VmoAjx, Smoxm, 1.0);
 Wjm_I->gemm(false, false, -1.0, VmoAjy, Smoym, 1.0);

 // ==> Vin_II and Wjm_II <== //

 // Vin_II = Vin_I + Sij*(Tnj - Sxn*Txj - Syn*Tyj) //
 Tmonj->gemm(true, false, -1.0, Smoxn, Tmoxj, 1.0);
 Tmonj->gemm(true, false, -1.0, Smoyn, Tmoyj, 1.0);
 psi::SharedMatrix Vin_II = psi::Matrix::doublet(Smoij, Tmonj, false, true);
 Vin_II->add(Vin_I);

 // Wjm_II = Wjm_I + Sij*(Tmi - Sxm*Txi - Sym*Tyi)
 Tmomi->gemm(true, false, -1.0, Smoxm, Tmoxi, 1.0);
 Tmomi->gemm(true, false, -1.0, Smoym, Tmoyi, 1.0);
 psi::SharedMatrix Wjm_II = psi::Matrix::doublet(Smoij, Tmomi, true, true);
 Wjm_II->add(Wjm_I);


 // ==> Final calculation of CT energy <== //

 double e_ct_AB = 0.0; // A ---> B
 for (int i=0; i<ndocc_1; ++i)
      {
      for (int n=0; n<nvir_2; ++n)
           {
           // Vin_I/(F_ii - T_nn)
           double value = Vin_I->get(i, n) / (Fmoii->get(i, i) - Tmonn->get(n, n));

           // (1 - S2_mn) -> normalizing factor
           double s2 = 0.0;
           for (int x=0; x<ndocc_1; ++x)
                {
                double s = Smoxn->get(x, n);
                s2 += s*s;
                }
           for (int y=0; y<nvir_1; ++y)
                {
                double s = Smoyn->get(y, n);
                s2 += s*s;
                }
           value /= (1.0 - s2);

           // Vin_I * Vin_II/(F_ii-T_nn)/(1-S2mn)
           value *= Vin_II->get(i, n);

           e_ct_AB += value;
           }
       }

 double e_ct_BA = 0.0; // B ---> A
 for (int j=0; j<ndocc_2; ++j)
      {
      for (int m=0; m<nvir_1; ++m)
           {
           // Wjm_I/(F_jj - T_mm)
           double value = Wjm_I->get(j, m) / (Fmojj->get(j, j) - Tmomm->get(m, m));

           // 1 - S2mn -> normalizing factor
           double s2 = 0.0;
           for (int x=0; x<ndocc_2; ++x)
                {
                double s = Smoxm->get(x, m);
                s2 += s*s;
                }
           for (int y=0; y<nvir_2; ++y)
                {
                double s = Smoym->get(y, m);
                s2 += s*s;
                }
           value /= (1.0 - s2);

           // Wjm_I * Wjm_II/(F_jj-T_mm)/(1-S2mn)
           value *= Wjm_II->get(j, m);

           e_ct_BA += value;
           }
       }

 // --> Compute total CT energy <-- //
 e_ct_AB *= 2.0;
 e_ct_BA *= 2.0;
 e_ct = e_ct_AB + e_ct_BA;
//cout << e_ct_AB << " " << e_ct_BA << endl;

 psi::Process::environment.globals["EINT CT EFP2 KCAL"] = e_ct * OEPDEV_AU_KcalPerMole;

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME EFP2:CT    " << ((double)t_time/CLOCKS_PER_SEC) << endl;

 return e_ct;
}
double oepdev::GenEffFrag::compute_pairwise_energy_efp2_disp(std::shared_ptr<GenEffFrag> other) const {
 // ---> Timer-on <--- //
 //clock_t t_time = -clock();

 // ---> Timer-off <--- //
 //t_time += clock(); // Clock END
 //cout << " o TIME EFP2:DISP  " << ((double)t_time/CLOCKS_PER_SEC) << endl;
 //TODO
 return 0.0;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_exrep(std::shared_ptr<GenEffFrag> other) const {

 /* Note:
  *
  * This implements repulsion energy from OEP-Murrell-etal by using localized molecular orbitals for S-1 term 
  * as well as a conventional EFP2 S-2 term by using the same orbitals
  * This saves a bit of time because the same orbital sets are used for S-1 and S-2. 
  * Theoretically, canonical orbitals could be also used for the S-1 term here with exactly the same result, 
  * but then more calculations and memory is required.
  */

 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 // Initialize results
 double e_s1 = 0.0, e_s2 = 0.0, e_ex = 0.0, e_tot = 0.0;

 // Extract EFP2 parameters from fragments
 const int nbf_p1 =  this->nbf();
 const int nbf_p2 = other->nbf();
 const int nat_1 =  this->natom();
 const int nat_2 = other->natom();
 const int ndocc_1=  this->ndocc();
 const int ndocc_2= other->ndocc();
 psi::SharedMolecule mol_1 =  this->molecule();
 psi::SharedMolecule mol_2 = other->molecule();
 psi::SharedBasisSet primary_1 =  this->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet primary_2 = other->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet auxiliary_1 =  this->parameters.at("efp2")->basisset("auxiliary");
 psi::SharedBasisSet auxiliary_2 = other->parameters.at("efp2")->basisset("auxiliary");
 psi::SharedMatrix lmoo_1 =  this->parameters.at("efp2")->matrix("lmoo-oep");
 psi::SharedMatrix lmoo_2 = other->parameters.at("efp2")->matrix("lmoo-oep");
 psi::SharedMatrix lmoc_1 =  this->parameters.at("efp2")->matrix("lmoc-oep");
 psi::SharedMatrix lmoc_2 = other->parameters.at("efp2")->matrix("lmoc-oep");
 const int nbf_a1 =auxiliary_1->nbf();
 const int nbf_a2 =auxiliary_2->nbf();

 // ===> Compute S^-1 term <=== //
 oepdev::SharedOEPotential oep_1 =  this->parameters.at("efp2")->oep("rep");
 oepdev::SharedOEPotential oep_2 = other->parameters.at("efp2")->oep("rep");

 psi::SharedMatrix Sao_1p2p     = std::make_shared<psi::Matrix>("Sao 1p2p", nbf_p1, nbf_p2);
 psi::SharedMatrix Sao_1a2p     = std::make_shared<psi::Matrix>("Sao 1a2p", nbf_a1, nbf_p2);
 psi::SharedMatrix Sao_1p2a     = std::make_shared<psi::Matrix>("Sao 1p2a", nbf_p1, nbf_a2);

 t_time += clock();
 psi::IntegralFactory fact_1p2p(  primary_1,   primary_2,   primary_1,   primary_2);
 psi::IntegralFactory fact_1a2p(auxiliary_1,   primary_2, auxiliary_1,   primary_2);
 psi::IntegralFactory fact_1p2a(  primary_1, auxiliary_2,   primary_1, auxiliary_2);
 t_time -= clock();

 std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2p(fact_1p2p.ao_overlap());
 std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
 std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

 ovlInt_1p2p->compute(Sao_1p2p);
 ovlInt_1a2p->compute(Sao_1a2p);
 ovlInt_1p2a->compute(Sao_1p2a);

 psi::SharedMatrix Smo = psi::Matrix::triplet(lmoo_1, Sao_1p2p, lmoo_2, true, false, false);
 psi::SharedMatrix Sba = psi::Matrix::doublet(lmoo_2, Sao_1a2p, true, true);
 psi::SharedMatrix Sab = psi::Matrix::doublet(lmoo_1, Sao_1p2a, true, false);

 psi::SharedMatrix SSG1= psi::Matrix::triplet(Smo, Sba, oep_1->matrix("Murrell-etal.S1"), false, false, false);
 psi::SharedMatrix SSG2= psi::Matrix::triplet(Smo, Sab, oep_2->matrix("Murrell-etal.S1"), true, false, false);

 e_s1  = SSG1->trace() + SSG2->trace();
 e_s1 *= -2.0;

 // ===> Compute S^-2 term: EFP2 model <=== //
 double** S = Smo->pointer();
 psi::SharedVector ss1 = std::make_shared<psi::Vector>("s^2 (1)", ndocc_1);
 psi::SharedVector ss2 = std::make_shared<psi::Vector>("s^2 (2)", ndocc_2);
 double* s1 = ss1->pointer();
 double* s2 = ss2->pointer();
 double val, sab;
 for (int a = 0; a<ndocc_1; ++a) {
      val = 0.0;
      for (int b = 0; b<ndocc_2; ++b) {
           sab = S[a][b];
           val += sab * sab;
      }
      s1[a] = val;
 }
 for (int b = 0; b<ndocc_2; ++b) {
      val = 0.0;
      for (int a = 0; a<ndocc_1; ++a) {
           sab = S[a][b];
           val += sab * sab;
      }
      s2[b] = val;
 }

 psi::SharedVector ZR1 = std::make_shared<psi::Vector>("ZR1", ndocc_1);
 psi::SharedVector ZR2 = std::make_shared<psi::Vector>("ZR2", ndocc_2);
 psi::SharedVector IR1 = std::make_shared<psi::Vector>("IR1", ndocc_1);
 psi::SharedVector IR2 = std::make_shared<psi::Vector>("IR2", ndocc_2);

 double* pZR1 = ZR1->pointer();
 double* pZR2 = ZR2->pointer();
 double* pIR1 = IR1->pointer();
 double* pIR2 = IR2->pointer();
 
 double rbx, ray, rad, rbc, rab;

 for (int b=0; b<ndocc_2; ++b){
      val = 0.0;
      for (int x=0; x<nat_1; ++x) { 
           rbx  = sqrt(pow(mol_1->x(x)-lmoc_2->get(b,0), 2.0) + 
                       pow(mol_1->y(x)-lmoc_2->get(b,1), 2.0) + 
                       pow(mol_1->z(x)-lmoc_2->get(b,2), 2.0) );
           val -= (double)mol_1->Z(x) / rbx;
      }
      pZR2[b] = val;
 }
 for (int a=0; a<ndocc_1; ++a){
      val = 0.0;
      for (int y=0; y<nat_2; ++y) { 
           ray  = sqrt(pow(mol_2->x(y)-lmoc_1->get(a,0), 2.0) + 
                       pow(mol_2->y(y)-lmoc_1->get(a,1), 2.0) + 
                       pow(mol_2->z(y)-lmoc_1->get(a,2), 2.0) );
           val -= (double)mol_2->Z(y) / ray;
      }
      pZR1[a] = val;
 }
 for (int a=0; a<ndocc_1; ++a) {
      val = 0.0;
      for (int b=0; b<ndocc_2; ++b) {
           rad  = sqrt(pow(lmoc_1->get(a,0)-lmoc_2->get(b,0), 2.0) + 
                       pow(lmoc_1->get(a,1)-lmoc_2->get(b,1), 2.0) + 
                       pow(lmoc_1->get(a,2)-lmoc_2->get(b,2), 2.0) );
           val += 2.0 / rad;
      }
      pIR1[a] = val;
 }
 for (int b=0; b<ndocc_2; ++b) {
      val = 0.0;
      for (int a=0; a<ndocc_1; ++a) {
           rbc  = sqrt(pow(lmoc_1->get(a,0)-lmoc_2->get(b,0), 2.0) + 
                       pow(lmoc_1->get(a,1)-lmoc_2->get(b,1), 2.0) + 
                       pow(lmoc_1->get(a,2)-lmoc_2->get(b,2), 2.0) );
           val += 2.0 / rbc;
      }
      pIR2[b] = val;
 }
 ZR1->add(IR1);
 ZR2->add(IR2);

 e_s2 = ss1->vector_dot(ZR1) + ss2->vector_dot(ZR2);

 for (int a=0; a<ndocc_1; ++a) {
      for (int b=0; b<ndocc_2; ++b) { 
           sab = S[a][b];
           rab = sqrt(pow(lmoc_1->get(a,0) - lmoc_2->get(b,0), 2.0) +
                      pow(lmoc_1->get(a,1) - lmoc_2->get(b,1), 2.0) +
                      pow(lmoc_1->get(a,2) - lmoc_2->get(b,2), 2.0) );
           e_s2 -= sab * sab / rab;
      }
 }
 e_s2 *= 2.0;


 // ===> Compute S_exc term: EFP2-SGO model <=== //
 for (int a=0; a<ndocc_1; ++a) {
      for (int b=0; b<ndocc_2; ++b) {
           sab = S[a][b];
           rab = sqrt(pow(lmoc_1->get(a,0) - lmoc_2->get(b,0), 2.0) + 
                      pow(lmoc_1->get(a,1) - lmoc_2->get(b,1), 2.0) + 
                      pow(lmoc_1->get(a,2) - lmoc_2->get(b,2), 2.0) );
           e_ex += sqrt(-2.0*log(abs(sab)) / M_PI) * sab * sab / rab;
      }
 }
 e_ex *= -4.0;

 // Save
 double v_e_rep = psi::Process::environment.globals["EINT REP OEP EFP2 KCAL"];   
 double v_e_exc = psi::Process::environment.globals["EINT EXC OEP EFP2 KCAL"];   
 double v_e_exr = psi::Process::environment.globals["EINT EXR OEP EFP2 KCAL"];   
 double v_e_rs1 = psi::Process::environment.globals["EINT REP OEP EFP2:S1 KCAL"];
 double v_e_rs2 = psi::Process::environment.globals["EINT REP OEP EFP2:S2 KCAL"];
 v_e_rep +=(e_s1 + e_s2        )*OEPDEV_AU_KcalPerMole;
 v_e_exc +=(e_ex               )*OEPDEV_AU_KcalPerMole;
 v_e_exr +=(e_s1 + e_s2 + e_ex )*OEPDEV_AU_KcalPerMole;
 v_e_rs1 +=(e_s1               )*OEPDEV_AU_KcalPerMole;
 v_e_rs2 +=(e_s2               )*OEPDEV_AU_KcalPerMole;

 psi::Process::environment.globals["EINT REP OEP EFP2 KCAL"]    = v_e_rep;
 psi::Process::environment.globals["EINT EXC OEP EFP2 KCAL"]    = v_e_exc;
 psi::Process::environment.globals["EINT EXR OEP EFP2 KCAL"]    = v_e_exr;
 psi::Process::environment.globals["EINT REP OEP EFP2:S1 KCAL"] = v_e_rs1;
 psi::Process::environment.globals["EINT REP OEP EFP2:S2 KCAL"] = v_e_rs2;
 
 e_tot = e_s1 + e_s2 + e_ex;
//cout << e_s1 << " " << e_s2 << " " << e_ex << endl;

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME OEP :EXREP " << ((double)t_time/CLOCKS_PER_SEC) << endl;
 return e_tot;
}
double oepdev::GenEffFrag::compute_pairwise_energy_oep_efp2_ct(std::shared_ptr<GenEffFrag> other) const {

 // ---> Timer-on <--- //
 clock_t t_time = -clock();

 // Initialize results
 double e_ab = 0.0, e_ba = 0.0, e_tot = 0.0;

 // Extract EFP2 parameters from fragments
 const int nbf_p1 =  this->nbf();
 const int nbf_p2 = other->nbf();
 const int nat_1 =  this->natom();
 const int nat_2 = other->natom();
 const int nocc_1=  this->ndocc();
 const int nocc_2= other->ndocc();
 oepdev::SharedOEPotential oep_1 =  this->parameters.at("efp2")->oep("ct");
 oepdev::SharedOEPotential oep_2 = other->parameters.at("efp2")->oep("ct");
 psi::SharedMolecule mol_1 =  this->molecule();
 psi::SharedMolecule mol_2 = other->molecule();
 psi::SharedBasisSet primary_1 =  this->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet primary_2 = other->parameters.at("efp2")->basisset("primary");
 psi::SharedBasisSet auxiliary_1 =  this->parameters.at("efp2")->basisset("auxiliary");
 psi::SharedBasisSet auxiliary_2 = other->parameters.at("efp2")->basisset("auxiliary");
 psi::SharedMatrix cmoo_1 =  this->parameters.at("efp2")->matrix("cmoo");
 psi::SharedMatrix cmoo_2 = other->parameters.at("efp2")->matrix("cmoo");
 psi::SharedMatrix cmov_1 =  this->parameters.at("efp2")->matrix("cmov");
 psi::SharedMatrix cmov_2 = other->parameters.at("efp2")->matrix("cmov");
 psi::SharedMatrix lmoc_1 =  this->parameters.at("efp2")->matrix("lmoc-oep");
 psi::SharedMatrix lmoc_2 = other->parameters.at("efp2")->matrix("lmoc-oep");
 psi::SharedMatrix lmoo_1 =  this->parameters.at("efp2")->matrix("lmoo-oep");
 psi::SharedMatrix lmoo_2 = other->parameters.at("efp2")->matrix("lmoo-oep");
 psi::SharedMatrix U_1 =  this->parameters.at("efp2")->matrix("cmoo2lmoo");
 psi::SharedMatrix U_2 = other->parameters.at("efp2")->matrix("cmoo2lmoo");
 psi::SharedVector e_occ_1 =  this->parameters.at("efp2")->vector("epso");
 psi::SharedVector e_occ_2 = other->parameters.at("efp2")->vector("epso"); 
 psi::SharedVector e_vir_1 =  this->parameters.at("efp2")->vector("epsv"); 
 psi::SharedVector e_vir_2 = other->parameters.at("efp2")->vector("epsv"); 

 const int nbf_a1 =auxiliary_1->nbf();
 const int nbf_a2 =auxiliary_2->nbf();
 const int nvir_1= cmov_1->ncol();
 const int nvir_2= cmov_2->ncol();

 psi::SharedMatrix Sao_1p2p     = std::make_shared<psi::Matrix>("Sao 1p2p", nbf_p1, nbf_p2);
 psi::SharedMatrix Sao_1a2p     = std::make_shared<psi::Matrix>("Sao 1a2p", nbf_a1, nbf_p2);
 psi::SharedMatrix Sao_1p2a     = std::make_shared<psi::Matrix>("Sao 1p2a", nbf_p1, nbf_a2);

 t_time += clock();
  psi::IntegralFactory fact_1p2p(primary_1, primary_2, primary_1, primary_2);
  psi::IntegralFactory fact_1a2p(auxiliary_1, primary_2, auxiliary_1, primary_2);
  psi::IntegralFactory fact_1p2a(primary_1, auxiliary_2, primary_1, auxiliary_2);
 t_time -= clock();

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2p(fact_1p2p.ao_overlap());
  ovlInt_1p2p->compute(Sao_1p2p);

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

  ovlInt_1a2p->compute(Sao_1a2p);
  ovlInt_1p2a->compute(Sao_1p2a);

  psi::SharedMatrix S1 = psi::Matrix::doublet(cmoo_1, Sao_1p2a, true, false); // OCC(A) x AUX(B)
  psi::SharedMatrix S2 = psi::Matrix::doublet(cmoo_2, Sao_1a2p, true, true ); // OCC(B) x AUX(A)


  // ---> Localized occupied orbitals <--- //
  psi::SharedMatrix S12= psi::Matrix::triplet(lmoo_1, Sao_1p2p, lmoo_2, true, false, false); // LOCC(A) x LOCC(B)
  psi::SharedMatrix S1Y= psi::Matrix::triplet(lmoo_1, Sao_1p2p, cmov_2, true, false, false); // LOCC(A) x  VIR(B)
  psi::SharedMatrix S2X= psi::Matrix::triplet(lmoo_2, Sao_1p2p, cmov_1, true, true , false); // LOCC(B) x  VIR(A)
  

  // ---> Compute auxiliary tensors u and w <--- //
  psi::SharedVector u_1 = this->compute_u_vector(lmoc_1, lmoc_2, mol_2);
  psi::SharedVector u_2 = this->compute_u_vector(lmoc_2, lmoc_1, mol_1);
  psi::SharedMatrix w_1 = this->compute_w_matrix(mol_1, mol_2, lmoc_1);
  psi::SharedMatrix w_2 = this->compute_w_matrix(mol_2, mol_1, lmoc_2);
  if (psi::Process::environment.options.get_int("PRINT") > 2) {
  u_1->set_name("Auxiliary calculables: Vector u_1"); u_1->print();
  w_1->set_name("Auxiliary calculables: Matrix w_1"); w_1->print();
  u_2->set_name("Auxiliary calculables: Vector u_2"); u_2->print();
  w_2->set_name("Auxiliary calculables: Matrix w_2"); w_2->print();
  }

  // ---> Get distributed effective charges <--- //
  std::vector<psi::SharedMatrix> q_1 = oep_1->oep("Otto-Ladik.V3.CAMM-nj").dmtp->charges();
  std::vector<psi::SharedMatrix> q_2 = oep_2->oep("Otto-Ladik.V3.CAMM-nj").dmtp->charges();
//cout << q_1.size() << endl;
//cout << q_1[0]->nrow() << " " << q_1[0]->ncol() << endl;


  // ===> Compute V1 term <=== //
  psi::SharedMatrix v_ab_v1 = psi::Matrix::doublet(S1, oep_2->matrix("Otto-Ladik.V1.GDF"), false, false);
  psi::SharedMatrix v_ba_v1 = psi::Matrix::doublet(S2, oep_1->matrix("Otto-Ladik.V1.GDF"), false, false);
  //const double sc = 0.25;
  //v_ab_v1->scale(sc);
  //v_ba_v1->scale(sc);

  // ===> Compute V2 term <=== //
  psi::SharedMatrix v_ab_v2 = std::make_shared<psi::Matrix>("", nocc_1, nvir_2);
  psi::SharedMatrix v_ba_v2 = std::make_shared<psi::Matrix>("", nocc_2, nvir_1);
  for (int i=0; i<nocc_1; ++i) {
       for (int n=0; n<nvir_2; ++n) {
            double v = S1Y->get(i,n) * u_1->get(i);
            v_ab_v2->set(i, n, v); 
       }
  }
  for (int i=0; i<nocc_2; ++i) {
       for (int n=0; n<nvir_1; ++n) {
            double v = S2X->get(i,n) * u_2->get(i);
            v_ba_v2->set(i, n, v); 
       }
  }
  psi::SharedMatrix v_ab_v2_copy = v_ab_v2->clone(); v_ab_v2->zero();
  psi::SharedMatrix v_ba_v2_copy = v_ba_v2->clone(); v_ba_v2->zero();
  v_ab_v2->gemm(false, false, 1.0, U_1, v_ab_v2_copy, 0.0);
  v_ba_v2->gemm(false, false, 1.0, U_2, v_ba_v2_copy, 0.0);

  // ===> Compute V3 term <=== //
  const bool old_implementation = true;
  psi::SharedMatrix v_ab_v3 = std::make_shared<psi::Matrix>("", nocc_1, nvir_2);
  psi::SharedMatrix v_ba_v3 = std::make_shared<psi::Matrix>("", nocc_2, nvir_1);

//cout << S12->ncol() << " " << S12->nrow() << endl;
//cout << w_1->ncol() << " " << w_1->nrow() << endl;

//cout << nvir_1*nocc_1 + nvir_1 << " " << q_1.size() << endl;
//cout << nvir_2*nocc_2 + nvir_2 << " " << q_2.size() << endl;
//cout << nvir_1 << " " << nocc_1 << " a" << endl;

  if (old_implementation) {
      for (int i=0; i<nocc_1; ++i) {                                                        
           for (int n=0; n<nvir_2; ++n) {
                double v = 0;
                for (int j=0; j<nocc_2; ++j) {              
                     for (int y=0; y<mol_2->natom(); ++y) {
                          v += S12->get(i,j) * w_1->get(i, y) * q_2[nvir_2*j+n]->get(y, 0);
                     }
                }
                v_ab_v3->set(i, n, v);
           }
      }
      //
      for (int i=0; i<nocc_2; ++i) {
           for (int n=0; n<nvir_1; ++n) {
                double v = 0;
                for (int j=0; j<nocc_1; ++j) {              
                     for (int y=0; y<mol_1->natom(); ++y) {
                          v += S12->get(j,i) * w_2->get(i, y) * q_1[nvir_1*j+n]->get(y, 0);
                     }
                }
                v_ba_v3->set(i, n, v);
           }
      }
  } else {
      double** ss = S12->pointer();
      double** w1p= w_1->pointer();
      double** w2p= w_2->pointer();
      for (int i=0; i<nocc_1; ++i) {                                                        
           for (int n=0; n<nvir_2; ++n) {
                double v = 0;
                for (int j=0; j<nocc_2; ++j) {                                              
                     double Q2 = q_2[nvir_2*j+n]->sum_of_squares();
                     if (Q2 > 0.0001) {
                         double** qq = q_2[nvir_2*j+n]->pointer();
                         for (int y=0; y<mol_2->natom(); ++y) {
                              v += ss[i][j] * w1p[i][y] * qq[y][0];
                         }
                     }
                }
                v_ab_v3->set(i, n, v);
           }
      }
      //
      for (int i=0; i<nocc_2; ++i) {
           for (int n=0; n<nvir_1; ++n) {
                double v = 0;
                for (int j=0; j<nocc_1; ++j) {                                    
                     double Q2 = q_1[nvir_1*j+n]->sum_of_squares();
                     if (Q2 > 0.0001) {
                         double** qq = q_1[nvir_1*j+n]->pointer();
                         for (int y=0; y<mol_1->natom(); ++y) {
                              v += ss[j][i] * w2p[i][y] * qq[y][0];
                         }
                     }
                }
                v_ba_v3->set(i, n, v);
           }
      }

  }

  psi::SharedMatrix v_ab_v3_copy = v_ab_v3->clone(); v_ab_v3->zero();
  psi::SharedMatrix v_ba_v3_copy = v_ba_v3->clone(); v_ba_v3->zero();
  v_ab_v3->gemm(false, false, 1.0, U_1, v_ab_v3_copy, 0.0);
  v_ba_v3->gemm(false, false, 1.0, U_2, v_ba_v3_copy, 0.0);


  // ---> Add coupling constant contributions <--- //
//std::shared_ptr<psi::Matrix> v_ab_v12 = v_ab_v1->clone(); v_ab_v12->add(v_ab_v2);
//std::shared_ptr<psi::Matrix> v_ba_v12 = v_ba_v1->clone(); v_ba_v12->add(v_ba_v2);
//std::shared_ptr<psi::Matrix> v_ab_v13 = v_ab_v1->clone(); v_ab_v13->add(v_ab_v3);
//std::shared_ptr<psi::Matrix> v_ba_v13 = v_ba_v1->clone(); v_ba_v13->add(v_ba_v3);
//std::shared_ptr<psi::Matrix> v_ab_v23 = v_ab_v2->clone(); v_ab_v23->add(v_ab_v3);
//std::shared_ptr<psi::Matrix> v_ba_v23 = v_ba_v2->clone(); v_ba_v23->add(v_ba_v3);
//std::shared_ptr<psi::Matrix> v_ab_v123= v_ab_v12->clone(); v_ab_v123->add(v_ab_v3);
//std::shared_ptr<psi::Matrix> v_ba_v123= v_ba_v12->clone(); v_ba_v123->add(v_ba_v3);

  v_ab_v1->add(v_ab_v2);
  v_ab_v1->add(v_ab_v3);

  v_ba_v1->add(v_ba_v2);
  v_ba_v1->add(v_ba_v3);


  // ===> Compute CT Energy <=== //
//double e_ab_v1  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v1);
//double e_ba_v1  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v1);
//double e_ab_v2  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v2);
//double e_ba_v2  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v2);
//double e_ab_v3  = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v3);
//double e_ba_v3  = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v3);
//double e_ab_v12 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v12);
//double e_ba_v12 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v12);
//double e_ab_v13 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v13);
//double e_ba_v13 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v13);
//double e_ab_v23 = this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v23);
//double e_ba_v23 = this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v23);
//double e_ab_v123= this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v123);
//double e_ba_v123= this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v123);

  double e_ab_v123= this->compute_ct_component(e_occ_1, e_vir_2, v_ab_v1  );
  double e_ba_v123= this->compute_ct_component(e_occ_2, e_vir_1, v_ba_v1  );

  e_ab = e_ab_v123;
  e_ba = e_ba_v123;

  e_tot = e_ab + e_ba;

 // ---> Timer-off <--- //
 t_time += clock(); // Clock END
 cout << " o TIME OEP :CT    " << ((double)t_time/CLOCKS_PER_SEC) << endl;

 return e_tot;
}

// --> quick helpers
psi::SharedVector oepdev::GenEffFrag::extract_xyz(psi::SharedMolecule mol) const
{
  auto xyz = std::make_shared<psi::Vector>(3 * mol->natom());
  double* xyz_p = xyz->pointer();
  for (int i = 0; i < mol->natom(); ++i) {
       *xyz_p++ = mol->x(i);
       *xyz_p++ = mol->y(i);
       *xyz_p++ = mol->z(i); 
  }
  return xyz;
}
//
psi::SharedVector oepdev::GenEffFrag::extract_dmtp(std::shared_ptr<oepdev::DMTPole> camm) const
{
  auto mult= std::make_shared<psi::Vector>((1 + 3 + 6 + 10) * camm->n_sites());
  psi::SharedMatrix m_0 = camm->charges(0)->clone();
  psi::SharedMatrix m_1 = camm->dipoles(0)->clone();
  psi::SharedMatrix m_2 = camm->quadrupoles(0)->clone();
  psi::SharedMatrix m_3 = camm->octupoles(0)->clone();
  //m_1->zero();
  //m_2->zero(); 
  if (psi::Process::environment.options.get_bool("EFP2_CT_NO_OCTUPOLES")) m_3->zero();

  double* mult_p = mult->pointer();
  double** p_0 = m_0->pointer();
  double** p_1 = m_1->pointer();
  double** p_2 = m_2->pointer();
  double** p_3 = m_3->pointer();

  for (int i = 0; i < camm->n_sites(); ++i) {
       *mult_p++ = p_0[i][0];   // 0

       *mult_p++ = p_1[i][0];   // X
       *mult_p++ = p_1[i][1];   // Y
       *mult_p++ = p_1[i][2];   // Z

        double t = 0.5 * (p_2[i][0] + p_2[i][3] + p_2[i][5]);
       *mult_p++ = p_2[i][0] * 1.5 - t;   // XX
       *mult_p++ = p_2[i][3] * 1.5 - t;   // YY
       *mult_p++ = p_2[i][5] * 1.5 - t;   // ZZ
       *mult_p++ = p_2[i][1] * 1.5;   // XY
       *mult_p++ = p_2[i][2] * 1.5;   // XZ
       *mult_p++ = p_2[i][4] * 1.5;   // YZ

        double tx = 0.5 * (p_3[i][0] + p_3[i][3] + p_3[i][5]);
        double ty = 0.5 * (p_3[i][6] + p_3[i][1] + p_3[i][8]);
        double tz = 0.5 * (p_3[i][9] + p_3[i][2] + p_3[i][7]);
       *mult_p++ = p_3[i][0] * 2.5 - 3.0 * tx;   // XXX
       *mult_p++ = p_3[i][6] * 2.5 - 3.0 * ty;   // YYY
       *mult_p++ = p_3[i][9] * 2.5 - 3.0 * tz;   // ZZZ
       *mult_p++ = p_3[i][1] * 2.5 -       ty;   // XXY
       *mult_p++ = p_3[i][2] * 2.5 -       tz;   // XXZ
       *mult_p++ = p_3[i][3] * 2.5 -       tx;   // XYY
       *mult_p++ = p_3[i][7] * 2.5 -       tz;   // YYZ
       *mult_p++ = p_3[i][5] * 2.5 -       tx;   // XZZ
       *mult_p++ = p_3[i][8] * 2.5 -       ty;   // YZZ
       *mult_p++ = p_3[i][4] * 2.5           ;   // XYZ
  }
  return mult;
}
psi::SharedVector oepdev::GenEffFrag::compute_u_vector(
  psi::SharedMatrix rmo_1, 
  psi::SharedMatrix rmo_2, 
  psi::SharedMolecule mol_2
) const
{
  const int nocc_1 = rmo_1->nrow();
  const int nocc_2 = rmo_2->nrow();
  const int nat_2  = mol_2->natom();

  psi::SharedVector u = std::make_shared<psi::Vector>("", nocc_1);

  for (int i=0; i<nocc_1; ++i) {
       double v = 0.0;
       //
       for (int y=0; y<nat_2; ++y) {
            double ryi = sqrt(pow(mol_2->x(y) - rmo_1->get(i,0), 2.0) +
                              pow(mol_2->y(y) - rmo_1->get(i,1), 2.0) +
                              pow(mol_2->z(y) - rmo_1->get(i,2), 2.0) );
            v += (double)mol_2->Z(y) / ryi;
       }
       //
       for (int j=0; j<nocc_2; ++j) {
            double rji = sqrt(pow(rmo_2->get(j,0) - rmo_1->get(i,0), 2.0) +
                              pow(rmo_2->get(j,1) - rmo_1->get(i,1), 2.0) +
                              pow(rmo_2->get(j,2) - rmo_1->get(i,2), 2.0) );
            v -= 2.0 / rji;
       }
       u->set(i, v);
  }
  return u;
}
std::shared_ptr<psi::Matrix> oepdev::GenEffFrag::compute_w_matrix(
  psi::SharedMolecule mol_1,
  psi::SharedMolecule mol_2,
  psi::SharedMatrix rmo_1
) const
{
  const int nat_1 = mol_1->natom();
  const int nat_2 = mol_2->natom();
  const int nocc_1= rmo_1->nrow();

  psi::SharedMatrix w = std::make_shared<psi::Matrix>("", nocc_1, nat_2);

  for (int y=0; y<nat_2; ++y) {
       double vy = 0.0;
       //
       for (int x=0; x<nat_1; ++x) { 
            double rxy = sqrt(pow(mol_1->x(x) - mol_2->x(y), 2.0) +
                              pow(mol_1->y(x) - mol_2->y(y), 2.0) +
                              pow(mol_1->z(x) - mol_2->z(y), 2.0) );
            vy += (double)mol_1->Z(x) / rxy;
       }
       //
       for (int k=0; k<nocc_1; ++k) {
            double rky = sqrt(pow(rmo_1->get(k,0) - mol_2->x(y), 2.0) +
                              pow(rmo_1->get(k,1) - mol_2->y(y), 2.0) +
                              pow(rmo_1->get(k,2) - mol_2->z(y), 2.0) );
            vy-= 2.0 / rky;
       }
       //
       for (int i=0; i<nocc_1; ++i) {
            //
            double riy = sqrt(pow(rmo_1->get(i,0) - mol_2->x(y), 2.0) +
                              pow(rmo_1->get(i,1) - mol_2->y(y), 2.0) +
                              pow(rmo_1->get(i,2) - mol_2->z(y), 2.0) );
            w->set(i, y, vy + 2.0/riy );
       }
  }
  w->scale(-1.0);
  return w;
}
double oepdev::GenEffFrag::compute_ct_component(
 psi::SharedVector eps_occ_X, psi::SharedVector eps_vir_Y, psi::SharedMatrix V) const
{
   // requires matrix elements in CMO (canonical SCF) basis
   const int nocc_X = eps_occ_X->dim();
   const int nvir_Y = eps_vir_Y->dim();
   double e_XY = 0.0;
   for (int i=0; i<nocc_X; ++i) { /* X-->Y term */
       for (int n=0; n<nvir_Y; ++n) {
            double vin = V->get(i, n); /* V: OCC_X x VIR_Y */
            e_XY += (vin * vin) / (eps_occ_X->get(i) - eps_vir_Y->get(n));
       }
  }
  e_XY *= 2.0;
  return e_XY;
}
