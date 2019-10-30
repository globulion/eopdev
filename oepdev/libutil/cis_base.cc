#include "cis.h"
#include<iostream>
#include "util.h"
#include "psi4/libmints/mintshelper.h"
#include "../../include/oepdev_files.h"

//#include "psi4/libpsi4util/PsiOutStream.h"

namespace oepdev{

const std::vector<std::string> CISComputer::reference_types = {"RHF", "UHF"};

CISComputer::CISComputer(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt, 
                         psi::IntegralTransform::TransformationType trans_type):
          DavidsonLiu(opt),
          ref_wfn_(wfn),
          //options_(opt),
          //Fa_oo_(nullptr),
          //Fa_vv_(nullptr),
          //Fb_oo_(nullptr),
          //Fb_vv_(nullptr),
          eps_a_o_(nullptr),
          eps_a_v_(nullptr),
          eps_b_o_(nullptr),
          eps_b_v_(nullptr),
          //E_(nullptr),
          //U_(nullptr),
          H_(nullptr),
          nmo_(ref_wfn_->nmo()),
          naocc_(ref_wfn_->nalpha()),
          nbocc_(ref_wfn_->nbeta()),
          navir_(nmo_ - naocc_),
          nbvir_(nmo_ - nbocc_),
          transformation_type_(trans_type),
          inttrans_(nullptr),
          jk_(nullptr)
{
  this->common_init(); 
}

CISComputer::~CISComputer() {}

void CISComputer::compute(void) {
 this->prepare_for_cis_();
 this->build_hamiltonian_();
 this->diagonalize_hamiltonian_(); 

 // Clear memory
 this->jk_->finalize();
}

void CISComputer::allocate_hamiltonian_(void) {
 H_ = std::make_shared<psi::Matrix>("CIS Excited State Hamiltonian", ndets_, ndets_);
}

void CISComputer::allocate_memory(void) {
 U_ = std::make_shared<psi::Matrix>("CIS Eigenvectors", ndets_, nstates_);
 E_ = std::make_shared<psi::Vector>("CIS Eigenvalues", nstates_);
 eps_a_o_ = ref_wfn_->epsilon_a_subset("MO", "OCC");
 eps_a_v_ = ref_wfn_->epsilon_a_subset("MO", "VIR");
 this->allocate_hamiltonian_();
}

void CISComputer::prepare_for_cis_(void) {
// Fa_oo_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","OCC"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","OCC"), true, false, false);
// Fa_vv_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","VIR"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","VIR"), true, false, false);
 this->allocate_memory();
 this->set_beta_();
 this->transform_integrals_();
}

void CISComputer::clear_dpd(void) {
  // Destruct the IntegralTransform object
  if (inttrans_) this->inttrans_.reset();
}

void CISComputer::transform_integrals_(void) {
 SharedMOSpaceVector spaces;
 spaces.push_back(psi::MOSpace::occ);
 spaces.push_back(psi::MOSpace::vir);
 inttrans_ = std::make_shared<psi::IntegralTransform>(ref_wfn_, spaces, this->transformation_type_, 
                                                   psi::IntegralTransform::OutputType::DPDOnly,
                                                   psi::IntegralTransform::MOOrdering::QTOrder,
                                                   psi::IntegralTransform::FrozenOrbitals::None); 
 inttrans_->set_keep_dpd_so_ints(true);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::occ);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::occ, psi::MOSpace::vir, psi::MOSpace::vir);
 inttrans_->transform_tei(psi::MOSpace::occ, psi::MOSpace::vir, psi::MOSpace::occ, psi::MOSpace::vir);
}

void CISComputer::diagonalize_hamiltonian_(void) {
 H_->diagonalize(U_, E_);
 //E_->scale(OEPDEV_AU_EV);
 if (options_.get_int("PRINT")>3) {
    E_->print_out();
    H_->print_out();
 }
}

std::shared_ptr<CISComputer> CISComputer::build(const std::string& type, 
                                                std::shared_ptr<psi::Wavefunction> ref_wfn, 
                                                psi::Options& opt, const std::string& reference) {

  // Determine reference if not specified
  std::string ref = reference;
  if (ref.empty()) {
      if (!ref_wfn->same_a_b_orbs() && !ref_wfn->same_a_b_dens()) {ref += "UHF";}
      else { ref += "RHF";}
  }

  // Sanity checks
  bool b = false;
  for (auto &refc : CISComputer::reference_types) {
       if (ref == refc) {b = true; break;}
  }
  if (!b) {throw psi::PSIEXCEPTION("Incorrect reference wavefunction type chosen. Only RHF and UHF are available");}

  if (ref =="RHF" and ref_wfn->molecule()->multiplicity() != 1)
   throw psi::PSIEXCEPTION("RHF reference cannot be set for open-shell system!");

  // Create
  std::shared_ptr<CISComputer> cis;
  std::string cis_type = opt.get_str("CIS_TYPE");

  if ((ref_wfn->molecule()->multiplicity() != 1) || (ref == "UHF")) { 
	 cis = std::make_shared<U_CISComputer>(ref_wfn, opt);
  }
  else {
     if (cis_type == "DIRECT") { 
	 cis = std::make_shared<R_CISComputer_Direct>(ref_wfn, opt); 
     } else 
     if (cis_type == "DAVIDSON_LIU") {
         cis = std::make_shared<R_CISComputer_DL>(ref_wfn, opt);
     } else { // Explicit CIS
	 cis = std::make_shared<R_CISComputer>(ref_wfn, opt);
     }
  }
  
  // Return 
  return cis;
}

void CISComputer::set_beta_(void) {}

void CISComputer::common_init(void) {
 ndets_ = naocc_ * navir_ + nbocc_ * nbvir_;
 if (true) {
     std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(ref_wfn_->basisset());
     mints->integrals();
 }
 nstates_ = ndets_; // Assumes Explicit CIS (diagonalization of entire Hamiltonian)

 // Construct the JK object
 jk_ = psi::JK::build_JK(ref_wfn_->basisset(), psi::BasisSet::zero_ao_basis_set(), options_);
 jk_->set_memory((options_.get_double("SCF_MEM_SAFETY_FACTOR")*(psi::Process::environment.get_memory() / 8L)));
 jk_->initialize();
 jk_->print_header();
}

std::pair<double,double> CISComputer::U_homo_lumo(int I, int h, int l) const {
  int i  = naocc_-1-h;
  int a  = 0+l;
  int ia = navir_*i + a;
  //
  int j  = nbocc_-1-h;
  int b  = 0+l;
  int jb = nbvir_*j + b + naocc_*navir_;
  //
  std::pair<double,double> t(U_->get(ia,I), U_->get(jb,I));
  return t;
}

psi::SharedMatrix CISComputer::Da_ao(int I) const {
  psi::SharedMatrix da_mo = this->Da_mo(I);
  psi::SharedMatrix C = ref_wfn_->Ca_subset("AO","ALL");
  psi::SharedMatrix d = psi::Matrix::triplet(C, da_mo, C, false, false, true);
  return d;
}

psi::SharedMatrix CISComputer::Db_ao(int I) const {
  psi::SharedMatrix db_mo = this->Db_mo(I);
  psi::SharedMatrix C = ref_wfn_->Cb_subset("AO","ALL");
  psi::SharedMatrix d = psi::Matrix::triplet(C, db_mo, C, false, false, true);
  return d;
}

psi::SharedMatrix CISComputer::Ta_ao(int J) const {
  psi::SharedMatrix Co = ref_wfn_->Ca_subset("AO","OCC");
  psi::SharedMatrix Cv = ref_wfn_->Ca_subset("AO","VIR");
  psi::SharedMatrix D = std::make_shared<psi::Matrix>(naocc_, navir_); 
  for (int i=0; i<naocc_; ++i) {
  for (int a=0; a<navir_; ++a) {
       int ia = navir_*i + a;
       D->set(i, a, U_->get(ia, J));
  }
  }
  psi::SharedMatrix d = psi::Matrix::triplet(Co, D, Cv, false, false, true);
  return d;
}

psi::SharedMatrix CISComputer::Tb_ao(int J) const {
  psi::SharedMatrix Co = ref_wfn_->Cb_subset("AO","OCC");
  psi::SharedMatrix Cv = ref_wfn_->Cb_subset("AO","VIR");
  psi::SharedMatrix D = std::make_shared<psi::Matrix>(nbocc_, nbvir_); 
  const int off = naocc_*navir_;
  for (int i=0; i<nbocc_; ++i) {
  for (int a=0; a<nbvir_; ++a) {
       int ia = nbvir_*i + a + off;
       D->set(i, a, U_->get(ia, J));
  }
  }
  psi::SharedMatrix d = psi::Matrix::triplet(Co, D, Cv, false, false, true);
  return d;
}

psi::SharedMatrix CISComputer::Ta_ao(int I, int J) const {
  throw psi::PSIEXCEPTION("Transition densities between excited states are not implemented yet!");
  // the below code is not correct
  psi::SharedMatrix T_0I = this->Ta_ao(I);
  psi::SharedMatrix T_0J = this->Ta_ao(J);
  T_0J->subtract(T_0I);
  return T_0J;
}

psi::SharedMatrix CISComputer::Tb_ao(int I, int J) const {
  throw psi::PSIEXCEPTION("Transition densities between excited states are not implemented yet!");
  // the below code is not correct
  psi::SharedMatrix T_0I = this->Tb_ao(I);
  psi::SharedMatrix T_0J = this->Tb_ao(J);
  T_0J->subtract(T_0I);
  return T_0J;
}

psi::SharedMatrix CISComputer::Da_mo(int I) const {
 psi::SharedMatrix D = std::make_shared<psi::Matrix>(string_sprintf("Alpha Density Matrix (MO) I=%3d", I), nmo_, nmo_);
 for (int i=0; i<naocc_; ++i) {
   D->set(i, i, 1.0);
 }
 double** U = U_->pointer();    
 double** d = D->pointer();    

 // OO block
 for (int p=0; p<naocc_; ++p) {
 for (int q=0; q<naocc_; ++q) {
      double v = 0.0;
      for (int a=0; a<navir_; ++a) {
           int pa = navir_*p + a;
           int qa = navir_*q + a;
           v += U[pa][I] * U[qa][I];
      }
      d[p][q] -= v;
 }
 } 

 // VV block
 for (int p=0; p<navir_; ++p) {
 for (int q=0; q<navir_; ++q) {
      double v = 0.0;
      for (int i=0; i<naocc_; ++i) {
           int pi = navir_*i + p;
           int qi = navir_*i + q;
           v += U[pi][I] * U[qi][I];
      }
      d[naocc_ + p][naocc_ + q] += v;
 }
 } 

 // VO and OV blocks: zero (for unrelaxed density which is the case)
 return D;
}

psi::SharedMatrix CISComputer::Db_mo(int I) const {
 psi::SharedMatrix D = std::make_shared<psi::Matrix>(string_sprintf("Beta Density Matrix (MO) I=%3d", I), nmo_, nmo_);
 for (int i=0; i<nbocc_; ++i) {
   D->set(i, i, 1.0);
 }
 const int off = navir_ * naocc_;
 double** U = U_->pointer();    
 double** d = D->pointer();    

 // OO block
 for (int p=0; p<nbocc_; ++p) {
 for (int q=0; q<nbocc_; ++q) {
      double v = 0.0;
      for (int a=0; a<nbvir_; ++a) {
           int pa = nbvir_*p + a;
           int qa = nbvir_*q + a;
           v += U[pa+off][I] * U[qa+off][I];
      }
      d[p][q] -= v;
 }
 } 

 // VV block
 for (int p=0; p<nbvir_; ++p) {
 for (int q=0; q<nbvir_; ++q) {
      double v = 0.0;
      for (int i=0; i<nbocc_; ++i) {
           int pi = nbvir_*i + p;
           int qi = nbvir_*i + q;
           v += U[pi+off][I] * U[qi+off][I];
      }
      d[naocc_ + p][naocc_ + q] += v;
 }
 } 

 // VO and OV blocks: zero (for unrelaxed density which is the case)
 return D;
}

SharedDMTPole CISComputer::camm(int j, bool symmetrize) const {
 psi::SharedMatrix T = this->Da_ao(j);
 T->add(this->Db_ao(j));
 if (symmetrize) {
   psi::SharedMatrix Tt = T->transpose();
   T->add(Tt);
   T->scale(0.5000000);
 }
 std::vector<psi::SharedMatrix> Tvec; Tvec.push_back(T);
 std::vector<bool> Bvec; Bvec.push_back(false);
 
 SharedDMTPole camm = oepdev::DMTPole::build("CAMM", ref_wfn_, 1);
 camm->compute(Tvec, Bvec);
 return camm;
}


SharedDMTPole CISComputer::trcamm(int j, bool symmetrize) const {
 psi::SharedMatrix T = this->Ta_ao(j);
 T->add(this->Tb_ao(j));
 if (symmetrize) {
   psi::SharedMatrix Tt = T->transpose();
   T->add(Tt);
   T->scale(0.5000000);
 }
 std::vector<psi::SharedMatrix> Tvec; Tvec.push_back(T);
 std::vector<bool> Bvec; Bvec.push_back(true);
 
 SharedDMTPole camm = oepdev::DMTPole::build("CAMM", ref_wfn_, 1);
 camm->compute(Tvec, Bvec);
 return camm;
}

SharedDMTPole CISComputer::trcamm(int i, int j, bool symmetrize) const {
 psi::SharedMatrix T = this->Ta_ao(i, j);
 T->add(this->Tb_ao(i, j));
 if (symmetrize) {
   psi::SharedMatrix Tt = T->transpose();
   T->add(Tt);
   T->scale(0.5000000);
 }
 std::vector<psi::SharedMatrix> Tvec; Tvec.push_back(T);
 std::vector<bool> Bvec; Bvec.push_back(true);
 
 SharedDMTPole camm = oepdev::DMTPole::build("CAMM", ref_wfn_, 1);
 camm->compute(Tvec, Bvec);
 return camm;
}

SharedVector CISComputer::transition_dipole(int j) const {
 psi::SharedMatrix T = this->Ta_ao(j);
 T->add(this->Tb_ao(j));

 std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(ref_wfn_->basisset());
 std::vector<psi::SharedMatrix> Dint = mints->ao_dipole();
 SharedVector tr = std::make_shared<psi::Vector>("Transition dipole", 3);
 for (int z=0; z<3; ++z) {
  tr->set(z, psi::Matrix::doublet(Dint[z], T)->trace());
 }
 return tr;
}

SharedVector CISComputer::transition_dipole(int i, int j) const {
 psi::SharedMatrix T = this->Ta_ao(i, j);
 T->add(this->Tb_ao(i, j));

 std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(ref_wfn_->basisset());
 std::vector<psi::SharedMatrix> Dint = mints->ao_dipole();
 SharedVector tr = std::make_shared<psi::Vector>("Transition dipole", 3);
 for (int z=0; z<3; ++z) {
  tr->set(z, psi::Matrix::doublet(Dint[z], T)->trace());
 }
 return tr;
}

double CISComputer::oscillator_strength(int j) const {
 double d = this->transition_dipole(j)->sum_of_squares();
 return (2.0/3.0) * d * E_->get(j);
}

double CISComputer::oscillator_strength(int i, int j) const {
 double d = this->transition_dipole(i, j)->sum_of_squares();
 return (2.0/3.0) * d * (E_->get(j) - E_->get(i));
}

void CISComputer::determine_electronic_state(int& I) {
 if (I<1) {
   int count = 1;
   const double ft = options_.get_double("OSCILLATOR_STRENGTH_THRESHOLD");
   for (int i=0; i<this->nstates(); ++i) {
        if (this->oscillator_strength(i) > ft) {
            if (count == -I) {
                I = i+1;
                break;
            }
            count += 1;
        }
   }
 } 
 I -= 1; // transform to C++ indexing (from 0)
}

std::shared_ptr<CISData> CISComputer::data(int I, int h, int l, bool symm) {
  // Only for closed shells as for now
  //
  double E_ex = this->eigenvalues()->get(I);                        // Excitation energy wrt ground state
  double t = this->U_homo_lumo(I, h, l).first;                      // CIS amplitude (HOMO-h)-->(LUMO+l)
  //
  SharedMatrix Pe  = this->Da_ao(I); Pe ->add(this->Db_ao(I));      // Excited state bond order matrices of monomers
  SharedMatrix Peg = this->Ta_ao(I); Peg->add(this->Tb_ao(I));      // Transition density matrices of monomers
  //
  SharedDMTPole trcamm = this->trcamm(I, symm);                     // TrCAMM moments
  //
  int nbf = this->ref_wfn_->basisset()->nbf();
  SharedVector ch = this->ref_wfn_->Ca_subset("AO","OCC")->get_column(0, this->ref_wfn_->nalpha() - 1 - h);
  SharedVector cl = this->ref_wfn_->Ca_subset("AO","VIR")->get_column(0, l);
  SharedMatrix D_homo = std::make_shared<psi::Matrix>("", nbf, nbf);
  SharedMatrix D_lumo = std::make_shared<psi::Matrix>("", nbf, nbf);
  for (int i=0; i<nbf; ++i) {
       for (int j=0; j<nbf; ++j) {
            D_homo->set(i, j, ch->get(i) * ch->get(j));
            D_lumo->set(i, j, cl->get(i) * cl->get(j));
       }
  }
  std::vector<psi::SharedMatrix> Hvec; Hvec.push_back(D_homo);
  std::vector<psi::SharedMatrix> Lvec; Lvec.push_back(D_lumo);
  std::vector<bool> Bvec; Bvec.push_back(true);

  SharedDMTPole camm_homo = DMTPole::build("CAMM", ref_wfn_, 1);    // CAMM of (HOMO-h) orbital
  camm_homo->compute(Hvec, Bvec);

  SharedDMTPole camm_lumo = DMTPole::build("CAMM", ref_wfn_, 1);    // CAMM of (LUMO+l) orbital
  camm_lumo->compute(Lvec, Bvec);

      psi::outfile->Printf("     State= %2d, f= %9.6f [a.u.] E= %9.3f [EV] t(H-%d->L+%d)= %9.6f\n", 
                                 I+1, this->oscillator_strength(I), E_ex*OEPDEV_AU_EV, h, l, t);

  std::shared_ptr<CISData> cis_data = std::make_shared<CISData>();
  cis_data->E_ex = E_ex;
  cis_data->t_homo_lumo = t;
  cis_data->Pe = Pe;
  cis_data->Peg= Peg;
  cis_data->trcamm = trcamm;
  cis_data->camm_homo = camm_homo;
  cis_data->camm_lumo = camm_lumo;
  return cis_data;
}




} // EndNameSpace oepdev
