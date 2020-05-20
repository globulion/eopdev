#include "cis.h"
#include<iostream>
#include "util.h"
#include "psi4/libmints/mintshelper.h"
#include "../../include/oepdev_files.h"


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
          E_(nullptr),
          U_(nullptr),
          H_(nullptr),
          nmo_(ref_wfn_->nmo()-ref_wfn_->nfrzc()),
          naocc_(ref_wfn_->nalpha()-ref_wfn_->nfrzc()),
          nbocc_(ref_wfn_->nbeta()-ref_wfn_->nfrzc()),
          navir_(nmo_ - naocc_),
          nbvir_(nmo_ - nbocc_),
          transformation_type_(trans_type),
          inttrans_(nullptr),
          jk_(nullptr)
{
  this->common_init(); 
}

CISComputer::~CISComputer() {}

void CISComputer::common_init(void) {
 ndets_ = naocc_ * navir_ + nbocc_ * nbvir_;
 if (true) {
     std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(ref_wfn_->basisset());
     mints->integrals();
 }
 // nstates_ is not set here but during computation since it depends on the type of CIS calculation

 // Construct the JK object
 if (ref_wfn_->basisset_exists("BASIS_DF_SCF")) {
     jk_ = psi::JK::build_JK(ref_wfn_->basisset(), ref_wfn_->get_basisset("BASIS_DF_SCF"), options_);
 } else {
     jk_ = psi::JK::build_JK(ref_wfn_->basisset(), BasisSet::zero_ao_basis_set(), options_);
 }
 jk_->set_memory((options_.get_double("SCF_MEM_SAFETY_FACTOR")*(psi::Process::environment.get_memory() / 8L)));
 jk_->initialize();
 jk_->print_header();
}

void CISComputer::print_header_(void) {
 std::string cis_algorithm = options_.get_str("CIS_TYPE");
 std::string cis_type;
 if (options_.get_str("REFERENCE") == "RHF") {cis_type = "R";}
 else {cis_type = "U";}

 psi::outfile->Printf("\n ===> Starting %s-CIS Calculation <===\n\n", cis_type.c_str());
 psi::outfile->Printf("  => Setup <=\n\n");
 psi::outfile->Printf("   Memory       = %11zu MB\n" , psi::Process::environment.get_memory() / (1024L*1024L));
 psi::outfile->Printf("   Algorithm    = %s\n"       , cis_algorithm.c_str());
 if (cis_algorithm == "DAVIDSON_LIU") {
 psi::outfile->Printf("     Conver     = %11.3E\n"   , options_.get_double("DAVIDSON_LIU_CONVER"));
 psi::outfile->Printf("     Maxiter    = %11d\n"     , options_.get_int("DAVIDSON_LIU_MAXITER"));
 psi::outfile->Printf("     BVecGuess  = %s\n"       , options_.get_str("DAVIDSON_LIU_GUESS").c_str());
 psi::outfile->Printf("     Thr_large  = %11.3E\n"   , options_.get_double("DAVIDSON_LIU_THRESH_LARGE"));
 psi::outfile->Printf("     Thr_small  = %11.3E\n"   , options_.get_double("DAVIDSON_LIU_THRESH_SMALL"));
 psi::outfile->Printf("     Nroots     = %11d\n"     , options_.get_int("DAVIDSON_LIU_NROOTS"));
 psi::outfile->Printf("     SpStart    = %11d\n"     , options_.get_int("DAVIDSON_LIU_SPACE_START"));
 psi::outfile->Printf("     SpMax      = %11d\n"     , options_.get_int("DAVIDSON_LIU_SPACE_MAX"));
 }
 psi::outfile->Printf("   SchwartzCut  = %11.3E\n"   , options_.get_double("CIS_SCHWARTZ_CUTOFF"));
 psi::outfile->Printf("   PrintAmpl    = %11.3f\n"   , options_.get_double("OEPDEV_AMPLITUDE_PRINT_THRESHOLD"));
 psi::outfile->Printf("   StandAmpl    = %s\n"       , options_.get_bool("OEPDEV_STANDARDIZE_AMPLITUDES")?"true":"false");
 psi::outfile->Printf("   Nbasis       = %11d\n"     , ref_wfn_->basisset()->nbf());
 psi::outfile->Printf("   Nsingles     = %11d\n"     , naocc_*navir_ + nbocc_*nbvir_);
 psi::outfile->Printf("   Naocc(Active)= %11d\n"     , naocc_);
 psi::outfile->Printf("   Navir(Active)= %11d\n"     , navir_);
 psi::outfile->Printf("   Nocc(Frozen) = %11d\n"     , ref_wfn_->nfrzc());
 psi::outfile->Printf("   Nvir(Frozen) = %11d\n"     , ref_wfn_->frzvpi()[0]);
 psi::outfile->Printf("   <S2>(Ref)    = %11.3f\n"   , this->compute_s2_reference());
 psi::outfile->Printf("   <S2>(Exact)  = %11.3f\n"   , this->compute_s2_exact());
 psi::outfile->Printf("\n");
 psi::outfile->Printf("  => Molecule <=\n\n");
 ref_wfn_->molecule()->print();
 psi::outfile->Printf("  => Primary Basis <=\n\n");
 ref_wfn_->basisset()->print();
}

double CISComputer::compute_s2_reference(void) {
 double s2 = this->compute_s2_exact() + (double)ref_wfn_->nbeta();

 psi::SharedMatrix S  = this->ref_wfn_->S();
 psi::SharedMatrix Ca = this->ref_wfn_->Ca_subset("AO","OCC");
 psi::SharedMatrix Cb = this->ref_wfn_->Cb_subset("AO","OCC");
 psi::SharedMatrix Sab = psi::Matrix::triplet(Ca, S, Cb, true, false, false);
 for (int a=0; a<ref_wfn_->nalpha(); ++a) {
 for (int b=0; b<ref_wfn_->nbeta (); ++b) {
      double s = Sab->get(a,b);
      s2 -= s*s;
 }
 }
 return s2;
}

double CISComputer::compute_s2_exact(void) {
 double ns = (double)ref_wfn_->nalpha() - (double)ref_wfn_->nbeta();
 return (ns*(ns+2.0)/4.0);
}

//double CISComputer::compute_overlap_between_singles_determinants(int i, int a, int j, int b, int s1, int s2,
//      psi::SharedMatrix& U1, psi::SharedMatrix& V1, psi::SharedVector& sigma1,
//      psi::SharedMatrix& U2, psi::SharedMatrix& V2, psi::SharedVector& sigma2,
//      psi::SharedMatrix& Ca_occ, psi::SharedMatrix& Cb_occ, 
//      psi::SharedMatrix& Ca_vir, psi::SharedMatrix& Cb_vir) {
//
//  // overlap = < D_ia(s1) | D_jb(s2) >
//
//  psi::SharedMatrix C1_a, C1_b, C2_a, C2_b; 
//  const int nfrzc = ref_wfn_->nfrzc();
//  // D_ia(s1)
//  if (s1 == 0) {
//      C1_a = Ca_occ->clone(); // substitute i-->a
//      C1_b = Cb_occ->clone(); // retain from reference
//
//      psi::SharedVector c_a = Ca_vir->get_column(0,a);
//      for (int k=0; k<c_a->dim(); ++k) {
//           C1_a->set(i+nfrzc, k, c_a->get(k));
//      }
//
//  } else {
//      C1_b = Cb_occ->clone(); // substitute i-->a
//      C1_a = Ca_occ->clone(); // retain from reference
//
//      psi::SharedVector c_a = Cb_vir->get_column(0,a);
//      for (int k=0; k<c_a->dim(); ++k) {
//           C1_b->set(i+nfrzc, k, c_a->get(k));
//      }
//  }
//  // D_jb(s2)
//  if (s1 == 0) {
//      C2_a = Ca_occ->clone(); // substitute j-->b
//      C2_b = Cb_occ->clone(); // retain from reference
//
//      psi::SharedVector c_b = Ca_vir->get_column(0,b);
//      for (int k=0; k<c_b->dim(); ++k) {
//           C2_a->set(j+nfrzc, k, c_b->get(k));
//      }
//
//  } else {
//      C2_b = Cb_occ->clone(); // substitute j-->b
//      C2_a = Ca_occ->clone(); // retain from reference
//
//      psi::SharedVector c_b = Cb_vir->get_column(0,b);
//      for (int k=0; k<c_b->dim(); ++k) {
//           C2_b->set(j+nfrzc, k, c_b->get(k));
//      }
//  }
//
//  psi::SharedMatrix Sab_a = psi::Matrix::triplet(C1_a, this->ref_wfn_->S(), C2_a, true, false, false);
//  psi::SharedMatrix Sab_b = psi::Matrix::triplet(C1_b, this->ref_wfn_->S(), C2_b, true, false, false);
//
//  U1->zero(); V1->zero(); sigma1->zero();
//  U2->zero(); V2->zero(); sigma2->zero();
//  Sab_a->svd(U1, sigma1, V1);
//  Sab_b->svd(U2, sigma2, V2);
//  double s=1.0;
//  for (int k=0; k<sigma1->dim(); ++k) s*= sigma1->get(k);
//  for (int k=0; k<sigma2->dim(); ++k) s*= sigma2->get(k);
//
//  return s;
//}

double CISComputer::s2(int I) const {
  // all reference molecular orbitals
  psi::SharedMatrix Ca_occ = this->ref_wfn_->Ca_subset("AO","OCC");
  psi::SharedMatrix Cb_occ = this->ref_wfn_->Cb_subset("AO","OCC");
  psi::SharedMatrix Ca_vir = this->ref_wfn_->Ca_subset("AO","VIR");
  psi::SharedMatrix Cb_vir = this->ref_wfn_->Cb_subset("AO","VIR");

  const int na = ref_wfn_->nalpha();
  const int nb = ref_wfn_->nbeta ();

  const double ns = (double)na - (double)nb;

  psi::SharedMatrix Dij = psi::Matrix::triplet(Ca_occ, this->ref_wfn_->S(), Cb_occ, true, false, false);
  psi::SharedMatrix Dab = psi::Matrix::triplet(Ca_vir, this->ref_wfn_->S(), Cb_vir, true, false, false);
  psi::SharedMatrix Dia = psi::Matrix::triplet(Ca_occ, this->ref_wfn_->S(), Cb_vir, true, false, false);
  psi::SharedMatrix Dai = psi::Matrix::triplet(Ca_vir, this->ref_wfn_->S(), Cb_occ, true, false, false);
  psi::SharedMatrix Pij_A = std::make_shared<psi::Matrix>("", na, na); 
  psi::SharedMatrix Pab_A = std::make_shared<psi::Matrix>("", navir_, navir_); 
  psi::SharedMatrix Pij_B = std::make_shared<psi::Matrix>("", nb, nb); 
  psi::SharedMatrix Pab_B = std::make_shared<psi::Matrix>("", nbvir_, nbvir_); 

  psi::SharedMatrix Qab_A   = psi::Matrix::doublet(Dai,Dai,false,true);
  psi::SharedMatrix Qij_A   = psi::Matrix::doublet(Dij,Dij,false,true);

  psi::SharedMatrix Qab_B   = psi::Matrix::doublet(Dia,Dia,true,false);
  psi::SharedMatrix Qij_B   = psi::Matrix::doublet(Dij,Dij,true,false);

  const int off = naocc_*navir_;
  const int nfrzc = ref_wfn_->nfrzc();

  // <S2>_UHF
  double s2 = ns*0.5*(ns*0.5 + 1.0) + (double)nb - Dij->vector_dot(Dij);
  double ds2= 0.0;

  // Pab_A
  for (int a=0; a<this->navir_; ++a) {
  for (int b=0; b<this->navir_; ++b) {
       double v = 0.0;
       for (int i=0; i<this->naocc_; ++i) {
            int ia = navir_*i + a;
            int ib = navir_*i + b;
            v += U_->get(ia,I) * U_->get(ib,I);
       }
       Pab_A->set(a,b,v);
  }
  }
  // Pab_B
  for (int a=0; a<this->nbvir_; ++a) {
  for (int b=0; b<this->nbvir_; ++b) {
       double v = 0.0;
       for (int i=0; i<this->nbocc_; ++i) {
            int ia = nbvir_*i + a;
            int ib = nbvir_*i + b;
            v += U_->get(ia+off,I) * U_->get(ib+off,I);
       }
       Pab_B->set(a,b,v);
  }
  }

  // Pij_A
  for (int i=0; i<this->naocc_; ++i) {
  for (int j=0; j<this->naocc_; ++j) {
       double v = 0.0;
       for (int a=0; a<this->navir_; ++a) {
            int ia = navir_*i + a;
            int ja = navir_*j + a;
            v -= U_->get(ia,I) * U_->get(ja,I);
       }
       Pij_A->set(i+nfrzc,j+nfrzc,v);
  }
  }
  // Pij_B
  for (int i=0; i<this->nbocc_; ++i) {
  for (int j=0; j<this->nbocc_; ++j) {
       double v = 0.0;
       for (int a=0; a<this->nbvir_; ++a) {
            int ia = nbvir_*i + a;
            int ja = nbvir_*j + a;
            v -= U_->get(ia+off,I) * U_->get(ja+off,I);
       }
       Pij_B->set(i+nfrzc,j+nfrzc,v);
  }
  }

  // UCIS contribution
  ds2 -= Qab_A->vector_dot(Pab_A);
  ds2 -= Qij_A->vector_dot(Pij_A);

  ds2 -= Qab_B->vector_dot(Pab_B);
  ds2 -= Qij_B->vector_dot(Pij_B);

  psi::SharedMatrix t_ia = std::make_shared<psi::Matrix>("", naocc_+nfrzc, navir_);
  psi::SharedMatrix t_jb = std::make_shared<psi::Matrix>("", nbocc_+nfrzc, nbvir_);
  for (int i=0; i<this->naocc_; ++i) {
  for (int a=0; a<this->navir_; ++a) {
       int ia = navir_*i + a;
       t_ia->set(i+nfrzc,a,U_->get(ia,I));
  }
  }
  for (int j=0; j<this->nbocc_; ++j) {
  for (int b=0; b<this->nbvir_; ++b) {
       int jb = nbvir_*j + b;
       t_jb->set(j+nfrzc,b,U_->get(jb+off,I));
  }
  }
  psi::SharedMatrix W = psi::Matrix::triplet(Dij, t_ia, Dab, true, false, false);
  ds2 -= 2.0 * W->vector_dot(t_jb);
 
  // slower code for quadruple sum over ia,jb
  //for (int i=0; i<this->naocc_; ++i) {
  //for (int a=0; a<this->navir_; ++a) {
  //     int ia = navir_*i + a;
  //     for (int j=0; j<this->nbocc_; ++j) {
  //     for (int b=0; b<this->nbvir_; ++b) {
  //          int jb = nbvir_*j + b;
  //          ds2 -= 2.0 * Dij->get(i+nfrzc,j+nfrzc) * Dab->get(a,b) * U_->get(ia,I) * U_->get(jb+off,I);
  //     }
  //     }
  //}
  //}
  s2 += ds2;

  return s2;
}


void CISComputer::print_excited_states_(void) {
 std::map<int, std::string> state_multiplicity;
 state_multiplicity[ 1] = "Singlet";
 state_multiplicity[ 2] = "Doublet";
 state_multiplicity[ 3] = "Triplet";
 state_multiplicity[ 4] = "Quadrup";
 state_multiplicity[ 5] = "Quintup";
 state_multiplicity[ 6] = "Quintup";
 state_multiplicity[ 7] = "Hextupl";
 state_multiplicity[ 8] = "Heptupl";
 state_multiplicity[ 9] = "Octuplt";
 state_multiplicity[10] = "Nonuplt";

 psi::outfile->Printf("\n ===> Excited States <===\n\n");
 for (int I=0; I<this->nstates_; ++I) {
       double E_ex = this->E_->get(I);
       psi::SharedVector tj = this->transition_dipole(I);
       double tjx = tj->get(0);
       double tjy = tj->get(1);
       double tjz = tj->get(2);

       double s2 = this->s2(I);
       double f  = this->oscillator_strength(I);

       double s = round(sqrt(4.0*s2+1.0));
       std::string mult = state_multiplicity.at((int)s);

       psi::outfile->Printf("     State= %2d [%6s], <S2>= %6.3f, f= %9.6f [a.u.] E= %9.3f [EV] TrMu= (%9.3f, %9.3f, %9.3f) [A.U.]\n", 
                                 I+1, mult.c_str(), s2, f, E_ex*OEPDEV_AU_EV, tjx, tjy, tjz);

       this->print_excited_state_character_(I);
 }
}

// Abstract
void CISComputer::print_excited_state_character_(int I) {}

void CISComputer::set_nstates_() {
  this->nstates_ = this->ndets_;
}
void CISComputer::set_beta_(void) {}

void CISComputer::allocate_memory(void) {
 eps_a_o_ = ref_wfn_->epsilon_a_subset("MO", "ACTIVE_OCC");
 eps_a_v_ = ref_wfn_->epsilon_a_subset("MO", "ACTIVE_VIR");
 this->allocate_hamiltonian_();
}

void CISComputer::allocate_hamiltonian_(void) {
 U_ = std::make_shared<psi::Matrix>("CIS Eigenvectors", ndets_, nstates_);
 E_ = std::make_shared<psi::Vector>("CIS Eigenvalues", nstates_);
 H_ = std::make_shared<psi::Matrix>("CIS Excited State Hamiltonian", ndets_, ndets_);
}

void CISComputer::compute(void) {
 this->print_header_();
 this->set_nstates_();    // Maybe better to move it to constructor? 
                          // (then Davidson-Liu must be re-adapted) - this needs to be now here
 this->allocate_memory(); // moved here to allocate E_ and U_ before Davidson-Liu needs to be initialized
 this->prepare_for_cis_();
 this->build_hamiltonian_();
 this->diagonalize_hamiltonian_(); 
 this->standardize_amplitudes_();
 this->print_excited_states_();

 // Clear memory
 this->jk_->finalize();
 psi::outfile->Printf("\n @CIS: Done.\n\n");
}

void CISComputer::prepare_for_cis_(void) {
// Fa_oo_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","OCC"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","OCC"), true, false, false);
// Fa_vv_ = psi::Matrix::triplet(ref_wfn_->Ca_subset("AO","VIR"), ref_wfn_->Fa(), ref_wfn_->Ca_subset("AO","VIR"), true, false, false);
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
                                                   psi::IntegralTransform::FrozenOrbitals::OccOnly); 
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

void CISComputer::standardize_amplitudes_(void) {
 if (options_.get_bool("CIS_STANDARDIZE_AMPLITUDES")==true) {

     for (int I=0; I<nstates_; ++I) {
          double dominant_amplitude = 0.0;
          for (int k=0; k<U_->nrow(); ++k) {
               double u = U_->get(k,I);
               if (std::abs(u) > std::abs(dominant_amplitude)) dominant_amplitude = u;
          }
          if (dominant_amplitude < 0.0) U_->scale_column(0,I,-1.0);
     }

 }
}

void CISComputer::davidson_liu_compute_sigma(void) {}
void CISComputer::davidson_liu_compute_diagonal_hamiltonian(void) {}

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
  const int nfrzc = ref_wfn_->nfrzc();
  psi::SharedMatrix Co = ref_wfn_->Ca_subset("AO","OCC");
  psi::SharedMatrix Cv = ref_wfn_->Ca_subset("AO","VIR");
  psi::SharedMatrix D = std::make_shared<psi::Matrix>(naocc_+nfrzc, navir_); 
  for (int i=0; i<naocc_; ++i) {
  for (int a=0; a<navir_; ++a) {
       int ia = navir_*i + a;
       D->set(i+nfrzc, a, U_->get(ia, J));
  }
  }
  psi::SharedMatrix d = psi::Matrix::triplet(Co, D, Cv, false, false, true);
  return d;
}

psi::SharedMatrix CISComputer::Tb_ao(int J) const {
  const int nfrzc = ref_wfn_->nfrzc();
  psi::SharedMatrix Co = ref_wfn_->Cb_subset("AO","OCC");
  psi::SharedMatrix Cv = ref_wfn_->Cb_subset("AO","VIR");
  psi::SharedMatrix D = std::make_shared<psi::Matrix>(nbocc_+nfrzc, nbvir_); 
  const int off = naocc_*navir_;
  for (int i=0; i<nbocc_; ++i) {
  for (int a=0; a<nbvir_; ++a) {
       int ia = nbvir_*i + a + off;
       D->set(i+nfrzc, a, U_->get(ia, J));
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
 const int nmo_all = ref_wfn_->nmo();
 const int nfrzc   = ref_wfn_->nfrzc();

 psi::SharedMatrix D = std::make_shared<psi::Matrix>(string_sprintf("Alpha Density Matrix (MO) I=%3d", I), nmo_all, nmo_all);
 for (int i=0; i<nfrzc; ++i) {
   D->set(i, i, 1.0);
 }
 for (int i=0; i<naocc_; ++i) {
   D->set(i+nfrzc, i+nfrzc, 1.0);
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
      d[p+nfrzc][q+nfrzc] -= v;
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
      d[naocc_ + p + nfrzc][naocc_ + q + nfrzc] += v;
 }
 } 

 // VO and OV blocks: zero (for unrelaxed density which is the case)
 return D;
}

psi::SharedMatrix CISComputer::Db_mo(int I) const {
 const int nmo_all = ref_wfn_->nmo();
 const int nfrzc   = ref_wfn_->nfrzc();

 psi::SharedMatrix D = std::make_shared<psi::Matrix>(string_sprintf("Beta Density Matrix (MO) I=%3d", I), nmo_all, nmo_all);
 for (int i=0; i<nfrzc; ++i) {
   D->set(i, i, 1.0);
 }
 for (int i=0; i<nbocc_; ++i) {
   D->set(i+nfrzc, i+nfrzc, 1.0);
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
      d[p+nfrzc][q+nfrzc] -= v;
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
      d[naocc_ + p + nfrzc][naocc_ + q + nfrzc] += v;
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

// Build routine
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
     if (cis_type == "DAVIDSON_LIU") {
	 cis = std::make_shared<U_CISComputer_DL>(ref_wfn, opt);
     } else {
	 cis = std::make_shared<U_CISComputer_Explicit>(ref_wfn, opt);
     }
  }
  else {
     if (cis_type == "DIRECT_EXPLICIT") { 
	 cis = std::make_shared<R_CISComputer_Direct>(ref_wfn, opt); 
     } else 
     if (cis_type == "DAVIDSON_LIU") {
         cis = std::make_shared<R_CISComputer_DL>(ref_wfn, opt);
     } else { // Explicit CIS
	 cis = std::make_shared<R_CISComputer_Explicit>(ref_wfn, opt);
     }
  }
  
  // Return 
  return cis;
}


} // EndNameSpace oepdev
