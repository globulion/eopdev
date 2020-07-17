#include "quambo.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/mintshelper.h"
#include "util.h"
#include "psi4/libscf_solver/rohf.h"


namespace psi {
 using SharedBasisSet = std::shared_ptr<BasisSet>;
}

namespace oepdev{
 using SharedQUAMBOData = std::shared_ptr<QUAMBOData>;

QUAMBO::QUAMBO(psi::SharedWavefunction wfn, bool acbs) : 
 wfn_(wfn), acbs(acbs), options_(wfn->options()),
 mol_(wfn->molecule()),
 nbas_atom_mini_({}),
 unpe_atom_({}),
 Sao_(nullptr),
 quambo_a_nonorthogonal_(nullptr),
 quambo_a_orthogonal_(nullptr),
 quambo_b_nonorthogonal_(nullptr),
 quambo_b_orthogonal_(nullptr),
 c_a_mini_vir_(nullptr),
 c_b_mini_vir_(nullptr),
 e_a_mini_vir_(nullptr),
 e_b_mini_vir_(nullptr),
 c_a_mini_(nullptr),
 c_b_mini_(nullptr),
 e_a_mini_(nullptr),
 e_b_mini_(nullptr),
 nbas_mini_(-1),
 naocc_mini_(-1),
 nbocc_mini_(-1),
 navir_mini_(-1),
 nbvir_mini_(-1),
 nbf_(-1)
{
 // Sanity checks
 //if (!this->acbs) throw psi::PSIEXCEPTION("MCBS is not possible to implement for QUAMBO routine yet!");

 // Populate the maps
 this->nbas_atom_mini_[ "H"] = 1;    this->unpe_atom_[ "H"] = 1;   
 this->nbas_atom_mini_["He"] = 1;    this->unpe_atom_["He"] = 0;
 //                                    
 this->nbas_atom_mini_["Li"] = 5;    this->unpe_atom_["Li"] = 1;
 this->nbas_atom_mini_["Be"] = 5;    this->unpe_atom_["Be"] = 0;
 this->nbas_atom_mini_[ "B"] = 5;    this->unpe_atom_[ "B"] = 1;
 this->nbas_atom_mini_[ "C"] = 5;    this->unpe_atom_[ "C"] = 4; // Assumed 1s2 2s1 2px 1 2py 1 2pz1
 this->nbas_atom_mini_[ "N"] = 5;    this->unpe_atom_[ "N"] = 3;
 this->nbas_atom_mini_[ "O"] = 5;    this->unpe_atom_[ "O"] = 2;
 this->nbas_atom_mini_[ "F"] = 5;    this->unpe_atom_[ "F"] = 1;
 this->nbas_atom_mini_["Ne"] = 5;    this->unpe_atom_["Ne"] = 0;
 //                                    
 this->nbas_atom_mini_["Na"] = 9;    this->unpe_atom_["Na"] = 1;
 this->nbas_atom_mini_["Mg"] = 9;    this->unpe_atom_["Mg"] = 0;
 this->nbas_atom_mini_["Al"] = 9;    this->unpe_atom_["Al"] = 1;
 this->nbas_atom_mini_["Si"] = 9;    this->unpe_atom_["Si"] = 4; // Assumed 1s2 2s1 2p6 3s1 3px 1 3py 1 3pz1
 this->nbas_atom_mini_[ "P"] = 9;    this->unpe_atom_[ "P"] = 3;
 this->nbas_atom_mini_[ "S"] = 9;    this->unpe_atom_[ "S"] = 2;
 this->nbas_atom_mini_["Cl"] = 9;    this->unpe_atom_["Cl"] = 1;
 this->nbas_atom_mini_["Ar"] = 9;    this->unpe_atom_["Ar"] = 0;
 //                                    
 this->nbas_atom_mini_[ "K"] =13;    this->unpe_atom_[ "K"] = 1;
 this->nbas_atom_mini_["Ca"] =13;    this->unpe_atom_["Ca"] = 0;
 //TODO Add more atoms here!

}
QUAMBO::~QUAMBO() {}
void QUAMBO::compute(void) {
 psi::outfile->Printf("\n ===> Computing QUAMBOs <===\n\n");

 // [1] Read the orbitals and Fock matrix
 psi::SharedMatrix ca_occ = this->wfn_->Ca_subset("AO","OCC");
 psi::SharedMatrix ca_vir = this->wfn_->Ca_subset("AO","VIR");
 psi::SharedMatrix cb_occ = this->wfn_->Cb_subset("AO","OCC");
 psi::SharedMatrix cb_vir = this->wfn_->Cb_subset("AO","VIR");

 psi::SharedBasisSet bfs = this->wfn_->basisset();

 psi::SharedVector eps_a_occ = this->wfn_->epsilon_a_subset("MO", "OCC");
 psi::SharedVector eps_a_vir = this->wfn_->epsilon_a_subset("MO", "VIR");
 psi::SharedVector eps_a_all = this->wfn_->epsilon_a_subset("MO", "ALL");
 psi::SharedVector eps_b_occ = this->wfn_->epsilon_b_subset("MO", "OCC");
 psi::SharedVector eps_b_vir = this->wfn_->epsilon_b_subset("MO", "VIR");
 psi::SharedVector eps_b_all = this->wfn_->epsilon_b_subset("MO", "ALL");
 
 psi::SharedMatrix Fa = this->wfn_->Fa();
 psi::SharedMatrix Fb = this->wfn_->Fb();

 this->Sao_ = this->wfn_->S ();

 // sizing
 const int ndocc = this->wfn_->doccpi()[0];
 const int nsocc = this->wfn_->soccpi()[0];

 const int naocc = ca_occ->ncol();
 const int navir = ca_vir->ncol();
 const int nbocc = cb_occ->ncol();
 const int nbvir = cb_vir->ncol();

 const int naocc_mini= this->wfn_->nalpha();       
 const int nbocc_mini= this->wfn_->nbeta();        
 const int nbas_mini = this->calculate_nbas_mini_();
 const int navir_mini= nbas_mini - naocc_mini;     
 const int nbvir_mini= nbas_mini - nbocc_mini;     

 psi::outfile->Printf("\n QUAMBO orbital analysis:\n");
 psi::outfile->Printf  (" Na_occ = %5d    Na_occ_mini = %5d\n", naocc, naocc_mini);
 psi::outfile->Printf  (" Nb_occ = %5d    Nb_occ_mini = %5d\n", nbocc, nbocc_mini);
 psi::outfile->Printf  (" Na_vir = %5d    Na_vir_mini = %5d\n", navir, navir_mini);
 psi::outfile->Printf  (" Nb_vir = %5d    Nb_vir_mini = %5d\n", nbvir, nbvir_mini);
 psi::outfile->Printf  (" There are %5d minimal basis molecular orbitals (QUAMBOs) for alpha and beta spin\n", nbas_mini);
 psi::outfile->Printf  (" Looking for %5d virtual ALPHA valence orbitals (VVOs)\n"  , navir_mini);
 psi::outfile->Printf  (" Looking for %5d virtual BETA  valence orbitals (VVOs)\n\n", nbvir_mini);

 // [2] Run ROHF for all free atoms
 std::vector<psi::SharedMatrix> a_occ_a_list, a_vir_a_list, a_occ_b_list, a_vir_b_list;

 std::vector<psi::SharedMolecule> atoms = this->atomize_();

 std::string reference_stash = this->options_.get_str("REFERENCE");
 double e_conver_stash = this->options_.get_double("E_CONVERGENCE");
 double d_conver_stash = this->options_.get_double("D_CONVERGENCE");
 int maxiter_stash = this->options_.get_int("MAXITER");
 this->options_.set_global_str("REFERENCE", "ROHF");
 this->options_.set_global_double("E_CONVERGENCE", 1.0E-05);
 this->options_.set_global_double("D_CONVERGENCE", 1.0E-05);
 this->options_.set_global_int("MAXITER", 200);

 for (int i=0; i<this->mol_->natom(); ++i) {
      // Solve ROHF for free atom
      psi::SharedMolecule free_atom = atoms[i];
    //psi::PSIOManager::shared_object()->psiclean();//TODO

      // Extract A* orbitals (free-atom minimal basis occupied valence+core orbitals)
      psi::SharedBasisSet atom_primary, atom_auxiliary;
      if (this->acbs) {
          atom_primary = oepdev::create_basisset_by_copy(this->wfn_->basisset(), free_atom);
          atom_auxiliary = oepdev::create_basisset_by_copy(this->wfn_->get_basisset("DF_BASIS_SCF"), free_atom);
      } else {
          atom_primary = oepdev::create_atom_basisset_by_copy(this->wfn_->basisset(), free_atom, i);
          atom_auxiliary = oepdev::create_atom_basisset_by_copy(this->wfn_->get_basisset("DF_BASIS_SCF"), free_atom, i);
      }

    //bool compute_mints = false;
    //if (compute_mints) {
    //   std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(atom_primary);
    //   mints->integrals();
    // }

      psi::SharedWavefunction scf_base = std::make_shared<psi::Wavefunction>(free_atom, atom_primary, this->options_);
      scf_base->set_basisset("DF_BASIS_SCF", atom_auxiliary);
      std::shared_ptr<psi::SuperFunctional> functional = oepdev::create_superfunctional("HF", this->options_);
      std::shared_ptr<psi::scf::ROHF> scf = std::make_shared<psi::scf::ROHF>(scf_base, functional);
      scf->compute_energy();
      psi::SharedMatrix ca_a = scf->Ca_subset("AO", "OCC");
    //psi::PSIOManager::shared_object()->psiclean();//TODO

      // Compute a* coefficients (projections of A* onto molecule's occupied and virtual orbitals)
      psi::SharedMatrix a_occ_a;
      psi::SharedMatrix a_vir_a;
      psi::SharedMatrix a_occ_b;
      psi::SharedMatrix a_vir_b;

      // ACBS mode
      if (this->acbs) {
          a_occ_a = psi::Matrix::triplet(ca_occ, this->Sao_, ca_a, true, false, false); 
          a_vir_a = psi::Matrix::triplet(ca_vir, this->Sao_, ca_a, true, false, false);
          a_occ_b = psi::Matrix::triplet(cb_occ, this->Sao_, ca_a, true, false, false);
          a_vir_b = psi::Matrix::triplet(cb_vir, this->Sao_, ca_a, true, false, false);
      // MCBS mode
      } else {
          std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(this->wfn_->basisset());
          psi::SharedMatrix S_a = mints->ao_overlap(this->wfn_->basisset(), atom_primary);
          a_occ_a = psi::Matrix::triplet(ca_occ, S_a, ca_a, true, false, false); 
          a_vir_a = psi::Matrix::triplet(ca_vir, S_a, ca_a, true, false, false);
          a_occ_b = psi::Matrix::triplet(cb_occ, S_a, ca_a, true, false, false);
          a_vir_b = psi::Matrix::triplet(cb_vir, S_a, ca_a, true, false, false);
      }

      // Append to the lists
      a_occ_a_list.push_back(a_occ_a);
      a_vir_a_list.push_back(a_vir_a);
      a_occ_b_list.push_back(a_occ_b);
      a_vir_b_list.push_back(a_vir_b);
 }
 psi::SharedMatrix a_occ_a_set = psi::Matrix::horzcat(a_occ_a_list);
 psi::SharedMatrix a_vir_a_set = psi::Matrix::horzcat(a_vir_a_list);
 psi::SharedMatrix a_occ_b_set = psi::Matrix::horzcat(a_occ_b_list);
 psi::SharedMatrix a_vir_b_set = psi::Matrix::horzcat(a_vir_b_list);

 // [3] Compute QUAMBOs 
 SharedQUAMBOData q_a = 
 this->compute_quambo_data_(ca_occ, ca_vir, eps_a_occ, eps_a_vir, Fa, a_occ_a_set, a_vir_a_set, naocc_mini, navir_mini, "ALPHA");

 SharedQUAMBOData q_b = 
 this->compute_quambo_data_(cb_occ, cb_vir, eps_b_occ, eps_b_vir, Fb, a_occ_b_set, a_vir_b_set, nbocc_mini, nbvir_mini, "BETA");

 // [4] Save
 this->quambo_a_nonorthogonal_ = q_a->quambo_nonorthogonal;
 this->quambo_b_nonorthogonal_ = q_b->quambo_nonorthogonal;
 this->quambo_a_orthogonal_    = q_a->quambo_orthogonal; 
 this->quambo_b_orthogonal_    = q_b->quambo_orthogonal;

 this->c_a_mini_vir_ = q_a->c_mini_vir;
 this->c_b_mini_vir_ = q_b->c_mini_vir;

 this->e_a_mini_vir_ = q_a->e_mini_vir;
 this->e_b_mini_vir_ = q_b->e_mini_vir;

 this->c_a_mini_ = q_a->c_mini;
 this->c_b_mini_ = q_b->c_mini;

 this->e_a_mini_ = q_a->e_mini;
 this->e_b_mini_ = q_b->e_mini;

 this->nbas_mini_ = nbas_mini;
 this->naocc_mini_= naocc_mini;
 this->nbocc_mini_= nbocc_mini;
 this->navir_mini_= navir_mini;
 this->nbvir_mini_= nbvir_mini;
 this->nbf_       = this->wfn_->basisset()->nbf();

 // [5] Restore initial options
 this->options_.set_global_str("REFERENCE", reference_stash);
 this->options_.set_global_double("E_CONVERGENCE", e_conver_stash);
 this->options_.set_global_double("D_CONVERGENCE", d_conver_stash);
 this->options_.set_global_int("MAXITER", maxiter_stash);

 psi::outfile->Printf(" @QUAMBO: Done.\n");

}

int QUAMBO::calculate_nbas_mini_(void) {
  int nbf = 0;
  for (int i=0; i<this->mol_->natom(); ++i) {
       std::string atom = this->mol_->symbol(i);
       int nb = this->nbas_atom_mini_.at(atom);
       nbf += nb;
  }
  return nbf;
}

double QUAMBO::compute_error_between_two_vectors_(psi::SharedVector a, psi::SharedVector b) {
 double error = 0.0;
 psi::SharedVector e = std::make_shared<psi::Vector>(*a);
 e->subtract(b);
 error = e->rms();
 return error;
}

std::vector<psi::SharedMolecule> QUAMBO::atomize_(void) {
 std::vector<psi::SharedMolecule> atoms;

 std::string log;
 std::string s;
 std::string ghost_prefix = "@";
 std::string endline = "symmetry c1\nunits bohr\nno_reorient\nno_com\n";

 int m;

 if (this->acbs) {
 
     int I = 0;                                                                                       
     for (int i=0; i<this->mol_->natom(); ++i) {
          std::string coord_log = "";
          for (int j=0; j<this->mol_->natom(); ++j) {
               double x = this->mol_->x(j);
               double y = this->mol_->y(j);
               double z = this->mol_->z(j);
               if (i==j) { 
                   s = this->mol_->symbol(j);
                   m = this->unpe_atom_.at(s);
               }
               else {
                   s = ghost_prefix + this->mol_->symbol(j);
               }
               coord_log += oepdev::string_sprintf("%4s %14.8f %14.8f %14.8f\n", s.c_str(), x, y, z);
          }
          log  = oepdev::string_sprintf("\n0 %i\n", m+1);
          log += coord_log + endline;
                                                                                                      
          cout << log;
          psi::SharedMolecule atom = psi::Molecule::create_molecule_from_string(log);
          atom->set_molecular_charge(0);
          atom->set_multiplicity(m+1);
          atom->update_geometry();
          atoms.push_back(atom);
                                                                                                      
          I += 1;
     }
 } else {
     // MCBS mode
     for (int i=0; i<this->mol_->natom(); ++i) {
          double x = this->mol_->x(i);
          double y = this->mol_->y(i);
          double z = this->mol_->z(i);
          s = this->mol_->symbol(i);
          m = this->unpe_atom_.at(s);

          log  = oepdev::string_sprintf("\n0 %i\n", m+1);
          log += oepdev::string_sprintf("%4s %14.8f %14.8f %14.8f\n", s.c_str(), x, y, z);
          log += endline;

          cout << log;

          psi::SharedMolecule atom = psi::Molecule::create_molecule_from_string(log);
          atom->set_molecular_charge(0);
          atom->set_multiplicity(m+1);
          atom->update_geometry();
          atoms.push_back(atom);
     }

 }
 return atoms;
}

SharedQUAMBOData QUAMBO::compute_quambo_data_(
     psi::SharedMatrix c_occ,     // Occupied MOs of the molecule
     psi::SharedMatrix c_vir,     // Virtual MOs of the molecule
     psi::SharedVector eps_occ,   // Occupied MO energies of the molecule
     psi::SharedVector eps_vir,   // Virtual MO energies of the molecule
     psi::SharedMatrix F,         // Fock matrix in non-orthogonal AO basis of the molecule
     psi::SharedMatrix a_occ_set, // The a_{nj}^* coefficients from Ref. [1] (projections of free-atom orbitals on c_occ)
     psi::SharedMatrix a_vir_set, // The a_{vj}^* coefficients from Ref. [1] (projections of free-atom orbitals on c_vir)
     int nocc_mini,               // Number of occupied orbitals in the minimal QUAMBO basis
     int nvir_mini,               // Number of virtual orbitals (VVOs) in the minimal QUAMBO basis
     std::string label            // Label of the orbitals
    )
{
  // [3] Compute B matrix
  psi::SharedMatrix B = psi::Matrix::doublet(a_vir_set, a_vir_set, false, true);
  const int dim = B->ncol(); // number of all Virtual MOs
  const int nquambo= a_vir_set->ncol(); // number of QUAMBOs
  const int nbf = c_occ->nrow(); // number of AOs
 
  // [4] Diagonalize to find transformation matrix T 
  psi::SharedMatrix U = std::make_shared<psi::Matrix>("", dim, dim);
  psi::SharedVector E = std::make_shared<psi::Vector>("", dim);
  B->diagonalize(U, E, psi::diagonalize_order::descending);
  psi::SharedMatrix T = std::make_shared<psi::Matrix>("", dim, nvir_mini);
  for (int i=0; i<nvir_mini; ++i) {
    psi::SharedVector col = U->get_column(0, i);
    T->set_column(0, i, col);
  }
  
  // [5] Calculate R matrix 
  psi::SharedMatrix R = psi::Matrix::doublet(T, T, false, true);

  // [6] Calculate normalization
  psi::SharedVector Dj = std::make_shared<psi::Vector>("", nquambo);
  psi::SharedMatrix r = psi::Matrix::doublet(R, a_vir_set, false, false);
  for (int i=0; i<nquambo; ++i) {
       psi::SharedVector aocc_i = a_occ_set->get_column(0, i);
       psi::SharedVector avir_i = a_vir_set->get_column(0, i);
       psi::SharedVector     ri =         r->get_column(0, i);
       double dj = aocc_i->vector_dot(aocc_i) + avir_i->vector_dot(ri);
       if (dj < 0.0) dj = 1.0;
       dj = 1.0/sqrt(dj);
       Dj->set(i, dj);
  }
  
  // [7] Calculate non-orthogonal QUAMBOs 
  psi::SharedMatrix a_occ_set_nonorthogonal = a_occ_set->clone();
  psi::SharedMatrix a_vir_set_nonorthogonal =         r->clone();
  for (int i=0; i<nquambo; ++i) {
       a_occ_set_nonorthogonal->scale_column(0, i, Dj->get(i)); // a_nj from Ref. [1]
       a_vir_set_nonorthogonal->scale_column(0, i, Dj->get(i)); // a_vj from Ref. [1]
  }
  psi::SharedMatrix quambo_nonorthogonal = std::make_shared<psi::Matrix>("", nbf, nquambo);
  quambo_nonorthogonal->gemm(false, false, 1.0, c_occ, a_occ_set_nonorthogonal, 1.0);
  quambo_nonorthogonal->gemm(false, false, 1.0, c_vir, a_vir_set_nonorthogonal, 1.0);

  // [8] Calculate orthogonal QUAMBOs
  psi::SharedMatrix S = std::make_shared<psi::Matrix>("", nquambo, nquambo);
  S->gemm(true, false, 1.0, a_occ_set_nonorthogonal, a_occ_set_nonorthogonal, 1.0);
  S->gemm(true, false, 1.0, a_vir_set_nonorthogonal, a_vir_set_nonorthogonal, 1.0);
  psi::SharedMatrix Sm12 = S->clone(); Sm12->power(-0.5);
  psi::SharedMatrix a_occ_set_orthogonal = psi::Matrix::doublet(a_occ_set_nonorthogonal, Sm12, false, false);// a_nj' from Ref. [1]
  psi::SharedMatrix a_vir_set_orthogonal = psi::Matrix::doublet(a_vir_set_nonorthogonal, Sm12, false, false);// a_vj' from Ref. [1]
  psi::SharedMatrix quambo_orthogonal = std::make_shared<psi::Matrix>("", nbf, nquambo);
  quambo_orthogonal->gemm(false, false, 1.0, c_occ, a_occ_set_orthogonal, 1.0);
  quambo_orthogonal->gemm(false, false, 1.0, c_vir, a_vir_set_orthogonal, 1.0);

  // [9] Diagonalize Fock matrix in QUAMBO basis
  psi::SharedMatrix quambo = quambo_orthogonal;
  psi::SharedMatrix F_quambo = psi::Matrix::triplet(quambo, F, quambo, true, false, false);
  psi::SharedVector e_quambo = std::make_shared<psi::Vector>("", nquambo);
  psi::SharedMatrix c_quambo = std::make_shared<psi::Matrix>("", nquambo, nquambo);
  F_quambo->diagonalize(c_quambo, e_quambo, psi::diagonalize_order::ascending);

  // [9a] Test whether occupied orbital energies are correctly obtained
  psi::SharedVector e_quambo_occ = std::make_shared<psi::Vector>("", nocc_mini); 
  for (int i=0; i<nocc_mini; ++i) {
    e_quambo_occ->set(i, e_quambo->get(i));
  }
  double error = this->compute_error_between_two_vectors_(e_quambo_occ, eps_occ);
  if (error > this->options_.get_double("QUAMBO_EPS_THRESHOLD_CHECK")) { 
     psi::outfile->Printf("Error in QUAMBO calculations = %14.5f", error);
     std::cout << "Error in QUAMBO calculations = " << error << std::endl;
     throw psi::PSIEXCEPTION("Error in QUAMBO calculations!"); 
  }

  // [10] Compute virtual minimal MOs (VMO)
  psi::SharedVector e_quambo_vir = std::make_shared<psi::Vector>("", nvir_mini); 
  psi::SharedMatrix vir_mini_mo = std::make_shared<psi::Matrix>("", nquambo, nvir_mini); 
  for (int i=0; i<nvir_mini; ++i) {
    e_quambo_vir->set(i, e_quambo->get(nocc_mini+i));
    psi::SharedVector col = c_quambo->get_column(0, nocc_mini+i);
    vir_mini_mo->set_column(0, i, col);
  }
  psi::SharedMatrix vir_mini = psi::Matrix::doublet(quambo, vir_mini_mo, false, false);
  psi::SharedMatrix vir      = psi::Matrix::doublet(quambo, c_quambo, false, false);

  psi::outfile->Printf("\n Eigenvalues of %6s Fock Operator in QUAMBO Basis\n\n", label.c_str());
  psi::outfile->Printf  (" Occupied\n");
  psi::outfile->Printf  (" --------\n\n");
  psi::outfile->Printf  ("          Full HF       QUAMBO        Error\n");
  for (int i=0; i<nocc_mini; ++i) {
      psi::outfile->Printf(" %3d %13.6f %13.6f %13.6f\n", i+1, 
         eps_occ->get(i), e_quambo_occ->get(i), e_quambo_occ->get(i)-eps_occ->get(i));
  }
  psi::outfile->Printf("\n Virtual\n");
  psi::outfile->Printf  (" -------\n\n");
  psi::outfile->Printf  ("          QUAMBO\n");
  for (int i=0; i<nvir_mini; ++i) {
      psi::outfile->Printf(" %3d %13.6f\n", i+1, e_quambo_vir->get(i));
  }
  psi::outfile->Printf("\n");


  // Save And Return
  SharedQUAMBOData data = std::make_shared<QUAMBOData>();
  data->quambo_nonorthogonal = quambo_nonorthogonal;
  data->quambo_orthogonal = quambo_orthogonal;
  data->c_mini_vir = vir_mini;
  data->e_mini_vir = e_quambo_vir;
  data->c_mini = vir;
  data->e_mini = e_quambo;

  return data;
}

psi::SharedMatrix QUAMBO::quambo(const std::string& spin, const std::string& type) {
  psi::SharedMatrix q;
  std::string name;
  if (spin == "ALPHA") {
      if (type == "ORTHOGONAL") { 
          name = "QUAMBOs Alpha, orthogonal";
          q = this->quambo_a_orthogonal_->clone();
      } else
      if (type == "NONORTHOGONAL") {
          name = "QUAMBOs Alpha, non-orthogonal";
          q = this->quambo_a_nonorthogonal_->clone();
      } else {throw psi::PSIEXCEPTION("Unknown QUAMBO type chosen! Only ORTHOGONAL and NONORTHOGONAL keywords work.");} 
  } else 
  if (spin == "BETA") {
      if (type == "ORTHOGONAL") { 
          name = "QUAMBOs Beta, orthogonal";
          q = this->quambo_b_orthogonal_->clone();
      } else
      if (type == "NONORTHOGONAL") {
          name = "QUAMBOs Beta, non-orthogonal";
          q = this->quambo_b_nonorthogonal_->clone();
      } else {throw psi::PSIEXCEPTION("Unknown QUAMBO type chosen! Only ORTHOGONAL and NONORTHOGONAL keywords work.");} 
  } else {throw psi::PSIEXCEPTION("Unknown electron spin chosen for QUAMBO! Only ALPHA and BETA keywords work.");}
  q->set_name(name);
  return q;
}

psi::SharedVector QUAMBO::epsilon_a_subset(const std::string& space, const std::string& subset) {
  psi::SharedVector C = this->epsilon_subset_helper_(this->e_a_mini_, "Alpha", this->naocc_mini_, space, subset);
  return C;
}

psi::SharedVector QUAMBO::epsilon_b_subset(const std::string& space, const std::string& subset) {
  psi::SharedVector C = this->epsilon_subset_helper_(this->e_b_mini_, "Beta", this->nbocc_mini_, space, subset);
  return C;
}

psi::SharedMatrix QUAMBO::Ca_subset(const std::string& space, const std::string& subset) {
  psi::SharedMatrix C = this->C_subset_helper_(this->c_a_mini_, "Alpha", this->naocc_mini_, space, subset);
  return C;
}

psi::SharedMatrix QUAMBO::Cb_subset(const std::string& space, const std::string& subset) {
  psi::SharedMatrix C = this->C_subset_helper_(this->c_b_mini_, "Beta", this->nbocc_mini_, space, subset);
  return C;
}

psi::SharedVector QUAMBO::epsilon_subset_helper_(
  psi::SharedVector full_eps, const std::string& label, const int& nocc, const std::string& space, const std::string& subset)
{
  if (space != "MO") throw psi::PSIEXCEPTION("Only MO space is now supported!");

  psi::SharedVector C;
  std::string name;
  if (subset == "OCC") {
      name = oepdev::string_sprintf("Occupied %6s MO energies from diagonalization of Fock matrix in QUAMBO basis", label.c_str() );
      C = std::make_shared<psi::Vector>(label, nocc);
      for (int i=0; i<nocc; ++i) {
           C->set(i, full_eps->get(i));
      }
  } else 
  if (subset == "VIR") {
      name = oepdev::string_sprintf("Virtual %6s MO energies from diagonalization of Fock matrix in QUAMBO basis", label.c_str() );
      C = std::make_shared<psi::Vector>(label, this->nbas_mini_-nocc);
      for (int i=0; i<this->nbas_mini_-nocc; ++i) {
           C->set(i, full_eps->get(nocc+i));
      }
  } else 
  if (subset == "ALL") {
    C = std::make_shared<psi::Vector>(*full_eps);
  } else {throw psi::PSIEXCEPTION("Wrong subset chosen!");}  
  C->set_name(name);
  return C;
}
psi::SharedMatrix QUAMBO::C_subset_helper_(
  psi::SharedMatrix full_C, const std::string& label, const int& nocc, const std::string& space, const std::string& subset) 
{
  if (space != "AO") throw psi::PSIEXCEPTION("Only AO space is now supported!");

  psi::SharedMatrix C;
  std::string name;
  if (subset == "OCC") {
      name = oepdev::string_sprintf("Occupied %6s MOs from diagonalization of Fock matrix in QUAMBO basis", label.c_str() );
      C = std::make_shared<psi::Matrix>(label, this->nbf_, nocc);
      for (int i=0; i<nocc; ++i) {
           psi::SharedVector ci = full_C->get_column(0, i);
           C->set_column(0, i, ci);
      }
  } else 
  if (subset == "VIR") {
      name = oepdev::string_sprintf("Virtual %6s MOs from diagonalization of Fock matrix in QUAMBO basis", label.c_str() );
      C = std::make_shared<psi::Matrix>(label, this->nbf_, this->nbas_mini_-nocc);
      for (int i=0; i<this->nbas_mini_-nocc; ++i) {
           psi::SharedVector ci = full_C->get_column(0, nocc+i);
           C->set_column(0, i, ci);
      }
  } else 
  if (subset == "ALL") {
    C = full_C->clone();
  } else {throw psi::PSIEXCEPTION("Wrong subset chosen!");}  
  C->set_name(name);
  return C;
}



} // EndNameSpace oepdev
