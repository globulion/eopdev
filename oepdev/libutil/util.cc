#include "util.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/mintshelper.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libiwl/iwl.hpp"


namespace oepdev{

extern "C" PSI_API
void preambule(void) {
      outfile->Printf("                                                                             \n");
      outfile->Printf("    -------------------------------------------------------------------------\n");
      outfile->Printf("          OepDev: One-Electron Effective Potentials Development Routines     \n");
      outfile->Printf("                               OepDev 0.2 (pre-release)                      \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("                         Git: Rev {203ae6c}                                  \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("    Bartosz Błasiak                                                          \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("    References:                                                              \n");
      outfile->Printf("     [1] B. Błasiak, J. Chem. Phys. 149, 16 (2018), 164115                   \n");
      outfile->Printf("     [2] B. Błasiak, J. D. Bednarska, M. Chołuj, W. Bartkowiak, (2020),      \n");
      outfile->Printf("                   arXiv:2002.00766                                          \n");
      outfile->Printf("     [3] B. Błasiak, R. W. Góra, W. Bartkowiak, (2020), arXiv:2002.00778     \n");
      outfile->Printf("    -------------------------------------------------------------------------\n");
      outfile->Printf("                                                                             \n");
}

extern "C" PSI_API
std::shared_ptr<SuperFunctional> 
create_superfunctional(std::string name, Options& options) {
  std::shared_ptr<SuperFunctional> functional = SuperFunctional::blank();
  if (name==(std::string)"HF") {
      int npoints = options.get_int("DFT_BLOCK_MAX_POINTS");                   
      int deriv   = 1;
      functional->set_max_points(npoints);
      functional->set_deriv(deriv);
      functional->set_name("HF");
      functional->set_description("    Equivalent to the Hartree-Fock-Roothaan-Hall theory.)\n");
      functional->set_citation("   \n");
      functional->set_x_alpha(1.0);
      functional->set_c_omega(0.0);
      functional->set_x_omega(0.0);
      functional->set_c_alpha(0.0);
  }
  else {
    throw PSIEXCEPTION("Only HF wavefunctions are now supported for OEPDEV plugin!");
  }
  return functional;
}

extern "C" PSI_API
SharedBasisSet
create_basisset_by_copy(SharedBasisSet basis_ref, SharedMolecule molecule_target) {

  std::map<std::string, std::map<std::string, std::vector<psi::ShellInfo>>> shellmap;
  std::map<std::string, std::map<std::string, std::vector<psi::ShellInfo>>> ecp_shellmap;

  std::string name = basis_ref->name();
  std::string key  = basis_ref->key();
  std::string label= name; // ? this is just assumed, but should be a 'blend'

  molecule_target->set_basis_all_atoms(name, key);

  //cout << molecule_target->natom() << " " << basis_ref->molecule()->natom() << " " << endl;
  if (molecule_target->nallatom() != basis_ref->molecule()->natom()) throw psi::PSIEXCEPTION("ERROR: Copy basis has different number of atoms than target molecule!!!");
  for (int a=0; a<molecule_target->nallatom(); ++a) {
       std::vector<psi::ShellInfo> shellinfos;
       int is = basis_ref->shell_on_center(a, 0);
       int in = basis_ref->nshell_on_center(a);
       for (int iis = 0; iis < in; ++iis) {
            int si = is + iis;
            const psi::GaussianShell& s = basis_ref->shell(si);

            std::vector<double> c, e;
            for (int pi=0; pi<s.nprimitive(); ++pi) {
                 c.push_back(s.original_coef(pi));
                 e.push_back(s.exp (pi));
            }

            psi::ShellInfo infos(s.am(), c, e, psi::GaussianType::Cartesian, psi::PrimitiveType::Unnormalized);
            shellinfos.push_back(infos);
       }

       std::string atomlabel = molecule_target->fsymbol(a);
       std::string hash = ""; // ? it is a comment of an atom I guess...

       molecule_target->set_shell_by_label(atomlabel, hash, key);

       shellmap[name][atomlabel] = shellinfos;
  }

  molecule_target->update_geometry();

  SharedBasisSet basis_target = std::make_shared<BasisSet>(key, molecule_target, shellmap, ecp_shellmap);

  basis_target->set_name(name);
  basis_target->set_key(key);
  basis_target->set_target(label);

  return basis_target;
}



extern "C" PSI_API
std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id) {
    std::vector<int> real_list; real_list.push_back(id-1);
    std::vector<int> ghost_list;
    return molecule_dimer->extract_subsets(real_list, ghost_list);
}

extern "C" PSI_API
double compute_distance(psi::SharedVector v1,
                        psi::SharedVector v2)
{
  if (v1->dim() != v2->dim()) throw psi::PSIEXCEPTION("Vectors have different lengths!");
  double r = 0.0;
  for (int i=0; i<v1->dim(); ++i) r += pow(v2->get(i) - v1->get(i), 2.0);
  r = sqrt(r);
  return r;
}

extern "C" PSI_API
std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
	  std::shared_ptr<BasisSet> auxiliary,
	  std::shared_ptr<BasisSet> guess,
          std::shared_ptr<SuperFunctional> functional,
          Options& options,
          std::shared_ptr<PSIO> psio, 
	  bool compute_mints) 
{
   /* 
      Important Note:

         compute_mints temporarily denotes psiclean routine to be performed between basis_guess and primary scf runs
         to solve SCF. This is used for WavefunctionUnions created via second constructor (by providing externally 
         wavefunctions for monomers). This is only a temporary solution, because this issue needs to be resolved in
         the future at the psi4 scratch file management level.

    */

    psi::SharedWavefunction scf;
    //psi::PSIOManager::shared_object()->psiclean();//TODO

    // Guess: Hcore in guess basis --> projection to primary basis
    if (options.get_bool("OEPDEV_BASIS_GUESS")==true) {

       // ===> Step 1: Guess SCF <=== //
       outfile->Printf("\n @solve_scf: Starting SCF in Guess Basis... \n\n");

       //compute_mints = true; 
       //if (compute_mints) {
       //  //psi::PSIOManager::shared_object()->psiclean();//TODO
       //  std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(guess);
       //  mints->integrals();
       //}
       SharedWavefunction scf_base_guess(new Wavefunction(molecule, guess, options));                     
       //bool opt_stash = options.get_bool("DF_SCF_GUESS");
       //options.set_bool("SCF", "DF_SCF_GUESS", false);
       scf_base_guess->set_basisset("DF_BASIS_SCF", auxiliary);//TODO -> perhaps better use smaller basis adequate for 3-21G
      
       SharedWavefunction scf_guess = std::make_shared<psi::scf::RHF>(scf_base_guess, functional, options, psio);
       scf_guess->compute_energy();

       // ===> Step 2: Basis set projection
       outfile->Printf("\n @solve_scf: Basis Set Projection to Target Basis... \n\n");
       psi::SharedMatrix pCa = scf_guess->basis_projection(scf_guess->Ca_subset("AO","OCC"),scf_guess->nalphapi(),
                                                           guess, primary);
       psi::SharedMatrix pCb = scf_guess->basis_projection(scf_guess->Cb_subset("AO","OCC"),scf_guess->nbetapi(),
                                                           guess, primary);
       //options.set_bool("SCF", "DF_SCF_GUESS", opt_stash);
     //psi::PSIOManager::shared_object()->print();
       if (compute_mints) psi::PSIOManager::shared_object()->psiclean();//TODO --> this needs to be set when calling from Python level! ??? --> i.e., when using second constructor of WavefunctionUnion
     //std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(primary); mints->one_electron_integrals();
    //IWL ERIOUT(psio.get(), PSIF_SO_TEI, 0.0, 0, 0);
    //ERIOUT.flush(1);
    //ERIOUT.set_keep_flag(true);
    //ERIOUT.close();

       // ===> Step 3: Target SCF <=== //
       outfile->Printf("\n @solve_scf: Starting SCF in Target Basis... \n\n");

       //if (compute_mints) {
       //  //psi::PSIOManager::shared_object()->psiclean();//TODO
       //  std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(primary);
       //  mints->integrals();
       //}
       SharedWavefunction scf_base(new Wavefunction(molecule, primary, options));                     
       scf_base->set_basisset("DF_BASIS_SCF", auxiliary);
                                                                                                     
       std::shared_ptr<psi::scf::RHF> scf_ = std::make_shared<psi::scf::RHF>(scf_base, functional, options, psio);
       scf_->guess_Ca(pCa);
       scf_->guess_Cb(pCb);
       scf_->compute_energy();
       scf = scf_;


    // Guess: Hcore in primary basis
    } else {
       outfile->Printf("\n @solve_scf: Starting SCF in Target Basis... \n\n");
       SharedWavefunction scf_base(new Wavefunction(molecule, primary, options));                     
       scf_base->set_basisset("DF_BASIS_SCF", auxiliary);
                                                                                                      
       // Compute integrals (write IWL entry to PSIO)
       //if (compute_mints) {
       //  std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(primary);
       //  mints->integrals();
       //}
       scf = std::make_shared<psi::scf::RHF>(scf_base, functional, options, psio);
       scf->compute_energy();


    }
    outfile->Printf("\n @solve_scf: Done. \n\n");
    //psi::PSIOManager::shared_object()->psiclean();//TODO

    return scf;
}

extern "C" PSI_API
std::shared_ptr<Wavefunction>
solve_scf_sad(std::shared_ptr<Molecule> molecule, 
              std::shared_ptr<BasisSet> primary,              
	      std::shared_ptr<BasisSet> auxiliary,
              std::vector<std::shared_ptr<BasisSet>> sad,
              std::vector<std::shared_ptr<BasisSet>> sad_fit,
              std::shared_ptr<SuperFunctional> functional,
              Options& options,
              std::shared_ptr<PSIO> psio,
	      bool compute_mints) 
{
    SharedWavefunction scf_base(new Wavefunction(molecule, primary, options));
    scf_base->set_basisset("DF_BASIS_SCF", auxiliary);

    // Compute integrals (write IWL entry to PSIO)
    if (compute_mints) {
      std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(primary);
      mints->integrals();
    }
    std::shared_ptr<psi::scf::RHF> hf = std::make_shared<psi::scf::RHF>(scf_base, functional, options, psio);
    if (!sad    .empty()) hf->set_sad_basissets(sad);
    if (!sad_fit.empty()) hf->set_sad_fitting_basissets(sad_fit);
    hf->compute_energy();
    SharedWavefunction scf = hf;
    return scf;
}


extern "C" PSI_API
double average_moment(std::shared_ptr<psi::Vector> moment)
{
  const int l = moment->dim();
  if (l == 3) { // Dipole
     double da_x = moment->get(0);
     double da_y = moment->get(1);
     double da_z = moment->get(2);
     return sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
  } else if (l == 6) { // Quadrupole
     double qa_xx= moment->get(0);
     double qa_xy= moment->get(1);
     double qa_xz= moment->get(2);
     double qa_yy= moment->get(3);
     double qa_yz= moment->get(4);
     double qa_zz= moment->get(5);
     double ta = (1.0 / 3.0) * (qa_xx + qa_yy + qa_zz);
     return sqrt((qa_zz-ta)*(qa_zz-ta) + (4.0/3.0)*(qa_xy*qa_xy+qa_xz*qa_xz+qa_yz*qa_yz) 
                                                  + (1.0/3.0)*(pow((qa_xx-ta) - (qa_yy-ta), 2.0)));
  } else {throw PSIEXCEPTION("Wrong size of multipole moment vector!");}
}

extern "C" PSI_API
std::vector<std::shared_ptr<psi::Matrix>> calculate_JK(std::shared_ptr<psi::Wavefunction> wfn, std::shared_ptr<psi::Matrix> C){

  // Initialize the J_ij and K_ij matrix
  int n = C->ncol();
  std::shared_ptr<psi::Matrix> Jij = std::make_shared<psi::Matrix>("Coulomb Integrals in MO basis", n, n);
  std::shared_ptr<psi::Matrix> Kij = std::make_shared<psi::Matrix>("Exchange Integrals in MO basis", n, n);
  double** pJij = Jij->pointer();
  double** pKij = Kij->pointer();

  std::vector<std::shared_ptr<psi::Matrix>> JK;
  JK.push_back(Jij);
  JK.push_back(Kij);

  // Compute ERI's
  std::shared_ptr<psi::MintsHelper> mints = std::make_shared<psi::MintsHelper>(wfn->basisset());
  mints->integrals();

  // Transform ERI's
  std::shared_ptr<psi::MOSpace> space = psi::MOSpace::all;
  std::vector<std::shared_ptr<psi::MOSpace>> spaces; spaces.push_back(space);

  psi::IntegralTransform tr(wfn, spaces,
                  psi::IntegralTransform::TransformationType::Restricted,
                  psi::IntegralTransform::OutputType::DPDOnly,
                  psi::IntegralTransform::MOOrdering::QTOrder,
                  psi::IntegralTransform::FrozenOrbitals::None);

  tr.set_orbitals(C);
  tr.transform_tei(space, space, space, space);

  // Read integrals and save
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();

  dpd_set_default(tr.get_dpd_id());
  dpdbuf4 buf;
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                         tr.DPD_ID("[A,A]"  ), tr.DPD_ID("[A,A]"  ),
                         tr.DPD_ID("[A>=A]+"), tr.DPD_ID("[A>=A]+"  ), 0, "MO Ints (AA|AA)");

  for (int h = 0; h < wfn->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
       for (int pq = 0; pq < buf.params->rowtot[h]; ++pq) {
            int p = buf.params->roworb[h][pq][0];
            int q = buf.params->roworb[h][pq][1];
            for (int rs = 0; rs < buf.params->coltot[h]; ++rs) {
         	   int r = buf.params->colorb[h][rs][0];
        	   int s = buf.params->colorb[h][rs][1];
		   /* J */
		   if ((p==q) && (r==s)) {
		       pJij[p][r] = buf.matrix[h][pq][rs];
	           }
		   /* K */
        	   if ((p==r) && (q==s)) {
        	       pKij[p][q] = buf.matrix[h][pq][rs];
      	           }
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  //Kij->print();
  return JK;
}

extern "C" PSI_API
std::vector<std::shared_ptr<psi::Matrix>> calculate_JK_r(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, std::shared_ptr<psi::Matrix> Dij){

  // Initialize the J_ij and K_ij matrix
  int n = Dij->ncol();
  std::shared_ptr<psi::Matrix> Jij = std::make_shared<psi::Matrix>("Coulomb Integrals in MO basis", n, n);
  std::shared_ptr<psi::Matrix> Kij = std::make_shared<psi::Matrix>("Exchange Integrals in MO basis", n, n);
  double** pJij = Jij->pointer();
  double** pKij = Kij->pointer();
  double** pD   = Dij->pointer();

  std::vector<std::shared_ptr<psi::Matrix>> JK;
  JK.push_back(Jij);
  JK.push_back(Kij);

  // Read integrals and save
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();

  dpd_set_default(tr->get_dpd_id());
  dpdbuf4 buf;
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                         tr->DPD_ID("[A,A]"  ), tr->DPD_ID("[A,A]"  ),
                         tr->DPD_ID("[A>=A]+"), tr->DPD_ID("[A>=A]+"  ), 0, "MO Ints (AA|AA)");

  // J and K
  for (int h = 0; h < wfn->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
       for (int pq = 0; pq < buf.params->rowtot[h]; ++pq) {
            int p = buf.params->roworb[h][pq][0];
            int q = buf.params->roworb[h][pq][1];
	    double vj_pq = 0.0;
            for (int rs = 0; rs < buf.params->coltot[h]; ++rs) {
         	   int r = buf.params->colorb[h][rs][0];
        	   int s = buf.params->colorb[h][rs][1];
		   vj_pq      += buf.matrix[h][pq][rs] * pD[r][s];
		   pKij[p][r] += buf.matrix[h][pq][rs] * pD[q][s];
            }
	    pJij[p][q] = vj_pq;
       }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  //Kij->print();
  return JK;
}

extern "C" PSI_API
std::shared_ptr<psi::Matrix>
calculate_der_D(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, 
		std::shared_ptr<psi::Matrix> C,
		std::vector<std::shared_ptr<psi::Matrix>> A) {
  
  // Initialize derivatives
  int N = C->ncol(); /* MO-A */
  int M = C->nrow(); /* MO-B */
  //std::shared_ptr<psi::Matrix> Deriv = std::make_shared<psi::Matrix>("Derivative of E_XC", N, N);
  std::shared_ptr<psi::Matrix> T     = std::make_shared<psi::Matrix>("Temporary"         , N, M);
  double** pT   = T->pointer();
  double** pC   = C->pointer();
  std::vector<double**> vA;
  for (int m=0; m<M; ++m) {
       vA.push_back(A[m]->pointer());
  }

  // Read integrals and save
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();

  dpd_set_default(tr->get_dpd_id());
  dpdbuf4 buf;
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                         tr->DPD_ID("[A,A]"  ), tr->DPD_ID("[A,A]"  ),
                         tr->DPD_ID("[A>=A]+"), tr->DPD_ID("[A>=A]+"  ), 0, "MO Ints (AA|AA)");

  for (int h = 0; h < wfn->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
	for (int sb = 0; sb < buf.params->rowtot[h]; ++sb) {
	     int s = buf.params->roworb[h][sb][0];
	     int b = buf.params->roworb[h][sb][1];
	     for (int cd = 0; cd < buf.params->coltot[h]; ++cd) {
	          int c = buf.params->colorb[h][cd][0];
	          int d = buf.params->colorb[h][cd][1];
	          double sbcd = buf.matrix[h][sb][cd];
                  for (int m = 0; m < M; ++m) {
	              pT[s][m] += sbcd * pC[m][c] * vA[m][b][d];
	          }
             }
        }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  // Symmetrize and scale
  std::shared_ptr<psi::Matrix> Deriv = psi::Matrix::doublet(T, C, false, false);
  //std::shared_ptr<psi::Matrix> Deriv_= psi::Matrix::doublet(C, T, true, true);
  //Deriv->add(Deriv_);
  Deriv->add(Deriv->clone()->transpose());
  Deriv->scale(-1.0);
  return Deriv;
}

extern "C" PSI_API
double calculate_e_xc(std::shared_ptr<psi::Wavefunction> wfn, 
		std::shared_ptr<psi::IntegralTransform> tr, 
		std::shared_ptr<psi::Matrix> f,
		std::shared_ptr<psi::Matrix> C){


  // Initialize derivatives
  double E = 0.0;
  int M = C->nrow(); /* MO-SCF */
  int N = C->ncol(); /* MO-NEW */
  double** pf   = f->pointer();
  double** pC   = C->pointer();

  // Read integrals and save
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();

  dpd_set_default(tr->get_dpd_id());
  dpdbuf4 buf;
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                         tr->DPD_ID("[A,A]"  ), tr->DPD_ID("[A,A]"  ),
                         tr->DPD_ID("[A>=A]+"), tr->DPD_ID("[A>=A]+"  ), 0, "MO Ints (AA|AA)");

  for (int h = 0; h < wfn->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
	for (int ab = 0; ab < buf.params->rowtot[h]; ++ab) {
	     int a = buf.params->roworb[h][ab][0];
	     int b = buf.params->roworb[h][ab][1];
	     for (int cd = 0; cd < buf.params->coltot[h]; ++cd) {
	          int c = buf.params->colorb[h][cd][0];
	          int d = buf.params->colorb[h][cd][1];
	          double abcd = buf.matrix[h][ab][cd];
                  for (int i = 0; i < N; ++i) {
		  for (int j = 0; j < N; ++j) {
		       E += abcd * pC[a][i] * pC[b][j] * pC[c][i] * pC[d][j] * pf[i][j];
	          }
		  }
             }
        }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  return -E;
}

extern "C" PSI_API
std::shared_ptr<psi::Matrix>
matrix_power_derivative(std::shared_ptr<psi::Matrix> A,
              		double g, double step){

  double hi = 1.0/ (step * 2.0);
  //int n = A->ncol();
  //std::shared_ptr<psi::Matrix> D = std::make_shared<psi::Matrix>("",n,n);
  //double** d = D->pointer();

  std::shared_ptr<psi::Matrix> Ag = A->clone();
  Ag->power(g);

  std::shared_ptr<psi::Matrix> Ag1 = A->clone();
  std::shared_ptr<psi::Matrix> I = A->clone(); I->identity();
  I->scale(step * 2.0);
  Ag1->add(I);
  Ag1->power(g);
  Ag1->subtract(Ag);
  Ag1->scale(hi);

  //for (int i=0; i<n; ++i) {
  //  for (int j=0; j<=i; ++j) {
  //       std::shared_ptr<psi::Matrix> Ag1 = A->clone();
  //       double** a = Ag1->pointer();
  //       a[i][j] += step;
  //       a[j][i] += step;
  //       Ag1->power(g);
  //       Ag1->subtract(Ag);
  //       Ag1->scale(hi);
  //       double v = Ag1->trace();
  //       d[i][j] = v;
  //       d[j][i] = v;
  //  }
  //}
  return Ag1;
}

extern "C"
std::shared_ptr<psi::Matrix> _calculate_DFI_Vel(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::IntegralFactory> f_abab,
                std::shared_ptr<psi::Matrix> db)
{
  std::shared_ptr<psi::Matrix> V = std::make_shared<psi::Matrix>("", db->nrow(), db->ncol());
  double** d = db->pointer();
  double** v = V ->pointer();

  // J contribution
  std::shared_ptr<oepdev::ShellCombinationsIterator> s_aabb = oepdev::ShellCombinationsIterator::build(f_aabb, "ALL");
  std::shared_ptr<psi::TwoBodyAOInt> t_aabb(f_aabb->eri()); const double* b_aabb = t_aabb->buffer();

  for (s_aabb->first(); s_aabb->is_done() == false; s_aabb->next()) {
    s_aabb->compute_shell(t_aabb);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_aabb->ao_iterator("ALL");
    for (ii->first(); ii->is_done() == false; ii->next()) {
         int i = ii->i(); // A
         int j = ii->j(); // A
         int k = ii->k(); // B
         int l = ii->l(); // B
         double eri = b_aabb[ii->index()];

         v[i][j] += 2.0 * d[k][l] * eri;
    }
  }

  // K contribution
  if (f_abab) {
     std::shared_ptr<oepdev::ShellCombinationsIterator> s_abab = oepdev::ShellCombinationsIterator::build(f_abab, "ALL");
     std::shared_ptr<psi::TwoBodyAOInt> t_abab(f_abab->eri()); const double* b_abab = t_abab->buffer();

     for (s_abab->first(); s_abab->is_done() == false; s_abab->next()) {             
       s_abab->compute_shell(t_abab);
       std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_abab->ao_iterator("ALL");
       for (ii->first(); ii->is_done() == false; ii->next()) {
            int i = ii->i(); // A
            int j = ii->j(); // B
            int k = ii->k(); // A
            int l = ii->l(); // B
            double eri = b_abab[ii->index()];
                                                                                     
            v[i][k] -= d[j][l] * eri;
       }
     }

  }
  return V;
}

extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_DFI_Vel_JK(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::IntegralFactory> f_abab,
                std::shared_ptr<psi::Matrix> db)
{
 return _calculate_DFI_Vel(f_aabb, f_abab, db);
}

extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_DFI_Vel_J(
                std::shared_ptr<psi::IntegralFactory> f_aabb,
                std::shared_ptr<psi::Matrix> db)
{
 return _calculate_DFI_Vel(f_aabb, nullptr, db);
}

extern "C" PSI_API
std::shared_ptr<psi::Matrix> calculate_OEP_basisopt_V(const int& nt,
                std::shared_ptr<psi::IntegralFactory> f_pppt,
                std::shared_ptr<psi::Matrix> ca, std::shared_ptr<psi::Matrix> da)
{
  std::shared_ptr<psi::Matrix> V = std::make_shared<psi::Matrix>("", nt, ca->ncol());
  double** c = ca->pointer();
  double** d = da->pointer();
  double** v = V ->pointer();
  const int nI = ca->ncol();

  std::shared_ptr<oepdev::ShellCombinationsIterator> s_pppt = oepdev::ShellCombinationsIterator::build(f_pppt, "ALL");
  std::shared_ptr<psi::TwoBodyAOInt> t_pppt(f_pppt->eri()); const double* b_pppt = t_pppt->buffer();

  for (s_pppt->first(); s_pppt->is_done() == false; s_pppt->next()) {

    s_pppt->compute_shell(t_pppt);
    std::shared_ptr<oepdev::AOIntegralsIterator> ii = s_pppt->ao_iterator("ALL");

    for (ii->first(); ii->is_done() == false; ii->next()) {

         double eri = b_pppt[ii->index()];

         if (std::abs(eri) > 0.0) {

             int i = ii->i(); // P 
             int j = ii->j(); // P
             int k = ii->k(); // P
             int l = ii->l(); // T

             double dij = d[i][j];
             double djk = d[j][k];

             for (int I=0; I<nI; ++I) v[l][I] += eri * (2.0 * c[k][I] * dij - c[i][I] * djk);

         }
    }
  }

  return V;
}



} // EndNameSpace oepdev
