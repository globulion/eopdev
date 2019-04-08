#include "util.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/mintshelper.h"
//#include "psi4/libdpd/dpd.h"

namespace oepdev{

extern "C" PSI_API
void preambule(void) {
      outfile->Printf("                                                                             \n");
      outfile->Printf("    -------------------------------------------------------------------------\n");
      outfile->Printf("          OepDev: One-Electron Effective Potentials Development Routine      \n");
      outfile->Printf("                               OepDev 0.1 (no release)                       \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("                         Git: Rev {}                                         \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("    Bartosz Błasiak                                                          \n");
      outfile->Printf("                                                                             \n");
      outfile->Printf("    References:                                                              \n");
      outfile->Printf("                                                                             \n");
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
std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id) {
    std::vector<int> real_list; real_list.push_back(id-1);
    std::vector<int> ghost_list;
    return molecule_dimer->extract_subsets(real_list, ghost_list);
}

extern "C" PSI_API
std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
	  std::shared_ptr<BasisSet> auxiliary,
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
    SharedWavefunction scf = std::make_shared<psi::scf::RHF>(scf_base, functional, options, psio);
    scf->compute_energy();
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
  psio->tocprint(PSIF_LIBTRANS_DPD);

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
  psio->tocprint(PSIF_LIBTRANS_DPD);

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




} // EndNameSpace oepdev
