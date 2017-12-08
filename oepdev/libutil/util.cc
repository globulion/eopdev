#include "util.h"

//namespace util{

extern "C" void preambule(void) {
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

extern "C" std::shared_ptr<SuperFunctional> 
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

extern "C" std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id) {
    std::vector<int> real_list; real_list.push_back(id-1);
    std::vector<int> ghost_list;
    return molecule_dimer->extract_subsets(real_list, ghost_list);
}

extern "C" std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
          std::shared_ptr<SuperFunctional> functional,
          Options& options,
          std::shared_ptr<PSIO> psio) 
{
    SharedWavefunction scf_base(new Wavefunction(molecule, primary, options));
    SharedWavefunction scf = SharedWavefunction(new scf::RHF(scf_base, functional, options, psio));
    scf->compute_energy();
    return scf;
}

//}
