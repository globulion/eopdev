#!/usr/bin/python
import oepdev, psi4

__all__ = ["Scanner", "make_union"]

class Scanner:
    """
 -----------------------------------------------------------------------------------------
 Helps in making distance- and angular-displacement scans for a dimer. 
 
 It is assumed that the monomer 1 (first molecule) is constant and does not
 change atomic coordinates. Monomer 2 (second molecule) is subject to user-defined
 modifications.

 Usage:
  scanner = Scanner(dimer)
  # make some changes in coordinates of dimer
  # ...
  # create Wfn union object that can be used for oepdev.OEPDevSolver
  un = scanner.make()
  solver = oepdev.OEPDevSolver.build("ELECTROSTATIC ENERGY", un)
  solver.compute_benchmark()
  un.clear_dpd() # necessary!

 Notes: 
   o Z-matrix does not work properly with fragments. Use Cartesian coordinates 
     for geometry specification.
 -----------------------------------------------------------------------------------------
                                                              Last Revision: 22 Nov 2018
 """
    def __init__(self, dimer):
        self.dimer = dimer
        self._initialize(dimer)
    def _initialize(self, dimer):
        """Grab first monomer (unchanging)"""
        dimer.update_geometry()
        #
        self.molecule_A = dimer.extract_subsets(1)
        # --- primary (SCF)
        self.basis_A = psi4.core.BasisSet.build(self.molecule_A, "BASIS", psi4.core.get_global_option("BASIS"), puream=False)
        # --- auxiliary (DF-SCF)
        self.basis_df_scf_A = psi4.core.BasisSet.build(self.molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
        # --- auxiliary (OEP)
        self.basis_df_oep_A = psi4.core.BasisSet.build(self.molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"), puream=False)
        # --- intermediate (OEP)
        self.basis_int_oep_A = psi4.core.BasisSet.build(self.molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_INT"), puream=False)
        self.en_1, self.wfn_1 = psi4.energy('hf', molecule=self.molecule_A, return_wfn=True)
    def get_basis(self):
        """Get the primary and auxiliary SCF-DF basis sets for a dimer"""
        self.dimer.update_geometry()
        basis        = psi4.core.BasisSet.build(self.dimer, "BASIS", psi4.core.get_global_option("BASIS"       ), puream=False)
        basis_df_scf = psi4.core.BasisSet.build(self.dimer, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
        return basis, basis_df_scf
    def make(self):
        """Creates union of wavefunctions"""
        self.basis, self.basis_df_scf = self.get_basis()
        self.molecule_B = self.dimer.extract_subsets(2)
        self.basis_B = psi4.core.BasisSet.build(self.molecule_B, "BASIS", psi4.core.get_global_option("BASIS"), puream=False)
        self.basis_df_scf_B = psi4.core.BasisSet.build(self.molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
        self.basis_df_oep_B = psi4.core.BasisSet.build(self.molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"), puream=False)
        self.basis_int_oep_B= psi4.core.BasisSet.build(self.molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_INT"), puream=False)
        self.en_2, self.wfn_2 = psi4.energy('hf', molecule=self.molecule_B, return_wfn=True)

        #print(self.molecule_A.geometry().array_interface())
        #print(self.molecule_B.geometry().array_interface())
        un = oepdev.WavefunctionUnion(self.dimer, self.basis, self.basis_df_scf,
                                  self.basis_A, self.basis_B,
                                  self.basis_df_oep_A, self.basis_df_oep_B,
                                  self.basis_df_scf_A, self.basis_df_scf_B,
                                  self.basis_int_oep_A, self.basis_int_oep_B,
                                  self.wfn_1, self.wfn_2,
                                  psi4.core.get_options())
        return un


def make_union(molecule):
    """Create WavefunctionUnion from the molecular dimer"""
    molecule.update_geometry()
    #
    molecule_A     = molecule.extract_subsets(1)                                                               
    molecule_B     = molecule.extract_subsets(2)

    # --- primary                                                                                                                 
    basis_A        = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("BASIS"), puream=False)
    basis_B        = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("BASIS"), puream=False)
    # --- auxiliary (DF-SCF)
    basis_df_scf_A = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
    basis_df_scf_B = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
    # --- auxiliary (OEP)
    basis_df_oep_A = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"), puream=False)
    basis_df_oep_B = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"), puream=False)
    # --- intermediate (OEP)
    basis_int_oep_A= psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_INT"), puream=False)
    basis_int_oep_B= psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_INT"), puream=False)

    basis        = psi4.core.BasisSet.build(molecule, "BASIS", psi4.core.get_global_option("BASIS"       ), puream=False)
    basis_df_scf = psi4.core.BasisSet.build(molecule, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=False)
    en_1, wfn_1 = psi4.energy('hf', molecule=molecule_A, return_wfn=True)
    en_2, wfn_2 = psi4.energy('hf', molecule=molecule_B, return_wfn=True)

    #psi4_io.set_default_path('/home/globulion/scr-d')
    psi4.core.clean()

    un = oepdev.WavefunctionUnion(molecule, basis, basis_df_scf,
                                  basis_A, basis_B,
                                  basis_df_oep_A, basis_df_oep_B,
                                  basis_df_scf_A, basis_df_scf_B,
                                  basis_int_oep_A, basis_int_oep_B,
                                  wfn_1, wfn_2,
                                  psi4.core.get_options())
    return un

