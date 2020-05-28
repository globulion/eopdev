#ifndef _oepdev_libutil_wavefunction_union_h
#define _oepdev_libutil_wavefunction_union_h
/** @file wavefunction_union.h */

#include<cstdio>
#include<string>
#include<map>

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/local.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"

#include "psi4/libscf_solver/rhf.h"
#include "psi4/libdpd/dpd.h"


namespace oepdev {

using namespace psi;
using namespace std;

using SharedMolecule           = std::shared_ptr<Molecule>;        
using SharedSuperFunctional    = std::shared_ptr<SuperFunctional>;
using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedVector             = std::shared_ptr<Vector>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedBasisSet           = std::shared_ptr<BasisSet>;
using SharedMOSpace            = std::shared_ptr<MOSpace>;
using SharedMOSpaceVector      = std::vector<std::shared_ptr<MOSpace>>;
using SharedIntegralTransform  = std::shared_ptr<IntegralTransform>;
using SharedLocalizer          = std::shared_ptr<Localizer>;
/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief Union of two Wavefunction objects.
 *  
 *  The WavefunctionUnion is the union of two unperturbed Wavefunctions.
 *
 *  __Notes:__ 
 *    1. Works only for C1 symmetry! Therefore `this->nirrep() = 1`.
 *    2. Does not set `reference_wavefunction_`
 *    3. Sets `oeprop_` for the union of uncoupled molecules
 *    3. Performs Hadamard sums on `H_`, `Fa_`, `Da_`, `Ca_` and `S_` 
 *       based on uncoupled wavefunctions. 
 *    4. Since it is based on shallow copy of the original Wavefunction,
 *       it __changes__ contents of this wavefunction. Reallocate and copy if 
 *       you want to keep the original wavefunction.
 *
 *  __Warnings:__
 *    1. Gradients, Hessians and frequencies are not touched, hence they are __wrong__!
 *    2. Lagrangian (if present) is not touched, hence its __wrong__!
 *    3. Ca/Cb and epsilon subsets were reimplemented from psi::Wavefunction to remove sorting of orbitals.
 *       However, the corresponding member functions are not virtual in psi::Wavefunction.
 *       This could bring problems when upcasting.
 * 
 *  The following variables are _shallow_ copies of variables inside 
 *  the Wavefunction object, that is created for the _whole_ molecule cluster:
 *    - `basissets_` (DF/RI/F12/etc basis sets)_
 *    - `basisset_` (ORBITAL basis set)
 *    - `sobasisset_` (Primary basis set for SO integrals)
 *    - `AO2SO_` (AO2SO conversion matrix (AO in rows, SO in cols)
 *    - `molecule_` (Molecule that this wavefunction is run on)
 *    - `options_` (Options object)
 *    - `psio_` (PSI file access variables)
 *    - `integral_` (Integral factory)
 *    - `factory_` (Matrix factory for creating standard sized matrices)
 *    - `memory_` (How much memory you have access to)
 *    - `nalpha_`, `nbeta_` (Total alpha and beta electrons)
 *    - `nfrzc_` (Total frozen core orbitals)
 *    - `doccpi_` (Number of doubly occupied per irrep)
 *    - `soccpi_` (Number of singly occupied per irrep)
 *    - `frzcpi_` (Number of frozen core per irrep)
 *    - `frzvpi_` (Number of frozen virtuals per irrep)
 *    - `nalphapi_` (Number of alpha electrons per irrep)
 *    - `nbetapi_` (Number of beta electrons per irrep)
 *    - `nsopi_` (Number of so per irrep)
 *    - `nmopi_` (Number of mo per irrep)
 *    - `nso_` (Total number of SOs)
 *    - `nmo_` (Total number of MOs)
 *    - `nirrep_` (Number of irreps; must be equal to 1 due to symmetry reasons)
 *    - `same_a_b_dens_` and `same_a_b_orbs_`
 *  The rest is altered so that the Wavefunction parameters reflect
 *  a cluster of non-interacting (uncoupled, isolated, unrelaxed) 
 *  molecular electron densities.
 *  
 */
class WavefunctionUnion : public Wavefunction
{
  protected:

    /// Number of isolated molecules
    int nIsolatedMolecules_;

    /// The wavefunction for a dimer (electrons relaxed in the field of monomers)
    SharedWavefunction dimer_wavefunction_;

    /// Integral transform object (2- and 4-index transformations)
    SharedIntegralTransform integrals_;

    /// whether orbitals of the union were localized (or not)
    bool hasLocalizedOrbitals_;

    /// Dictionary of MO spaces for the entire union (`OCC` and `VIR`)
    std::map<const std::string, SharedMOSpace> mospacesUnion_;


    // ---> Monomer Lists <--- //

    /// List of molecules
    std::vector<SharedMolecule> l_molecule_;
    /// List of primary basis functions per molecule
    std::vector<SharedBasisSet> l_primary_;
    /// List of auxiliary basis functions per molecule
    std::vector<SharedBasisSet> l_auxiliary_;
    /// List of intermediate basis functions per molecule
    std::vector<SharedBasisSet> l_intermediate_;
    /// List of guess basis functions per molecule
    std::vector<SharedBasisSet> l_guess_;
    /// List of original isolated wavefunctions (electrons unrelaxed)
    std::vector<SharedWavefunction> l_wfn_;
    /// List of names of isolated wavefunctions
    std::vector<std::string> l_name_;
    /// List of basis function numbers per molecule
    std::vector<int> l_nbf_;
    /// List of numbers of molecular orbitals (MO's) per molecule
    std::vector<int> l_nmo_;
    /// List of numbers of SO's per molecule
    std::vector<int> l_nso_;
    /// List of numbers of doubly occupied orbitals per molecule
    std::vector<int> l_ndocc_;
    /// List of numbers of virtual orbitals per molecule
    std::vector<int> l_nvir_;
    /// List of basis set offsets per molecule
    std::vector<int> l_noffs_ao_;
    /// List of energies of isolated wavefunctions
    std::vector<double> l_energy_;
    /// List of frozen-core energies per isolated wavefunction
    std::vector<double> l_efzc_;
    /// List of information per wfn whether it was obtained using DF or not
    std::vector<bool> l_density_fitted_;
    /// List of numbers of alpha electrons per isolated wavefunction
    std::vector<int> l_nalpha_; 
    /// List of numbers of beta electrons per isolated wavefunction
    std::vector<int> l_nbeta_; 
    /// List of numbers of frozen-core orbitals per isolated molecule
    std::vector<int> l_nfrzc_;
    /// List of orbital localizers
    std::vector<SharedLocalizer> l_localizer_;
    /// List of dictionaries of MO spaces
    std::vector<std::map<const std::string, SharedMOSpace>> l_mospace_;

    /// One-Electron Property
    std::shared_ptr<psi::OEProp> oeprop_;


  public:
    /** \brief Constructor. 
     *
     *  Provide wavefunction with molecule containing at least 2 fragments.
     *  @param ref_wfn - reference wavefunction
     *  @param options - Psi4 options
     */
    WavefunctionUnion(SharedWavefunction ref_wfn, Options& options);

    /** \brief Constructor. 
     *
     *  Provide molecule dimer and all the required monomer basis sets and wavefunctions.
     *  @param dimer          - molecule object
     *  @param primary        - basis set object: dimer (primary)
     *  @param auxiliary_df   - basis set object: dimer (DF SCF)
     *  @param guess          - basis set object: dimer (guess)
     *  @param primary_1      - basis set object for 1st monomer
     *  @param primary_2      - basis set object for 2nd monomer
     *  @param auxiliary_1    - basis set object for 1st monomer
     *  @param auxiliary_2    - basis set object for 2nd monomer
     *  @param auxiliary_df_1 - basis set object for 1st monomer
     *  @param auxiliary_df_2 - basis set object for 2nd monomer
     *  @param intermediate_1 - basis set object for 1st monomer
     *  @param intermediate_2 - basis set object for 2nd monomer
     *  @param guess_1        - basis set object for 1st monomer
     *  @param guess_2        - basis set object for 2nd monomer
     *  @param wfn_1          - unperturbed wavefunction object
     *  @param wfn_2          - unperturbed wavefunction object
     *  @param options        - Psi4 options
     */
     WavefunctionUnion(
//		SharedWavefunction ref_wfn,
		SharedMolecule dimer,
		SharedBasisSet primary,
		SharedBasisSet auxiliary_df,
                SharedBasisSet guess,
		SharedBasisSet primary_1,
		SharedBasisSet primary_2,
		SharedBasisSet auxiliary_1,
		SharedBasisSet auxiliary_2,
		SharedBasisSet auxiliary_df_1,
		SharedBasisSet auxiliary_df_2,
		SharedBasisSet intermediate_1,
		SharedBasisSet intermediate_2,
		SharedBasisSet guess_1,
		SharedBasisSet guess_2,
		SharedWavefunction wfn_1,
		SharedWavefunction wfn_2,
		Options& options
		);


    /// Destructor
    virtual ~WavefunctionUnion();


    // ---> Computers <--- //

    /// Compute Energy (now blank)
    virtual double compute_energy();

    /// Compute Nuclear Repulsion Energy between unions
    virtual double nuclear_repulsion_interaction_energy();

    /// Localize Molecular Orbitals
    void localize_orbitals();

    /// Transform Integrals (2- and 4-index transformations)
    void transform_integrals(); 

    /// Close the DPD instance
    void clear_dpd();


    // ---> Accessors <--- //
    // ---> Properties of Particular Fragments <--- //

    /// Get number of molecular orbitals of the *n*th fragment
    int l_nmo (int n) const {return l_nmo_[n];}

    /// Get number of symmetry orbitals of the *n*th fragment
    int l_nso (int n) const {return l_nso_[n];}

    /// Get number of doubly occupied orbitals of the *n*th fragment
    int l_ndocc (int n) const {return l_ndocc_[n];}

    /// Get number of virtual orbitals of the *n*th fragment
    int l_nvir (int n) const {return l_nvir_[n];}

    /// Get the number of the alpha electrons of the *n*th fragment
    int l_nalpha (int n) const {return l_nalpha_[n];}

    /// Get the number of the beta electrons of the *n*th fragment
    int l_nbeta (int n) const {return l_nbeta_[n];}

    /// Get number of basis functions of the *n*th fragment
    int l_nbf (int n) const {return l_nbf_[n];}

    /// Get the basis set offset of the *n*th fragment
    int l_noffs_ao (int n) const {return l_noffs_ao_[n];}

    /// Get the reference energy of the *n*th fragment
    double l_energy (int n) const {return l_energy_[n];}

    /// Get the molecule object of the *n*th fragment
    SharedMolecule l_molecule (int n) const {return l_molecule_[n];}

    /// Get the primary basis set object of the *n*th fragment
    SharedBasisSet l_primary (int n) const {return l_primary_[n];}

    /// Get the auxiliary basis set object of the *n*th fragment
    SharedBasisSet l_auxiliary (int n) const {return l_auxiliary_[n];}

    /// Get the intermediate basis set object of the *n*th fragment
    SharedBasisSet l_intermediate (int n) const {return l_intermediate_[n];}

    /// Get the guess basis set object of the *n*th fragment
    SharedBasisSet l_guess (int n) const {return l_guess_[n];}

    /// Get the wavefunction object of the *n*th fragment
    SharedWavefunction l_wfn (int n) const {return l_wfn_[n];}

    /// Get the MO space named `label` (either `OCC` or `VIR`) of the *n*th fragment 
    SharedMOSpace l_mospace (int n, const std::string& label) const {return l_mospace_[n].at(label);}

    /// Get the orbital localizer object of the *n*th fragment
    SharedLocalizer l_localizer (int n) const;

    // ---> Properties of Entire Union <--- //

    /// Get the integral transform object of the entire union
    SharedIntegralTransform integrals (void ) const;

    /// If union got its molecular orbital localized or not
    bool has_localized_orbitals(void ) const {return hasLocalizedOrbitals_;}

    /// Get the primary basis set for the entire union
    SharedBasisSet primary (void ) const {return basisset_;}

    /// Get the MO space named `label` (either `OCC` or `VIR`) 
    SharedMOSpace mospace (const std::string& label) const {return mospacesUnion_.at(label);}

    /**
    * Return a subset of the Ca matrix in a desired basis
    * @param basis the symmetry basis to use
    *  AO, SO
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    * @return the matrix in Pitzer order in the desired basis
    **/
    SharedMatrix Ca_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /**
    * Return a subset of the Cb matrix in a desired basis
    * @param basis the symmetry basis to use
    *  AO, SO
    * @param subset the subset of orbitals to return
    *  ALL, ACTIVE, FROZEN, OCC, VIR, FROZEN_OCC, ACTIVE_OCC, ACTIVE_VIR, FROZEN_VIR
    * @return the matrix in Pitzer order in the desired basis
    **/
    SharedMatrix Cb_subset(const std::string& basis = "SO", const std::string& subset = "ALL");

    /// Helpers for Ca_ and Cb_ matrix transformers
    SharedMatrix C_subset_helper(SharedMatrix C, const Dimension& noccpi, SharedVector epsilon, 
                                 const std::string& basis, const std::string& subset);

    /// Helper for epsilon transformer
    SharedVector epsilon_subset_helper(SharedVector epsilon, const Dimension& noccpi, 
                                       const std::string& basis, const std::string& subset);


    // <--- Printers ---> //

    /// Print information about this wavefunction union
    void print_header(void);

    /// Print the MO ingegrals
    void print_mo_integrals(void);


  private:
    /// Finish initialising the object
    void common_init(
		SharedWavefunction ref_wfn,
		SharedMolecule molecule_1,
		SharedMolecule molecule_2,
		SharedBasisSet primary_1,
		SharedBasisSet primary_2,
		SharedBasisSet auxiliary_1,
		SharedBasisSet auxiliary_2,
		SharedBasisSet auxiliary_df_1,
		SharedBasisSet auxiliary_df_2,
		SharedBasisSet intermediate_1,
		SharedBasisSet intermediate_2,
		SharedBasisSet guess_1,
		SharedBasisSet guess_2,
		SharedWavefunction wfn_1,
		SharedWavefunction wfn_2
		);
};

/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_wavefunction_union_h
