#ifndef _oepdev_libutil_util_h
#define _oepdev_libutil_util_h

#include<cstdio>
#include<string>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libfunctional/superfunctional.h"

#include "psi4/libscf_solver/rhf.h"

namespace oepdev_libutil{

using namespace psi;
using namespace std;


using SharedMolecule        = std::shared_ptr<Molecule>;
using SharedSuperFunctional = std::shared_ptr<SuperFunctional>;
using SharedWavefunction    = std::shared_ptr<Wavefunction>;
using SharedVector          = std::shared_ptr<Vector>;
using SharedMatrix          = std::shared_ptr<Matrix>;
using SharedBasisSet        = std::shared_ptr<BasisSet>;

/** \brief Print preambule for module OEPDEV
 */
extern "C" void preambule(void);

/** \brief Set up DFT functional.
 *
 *  Now it accepts only pure HF functional.
 *  @param name name of the functional ("HF" is now only available)
 *  @param options psi::Options object
 *  @return psi::SharedSuperFunctional object with functional.
 */
extern "C" std::shared_ptr<SuperFunctional> 
create_superfunctional(std::string name, Options& options);

/** \brief Extract molecule from dimer.
 *
 *  @param molecule_dimer psi::SharedMolecule object with dimer
 *  @param id index of a molecule (starts from 1)
 *  @return psi::SharedMolecule object with indicated monomer
 */
extern "C" std::shared_ptr<Molecule>
extract_monomer(std::shared_ptr<const Molecule> molecule_dimer, int id);

/** \brief Solve RHF-SCF equations for a given molecule in a given basis set.
 *
 *  @param molecule psi::SharedMolecule object with molecule
 *  @param primary shared primary basis set 
 *  @param functional DFT functional
 *  @param options psi::Options object
 *  @param psio psi::PSIO object
 *  @return psi::SharedWavefunction SCF wavefunction of the molecule
 */
extern "C" std::shared_ptr<Wavefunction>
solve_scf(std::shared_ptr<Molecule> molecule, 
          std::shared_ptr<BasisSet> primary, 
          std::shared_ptr<SuperFunctional> functional,
          Options& options,
          std::shared_ptr<PSIO> psio);

/** \brief Union of two Wavefunction objects.
 *  
 *  The WavefunctionUnion is the exact union of two unperturbed Wavefunctions.
 *  __Note:__ Works only for C1 symmetry! Therefore `this->nirrep() = 1`.
 *  The following variables are shallow copies of variables inside Wavefunction object,
 *  that is created for the _whole_ molecule cluster:
 *    - basissets (DF/RI/F12/etc basis sets)_
 *    - basisset_ (ORBITAL basis set)
 *    - sobasisset_ (Primary basis set for SO integrals)
 *    - AO2SO_ (AO2SO conversion matrix (AO in rows, SO in cols)
 *    - molecule_ (Molecule that this wavefunction is run on)
 *    - options_ (Options object)
 *    - psio_ (PSI file access variables)
 *    - integral_ (Integral factory)
 *    - factory_ (Matrix factory for creating standard sized matrices)
 *    - memory_ (How much memory you have access to)
 *    - nalpha_, nbeta_ (Total alpha and beta electrons)
 *    - nfrzc_ (Total frozen core orbitals)
 *    - doccpi_ (Number of doubly occupied per irrep)
 *    - soccpi_ (Number of singly occupied per irrep)
 *    - frzcpi_ (Number of frozen core per irrep)
 *    - frzvpi_ (Number of frozen virtuals per irrep)
 *    - nalphapi_ (Number of alpha electrons per irrep)
 *    - nbetapi_ (Number of beta electrons per irrep)
 *    - nsopi_ (Number of so per irrep)
 *    - nmopi_ (Number of mo per irrep)
 *    - nso_ (Total number of SOs)
 *    - nmo_ (Total number of MOs)
 *    - nirrep_ (Number of irreps)
 *  The rest is altered so that the Wavefunction parameters reflect
 *  a cluster of non-interacting (uncoupled, isolated, unrelaxed) 
 *  molecular electron densities.
 */
class WavefunctionUnion : public Wavefunction//, public std::enable_shared_from_this<WavefunctionUnion>
{
  private:
    void common_init(SharedWavefunction ref_wfn);

  protected:
    /// Number of isolated molecules
    int nIsolatedMolecules_;
    /// List of molecules
    std::vector<SharedMolecule> l_molecule_;
    /// List of primary basis functions per molecule
    std::vector<SharedBasisSet> l_primary_;
    /// List of auxiliary basis functions per molecule
    std::vector<SharedBasisSet> l_auxiliary_;
    /// List of isolated wavefunctions (electrons unrelaxed)
    std::vector<SharedWavefunction> l_wfn_;
    /// List of names of isolated wavefunctions
    std::vector<std::string> l_name_;
    /// List of numbers of molecular orbitals (MO's) per molecule
    std::vector<int> l_nmo_;
    /// List of numbers of SO's per molecule
    std::vector<int> l_nso_;
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

    /// The wavefunction for a dimer (electrons relaxed in the field of monomers)
    SharedWavefunction dimer_wavefunction_;

  public:
    /* \brief Constructor. 
     *
     *  Provide wavefunction with molecule containing at least 2 fragments.
     *  @param ref_wfn - reference wavefunction
     *  @param options - Psi4 options
     */
    WavefunctionUnion(SharedWavefunction ref_wfn, Options& options);

    /// Destructor
    virtual ~WavefunctionUnion();

    /// Compute Energy (now blank)
    virtual double compute_energy();

    // <--- Getters ---> //
    int                  l_nmo       (int n) const {return l_nmo_          [n];}         
    int                  l_nso       (int n) const {return l_nso_          [n];}
    int                  l_nalpha    (int n) const {return l_nalpha_       [n];}
    int                  l_nbeta     (int n) const {return l_nbeta_        [n];}
    double               l_energy    (int n) const {return l_energy_       [n];}
    SharedMolecule       l_molecule  (int n) const {return l_molecule_     [n];}
    SharedBasisSet       l_primary   (int n) const {return l_primary_      [n];}
    SharedBasisSet       l_auxiliary (int n) const {return l_auxiliary_    [n];}
    SharedWavefunction   l_wfn       (int n) const {return l_wfn_          [n];}

};

} // EndNameSpace oepdev_libutil


#endif // _oepdev_libutil_util_h

