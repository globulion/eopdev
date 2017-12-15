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
 *  __Note:__ Works only for C1 symmetry!
 *  The following variables are exact deep copies of variables inside Wavefunction object,
 *  that is created for the _whole_ molecule cluster:
 *    - *basissets (DF/RI/F12/etc basis sets)_
 *    - *basisset_ (ORBITAL basis set)
 *    - *sobasisset_ (Primary basis set for SO integrals)
 *    - *AO2SO_ (AO2SO conversion matrix (AO in rows, SO in cols)
 *    - *molecule_ (Molecule that this wavefunction is run on)
 *    - *options_ (Options object)
 *    - *psio_ (PSI file access variables)
 *    - *integral_ (Integral factory)
 *    - *factory_ (Matrix factory for creating standard sized matrices)
 *    - memory_ (How much memory you have access to)
 *    - *nalpha_, nbeta_ (Total alpha and beta electrons)
 *    - nfrzc_ (Total frozen core orbitals)
 *    - *doccpi_ (Number of doubly occupied per irrep)
 *    - *soccpi_ (Number of singly occupied per irrep)
 *    - *frzcpi_ (Number of frozen core per irrep)
 *    - *frzvpi_ (Number of frozen virtuals per irrep)
 *    - *nalphapi_ (Number of alpha electrons per irrep)
 *    - *nbetapi_ (Number of beta electrons per irrep)
 *    - *nsopi_ (Number of so per irrep)
 *    - *nmopi_ (Number of mo per irrep)
 *    - *nso_ (Total number of SOs)
 *    - *nmo_ (Total number of MOs)
 *    - *nirrep_ (Number of irreps)
 *    - 
 */
class WavefunctionUnion : public Wavefunction
{
  private:
    void common_init(void);
  protected:
    /// Number of molecular orbitals in isolated monomer 1 and 2
    int nmo_1_, nmo_2_;
    /// Numberof SO's in isolated monomer 1 and 2
    int nso_1_, nso_2_;
    /// Number of alpha and beta electrons
    int nalpha_1_, nalpha_2_, nbeta_1_, nbeta_2_;
    /// Energies of isolated molecules
    double energy_1_, energy_2_;
    /// The wavefunction for a dimer (electrons relaxed in the field of monomers)
    SharedWavefunction dimer_wavefunction_;
    /// Unrelaxed wavefunctions of monomers
    SharedWavefunction wfn_1_, wfn_2_;
    /// Basis functions of monomenrs
    SharedBasisSet primary_1_, primary_2_;
    /// Molecules of monomers
    SharedMolecule molecule_1_, molecule_2_;

  public:
    WavefunctionUnion(SharedWavefunction ref_wfn, Options& options);
    virtual ~WavefunctionUnion();
    double compute_energy();

    // <--- Getters ---> //
    int nmo_1() const {return nmo_1_;}
    int nso_1() const {return nso_1_;}
    int nmo_2() const {return nmo_2_;}
    int nso_2() const {return nso_2_;}
    int nalpha_1() const {return nalpha_1_;}
    int nalpha_2() const {return nalpha_2_;}
    int nbeta_1() const {return nbeta_1_;}
    int nbeta_2() const {return nbeta_2_;}
    double energy_1() const {return energy_1_;}
    double energy_2() const {return energy_2_;}
    SharedMolecule molecule_1() const {return molecule_1_;}
    SharedMolecule molecule_2() const {return molecule_2_;}
    SharedWavefunction wfn_1() const {return wfn_1_;}
    SharedWavefunction wfn_2() const {return wfn_2_;}
    SharedBasisSet primary_1() const {return primary_1_;}
    SharedBasisSet primary_2() const {return primary_2_;}

};

} // EndNameSpace oepdev_libutil


#endif // _oepdev_libutil_util_h

