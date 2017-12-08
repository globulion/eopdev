/*
 * @BEGIN LICENSE
 *
 * oepdev by Bartosz Błasiak (email: blasiak.bartosz@gmail.com)
 * A plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#define DIIS_BB

#include <cstdlib>
#include <cstdio>
#include <string>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libfunctional/superfunctional.h"

#include "oepdev/libutil/util.h"
#include "oepdev/libutil/cphf.h"

namespace psi{ 
namespace oepdev{

/** \brief Options for the OEPDEV plugin.
 *
 *  @param name name of driver function
 *  @param options psi::Options object
 *  @return true
 */
extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "OEPDEV" || options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Basis set for dimer A -*/
        options.add_str("BASIS_A", "");
        /*- Basis set for dimer B -*/
        options.add_str("BASIS_B", "");
        /*- CPHF maximum iterations -*/
        options.add_int("CPHF_MAXITER", 50);
        /*- CPHF convergence -*/
        options.add_double("CPHF_CONVERGENCE", 1.0E-8);
        /*- whether use DIIS for CPHF -*/
        options.add_bool("CPHF_DIIS", false);
        /*- size of DIIS subspace for CPHF -*/
        options.add_int("CPHF_DIIS_DIM", 3);

    }
    return true;
}

/** \brief Main routine of the OEPDEV plugin.
 *
 *  Created with intention to test various models of the interaction energy 
 *  between two molecules, described by the Hartree-Fock-Roothaan-Hall theory 
 *  or the configuration interaction with singles theory. 
 *
 *  In particular, the plugin tests the models of:
 *
 *   1. the Pauli exchange-repulsion interaction energy    (Project II ) 
 *   2. the Induction interaction energy                   (Project III)
 *   3. the excitation energy transfer couplings           (Project I  )
 *
 *  against benchmarks (exact or reference solutions). Detailed list of models 
 *  is given below:
 *
 *  +----------------------------------------------------------------------------+
 *  |                              Interaction Property                          |
 *  +--------------------------+--------------------------+----------------------+
 *  | Pauli energy             | Induction energy         | EET Coupling         |
 *  +--------------------------+--------------------------+----------------------+
 *  |                                    Methods                                 |       
 *  +==========================+==========================+======================+
 *  | EFP2-Pauli               | EFP2-Induced Dipoles     | TrCAMM               |
 *  +--------------------------+--------------------------+----------------------+
 *  | Murrel et al.'s theory   | Density Susceptibility   | OEP-ET/HT            |
 *  +--------------------------+--------------------------+----------------------+
 *  | OEP-Murrel et al.'s      |                          | TDFI-TI              |
 *  +--------------------------+--------------------------+----------------------+ 
 *  |                          |                          | FED                  |
 *  +--------------------------+--------------------------+----------------------+ 
 *  | Exact (Stone's)          | Exact (incl. CT)         | Exact (ESD)          |
 *  +--------------------------+--------------------------+----------------------+
 *
 *  The target models introduced in the Project shall be tested against the
 *  following benchmarks and compared with the following state-of-the-art models:
 *
 *  +--------------------------+--------------------------+----------------------+
 *  | Target Model             | Benchmarks               | Competing Model      |
 *  +==========================+==========================+======================+
 *  | OEP-Murrel et al.'s      | Murrel et al.'s          | EFP2-Pauli           |
 *  |                          | Exact (Stone's)          |                      |
 *  +--------------------------+--------------------------+----------------------+
 *  | OEP-ET/HT + TrCAMM       | Exact (ESD)              | TDFI-TI              |
 *  |                          | FED                      | FED                  |
 *  |                          | TDFI-TI                  |                      |
 *  +--------------------------+--------------------------+----------------------+
 *  | Density Susceptibility   | Exact (incl. CT)         | EFP2-Induced Dipoles |
 *  +--------------------------+--------------------------+----------------------+
 *
 *  Parameters of the plugin driver:
 *
 *  @param ref_wfn shared wavefunction of a dimer
 *  @param options psi::Options object
 *  @return psi::SharedWavefunction (as for now the same as ref_wfn)
 */
extern "C"
SharedWavefunction oepdev(SharedWavefunction ref_wfn, Options& options)
{
    oepdev_libutil::preambule();

    int print = options.get_int("PRINT");

    // Parse molecules, fragments, basis sets and other primary informations
    std::shared_ptr<Wavefunction>    scf_A;
    std::shared_ptr<Wavefunction>    scf_B;
    std::shared_ptr<PSIO>            psio           = PSIO::shared_object();
    std::shared_ptr<SuperFunctional> functional     = oepdev_libutil::create_superfunctional("HF", options);
    std::shared_ptr<Molecule>        molecule_dimer = ref_wfn->molecule();
    std::shared_ptr<BasisSet>        primary        = ref_wfn->basisset();
    std::shared_ptr<BasisSet>        primary_A      = ref_wfn->get_basisset("BASIS_SCF_A");
    std::shared_ptr<BasisSet>        primary_B      = ref_wfn->get_basisset("BASIS_SCF_B");
    std::shared_ptr<Molecule>        molecule_A     = oepdev_libutil::extract_monomer(molecule_dimer, 1);
    std::shared_ptr<Molecule>        molecule_B     = oepdev_libutil::extract_monomer(molecule_dimer, 2);
    molecule_A->set_name("Monomer 1");
    molecule_B->set_name("Monomer 2");
    molecule_dimer->set_name("Aggregate (Dimer)");

    outfile->Printf("  ==> Molecules in the aggregate <==\n");
    molecule_A->print();
    molecule_B->print();

    //TRIAL create BasisSet objects
    //TRIAL std::shared_ptr<BasisSet> primary_A = BasisSet::pyconstruct_orbital(molecule_A, options, "BASIS");

    // solve SCF for each monomer   
    outfile->Printf("  ====> Computations for Monomer A <====\n");
    scf_A = oepdev_libutil::solve_scf(molecule_A, primary_A, functional, options, psio);
    outfile->Printf("  ====> Computations for Monomer B <====\n");
    scf_B = oepdev_libutil::solve_scf(molecule_B, primary_B, functional, options, psio);

    // solve CPHF equations for each monomer
    std::shared_ptr<oepdev_libutil::CPHF> cphf_A(new oepdev_libutil::CPHF(scf_A, options));
    cphf_A->compute();
    std::shared_ptr<Matrix> pol_A = cphf_A->get_molecular_polarizability();
    pol_A->print();

    std::shared_ptr<oepdev_libutil::CPHF> cphf_B(new oepdev_libutil::CPHF(scf_B, options));
    cphf_B->compute();
    std::shared_ptr<Matrix> pol_B = cphf_B->get_molecular_polarizability();
    pol_B->print();


    return ref_wfn;
}

} // EndNameSpace oepdev
} // EndNameSpace psi

