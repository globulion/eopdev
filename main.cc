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

#include <cstdlib>
#include <cstdio>
#include <string>

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"

#include "include/oepdev_files.h"

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

#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"

using SharedMolecule           = std::shared_ptr<Molecule>;                            
using SharedSuperFunctional    = std::shared_ptr<SuperFunctional>;
using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedVector             = std::shared_ptr<Vector>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedBasisSet           = std::shared_ptr<BasisSet>;
using SharedUnion              = std::shared_ptr<oepdev_libutil::WavefunctionUnion>; 
using SharedPSIO               = std::shared_ptr<PSIO>;
using SharedCPHF               = std::shared_ptr<oepdev_libutil::CPHF>;
using SharedMOSpace            = std::shared_ptr<MOSpace>;
using SharedIntegralTransform  = std::shared_ptr<IntegralTransform>;
using SharedIntegralFactory    = std::shared_ptr<IntegralFactory>;
using SharedTwoBodyAOInt       = std::shared_ptr<TwoBodyAOInt>;
using SharedMOSpaceVector      = std::vector<std::shared_ptr<MOSpace>>;
using intVector                = std::vector<int>;

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

    // Psi4 Input/Output Stream
    SharedPSIO              psio           = PSIO::shared_object();

    // Create Wavefunction Union of two monomers
    SharedUnion             wfn_union(     new oepdev_libutil::WavefunctionUnion(ref_wfn, options));
    SharedWavefunction      wfn_union_base = static_cast<SharedWavefunction>(wfn_union);

    // Parse molecules, fragments, basis sets and other primary informations
    SharedBasisSet          primary_1      = wfn_union->primary_1();    
    SharedBasisSet          primary_2      = wfn_union->primary_2();
    SharedMolecule          molecule_1     = wfn_union->molecule_1();
    SharedMolecule          molecule_2     = wfn_union->molecule_2();
    SharedMolecule          molecule       = wfn_union->molecule();
    SharedBasisSet          primary        = wfn_union->basisset();
    SharedWavefunction      scf_1          = wfn_union->wfn_1(); 
    SharedWavefunction      scf_2          = wfn_union->wfn_2();
    
    //TRIAL create BasisSet objects
    //TRIAL std::shared_ptr<BasisSet> primary_A = BasisSet::pyconstruct_orbital(molecule_A, options, "BASIS");

    // Solve CPHF equations for each monomer
    SharedCPHF cphf_1(new oepdev_libutil::CPHF(scf_1, options));
    SharedCPHF cphf_2(new oepdev_libutil::CPHF(scf_1, options));
    cphf_1->compute();
    cphf_2->compute();
    SharedMatrix pol_1 = cphf_1->get_molecular_polarizability();
    SharedMatrix pol_2 = cphf_2->get_molecular_polarizability();
    pol_1->print(); pol_2->print();

    // Wavefunction coefficients for isolated monomers
    SharedMatrix c_1 = scf_1->Ca_subset("AO","ALL");
    SharedMatrix c_2 = scf_2->Ca_subset("AO","ALL");

    // compute AO-ERI for the dimer
    if (false) {
    SharedIntegralFactory ints_dimer(new IntegralFactory(primary, primary, primary, primary));
    SharedTwoBodyAOInt    eri_dimer(ints_dimer->eri());
    AOShellCombinationsIterator shellIter(primary, primary, primary, primary);
    const double * buffer = eri_dimer->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
         eri_dimer->compute_shell(shellIter);
         outfile->Printf("( %d %d | %d %d )\n", 
                           shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
    }
    }

    // Perform the integral transformation to MO basis
    intVector orbitals_1, orbitals_2, indices;
    for (int i=0; i<scf_1->nmopi()[0]; i++) orbitals_1.push_back(i);
    for (int i=0; i<scf_2->nmopi()[0]; i++) orbitals_2.push_back(i);
    SharedMOSpace space1(new MOSpace('1', orbitals_1, indices));
    SharedMOSpace space2(new MOSpace('2', orbitals_2, indices));

    SharedMOSpaceVector spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    spaces.push_back(space1);
    spaces.push_back(space2);
    SharedIntegralTransform transform(new IntegralTransform(wfn_union_base, spaces, IntegralTransform::Restricted));

    transform->set_keep_dpd_so_ints(1);
    // Trans (AA|AA)
    timer_on("Trans (AA|AA)");
    transform->transform_tei(space1, space1, space1, space1, IntegralTransform::MakeAndKeep);
    timer_off("Trans (AA|AA)");

    // Trans (AA|AB)
    timer_on("Trans (AA|AB)");
    transform->transform_tei(space1, space1, space1, space2, IntegralTransform::ReadAndKeep);
    timer_off("Trans (AA|AB)");

    // Trans (AA|BB)
    timer_on("Trans (AA|BB)");
    transform->transform_tei(space1, space1, space2, space2, IntegralTransform::ReadAndNuke);
    timer_off("Trans (AA|BB)");

    // Trans (AB|AB)
    timer_on("Trans (AB|AB)");
    transform->transform_tei(space1, space2, space1, space2, IntegralTransform::MakeAndKeep);
    timer_off("Trans (AB|AB)");

    // Trans (AB|BB)
    timer_on("Trans (AB|BB)");
    transform->transform_tei(space1, space2, space2, space2, IntegralTransform::ReadAndNuke);
    timer_off("Trans (AB|BB)");

    // Trans (BB|BB)
    timer_on("Trans (BB|BB)");
    transform->transform_tei(space2, space2, space2, space2);
    timer_off("Trans (BB|BB)");


    // Read the integrals
    dpd_set_default(transform->get_dpd_id());
    dpdbuf4 ABBB, AAAB;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&ABBB, PSIF_LIBTRANS_DPD, 0, // (AB|BB)
                           transform->DPD_ID("[1,2]"), transform->DPD_ID("[2,2]"),
                           transform->DPD_ID("[1,2]"), transform->DPD_ID("[2,2]"), 0, "MO Ints (12|22)");
    global_dpd_->buf4_init(&AAAB, PSIF_LIBTRANS_DPD, 0, // (AA|AB)
                           transform->DPD_ID("[1,1]"), transform->DPD_ID("[1,2]"),
                           transform->DPD_ID("[1,1]"), transform->DPD_ID("[1,2]"), 0, "MO Ints (11|12)");

    //
    global_dpd_->buf4_close(&ABBB);
    global_dpd_->buf4_close(&AAAB);

    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    return ref_wfn;
}

} // EndNameSpace oepdev
} // EndNameSpace psi

