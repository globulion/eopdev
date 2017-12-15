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
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"

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

    // Wavefunction coefficients for isolated monomers
    std::shared_ptr<Matrix> c_A = scf_A->Ca_subset("AO","ALL");
    std::shared_ptr<Matrix> c_B = scf_B->Ca_subset("AO","ALL");

    // compute AO-ERI for the dimer
    std::shared_ptr<IntegralFactory> ints_dimer(new IntegralFactory(primary, primary, primary, primary));
    std::shared_ptr<TwoBodyAOInt>    eri_dimer(ints_dimer->eri());
    AOShellCombinationsIterator shellIter(primary, primary, primary, primary);
    const double * buffer = eri_dimer->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
         eri_dimer->compute_shell(shellIter);
         outfile->Printf("( %d %d | %d %d )\n", 
                           shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
    }

    // Towards making union of molecules
    std::vector<int> orbitals_A;
    std::vector<int> orbitals_B;
    std::vector<int> indices;
    for (int i=0; i<scf_A->nmopi()[0]; i++) orbitals_A.push_back(i);
    for (int i=0; i<scf_B->nmopi()[0]; i++) orbitals_B.push_back(i);
    std::shared_ptr<MOSpace> spaceA(new MOSpace('X', orbitals_A, indices));
    std::shared_ptr<MOSpace> spaceB(new MOSpace('Y', orbitals_B, indices));

    std::vector<std::shared_ptr<MOSpace>> spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    spaces.push_back(spaceA);
    spaces.push_back(spaceB);
    std::shared_ptr<IntegralTransform> transform(new IntegralTransform(ref_wfn, spaces, IntegralTransform::Restricted));
    transform->set_keep_dpd_so_ints(1);
    // Trans (AA|AA)
    timer_on("Trans (AA|AA)");
    transform->transform_tei(spaceA, spaceA, spaceA, spaceA, IntegralTransform::MakeAndKeep);
    timer_off("Trans (AA|AA)");

    // Trans (AA|AB)
    timer_on("Trans (AA|AB)");
    transform->transform_tei(spaceA, spaceA, spaceA, spaceB, IntegralTransform::ReadAndKeep);
    timer_off("Trans (AA|AB)");

    // Trans (AA|BB)
    timer_on("Trans (AA|BB)");
    transform->transform_tei(spaceA, spaceA, spaceB, spaceB, IntegralTransform::ReadAndNuke);
    timer_off("Trans (AA|BB)");

    // Trans (AB|AB)
    timer_on("Trans (AB|AB)");
    transform->transform_tei(spaceA, spaceB, spaceA, spaceB, IntegralTransform::MakeAndKeep);
    timer_off("Trans (AB|AB)");

    // Trans (AB|BB)
    timer_on("Trans (AB|BB)");
    transform->transform_tei(spaceA, spaceB, spaceB, spaceB, IntegralTransform::ReadAndNuke);
    timer_off("Trans (AB|BB)");

    // Trans (BB|BB)
    timer_on("Trans (BB|BB)");
    transform->transform_tei(spaceB, spaceB, spaceB, spaceB);
    timer_off("Trans (BB|BB)");


    // Read the integrals
    dpd_set_default(transform->get_dpd_id());
    dpdbuf4 ABBB, AAAB;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&ABBB, PSIF_LIBTRANS_DPD, 0, // (AB|BB)
                           transform->DPD_ID("[X,Y]"), transform->DPD_ID("[Y,Y]"),
                           transform->DPD_ID("[X,Y]"), transform->DPD_ID("[Y,Y]"), 0, "MO Ints (XY|YY)");
    global_dpd_->buf4_init(&AAAB, PSIF_LIBTRANS_DPD, 0, // (AA|AB)
                           transform->DPD_ID("[X,X]"), transform->DPD_ID("[X,Y]"),
                           transform->DPD_ID("[X,X]"), transform->DPD_ID("[X,Y]"), 0, "MO Ints (XX|XY)");

    //
    global_dpd_->buf4_close(&ABBB);
    global_dpd_->buf4_close(&AAAB);

    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);




    
   
    // compute ERI of the type (AB|BB), (AA|AB) and (AA|BB)
    //std::shared_ptr<IntegralFactory> factory_ABBB(new IntegralFactory(primary_A, primary_B, primary_B, primary_B));
    //std::shared_ptr<IntegralFactory> factory_AAAB(new IntegralFactory(primary_A, primary_A, primary_A, primary_B));
    //std::shared_ptr<IntegralFactory> factory_AABB(new IntegralFactory(primary_A, primary_A, primary_B, primary_B));

    //std::shared_ptr<Matrix> sao_AB(new Matrix("Overlap matrix AB", primary_A->nbf(), primary_B->nbf()));
    //std::shared_ptr<OneBodyAOInt> ints_AB  (factory_ABBB->ao_overlap());
    //std::shared_ptr<TwoBodyAOInt> ints_ABBB(factory_ABBB->eri());
    //ints_AB->compute(sao_AB);
    //sao_AB->print();


    //std::vector<std::shared_ptr<MOSpace> > spaces;
    //spaces.push_back(MOSpace::all);
    //IntegralTransform tran(c_A, c_B, c_B, c_B, spaces, IntegralTransform::Restricted);
    //tran.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

    return ref_wfn;
}

} // EndNameSpace oepdev
} // EndNameSpace psi

