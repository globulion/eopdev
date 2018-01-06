/*
 * @BEGIN LICENSE
 *
 * oepdev by Bartosz Błasiak (email: blasiak.bartosz@gmail.com)
 *
 * Copyright (c) 2017 Bartosz Błasiak.
 * NOTE: oepdev is not open-source and cannot be redistributed, copied and downloaded
 * without a written consent of the Project Administrator.
 * 
 * oepdev is a plugin to:
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
#include "psi4/libdpd/dpd.h"


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
using SharedLocalizer          = std::shared_ptr<Localizer>;

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
    // ==> Reference information <== //
    oepdev_libutil::preambule();

    // ==> Determine what to do <== //
    int print = options.get_int("PRINT");

    // ==> Psi4 Input/Output Stream <== //
    SharedPSIO psio = PSIO::shared_object();

    // ==> Create Wavefunction Union of two monomers <==
    SharedUnion wfn_union = std::make_shared<oepdev_libutil::WavefunctionUnion>(ref_wfn, options);

    // ==> Localize Molecular orbitals of the Union (optionally) <== //
    wfn_union->localize_orbitals();

    // ==> Perform the integral transformation to MO basis <== //
    wfn_union->transform_integrals();



    /* Below there are a few tests needed to develop Initial Version of the Plugin.
     * At the end of the process, all functionalities of the plugin should be using only:
     *  - the WavefunctionUnion object
     *  - the PSIO object and
     *  - the Options object.
     * Therefore, these tests shall be removed once particular functionalities
     * will be developed during the Project. This means that the includes and typedefs 
     * at the top of `main.cc` shall be removed.
     */

    SharedWavefunction      wfn_union_base = wfn_union;
    SharedIntegralTransform transform      = wfn_union->integrals();

    // Parse molecules, fragments, basis sets and other primary informations
    SharedBasisSet          primary_1      = wfn_union->l_primary(0);    
    SharedBasisSet          primary_2      = wfn_union->l_primary(1);
    SharedMolecule          molecule_1     = wfn_union->l_molecule(0);
    SharedMolecule          molecule_2     = wfn_union->l_molecule(1);
    SharedMolecule          molecule       = wfn_union->molecule();
    SharedBasisSet          primary        = wfn_union->basisset();
    SharedWavefunction      scf_1          = wfn_union->l_wfn(0); 
    SharedWavefunction      scf_2          = wfn_union->l_wfn(1);
    
    // Solve CPHF equations for each monomer
    SharedCPHF cphf_1(new oepdev_libutil::CPHF(scf_1, options));
    SharedCPHF cphf_2(new oepdev_libutil::CPHF(scf_2, options));
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
    SharedIntegralFactory ints_dimer = std::make_shared<IntegralFactory>(primary, primary, primary, primary);
    SharedTwoBodyAOInt    eri_dimer(ints_dimer->eri());
    AOShellCombinationsIterator shellIter(primary, primary, primary, primary);
    const double * buffer = eri_dimer->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
         eri_dimer->compute_shell(shellIter);
         outfile->Printf("( %d %d | %d %d )\n", 
                           shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
    }
    }

    // Read the MO integrals
    dpd_set_default(transform->get_dpd_id());
    dpdbuf4 buf_1122, buf_1222, buf_1112, buf_1212;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio->tocprint(PSIF_LIBTRANS_DPD);

    global_dpd_->buf4_init(&buf_1112, PSIF_LIBTRANS_DPD, 0, 
                           transform->DPD_ID("[1,1]"  ), transform->DPD_ID("[1,2]"  ),
                           transform->DPD_ID("[1>=1]+"), transform->DPD_ID("[1,2]"  ), 0, "MO Ints (11|12)");
    global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0, 
                           transform->DPD_ID("[1,1]"  ), transform->DPD_ID("[2,2]"  ),
                           transform->DPD_ID("[1>=1]+"), transform->DPD_ID("[2>=2]+"), 0, "MO Ints (11|22)");
    global_dpd_->buf4_init(&buf_1212, PSIF_LIBTRANS_DPD, 0, 
                           transform->DPD_ID("[1,2]"  ), transform->DPD_ID("[1,2]"  ),
                           transform->DPD_ID("[1,2]"  ), transform->DPD_ID("[1,2]"  ), 0, "MO Ints (12|12)");
    global_dpd_->buf4_init(&buf_1222, PSIF_LIBTRANS_DPD, 0,
                           transform->DPD_ID("[1,2]"  ), transform->DPD_ID("[2,2]"  ),
                           transform->DPD_ID("[1,2]"  ), transform->DPD_ID("[2>=2]+"), 0, "MO Ints (12|22)");

    outfile->Printf("\n <=== buf_1112 MO Integrals ===>\n\n");
    for (int h = 0; h < wfn_union->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_1112, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_1112, h);
         for (int pq = 0; pq < buf_1112.params->rowtot[h]; ++pq) {
              int p = buf_1112.params->roworb[h][pq][0];
              int q = buf_1112.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_1112.params->coltot[h]; ++rs) {
                   int r = buf_1112.params->colorb[h][rs][0];
                   int s = buf_1112.params->colorb[h][rs][1];
                   outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1112.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_1112, h);
    }

    outfile->Printf("\n <=== buf_1122 MO Integrals ===>\n\n");
    for (int h = 0; h < wfn_union->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
         for (int pq = 0; pq < buf_1122.params->rowtot[h]; ++pq) {
              int p = buf_1122.params->roworb[h][pq][0];
              int q = buf_1122.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_1122.params->coltot[h]; ++rs) {
                   int r = buf_1122.params->colorb[h][rs][0];
                   int s = buf_1122.params->colorb[h][rs][1];
                   outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1122.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
    }

    outfile->Printf("\n <=== buf_1222 MO Integrals ===>\n\n");
    for (int h = 0; h < wfn_union->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_1222, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_1222, h);
         for (int pq = 0; pq < buf_1222.params->rowtot[h]; ++pq) {
              int p = buf_1222.params->roworb[h][pq][0];
              int q = buf_1222.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_1222.params->coltot[h]; ++rs) {
                   int r = buf_1222.params->colorb[h][rs][0];
                   int s = buf_1222.params->colorb[h][rs][1];
                   outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1222.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_1222, h);
    }

    outfile->Printf("\n <=== buf_1222 MO Integrals ===>\n\n");
    for (int h = 0; h < wfn_union->nirrep(); ++h) {
         global_dpd_->buf4_mat_irrep_init(&buf_1222, h);
         global_dpd_->buf4_mat_irrep_rd(&buf_1222, h);
         for (int pq = 0; pq < buf_1222.params->rowtot[h]; ++pq) {
              int p = buf_1222.params->roworb[h][pq][0];
              int q = buf_1222.params->roworb[h][pq][1];
              for (int rs = 0; rs < buf_1222.params->coltot[h]; ++rs) {
                   int r = buf_1222.params->colorb[h][rs][0];
                   int s = buf_1222.params->colorb[h][rs][1];
                   outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1222.matrix[h][pq][rs]);
              }
         }
         global_dpd_->buf4_mat_irrep_close(&buf_1222, h);
    }

    //
    global_dpd_->buf4_close(&buf_1122);
    global_dpd_->buf4_close(&buf_1222);
    global_dpd_->buf4_close(&buf_1112);

    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    return ref_wfn;
}

} // EndNameSpace oepdev
} // EndNameSpace psi

