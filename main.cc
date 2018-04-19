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
/** @file main.cc */

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
#include "oepdev/libutil/wavefunction_union.h"
#include "oepdev/libutil/integrals_iter.h"
#include "oepdev/libutil/space3d.h"
#include "oepdev/libutil/esp.h"
#include "oepdev/libutil/solver.h"
#include "oepdev/liboep/oep.h"
#include "oepdev/libpsi/potential.h"

#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"

#include "oepdev/libints/eri.h"
#include "oepdev/libpsi/integral.h"

#include "oepdev/libtest/test.h"


using SharedMolecule           = std::shared_ptr<Molecule>;                            
using SharedSuperFunctional    = std::shared_ptr<SuperFunctional>;
using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedVector             = std::shared_ptr<Vector>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedBasisSet           = std::shared_ptr<BasisSet>;
using SharedUnion              = std::shared_ptr<oepdev::WavefunctionUnion>; 
using SharedPSIO               = std::shared_ptr<PSIO>;
using SharedCPHF               = std::shared_ptr<oepdev::CPHF>;
using SharedMOSpace            = std::shared_ptr<MOSpace>;
using SharedIntegralTransform  = std::shared_ptr<IntegralTransform>;
using SharedIntegralFactory    = std::shared_ptr<IntegralFactory>;
using SharedTwoBodyAOInt       = std::shared_ptr<TwoBodyAOInt>;
using SharedMOSpaceVector      = std::vector<std::shared_ptr<MOSpace>>;
using intVector                = std::vector<int>;
using SharedLocalizer          = std::shared_ptr<Localizer>;
using SharedOEPotential        = std::shared_ptr<oepdev::OEPotential>;
using SharedField3D            = std::shared_ptr<oepdev::ScalarField3D>;

namespace psi{ 


/** \brief Options for the OEPDev plugin.
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
        options.add_int    ("PRINT"                 , 1                          );
        /*- Basis set for OEP density fitting -*/                                   
        options.add_str    ("DF_BASIS_OEP"          , "sto-3g"                   );
        /*- CPHF maximum iterations -*/                                             
        options.add_int    ("CPHF_MAXITER"          , 50                         );
        /*- CPHF convergence -*/                                                    
        options.add_double ("CPHF_CONVER"           , 1.0E-8                     );
        /*- Whether use DIIS for CPHF -*/                                           
        options.add_bool   ("CPHF_DIIS"             , true                       );
        /*- Size of DIIS subspace for CPHF -*/                                      
        options.add_int    ("CPHF_DIIS_DIM"         , 3                          );
        /*- Localizer for CPHF routine -*/                                          
        options.add_bool   ("CPHF_LOCALIZE"         , true                       );
        /*- Localizer for CPHF routine -*/                                          
        options.add_str    ("CPHF_LOCALIZER"        , "BOYS"                     );
        /*- ESP: number of points per atom -*/                                      
        options.add_int    ("ESP_NPOINTS_PER_ATOM"  , 1500                       );
        /*- ESP: padding of a sphere enclosing the molecule -*/                     
        options.add_double ("ESP_PAD_SPHERE"        , 10.0                       );
        /*- Target: What to do? -*/
        options.add_str    ("OEPDEV_TARGET"         , "SOLVER"                   );
        /*- Whether localize MO's or not -*/
        options.add_bool   ("OEPDEV_LOCALIZE"       , false                      );
        /*- Whether enable trial tests in main.cc or not -*/
        options.add_bool   ("OEPDEV_ENABLE_TRIAL"   , false                      );
        /*- Which OEP to build? -*/
        options.add_str    ("OEPDEV_OEP_BUILD_TYPE" , "ELECTROSTATIC_ENERGY"     );
        /*- Which Solver to use? -*/
        options.add_str    ("OEPDEV_SOLVER_TYPE"    , "REPULSION_ENERGY"         );
        /*- OEPDev test to perform -*/
        options.add_str    ("OEPDEV_TEST_NAME"      , ""                         );
        /*- Number of random electric fields in DmatPol models parameterization -*/
        options.add_int    ("DMATPOL_NFIELDS"       , 30                         );
        /*- Electric field scale factor in DmatPol models parameterization -*/
        options.add_double ("DMATPOL_FIELD_SCALE"   , 0.01                       );
     }

    return true;
}

/** \brief Main routine of the OEPDev plugin.
 *
 *  Created with intention to test various models of the interaction energy 
 *  between two molecules, described by the Hartree-Fock-Roothaan-Hall theory 
 *  or the configuration interaction with singles theory. 
 *
 *  In particular, the plugin tests the models of:
 *
 *   1. the Pauli repulsion and CT interaction energy      (Project II ) 
 *   2. the Induction interaction energy                   (Project III)
 *   3. the excitation energy transfer couplings           (Project I  )
 *
 *  against benchmarks (exact or reference solutions). 
 *  The list of implemented models can be found in \ref pimplementedmodels .
 *
 *  @param ref_wfn shared wavefunction of a dimer
 *  @param options psi::Options object
 *  @return psi::SharedWavefunction (either ref_wfn or wavefunction union)
 */
extern "C"
SharedWavefunction oepdev(SharedWavefunction ref_wfn, Options& options)
{
    // ==> Reference information <== //
    oepdev::preambule();

    // ==> Determine what to do <== //
    std::string o_task          = options.get_str  ("OEPDEV_TARGET"        );
    std::string o_oep_build_type= options.get_str  ("OEPDEV_OEP_BUILD_TYPE");
    std::string o_solver_type   = options.get_str  ("OEPDEV_SOLVER_TYPE"   );
    bool        o_local         = options.get_bool ("OEPDEV_LOCALIZE"      );
    bool        o_enable_trial  = options.get_bool ("OEPDEV_ENABLE_TRIAL"  );
    int         o_print         = options.get_int  ("PRINT"                );

    // ==> Psi4 Input/Output Stream <== //
    SharedPSIO psio = PSIO::shared_object();

    // ==> Create OEP's <==
    if (o_task == "OEP_BUILD" || o_task == "OEP") 
    {

        SharedOEPotential oep;

        if      (o_oep_build_type == "ELECTROSTATIC_ENERGY") 
                   oep = oepdev::OEPotential::build("ELECTROSTATIC ENERGY"  , ref_wfn, options);        
        else if (o_oep_build_type == "REPULSION_ENERGY")
                   oep = oepdev::OEPotential::build("REPULSION ENERGY"      , ref_wfn, ref_wfn->basisset(), options);
        else if (o_oep_build_type == "CT_ENERGY")
                   oep = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY", ref_wfn, ref_wfn->basisset(), options);
        else if (o_oep_build_type == "EET_COUPLING")
                   oep = oepdev::OEPotential::build("EET COUPLING CONSTANT" , ref_wfn, options);
        else 
             throw psi::PSIEXCEPTION("OEPDEV: Invalid OEP build category selected!");

        oep->compute();

    }
    // ==> Create Wavefunction Union of two monomers <==
    else if (o_task == "SOLVER")
    {

        SharedUnion wfn_union = std::make_shared<oepdev::WavefunctionUnion>(ref_wfn, options); 
        wfn_union->print_header();

        // ==> Localize Molecular orbitals of the Union (optionally) <== //                                                 
        if (o_local) wfn_union->localize_orbitals();
                                                                                                                            
        // ==> Perform the integral transformation to MO basis <== //
        wfn_union->transform_integrals();
        if (o_print > 2) wfn_union->print_mo_integrals();
                                                                                                                            
        // ==> Perform the work <== //
        if (o_solver_type == "ELECTROSTATIC_ENERGY") 
        {
            std::shared_ptr<oepdev::OEPDevSolver> solver = oepdev::OEPDevSolver::build("ELECTROSTATIC ENERGY", wfn_union); 
                                                                                                                            
            double el1 = solver->compute_benchmark("AO_EXPANDED"    );
            double el2 = solver->compute_benchmark("MO_EXPANDED"    );
            double el3 = solver->compute_oep_based("ESP_SYMMETRIZED");
        }
        else if (o_solver_type == "REPULSION_ENERGY") 
        {
            std::shared_ptr<oepdev::OEPDevSolver> solver = oepdev::OEPDevSolver::build("REPULSION ENERGY", wfn_union);
                                                                                                                            
            double e_stone = solver->compute_benchmark("HAYES_STONE"  );
            double e_dens  = solver->compute_benchmark("DENSITY_BASED");
            double e_murr  = solver->compute_benchmark("MURRELL_ETAL" );
            double e_oep1  = solver->compute_oep_based("MURRELL_ETAL_MIX");
            double e_otla  = solver->compute_benchmark("OTTO_LADIK"   );
            double e_efp2  = solver->compute_benchmark("EFP2"         );        
        } 
        else 
           throw PSIEXCEPTION("Incorrect solver type chosen!\n");

    } else if (o_task == "TEST") {
            oepdev::test::Test test(ref_wfn, options);
            double result = test.run();
    } else 
           throw PSIEXCEPTION("Incorrect target for oepdev program!\n");



        /* Below there are a few tests needed to develop Initial Version of the Plugin.                                                      
         * At the end of the process, all functionalities of the plugin should be using only:
         *  - the WavefunctionUnion object
         *  - the Options object.
         * Therefore, these tests shall be removed once particular functionalities
         * will be developed during the Project. This means that the includes and typedefs 
         * at the top of `main.cc` shall be removed.
         */
        if (o_enable_trial) {


        //SharedWavefunction      wfn_union_base = wfn_union;
        //SharedIntegralTransform transform      = wfn_union->integrals();
        //                                                                                                                                     
        //// Parse molecules, fragments, basis sets and other primary informations
        //SharedBasisSet          primary_1      = wfn_union->l_primary(0);    
        //SharedBasisSet          primary_2      = wfn_union->l_primary(1);
        //SharedMolecule          molecule_1     = wfn_union->l_molecule(0);
        //SharedMolecule          molecule_2     = wfn_union->l_molecule(1);
        //SharedMolecule          molecule       = wfn_union->molecule();
        //SharedBasisSet          primary        = wfn_union->basisset();
        //SharedBasisSet          auxiliary_1    = wfn_union->l_auxiliary(0);
        //SharedBasisSet          auxiliary_2    = wfn_union->l_auxiliary(1);
        //SharedWavefunction      scf_1          = wfn_union->l_wfn(0); 
        //SharedWavefunction      scf_2          = wfn_union->l_wfn(1);
        //
        //// Solve CPHF equations for each monomer
        //SharedCPHF cphf_1(new oepdev::CPHF(scf_1, options));
        //SharedCPHF cphf_2(new oepdev::CPHF(scf_2, options));
        //cphf_1->compute();
        //cphf_2->compute();
        //SharedMatrix pol_1 = cphf_1->get_molecular_polarizability();
        //SharedMatrix pol_2 = cphf_2->get_molecular_polarizability();
        //pol_1->print(); pol_2->print();
        //                                                                                                                                     
        //// Wavefunction coefficients for isolated monomers
        //SharedMatrix c_1 = scf_1->Ca_subset("AO","ALL");
        //SharedMatrix c_2 = scf_2->Ca_subset("AO","ALL");
        //                                                                                                                                     
        //// Create some OEP's
        //SharedOEPotential oep_cou = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", scf_1, options);
        //SharedOEPotential oep_rep = oepdev::OEPotential::build("REPULSION ENERGY", scf_1, primary_1, options);
        //SharedOEPotential oep_eet = oepdev::OEPotential::build("EET COUPLING", scf_1, options);
        //                                                                                                                                     
        //oep_cou->write_cube("V", "oep");
        //                                                                                                                                     
        //// Compute potentials
        //if (false) {
        //SharedField3D potential_cube = oepdev::ScalarField3D::build("ELECTROSTATIC", 60, 60, 60, 10.0, 10.0, 10.0, scf_1, options);
        //SharedField3D potential_random = oepdev::ScalarField3D::build("ELECTROSTATIC", 50000, 10.0, scf_1, options);
        //                                                                                                                                     
        //potential_cube->print();
        //potential_random->print();
        //potential_random->compute();
        //                                                                                                                                     
        //oepdev::ESPSolver esp(potential_random);
        //esp.compute();
        //esp.charges()->print();
        ////std::shared_ptr<oepdev::CubeDistribution3D> di = std::dynamic_pointer_cast<oepdev::CubeDistribution3D>(potential->distribution());
        ////di->print_header();
        ////potential->write_cube_file("pot");
        //}
        } // End of trial

    return ref_wfn;
}

} // EndNameSpace psi

