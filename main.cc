/*
 * @BEGIN LICENSE
 *
 * oepdev by Bartosz Błasiak (email: blasiak.bartosz@gmail.com)
 *
 * Copyright (c) 2017 Bartosz Błasiak.
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

// Basic
#include <string>

// Defines
#include "include/oepdev_files.h"

// Psi4
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"

// OEPDev
#include "include/oepdev_options.h"
#include "oepdev/liboep/oep.h"
#include "oepdev/libgefp/gefp.h"
#include "oepdev/libsolver/solver.h"
#include "oepdev/libtest/test.h"

// PyBind11
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Typedefs
using SharedWavefunction       = std::shared_ptr<psi::Wavefunction>;
using SharedUnion              = std::shared_ptr<oepdev::WavefunctionUnion>; 
using SharedOEPotential        = std::shared_ptr<oepdev::OEPotential>;
using SharedGEFPFactory        = std::shared_ptr<oepdev::GenEffParFactory>;
using SharedGEFPParameters     = std::shared_ptr<oepdev::GenEffPar>;

namespace psi{ 

// Helpers for Python
void export_dmtp(py::module&);
void export_cphf(py::module&);
void export_solver(py::module&);
void export_util(py::module&);
void export_oep(py::module&);
void export_gefp(py::module&);

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
extern "C" PSI_API
SharedWavefunction oepdev(SharedWavefunction ref_wfn, Options& options)
{
    // ==> Reference information <== //
    oepdev::preambule();

    // ==> Determine what to do <== //
    std::string o_task          = options.get_str  ("OEPDEV_TARGET"        );
    bool        o_local         = options.get_bool ("OEPDEV_LOCALIZE"      );
    int         o_print         = options.get_int  ("PRINT"                );

    // ==> Create OEP's <==
    if (o_task == "OEP_BUILD" || o_task == "OEP") 
    {

        std::string o_oep_build_type= options.get_str  ("OEPDEV_OEP_BUILD_TYPE");

        SharedOEPotential oep;

        if      (o_oep_build_type == "ELECTROSTATIC_ENERGY") 
                   oep = oepdev::OEPotential::build("ELECTROSTATIC ENERGY"  , ref_wfn, options);        
        else if (o_oep_build_type == "REPULSION_ENERGY")
                   oep = oepdev::OEPotential::build("REPULSION ENERGY", ref_wfn, 
				   ref_wfn->get_basisset("BASIS_DF_OEP"), 
				   ref_wfn->get_basisset("BASIS_DF_INT"), 
				   options);
        else if (o_oep_build_type == "CT_ENERGY")
                   oep = oepdev::OEPotential::build("CHARGE TRANSFER ENERGY", ref_wfn, 
				   ref_wfn->get_basisset("BASIS_DF_OEP"), 
				   ref_wfn->get_basisset("BASIS_DF_INT"), 
				   options);
        else if (o_oep_build_type == "EET_COUPLING")
                   oep = oepdev::OEPotential::build("EET COUPLING CONSTANT" , ref_wfn, options);
        else 
             throw psi::PSIEXCEPTION("OEPDEV: Invalid OEP build category selected!");

      //oep->localize(); // add later option for it
        oep->compute();
	oep->print_header();
	oep->print();

    }
    // ===> Density Matrix Susceptibility Tensor Model <=== //
    else if (o_task == "DMATPOL") 
    {

        SharedGEFPFactory factory = oepdev::GenEffParFactory::build("POLARIZATION", ref_wfn, options);
        SharedGEFPParameters parameters = factory->compute();

    } 
    // ==> Create Wavefunction Union of two monomers <==
    else if (o_task == "SOLVER")
    {

        std::string o_solver_type   = options.get_str  ("OEPDEV_SOLVER_TYPE");

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
                                                                                                                            
            if (options.get_bool("OEPDEV_SOLVER_EINT_COUL_AO"  )) solver->compute_benchmark("AO_EXPANDED"    );
            if (options.get_bool("OEPDEV_SOLVER_EINT_COUL_MO"  )) solver->compute_benchmark("MO_EXPANDED"    );
            if (options.get_bool("OEPDEV_SOLVER_EINT_COUL_ESP" )) solver->compute_oep_based("ESP_SYMMETRIZED");
	    if (options.get_bool("OEPDEV_SOLVER_EINT_COUL_CAMM")) solver->compute_oep_based("CAMM"           );
            wfn_union->clear_dpd();
        }
        else if (o_solver_type == "REPULSION_ENERGY") 
        {
            std::shared_ptr<oepdev::OEPDevSolver> solver = oepdev::OEPDevSolver::build("REPULSION ENERGY", wfn_union);
                                                                                                                            
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_HS"  )) solver->compute_benchmark("HAYES_STONE"          );
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_DDS" )) solver->compute_benchmark("DDS"                  );
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_MRW" )) solver->compute_benchmark("MURRELL_ETAL"         );
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_OL"  )) solver->compute_benchmark("OTTO_LADIK"           );
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_OEP1")) solver->compute_oep_based("MURRELL_ETAL_GDF_ESP" );
	    if (options.get_bool("OEPDEV_SOLVER_EINT_REP_OEP2")) solver->compute_oep_based("MURRELL_ETAL_GDF_CAMM");
            if (options.get_bool("OEPDEV_SOLVER_EINT_REP_EFP2")) solver->compute_benchmark("EFP2"                 );        
            wfn_union->clear_dpd();
        } 
        else if (o_solver_type == "CHARGE_TRANSFER_ENERGY") 
        {
            std::shared_ptr<oepdev::OEPDevSolver> solver = oepdev::OEPDevSolver::build("CHARGE TRANSFER ENERGY", wfn_union);
                                                                                                                            
            if (options.get_bool("OEPDEV_SOLVER_EINT_CT_OL"  )) solver->compute_benchmark("OTTO_LADIK" );
            if (options.get_bool("OEPDEV_SOLVER_EINT_CT_OEP" )) solver->compute_oep_based("OTTO_LADIK" );
            if (options.get_bool("OEPDEV_SOLVER_EINT_CT_EFP2")) solver->compute_benchmark("EFP2"       );        
            wfn_union->clear_dpd();
        } 
        else if (o_solver_type == "EET_COUPLING")
        {
            std::shared_ptr<oepdev::OEPDevSolver> solver = oepdev::OEPDevSolver::build("EET COUPLING CONSTANT", wfn_union);
            if (options.get_bool("OEPDEV_SOLVER_EET_V_TDFI_TI")) solver->compute_benchmark("FUJIMOTO_TI_CIS");
        }
        else 
           throw PSIEXCEPTION("Incorrect solver type chosen!\n");

	// ==> Clear the DPD file for integral transformations <== //
	//wfn_union->clear_dpd();

    } else if (o_task == "TEST") {

            oepdev::test::Test test(ref_wfn, options);
            double result = test.run();

    } else 
           throw PSIEXCEPTION("Incorrect target for oepdev program!\n");

    return ref_wfn;
}


// Python interface
PYBIND11_MODULE(oepdev, m) {

  // Export OEPDev-Psi4 Interface
  m.doc() = "OEPDev: Quantum Chemistry for Extended Molecular Aggregates. Python Interface.";
  m.def("oepdev", &oepdev, "Run OEPDev plugin: Python interface.");
  m.def("read_options", &read_options, "Read options for OEPDev plugin. Python interface");

  // Export OEPDev Functionalities
  export_dmtp(m);
  // Export CHPF
  export_cphf(m);
  // Export OEPDevSolver and WavefunctionUnion
  export_solver(m);
  // Export Utilities
  export_util(m);
  // Export OEP's
  export_oep(m);
  // Export GEFP
  export_gefp(m);
}

} // EndNameSpace psi

