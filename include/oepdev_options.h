#ifndef oepdev_include_oepdev_options_h
#define oepdev_include_oepdev_options_h
/** @file oepdev_options.h */

namespace psi{

/** \brief Options for the OEPDev plugin.
 *
 *  @param name name of driver function
 *  @param options psi::Options object
 *  @return true
 */
extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "OEPDEV" || options.read_globals()) {

        // ---> Psi4: General Options  <--- //

        /*- The amount of information printed to the output file -*/
        options.add_int    ("PRINT"                 , 1                          );

        /*- Auxiliary basis set for OEP density fitting -*/                                   
        options.add_str    ("DF_BASIS_OEP"          , "sto-3g"                   );
        /*- Auxiliary basis set for OEP density fitting: Molecule A -*/                                   
        options.add_str    ("DF_BASIS_OEP_A"        , ""                         );
        /*- Auxiliary basis set for OEP density fitting: Molecule B -*/                                   
        options.add_str    ("DF_BASIS_OEP_B"        , ""                         );

        /*- Intermediate basis set for OEP density fitting -*/                                   
        options.add_str    ("DF_BASIS_INT"          , ""                         );

        /*- Type of Density Fitting -*/
        options.add_str    ("OEPDEV_DF_TYPE"        , "DOUBLE"                   );

        // ---> OEPDEV: Main Routine  <--- //

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

        /*- OEPDev test mode to perform -*/
        options.add_str    ("OEPDEV_TEST_MODE"      , "MONOMER"                  );


        // ---> OEPDevSolver Options  <--- //

        /*- MO Localizer for OEP-Based CT Solver -*/
        options.add_str    ("WFN_UNION_LOCALIZER"   , "BOYS"                     );
        options.add_str    ("SOLVER_CT_LOCALIZER"   , "BOYS"                     );

        /*- Type of EFP potential integrals used in CT energy calculations -*/
        options.add_str    ("EFP2_CT_POTENTIAL_INTS", "DMTP"                     );

        /*- Exclude octupoles in CT energy calculations? (relevant if above is "DMTP") -*/
        options.add_bool   ("EFP2_CT_NO_OCTUPOLES"  , true                       );

        /*- Which methods of interaction energy are to be used? -*/
        options.add_bool   ("OEPDEV_SOLVER_EINT_COUL_AO"  , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_COUL_MO"  , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_COUL_ESP" , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_COUL_CAMM", true                 ); 
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_HS"   , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_DDS"  , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_MRW"  , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_OL"   , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_OEP1" , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_OEP2" , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_REP_EFP2" , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_CT_OL"    , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_CT_OEP"   , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EINT_CT_EFP2"  , true                 );
        options.add_bool   ("OEPDEV_SOLVER_EET_V_TDFI_TI" , true                 );
        options.add_int    ("OEPDEV_SOLVER_EET_HOMO"      , 0                    );
        options.add_int    ("OEPDEV_SOLVER_EET_LUMO"      , 0                    );
        options.add_int    ("OEPDEV_SOLVER_EET_HOMO_A"    , 0                    );
        options.add_int    ("OEPDEV_SOLVER_EET_HOMO_B"    , 0                    );
        options.add_int    ("OEPDEV_SOLVER_EET_LUMO_A"    , 0                    );
        options.add_int    ("OEPDEV_SOLVER_EET_LUMO_B"    , 0                    );

        /*- Electronic excited states */
        options.add_double ("OEPDEV_AMPLITUDE_PRINT_THRESHOLD", 0.2              );
        options.add_int    ("EXCITED_STATE"               ,-1                    );
        options.add_int    ("EXCITED_STATE_A"             ,-1                    );
        options.add_int    ("EXCITED_STATE_B"             ,-1                    );                  
        options.add_double ("OSCILLATOR_STRENGTH_THRESHOLD", 0.01                );
        options.add_bool   ("TrCAMM_SYMMETRIZE"           , true                 );
        options.add_bool   ("TI_CIS_SCF_FOCK_MATRIX"      , false                );
        options.add_bool   ("TI_CIS_PRINT_FOCK_MATRIX"    , false                );
        options.add_int    ("PHASE_HOMO_A"                , 1                    );
        options.add_int    ("PHASE_HOMO_B"                , 1                    );
        options.add_int    ("PHASE_LUMO_A"                , 1                    );
        options.add_int    ("PHASE_LUMO_B"                , 1                    );
	options.add_str    ("CIS_TYPE"                    , "DAVIDSON_LIU"       );
        options.add_double ("CIS_SCHWARTZ_CUTOFF"         , 0.0                  );
        options.add_int    ("DAVIDSON_LIU_NROOTS"         , 1                    );
        options.add_double ("DAVIDSON_LIU_CONVER"         , 1.0E-10              );
        options.add_int    ("DAVIDSON_LIU_MAXITER"        , 500                  );
        options.add_str    ("DAVIDSON_LIU_GUESS"          , "RANDOM"             );
        options.add_double ("DAVIDSON_LIU_THRESH_LARGE"   , 1.0E-03              );
        options.add_double ("DAVIDSON_LIU_THRESH_SMALL"   , 1.0E-06              );
        options.add_int    ("DAVIDSON_LIU_SPACE_MAX"      , 200                  );
        options.add_int    ("DAVIDSON_LIU_SPACE_START"    ,-1                    );
        options.add_bool   ("DAVIDSON_LIU_STOP_WHEN_UNCONVERGED", true           );


        // ---> CPHF Global Options  <--- //

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


        // ---> ESP Global Options  <--- //

        /*- ESP: number of points per atom -*/                                      
        options.add_int    ("ESP_NPOINTS_PER_ATOM"  , 1500                       );

        /*- ESP: padding of a sphere enclosing the molecule -*/                     
        options.add_double ("ESP_PAD_SPHERE"        , 10.0                       );

        /*- ESP: vdW radii of atoms used for random point selections [A.U.] -*/            
        options.add_double ("ESP_VDW_RADIUS_C"      ,  3.0                       );
        options.add_double ("ESP_VDW_RADIUS_H"      ,  4.0                       );
        options.add_double ("ESP_VDW_RADIUS_N"      ,  2.4                       );
        options.add_double ("ESP_VDW_RADIUS_O"      ,  5.6                       );
        options.add_double ("ESP_VDW_RADIUS_F"      ,  2.3                       );
        options.add_double ("ESP_VDW_RADIUS_CL"     ,  2.9                       );
        options.add_double ("ESP_VDW_RADIUS_HE"     ,  4.0                       );
        options.add_double ("ESP_VDW_RADIUS_NE"     ,  4.0                       );
        options.add_double ("ESP_VDW_RADIUS_AR"     ,  4.0                       );
        options.add_double ("ESP_VDW_RADIUS_LI"     ,  3.0                       );


        // ---> DMTP Global Options  <--- //
        
	/*- DMTP: convergence level -*/
	options.add_str    ("DMTP_CONVER"           ,  "R5"                      );
	options.add_str    ("DMTP_CONVER_TI_CIS_E34",  "R5"                      );
	options.add_str    ("DMTP_CONVER_TI_CIS_CT" ,  "R5"                      );

        // ---> DmatPol Global Options  <--- //

        /*- DmatPol training mode -*/
        options.add_str    ("DMATPOL_TRAINING_MODE"  , "EFIELD"                  );

        /*- DmatPol field rank -*/
        options.add_int    ("DMATPOL_FIELD_RANK"     , 1                         );

        /*- DmatPol field gradient rank -*/
        options.add_int    ("DMATPOL_GRADIENT_RANK"  , 0                         );

        /*- Number of random samples (fields or sets of test charges) in DmatPol models parameterization -*/
        options.add_int    ("DMATPOL_NSAMPLES"      , 30                         );

        /*- Electric field scale factor in DmatPol models parameterization -*/
        options.add_double ("DMATPOL_FIELD_SCALE"   , 0.01                       );

        /*- Number of test charges per sample in DmatPol models parameterization -*/
        options.add_int    ("DMATPOL_NTEST_CHARGE"   , 1                         );

        /*- Test Charge value in DmatPol models parameterization -*/
        options.add_double ("DMATPOL_TEST_CHARGE"   , 0.001                      );

        /*- Level shifter for DMATPOL -*/
        options.add_double ("DMATPOL_SCALE_1"       , 1.0                        );

        /*- Test electric field for DMATPOL -*/
        options.add_double ("DMATPOL_TEST_FIELD_X"  , 0.000                      );
        options.add_double ("DMATPOL_TEST_FIELD_Y"  , 0.000                      );
        options.add_double ("DMATPOL_TEST_FIELD_Z"  , 0.008                      );

        /*- Output file name for statistical evaluation for DMATPOL -*/
        options.add_str    ("DMATPOL_OUT_STATS"     , "dmatpol.stats.dat"        );

        /*- Do the Ab Initio model calculations or not -*/
        options.add_bool   ("DMATPOL_DO_AB_INITIO"  , false                      );
        options.add_str    ("DMATPOL_OUT_STATS_AB_INITIO" , "dmatpol.stats.abinitio.dat"        );
        options.add_bool   ("DMATPOL_FF_AB_INITIO"  , false                      );
        options.add_double ("DMATPOL_EFIELD_FF_STEP", 0.0005                     );

    }

    return true;
};

} // EndNameSpace psi

#endif // oepdev_include_oepdev_options_h
