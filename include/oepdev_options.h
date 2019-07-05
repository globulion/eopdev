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

        /*- Intermediate basis set for OEP density fitting -*/                                   
        options.add_str    ("DF_BASIS_INT"          , "aug-cc-pVDZ-jkfit"        );

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
        options.add_str    ("SOLVER_CT_LOCALIZER"   , "BOYS"                     );

        /*- Type of EFP potential integrals used in CT energy calculations -*/
        options.add_str    ("EFP2_CT_POTENTIAL_INTS", "DMTP"                     );


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


        // ---> DMTP Global Options  <--- //
        
	/*- DMTP: convergence level -*/
	options.add_str    ("DMTP_CONVER"           ,  "R5"                      );

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
