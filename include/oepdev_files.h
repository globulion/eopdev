#ifndef oepdev_include_oepdev_files_h
#define oepdev_include_oepdev_files_h
/** @file oepdev_files.h */

/** \addtogroup OEPDEV_UTILITIES
 *@{
 */
/// Use DIIS from Psi4 (1) or OEPDev (0)?
#define OEPDEV_USE_PSI4_DIIS_MANAGER      0
/// L_max
#define OEPDEV_MAX_AM                     8
/// 2L_max+1
#define OEPDEV_N_MAX_AM                   17
/// ERI criterion for E12, E34, E123 and lambda*EXY coefficients
#define OEPDEV_CRIT_ERI                   1e-8
/// Size of R buffer (OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*3)
#define OEPDEV_SIZE_BUFFER_R              250563 
/// Size of D2 buffer (3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1)*OEPDEV_N_MAX_AM)
#define OEPDEV_SIZE_BUFFER_D2             3264

/** @}*/
#endif // oepdev_include_oepdev_files_h

/**
 * \namespace oepdev
 * \brief OEPDev module namespace.
 * 
 * Contains:
 */

/**
 * \namespace psi
 * \brief Psi4 package namespace.
 *
 * Contains all Psi4 functionalities.
 */

/**
 *  @defgroup OEPDEV_SOLVERS The OEPDev solver library
 *  Implementations various solvers for molecular properties
 *  as a functions of unperturbed monomeric wavefunctions.
 *  This is the place all target OEP-based models are implemented
 *  and compared with benchmark models.
 */

/** 
 *  @defgroup OEPDEV_OEPS The Generalized One-Electron Potentials library
 *  Implements the goal of this project: The Generalized One-Electron Potentials
 *  (OEP's). You will find here OEP's for computation of Pauli repulsion energy,
 *  charge-transfer energy and others.
 */

/**
 *  @defgroup OEPDEV_LIBINTS The Integral Package Library
 *  Implementations of various one- and two-body
 *  integrals via McMurchie-Davidson recurrence scheme.
 */

/**
 *  @defgroup OEPDEV_INTEGRAL_HELPERS The Integral Helper Library
 *  You will find here various iterators to go through shells
 *  while computing ERI, or iterators over ERI itself.
 */

/**
 *  @defgroup OEPDEV_3DFIELDS The Three-Dimensional Scalar Fields Library
 *  Handles all sorts of scalar distributions in 3D Euclidean space,
 *  such as potentials defined at particular collection of points.
 *  In this Module, you will also find handling both random and ordered
 *  points collections in a form of a G09 cube, as well as handling 
 *  G09 Cube files.
 */

/**
 *  @defgroup OEPDEV_ESP The Multipole Fitting Library
 *  Implements methods to fit the generalized multipole moments
 *  of a generalized density distribution based on the input generalized
 *  potential scalar field. Among others, you will find here the 
 *  electrostatic potential (ESP) fitting method.
 */

/**
 *  @defgroup OEPDEV_UTILITIES The OEPDev Utilities
 *  Contains utility functions such as printing OEPDev preambule 
 *  to the output file, class for wavefunction union etc.
 */


