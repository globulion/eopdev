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
#define OEPDEV_CRIT_ERI                   1e-9
/// Size of R buffer (OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*OEPDEV_N_MAX_AM*3)
#define OEPDEV_SIZE_BUFFER_R              250563 
/// Size of D2 buffer (3*(OEPDEV_MAX_AM+1)*(OEPDEV_MAX_AM+1)*OEPDEV_N_MAX_AM)
#define OEPDEV_SIZE_BUFFER_D2             3264
/// Energy converters
#define OEPDEV_AU_KcalPerMole             627.509
#define OEPDEV_AU_CMRec                   219474.63
#define OEPDEV_AU_EV                      27.21138

/** @}*/
#endif // oepdev_include_oepdev_files_h

// uncomment to disable assert()
// #define NDEBUG
