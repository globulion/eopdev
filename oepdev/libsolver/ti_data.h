#ifndef _oepdev_libutil_ti_data_h
#define _oepdev_libutil_ti_data_h
/** @file ti_data.h */

#include<cstdio>
#include<string>
#include<map>

#include "../lib3d/dmtp.h"

namespace oepdev{

using namespace std;

using SharedDMTPConvergence = std::shared_ptr<oepdev::MultipoleConvergence>;

/** \addtogroup OEPDEV_SOLVERS
 * @{
 */

/**\brief Transfer Integral EET Data.
 *
 * Container for storing and managing TI data for EET coupling calculations,
 * according to Fujimoto JCP 2012:
 *  * exciton Hamiltonian matrix elements
 *  * overlap integrals between basis states
 *  * TrCAMM EET coupling convergence object
 * 
 * Contains useful methods to process exciton Hamiltonian matrix elements:
 *  * compute direct and indirect EET coupling constants
 *  * compute overlap-corrected exciton Hamiltonian off-diagonal matrix elements
 *  * include or exclude environmental correction in the diagonal exciton Hamiltonian
 *  * activate TrCAMM approximation of \f$ V^{{\rm Coul},(0)} \f$
 *  * activate Mulliken approximation for \f$ V^{{\rm Exch},(0)} \f$ and \f$ V^{{\rm CT},(0)} \f$
 *
 * To activate/deactivate the various approximations and corrections listed above, set
 * the following attributes
 *  * `diagonal_correction`
 *  * `mulliken_approximation`
 *  * `overlap_correction`
 *  * `trcamm_approximation`
 * 
 * to `true`/`false`, accroding to your need.
 *
 * ## Example of usage.
 * 
 * ```c++
 *   // Set up exciton Hamiltonian
 *   TIData data = TIData();                                        
 *   data.set_s(S12, S13, S14, S32, S42, S34);
 *   data.set_e(E1, E2, E3, E4);
 *   data.set_de(E1 - E01, E2 - E02);
 *   data.v0["COUL"]= V0_Coul;
 *   data.v0["EXCH"]= V0_Exch;
 *   data.v0["ET1"] = V0_ET1;
 *   data.v0["ET2"] = V0_ET2;
 *   data.v0["HT1"] = V0_HT1;
 *   data.v0["HT2"] = V0_HT2;
 *   data.v0["CT" ] = V0_CT;
 *   data.v0["EXCH_M"]= V0_Exch_M;
 *   data.v0["CT_M"] = V0_CT_M;
 *                                                                  
 *   // Set up appriximations and corrections
 *   data.diagonal_correction = true;
 *   data.mulliken_approximation= false;
 *   data.trcamm_approximation = false;
 *   data.overlap_correction = true;
 *                                                                  
 *   // Compute overlap-corrected indirect coupling matrix elements
 *   double V_ET1 = data.overlap_corrected("ET1");
 *   double V_ET2 = data.overlap_corrected("ET2");
 *   double V_HT1 = data.overlap_corrected("HT1");
 *   double V_HT2 = data.overlap_corrected("HT2");
 *   double V_CT  = data.overlap_corrected("CT") ;
 *   double V_CT_M= data.overlap_corrected("CT_M");
 *                                                                  
 *   // Compute final coupling contributions
 *   double V_Coul = data.overlap_corrected("COUL");
 *   double V_Exch = data.overlap_corrected("EXCH");
 *   double V_Ovrl = data.overlap_corrected("OVRL");
 *                                                                  
 *   double V_Exch_M= data.overlap_corrected("EXCH_M");
 *                                                                  
 *   double V_TI_2 = data.coupling_indirect_ti2();
 *   double V_TI_3 = data.coupling_indirect_ti3();
 *                                                                  
 *   data.diagonal_correction = false;
 *   double V0_TI_2 = data.coupling_indirect_ti2();
 *   double V0_TI_3 = data.coupling_indirect_ti3();
 *                                                                  
 *   data.mulliken_approximation = true;
 *   double V0_TI_3_M = data.coupling_indirect_ti3();
 *   data.diagonal_correction = true;
 *   double V_TI_3_M = data.coupling_indirect_ti3();
 *                                                                  
 *   double V_direct = V_Coul + V_Exch + V_Ovrl;
 *   double V_indirect = V_TI_2 + V_TI_3;
 * ```
 *
 * \see oepdev::EETCouplingSolver
 */
class TIData 
{
 public:

   /// Constructor
   TIData();
   /// Destroctor
   virtual ~TIData();

   // <--- Modifiers ---> //

   /// Set the overlap integrals between basis states, \f$ S_{ij} \f$, for *ij*=12,13,14,23,24,34
   void set_s(double, double, double, double, double, double);
   /// Set the diagonal exciton Hamiltonian matrix elements \f$ E_n \f$ for *n*=1,2,3,4
   void set_e(double, double, double, double);
   /// Set environmental corrections \f$ \Delta E_1 \f$ and \f$ \Delta E_2 \f$
   void set_de(double, double);
   /// Set the convergence object for TrCAMM-based \f$ V^{{\rm Coul},(0)}\f$
   void set_trcamm_coupling(oepdev::SharedDMTPConvergence);

 //void set_output_coupling_units_converter(double c);


   // <--- Computers ---> //

   /** Compute Coulombic coupling approximated by TrCAMM.
    * @param rn - convergence of TrCAMM coupling. Can be from `R1` to `R5`, which corresponds
    *             to the \f$ R^{-n} \f$ series expansion of distributed multipoles.
    * @returns \f$ V^{{\rm Coul},(0)} \approx V^{{\rm TrCAMM},(0)} (R^{-n}) \f$
    */
   virtual double coupling_trcamm(const std::string& rn);
   /** Compute the direct EET coupling constant.
    * @returns \f$ V^{\rm Dir} = V^{{\rm Coul}} + V^{{\rm Exch}} + V^{\rm Ovrl} \f$
    *
    * Overlap and diagonal corrections as well as TrCAMM and Mulliken approximations 
    * for Coulomb and pure exchange parts can be set.
    */
   virtual double coupling_direct(void);
   /** Compute the direct EET coupling constant in Forster limit (Coulombic approximation)
    * @returns \f$ V^{\rm Dir} = V^{{\rm Coul}}\f$
    *
    * Overlap correction as well as TrCAMM approximation for Coulomb coupling can be set.
    */
   virtual double coupling_direct_coul(void);
   /** Compute the direct EET coupling constant due to pure exchange.
    * @returns \f$ V^{\rm Dir} = V^{{\rm Exch}} \f$
    *
    * Overlap correction as well Mulliken approximation for pure exchange coupling can be set.
    */
   virtual double coupling_direct_exch(void);
   /** Compute the indirect EET coupling constant.
    * @returns \f$ V^{\rm Indir} = V^{{\rm TI},(2)} + V^{{\rm TI},(3)} \f$
    *
    * Overlap and diagonal corrections as well as Mulliken approximations for \f$ V^{{\rm CT},(0)}\f$ can be set.
    */
   virtual double coupling_indirect(void);
   /** Compute the indirect EET coupling constant in second-order with respect to TI.
    * @returns \f$ V^{{\rm TI},(2)} = -\frac{V^{\rm ET1}V^{\rm HT2}}{E_3-E_1} -\frac{V^{\rm ET2}V^{\rm HT1}}{E_4-E_1} \f$
    *
    * Overlap and diagonal corrections can be set.
    */
   virtual double coupling_indirect_ti2(void);
   /** Compute the indirect EET coupling constant in third-order with respect to TI.
    * @returns \f$ V^{{\rm TI},(3)} = \frac{V^{\rm CT}\left( V^{\rm ET1}V^{\rm ET2} + V^{\rm HT1}V^{\rm HT2}\right)}{(E_3-E_1)(E_4-E_1)} \f$
    *
    * Overlap and diagonal corrections as well as Mulliken approximations for \f$ V^{{\rm CT},(0)}\f$ can be set.
    */
   virtual double coupling_indirect_ti3(void);
   /** Compute the *total* EET coupling constant.
    * @returns \f$ V^{\rm Total} = V^{\rm Dir} + V^{\rm Indir} \f$
    *
    * Overlap and diagonal corrections, TrCAMM approximation for Coulomb coupling, and Mulliken approximations for pure exchange 
    * coupling and \f$ V^{{\rm CT},(0)}\f$ can be set.
    */
   virtual double coupling_total(void);

   /** Compute overlap corrected matrix elements.
    * 
    * @param type - matrix element \f$ V^{{\rm type},(0)}\f$ subject to overlap correction,
    *               where type is one of the following: 
    *                * `COUL` - \f$ V^{{\rm Coul},(0)}\f$,
    *                * `EXCH` - \f$ V^{{\rm Exch},(0)}\f$,
    *                * `TrCAMM_R1` - \f$ V^{{\rm TrCAMM},(0)} (R^{-1})\f$,
    *                * `TrCAMM_R2` - \f$ V^{{\rm TrCAMM},(0)} (R^{-2})\f$,
    *                * `TrCAMM_R3` - \f$ V^{{\rm TrCAMM},(0)} (R^{-3})\f$,
    *                * `TrCAMM_R4` - \f$ V^{{\rm TrCAMM},(0)} (R^{-4})\f$,
    *                * `TrCAMM_R5` - \f$ V^{{\rm TrCAMM},(0)} (R^{-5})\f$,
    *                * `ET1` - \f$ V^{{\rm ET1},(0)}\f$,
    *                * `ET2` - \f$ V^{{\rm ET2},(0)}\f$,
    *                * `HT1` - \f$ V^{{\rm HT1},(0)}\f$,
    *                * `HT2` - \f$ V^{{\rm HT2},(0)}\f$,
    *                * `CT` - \f$ V^{{\rm CT},(0)}\f$,
    *                * `CT_M` - Mulliken-approximated \f$ V^{{\rm CT},(0)}\f$,
    *                * `EXCH_M` - Mulliken-approximated \f$ V^{{\rm Exch},(0)}\f$.
    * 
    *               If type = `OVRL`, the overlap-correction to the direct EET
    *               coupling constant is returned, \f$ V^{\rm Ovrl}\f$.
    * @returns overlap-corrected exciton Hamiltonian matrix element contribution of selected type
    *
    * Diagonal correction can be set.
    */
   virtual double overlap_corrected(const std::string& type);
   /** Compute overlap-corrected direct EET coupling constant.
    *
    * @returns \f$ V^{\rm Dir} = V^{{\rm Coul}} + V^{{\rm Exch}} + V^{\rm Ovrl} \f$
    *
    * Diagonal correction, TrCAMM approximation and Mulliken approximation can be set.
    */
   virtual double overlap_corrected_direct(void);
   /** Compute overlap-corrected direct EET coupling constant from value *v*.
    *
    * @returns \f$ \frac{v}{1 - S_{12}^2}  \f$
    */
   virtual double overlap_corrected_direct(double v);
   /** Compute overlap-corrected coupling constant from value *v* and associated overlap integral *s*.
    *
    * @returns \f$ \frac{1}{1 - s^2} \left( v - \frac{(E_1+E_2)s}{2}\right) \f$
    *
    * Diagonal correction can be set.
    */
   virtual double overlap_corrected_indirect(double v, double s);


 public:

   //@{
   /// Overlap matrix element between basis functions
   double s12, s13, s14, s23, s24, s34;
   //@}
   //@{
   /// Diagonal Hamiltonian matrix element
   double e1, e2, e3, e4;
   //@}
   //@{
   /// Environmental correction to the \f$ E_n \f$ for *n* =1,2
   double de1, de2;
   //@}
   /// Convergence object for Coulombic coupling under TrCAMM approximation
   oepdev::MultipoleConvergence::ConvergenceLevel trcamm_convergence;

   /// Environmental correction activated?
   bool diagonal_correction;
   /// Mulliken approximation activated?
   bool mulliken_approximation;
   /// Overlap correction acrivatved?
   bool overlap_correction;
   /// TrCAMM approximation activated?
   bool trcamm_approximation;


   /** Dictionary of all zeroth-order off-diagonal matrix elements.
    *  
    *  Use only the following keywords:
    *                * `COUL` - \f$ V^{{\rm Coul},(0)}\f$,
    *                * `EXCH` - \f$ V^{{\rm Exch},(0)}\f$,
    *                * `TrCAMM_R1` - \f$ V^{{\rm TrCAMM},(0)} (R^{-1})\f$,
    *                * `TrCAMM_R2` - \f$ V^{{\rm TrCAMM},(0)} (R^{-2})\f$,
    *                * `TrCAMM_R3` - \f$ V^{{\rm TrCAMM},(0)} (R^{-3})\f$,
    *                * `TrCAMM_R4` - \f$ V^{{\rm TrCAMM},(0)} (R^{-4})\f$,
    *                * `TrCAMM_R5` - \f$ V^{{\rm TrCAMM},(0)} (R^{-5})\f$,
    *                * `ET1` - \f$ V^{{\rm ET1},(0)}\f$,
    *                * `ET2` - \f$ V^{{\rm ET2},(0)}\f$,
    *                * `HT1` - \f$ V^{{\rm HT1},(0)}\f$,
    *                * `HT2` - \f$ V^{{\rm HT2},(0)}\f$,
    *                * `CT` - \f$ V^{{\rm CT},(0)}\f$,
    *                * `CT_M` - Mulliken-approximated \f$ V^{{\rm CT},(0)}\f$,
    *                * `EXCH_M` - Mulliken-approximated \f$ V^{{\rm Exch},(0)}\f$,
    *                * `OVRL` - \f$ V^{{\rm Ovrl}}\f$.
    */
   std::map<std::string, double> v0;

   /// V0_Coul multipole convergence
   oepdev::SharedDMTPConvergence v0_trcamm;

 protected:
   /// Conversion factor (unused now)
   double c_;
};




/** @}*/
}      // EndNameSpace oepdev
#endif //_oepdev_libutil_ti_data_h
