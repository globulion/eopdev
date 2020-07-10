#ifndef _oepdev_libutil_quambo_h
#define _oepdev_libutil_quambo_h
/** @file quambo.h */

#include <string>

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"


namespace oepdev{

/** \addtogroup OEPDEV_UTILITIES 
 * @{
 */

/** \brief Container to store the QUAMBO data
 *
 */
struct QUAMBOData {
 public:
    /// QUAMBO (non-orthogonal) 
    psi::SharedMatrix quambo_nonorthogonal;
    /// QUAMBO (orthogonal) 
    psi::SharedMatrix quambo_orthogonal;
    /// Virtual Valence Molecular Orbitals (VVO)
    psi::SharedMatrix c_mini_vir;
    /// VVO Energies
    psi::SharedVector e_mini_vir;
    /// All Molecular orbitals (OCC + VVO)
    psi::SharedMatrix c_mini;
    /// Energies of All Molecular Orbitals
    psi::SharedVector e_mini;
};

using SharedQUAMBOData = std::shared_ptr<QUAMBOData>;

/** \brief The Quasiatomic Minimal Basis Set Molecular Orbitals (QUAMBO)
 *
 *  TODO
 *  
 *  # Calculation Algorithm.
 *  TODO
 *
 *
 *  References: 
 * 
 *   [1] W. C. Lu, C. Z. Wang, M.W. Schmidt, L. Bytautas, K. M. Ho, K. Reudenberg, 
 *       J. Chem. Phys. 120, 2629 (2004) [original QUAMBO paper]
 *   [2] P. Xu, M. S. Gordon, J. Chem. Phys. 139, 194104 (2013) [application 
 *       of QUAMBO in EFP2 CT term]
 * 
 */
class QUAMBO
{
  public:
   /// Constructor
   QUAMBO(psi::SharedWavefunction wfn, bool acbs = true);
   /// Destructor
  ~QUAMBO();

   const bool acbs;
   void compute(void);

  protected:

    /// Psi4 options
    psi::Options& options_;
    /// Molecule
    psi::SharedMolecule mol_;
    /// Wavefunction
    psi::SharedWavefunction wfn_;
    /// numbers of minimal basis functions of free atoms
    std::map<std::string, int> nbas_atom_mini_;
    /// numbers of unpaired electrons in free atoms
    std::map<std::string, int> unpe_atom_;
    /// AO Overlap Matrix
    psi::SharedMatrix Sao_;
    /// QUAMBO (Alpha, non-orthogonal) 
    psi::SharedMatrix quambo_a_nonorthogonal_;
    /// QUAMBO (Alpha, orthogonal) 
    psi::SharedMatrix quambo_a_orthogonal_;
    /// QUAMBO (Beta, non-orthogonal) 
    psi::SharedMatrix quambo_b_nonorthogonal_;
    /// QUAMBO (Beta, orthogonal) 
    psi::SharedMatrix quambo_b_orthogonal_;
    /// Virtual Valence Molecular Orbitals (Alpha, VVO)
    psi::SharedMatrix c_a_mini_vir_;
    /// Virtual Valence Molecular Orbitals (Beta, VVO)
    psi::SharedMatrix c_b_mini_vir_;
    /// VVO Energies (Alpha)
    psi::SharedVector e_a_mini_vir_;
    /// VVO Energies (Beta)
    psi::SharedVector e_b_mini_vir_;
    /// All Molecular orbitals (Alpha, OCC + VVO)
    psi::SharedMatrix c_a_mini_;
    /// All Molecular orbitals (Beta, OCC + VVO)
    psi::SharedMatrix c_b_mini_;
    /// Energies of All Molecular Orbitals (Alpha)
    psi::SharedVector e_a_mini_;
    /// Energies of All Molecular Orbitals (Beta)
    psi::SharedVector e_b_mini_;

    double compute_error_between_two_vectors_(psi::SharedVector a, psi::SharedVector b);
    int calculate_nbas_mini_(void);
    psi::SharedMolecule atomize_(int i);
    SharedQUAMBOData compute_quambo_data_(psi::SharedMatrix, psi::SharedMatrix, psi::SharedVector, psi::SharedVector,
                                          psi::SharedMatrix, psi::SharedMatrix, psi::SharedMatrix, int, int, std::string);

};




/** @}*/


} // EndNameSpace oepdev


#endif // _oepdev_libutil_quambo_h
