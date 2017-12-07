oep-dev
=======

Generalized One-Electron Potentials: Development Platform

Description
-----------

Created with intention to test various models of the interaction energy 
between two molecules, described by the Hartree-Fock-Roothaan-Hall theory 
or the configuration interaction with singles theory. 

In particular, the plugin tests the models of:

*******
  1. the Pauli exchange-repulsion interaction energy    (Project II ) 
  2. the Induction interaction energy                   (Project III)
  3. the excitation energy transfer couplings           (Project I  )
*******

against benchmarks (exact or reference solutions). Detailed list of models 
is given below:

**Table 1.** Models subject to be implemented and analyzed within oep-dev.

 | Pauli energy             | Induction energy         | EET Coupling         |
 |--------------------------|--------------------------|----------------------|
 | EFP2-Pauli               | EFP2-Induced Dipoles     | TrCAMM               |
 | Murrel et al.'s theory   | Density Susceptibility   | OEP-ET/HT            |
 | OEP-Murrel et al.'s      |                          | TDFI-TI              |
 |                          |                          | FED                  |
 | Exact (Stone's)          | Exact (incl. CT)         | Exact (ESD)          |

