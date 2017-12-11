OEP Design
==========

OEP (One-Electron Potential) is associated with certain quantum one-electron operator 
<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{v}^A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{v}^A" title="\hat{v}^A" /></a>
that defines the ability of molecule *A* to interact in a particular way with other molecules. 
Technically, OEP can be understood as a **container object** (associated with the molecule in question)
that stores the information about the above mentioned quantum operator. 
Here, it is assumed that similar OEP
object is also defined for all other molecules in a molecular aggregate. 

In case of interaction between molecules *A* and *B*,
OEP object of molecule *A* interacts directly with wavefunction object
of the molecule *B*. Defining a 
 * Solver class that handles such interaction 
 * Wavefunction class and
 * OEP class
the universal design of OEP-based approaches can be established and developed.

> **Important:**
>  OEP and Wavefunction classes should not be restricted to Hartree-Fock; in generall any correlated 
>  wavefunction and derived OEP`s should be allowed to work with each other
>

OEP Classes
-----------

There are many types of OEP’s, but the underlying principle is the same and independent of the
type of intermolecular interaction. Therefore OEP’s should be implemented by using a multi-level class design.
In turn, this design depends on the way OEP`s enter the mathematical expressions, i.e., on the types
of matrix elements of the one-electron effective operator
<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{v}^A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{v}^A" title="\hat{v}^A" /></a>
.

### Structure of possible OEP-based expressions and their unification

Structure of OEP-based mathematical expressions is listed below:

| Type  | Matrix Element | 
|--------|----|
| Type 1: | <a href="https://www.codecogs.com/eqnedit.phplatex=\left(&space;I&space;\left|&space;\hat{v}^A&space;\right|&space;\right&space;K)&space;\quad&space;\text{where&space;}&space;I,K\in&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex \left(&space;I&space;\left|&space;\hat{v}^A&space;\right|&space;\right&space;K)&space;\quad&space;\text{where&space;}&space;I,K\in&space;A" title="\left( I \left| \hat{v}^A \right| \right K) \quad \text{where } I,K\in A" /></a>
|
| Type 2: | fhdfdhf |


 * Type 1: 
   <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;I&space;\left|&space;\hat{v}^A&space;\right|&space;\right&space;K)&space;\quad&space;\text{where&space;}&space;I,K\in&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;I&space;\left|&space;\hat{v}^A&space;\right|&space;\right&space;K)&space;\quad&space;\text{where&space;}&space;I,K\in&space;A" title="\left( I \left| \hat{v}^A \right| \right K) \quad \text{where } I,K\in A" /></a>
   Bulaaaa
 * Type 2:
   Bulaaa
