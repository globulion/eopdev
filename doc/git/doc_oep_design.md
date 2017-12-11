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
<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{v}^A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{v}^A" title="\hat{v}^A" /></a>.

### Structure of possible OEP-based expressions and their unification

Structure of OEP-based mathematical expressions is listed below:

| Type  | Matrix Element | Comment |
|--------|----|---|
| Type 1:| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;I&space;\vert&space;\hat{v}^A&space;\vert&space;\right&space;J)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;I&space;\vert&space;\hat{v}^A&space;\vert&space;\right&space;J)" title="\left( I \vert \hat{v}^A \vert \right J)" /></a>    | <a href="https://www.codecogs.com/eqnedit.php?latex=I\in&space;A,J\in&space;B" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I\in&space;A,J\in&space;B" title="I\in A,J\in B" /></a>  |
| Type 2:| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;J&space;\vert&space;\hat{v}^A&space;\vert&space;\right&space;L)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;J&space;\vert&space;\hat{v}^A&space;\vert&space;\right&space;L)" title="\left( J \vert \hat{v}^A \vert \right L)" /></a>    | <a href="https://www.codecogs.com/eqnedit.php?latex=J,L&space;\in&space;B" target="_blank"><img src="https://latex.codecogs.com/gif.latex?J,L&space;\in&space;B" title="J,L \in B" /></a>  |

In the above table, *I*, *J* and *K* indices correspond to basis functions or molecular orbitals. Basis functions can be primary or auxiliary ([OEP-specialized density-fitting](https://github.com/globulion/oepdev/blob/master/doc/git/doc_density_fitting.md)). Depending on the type of function and matrix element, there are many subtypes of resulting matrix elements that differ in their dimensionality. Examples are given below:

| Matrix Element | ESP-based | DF-based |
|----|---|---|
| | | |
| | | |
| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;j&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;l&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;j&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;l&space;\right)" title="\left( j \vert \hat{v}^{A[i]} \vert l \right)" /></a> | <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\iota\kappa\in&space;A}&space;{\color{DarkOrange}&space;S_{j\iota}&space;}&space;{\color{blue}&space;v_{\iota\kappa}^{A[i]}}&space;{\color{DarkOrange}&space;S_{l\kappa}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\iota\kappa\in&space;A}&space;{\color{DarkOrange}&space;S_{j\iota}&space;}&space;{\color{blue}&space;v_{\iota\kappa}^{A[i]}}&space;{\color{DarkOrange}&space;S_{l\kappa}}" title="\sum_{\iota\kappa\in A} {\color{DarkOrange} S_{j\iota} } {\color{blue} v_{\iota\kappa}^{A[i]}} {\color{DarkOrange} S_{l\kappa}}" /></a>| |
| | | |
