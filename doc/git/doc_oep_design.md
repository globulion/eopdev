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

There are many [types of OEP’s](https://github.com/globulion/oepdev/tree/master/doc/git/doc_oep_types.md), but the underlying principle is the same and independent of the
type of intermolecular interaction. Therefore, the OEP’s should be implemented by using a multi-level class design.
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

| Matrix Element | DF-based form | ESP-based form |
|----|---|---|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;\mu&space;\vert&space;\hat{v}^{A[\mu]}&space;\vert&space;\sigma&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;\mu&space;\vert&space;\hat{v}^{A[\mu]}&space;\vert&space;\sigma&space;\right)" title="\left( \mu \vert \hat{v}^{A[\mu]} \vert \sigma \right)" /></a> | <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\iota\in&space;A}&space;{\color{Blue}&space;v_{\mu\iota}^A}&space;{\color{DarkOrange}&space;S_{\iota\sigma}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\iota\in&space;A}&space;{\color{Blue}&space;v_{\mu\iota}^A}&space;{\color{DarkOrange}&space;S_{\iota\sigma}}" title="\sum_{\iota\in A} {\color{Blue} v_{\mu\iota}^A} {\color{DarkOrange} S_{\iota\sigma}}" /></a>| <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[\mu]}}&space;{\color{DarkOrange}&space;V_{\mu\sigma}^{(\alpha)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[\mu]}}&space;{\color{DarkOrange}&space;V_{\mu\sigma}^{(\alpha)}}" title="\sum_{\alpha\in A} {\color{Blue} q_\alpha^{A[\mu]}} {\color{DarkOrange} V_{\mu\sigma}^{(\alpha)}}" /></a>|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;i&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;j&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;i&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;j&space;\right)" title="\left( i \vert \hat{v}^{A[i]} \vert j \right)" /></a> | <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\iota\in&space;A}&space;{\color{Blue}&space;v_{i\iota}^A}&space;{\color{DarkOrange}&space;S_{\iota\sigma}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\iota\in&space;A}&space;{\color{Blue}&space;v_{i\iota}^A}&space;{\color{DarkOrange}&space;S_{\iota\sigma}}" title="\sum_{\iota\in A} {\color{Blue} v_{i\iota}^A} {\color{DarkOrange} S_{\iota\sigma}}" /></a>| <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[i]}}&space;{\color{DarkOrange}&space;V_{ij}^{(\alpha)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[i]}}&space;{\color{DarkOrange}&space;V_{ij}^{(\alpha)}}" title="\sum_{\alpha\in A} {\color{Blue} q_\alpha^{A[i]}} {\color{DarkOrange} V_{ij}^{(\alpha)}}" /></a>|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\left(&space;j&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;l&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(&space;j&space;\vert&space;\hat{v}^{A[i]}&space;\vert&space;l&space;\right)" title="\left( j \vert \hat{v}^{A[i]} \vert l \right)" /></a> | <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\iota\kappa\in&space;A}&space;{\color{DarkOrange}&space;S_{j\iota}&space;}&space;{\color{blue}&space;v_{\iota\kappa}^{A[i]}}&space;{\color{DarkOrange}&space;S_{\kappa&space;l}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\iota\kappa\in&space;A}&space;{\color{DarkOrange}&space;S_{j\iota}&space;}&space;{\color{blue}&space;v_{\iota\kappa}^{A[i]}}&space;{\color{DarkOrange}&space;S_{\kappa&space;l}}" title="\sum_{\iota\kappa\in A} {\color{DarkOrange} S_{j\iota} } {\color{blue} v_{\iota\kappa}^{A[i]}} {\color{DarkOrange} S_{\kappa l}}" /></a>| <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[i]}}&space;{\color{DarkOrange}&space;V_{jl}^{(\alpha)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{\alpha\in&space;A}&space;{\color{Blue}&space;q_\alpha^{A[i]}}&space;{\color{DarkOrange}&space;V_{jl}^{(\alpha)}}" title="\sum_{\alpha\in A} {\color{Blue} q_\alpha^{A[i]}} {\color{DarkOrange} V_{jl}^{(\alpha)}}" /></a> |

In the formulae above, the OEP-part (stored by OEP instances) is shown in blue whereas the Solver-part (to be computed by the Solver) is shown in brown. It is apparent that all OEP-parts have the form of 2nd- or 3rd-rank tensors with different class of axes (molecular orbitals, primary/auxiliary basis, atomic space). Therefore, they can be uniquely defined by a unified *tensor object* (storing double precision numbers) and unified *dimension object* storing the information of the axes classes.

In Psi4, a perfect candidate for the above is `psi4::Tensor` class declared in `psi4/libthce/thce.h`. Except from the numeric content its instances also store the information of the dimensions in a form of a vector of `psi4::Dimension` instances.

Another possibility is to use `psi::Matrix` objects, instead of `psi4::Tensor` objects, possibly putting them into a `std::vector` container in case there is more than two axes.




