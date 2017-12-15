Density-fitting specialized for OEP’s
=====================================

To get the ab-initio representation of a OEP, one can use a procedure similar to
the typical density fitting or resolution of identity, both of which are nowadays widely used 
to compute electron-repulsion integrals (ERI’s) more efficiently. 

An arbitrary one-electron potential of molecule *A* acting on any state vector associated with molecule *A* can be expanded in the *auxiliary space* centered on *A* as

<a href="https://www.codecogs.com/eqnedit.php?latex=v\vert&space;i)&space;=&space;\sum_{\xi}&space;v\vert&space;\xi)&space;(&space;\xi&space;\vert&space;i)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v\vert&space;i)&space;=&space;\sum_{\xi}&space;v\vert&space;\xi)&space;(&space;\xi&space;\vert&space;i)" title="v\vert i) = \sum_{\xi} v\vert \xi) ( \xi \vert i)" /></a>

under the necessary assumption that the auxiliary basis set is *complete*. In that case, formally one can write the following identity

<a href="https://www.codecogs.com/eqnedit.php?latex=(\eta\vert&space;v\vert&space;i)&space;=&space;\sum_{\xi}&space;(\eta\vert&space;v\vert&space;\xi)&space;S_{\xi&space;i}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(\eta\vert&space;v\vert&space;i)&space;=&space;\sum_{\xi}&space;(\eta\vert&space;v\vert&space;\xi)&space;S_{\xi&space;i}" title="(\eta\vert v\vert i) = \sum_{\xi} (\eta\vert v\vert \xi) S_{\xi i}" /></a>

The matrix elements of the OEP operator in auxiliary space can be computed in the same way as the matrix elements with any other basis function. In reality, it is almost impossible to reach the completness of the basis set, however, but it is possible to obtain the **effective** matrix elements of the OEP operator in auxiliary space, rather than compute them as they are in the above equation explicitly. We expand the LHS of the first equation on this page in a series of the auxiliary basis functions scaled by the undetermined expansion coefficients: 

<a href="https://www.codecogs.com/eqnedit.php?latex=v\vert&space;i)&space;=&space;\sum_{\xi}&space;{\color{Blue}&space;G_{i\xi}}&space;\vert&space;\xi)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v\vert&space;i)&space;=&space;\sum_{\xi}&space;{\color{Blue}&space;G_{i\xi}}&space;\vert&space;\xi)" title="v\vert i) = \sum_{\xi} {\color{Blue} G_{i\xi}} \vert \xi)" /></a>

The expansion coefficients (shown in blue) are the effective matrix elements of the OEP operator in auxiliary basis set. Now, multiplying both sides by another auxiliary basis function and subsequently inverting the equation one obtains the expansion coefficients:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boxed{&space;G_{i\eta}&space;=&space;\sum_\xi&space;(i&space;\vert&space;v&space;\vert&space;\eta)&space;\left[&space;{\bf&space;S}^{-1}&space;\right]_{\eta\xi}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boxed{&space;G_{i\eta}&space;=&space;\sum_\xi&space;(i&space;\vert&space;v&space;\vert&space;\eta)&space;\left[&space;{\bf&space;S}^{-1}&space;\right]_{\eta\xi}&space;}" title="\boxed{ G_{i\eta} = \sum_\xi (i \vert v \vert \eta) \left[ {\bf S}^{-1} \right]_{\eta\xi} }" /></a>

In this way, it is possible to approximately determine the matrix elements of the OEP operator with any other basis function in case the auxiliary basis set is not complete. **In particular**, when the other basis function does not belong to molecule *A* but to the (changing) environment, it is straightforward to compute the resulting matrix element, which is simply given as  

<a href="https://www.codecogs.com/eqnedit.php?latex=(j_{\in&space;B}&space;\vert&space;v^A&space;\vert&space;i_{\in&space;A})&space;=&space;\sum_\xi&space;{\color{DarkOrange}&space;S_{j\xi}}&space;{\color{Blue}&space;G_{i\xi}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(j_{\in&space;B}&space;\vert&space;v^A&space;\vert&space;i_{\in&space;A})&space;=&space;\sum_\xi&space;{\color{DarkOrange}&space;S_{j\xi}}&space;{\color{Blue}&space;G_{i\xi}}" title="(j_{\in B} \vert v^A \vert i_{\in A}) = \sum_\xi {\color{DarkOrange} S_{j\xi}} {\color{Blue} G_{i\xi}}" /></a>

where *j* denotes the other (environmental) basis function.

In the above equation, the OEP-part (fragment parameters for molecule *A* only) are shown in blue whereas the Solver-part (subject to be computed by solver on the fly) is shown in brown. This then forms a basis for fragment-based approach to solve Quantum Chemistry problems related to the extended molecular aggregates.
