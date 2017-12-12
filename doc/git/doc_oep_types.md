List of One-Electron Potentals
==============================

Here I provide the list of OEP’s that have been already derived within the scope of the OEPDev project. 

Pauli Repulsion OEP’s
---------------------

The following potentials are derived for the evaluation of the Pauli repulsion energy based on Murrel’s expressions.

### First-order contribution in overlap matrix expansion.

3D forms:

<a href="https://www.codecogs.com/eqnedit.php?latex=v({\bf&space;r})^{A[i]}_{S^{-1}}&space;=&space;-\sum_{x\in&space;A}&space;\frac{Z_x}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r}_x\vert}&space;&plus;&space;\sum_{\mu\nu\in&space;A}&space;\left\{&space;D_{\nu\mu}&space;-&space;C^{*}_{\mu&space;i}&space;C_{\nu&space;i}&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\mu({\bf&space;r'})&space;\phi_\nu({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v({\bf&space;r})^{A[i]}_{S^{-1}}&space;=&space;-\sum_{x\in&space;A}&space;\frac{Z_x}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r}_x\vert}&space;&plus;&space;\sum_{\mu\nu\in&space;A}&space;\left\{&space;D_{\nu\mu}&space;-&space;C^{*}_{\mu&space;i}&space;C_{\nu&space;i}&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\mu({\bf&space;r'})&space;\phi_\nu({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" title="v({\bf r})^{A[i]}_{S^{-1}} = -\sum_{x\in A} \frac{Z_x}{\vert {\bf r} - {\bf r}_x\vert} + \sum_{\mu\nu\in A} \left\{ D_{\nu\mu} - C^{*}_{\mu i} C_{\nu i} \right \} \int d{\bf r'} \frac{\phi^{*}_\mu({\bf r'}) \phi_\nu({\bf r'})}{\vert {\bf r} - {\bf r'}\vert}" /></a>

Matrix forms:

<a href="https://www.codecogs.com/eqnedit.php?latex=v_{\xi&space;i}({S^{-1}})&space;=&space;\sum_{\kappa\in&space;A}&space;C_{i\kappa}&space;\left\{&space;-\sum_{x\in&space;A}&space;V^{(x)}_{\kappa\xi}&space;&plus;&space;\sum_{\mu\nu\in&space;A}&space;\left\{&space;D_{\nu\mu}&space;-&space;C^{*}_{\mu&space;i}&space;C_{\nu&space;i}&space;\right&space;\}&space;\left(&space;\mu\nu&space;\vert&space;\xi\kappa&space;\right&space;)\right&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v_{\xi&space;i}({S^{-1}})&space;=&space;\sum_{\kappa\in&space;A}&space;C_{i\kappa}&space;\left\{&space;-\sum_{x\in&space;A}&space;V^{(x)}_{\kappa\xi}&space;&plus;&space;\sum_{\mu\nu\in&space;A}&space;\left\{&space;D_{\nu\mu}&space;-&space;C^{*}_{\mu&space;i}&space;C_{\nu&space;i}&space;\right&space;\}&space;\left(&space;\mu\nu&space;\vert&space;\xi\kappa&space;\right&space;)\right&space;\}" title="v_{\xi i}({S^{-1}}) = \sum_{\kappa\in A} C_{i\kappa} \left\{ -\sum_{x\in A} V^{(x)}_{\kappa\xi} + \sum_{\mu\nu\in A} \left\{ D_{\nu\mu} - C^{*}_{\mu i} C_{\nu i} \right \} \left( \mu\nu \vert \xi\kappa \right )\right \}" /></a>

### Second-order contribution in overlap matrix expansion.

Excitonic Energy Transfer OEP’s
-------------------------------

The following potentials are derived for the evaluation of the short-range EET couplings based on Fujimoto’s TDFI-TI method.

### ET contributions.

3D forms:

<a href="https://www.codecogs.com/eqnedit.php?latex=v({\bf&space;r})^{A[\mu]}_{1}&space;=&space;-C^*_{\mu&space;L}&space;\sum_{x\in&space;A}&space;\frac{Z_x}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r}_x\vert}&space;&plus;&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;C^*_{\mu&space;L}&space;D_{\nu\kappa}&space;-&space;\frac{1}{2}&space;C^{*}_{\nu&space;L}&space;D_{\mu\kappa}&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\nu({\bf&space;r'})&space;\phi_\kappa({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v({\bf&space;r})^{A[\mu]}_{1}&space;=&space;-C^*_{\mu&space;L}&space;\sum_{x\in&space;A}&space;\frac{Z_x}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r}_x\vert}&space;&plus;&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;C^*_{\mu&space;L}&space;D_{\nu\kappa}&space;-&space;\frac{1}{2}&space;C^{*}_{\nu&space;L}&space;D_{\mu\kappa}&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\nu({\bf&space;r'})&space;\phi_\kappa({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" title="v({\bf r})^{A[\mu]}_{1} = -C^*_{\mu L} \sum_{x\in A} \frac{Z_x}{\vert {\bf r} - {\bf r}_x\vert} + \sum_{\nu\kappa\in A} \left\{ C^*_{\mu L} D_{\nu\kappa} - \frac{1}{2} C^{*}_{\nu L} D_{\mu\kappa} \right \} \int d{\bf r'} \frac{\phi^{*}_\nu({\bf r'}) \phi_\kappa({\bf r'})}{\vert {\bf r} - {\bf r'}\vert}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=v({\bf&space;r})^{A[\mu]}_{2}&space;=&space;C_{\kappa&space;H}&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;2&space;C^*_{\nu&space;L}&space;C_{\mu&space;H}^*&space;-&space;C^{*}_{\nu&space;H}&space;C_{\mu&space;L}^*&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\nu({\bf&space;r'})&space;\phi_\kappa({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v({\bf&space;r})^{A[\mu]}_{2}&space;=&space;C_{\kappa&space;H}&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;2&space;C^*_{\nu&space;L}&space;C_{\mu&space;H}^*&space;-&space;C^{*}_{\nu&space;H}&space;C_{\mu&space;L}^*&space;\right&space;\}&space;\int&space;d{\bf&space;r'}&space;\frac{\phi^{*}_\nu({\bf&space;r'})&space;\phi_\kappa({\bf&space;r'})}{\vert&space;{\bf&space;r}&space;-&space;{\bf&space;r'}\vert}" title="v({\bf r})^{A[\mu]}_{2} = C_{\kappa H} \sum_{\nu\kappa\in A} \left\{ 2 C^*_{\nu L} C_{\mu H}^* - C^{*}_{\nu H} C_{\mu L}^* \right \} \int d{\bf r'} \frac{\phi^{*}_\nu({\bf r'}) \phi_\kappa({\bf r'})}{\vert {\bf r} - {\bf r'}\vert}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=v({\bf&space;r})^{A[\mu]}_{3}&space;=&space;v({\bf&space;r})^{A[\mu]}_{1}&space;&plus;&space;v({\bf&space;r})^{A[\mu]}_{1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v({\bf&space;r})^{A[\mu]}_{3}&space;=&space;v({\bf&space;r})^{A[\mu]}_{1}&space;&plus;&space;v({\bf&space;r})^{A[\mu]}_{1}" title="v({\bf r})^{A[\mu]}_{3} = v({\bf r})^{A[\mu]}_{1} + v({\bf r})^{A[\mu]}_{1}" /></a>

Matrix forms:

<a href="https://www.codecogs.com/eqnedit.php?latex=v_{\mu\xi}(1)&space;=&space;-C^*_{\mu&space;L}&space;\sum_{x\in&space;A}&space;V^{x}_{\mu\xi}&space;&plus;&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;C^*_{\mu&space;L}&space;D_{\nu\kappa}&space;-&space;\frac{1}{2}&space;C^{*}_{\nu&space;L}&space;D_{\mu\kappa}&space;\right&space;\}&space;\left(&space;\nu\kappa&space;\vert&space;\mu\xi&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v_{\mu\xi}(1)&space;=&space;-C^*_{\mu&space;L}&space;\sum_{x\in&space;A}&space;V^{x}_{\mu\xi}&space;&plus;&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;C^*_{\mu&space;L}&space;D_{\nu\kappa}&space;-&space;\frac{1}{2}&space;C^{*}_{\nu&space;L}&space;D_{\mu\kappa}&space;\right&space;\}&space;\left(&space;\nu\kappa&space;\vert&space;\mu\xi&space;\right&space;)" title="v_{\mu\xi}(1) = -C^*_{\mu L} \sum_{x\in A} V^{x}_{\mu\xi} + \sum_{\nu\kappa\in A} \left\{ C^*_{\mu L} D_{\nu\kappa} - \frac{1}{2} C^{*}_{\nu L} D_{\mu\kappa} \right \} \left( \nu\kappa \vert \mu\xi \right )" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=v_{\mu\xi}(2)=&space;C_{\kappa&space;H}&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;2&space;C^*_{\nu&space;L}&space;C_{\mu&space;H}^*&space;-&space;C^{*}_{\nu&space;H}&space;C_{\mu&space;L}^*&space;\right&space;\}&space;\left(&space;\nu\kappa&space;\vert&space;\mu\xi&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v_{\mu\xi}(2)=&space;C_{\kappa&space;H}&space;\sum_{\nu\kappa\in&space;A}&space;\left\{&space;2&space;C^*_{\nu&space;L}&space;C_{\mu&space;H}^*&space;-&space;C^{*}_{\nu&space;H}&space;C_{\mu&space;L}^*&space;\right&space;\}&space;\left(&space;\nu\kappa&space;\vert&space;\mu\xi&space;\right&space;)" title="v_{\mu\xi}(2)= C_{\kappa H} \sum_{\nu\kappa\in A} \left\{ 2 C^*_{\nu L} C_{\mu H}^* - C^{*}_{\nu H} C_{\mu L}^* \right \} \left( \nu\kappa \vert \mu\xi \right )" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=v_{\mu\xi}(3)=&space;v_{\mu\xi}(1)&space;&plus;&space;v_{\mu\xi}(2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v_{\mu\xi}(3)=&space;v_{\mu\xi}(1)&space;&plus;&space;v_{\mu\xi}(2)" title="v_{\mu\xi}(3)= v_{\mu\xi}(1) + v_{\mu\xi}(2)" /></a>

### HT contributions.

Do be derived.

### CT contributions.

To be derived.

Full HF Interaction OEP’s
-------------------------

The following potentials are derived for the evaluation of the full Hartree-Fock interaction energy based on the OEPDev equations.






****************

References
----------
