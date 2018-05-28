/** 
 *  @defgroup OEPDEV_OEPS The Generalized One-Electron Potentials Library
 *  @brief
 *  Implements the goal of this project: The Generalized One-Electron Potentials
 *  (OEP's). You will find here OEP's for computation of Pauli repulsion energy,
 *  charge-transfer energy and others.
 */

/**
 *  @defgroup OEPDEV_SOLVERS The OEPDev Solver Library
 *  @brief 
 *  Implementations of various solvers for molecular properties
 *  as a functions of unperturbed monomeric wavefunctions.
 *  This is the place all target OEP-based models are implemented
 *  and compared with benchmark models.
 */

/**
 *  @defgroup OEPDEV_GEFP The Generalized Effective Fragment Potentials Library
 *  @brief
 *  Implements the GEFP method, the far goal of the OEPDev project.
 *  Here you will find the containers for GEFP parameters, the density matrix
 *  susceptibility tensors and GEFP solvers.
 */

/**
 *  @defgroup OEPDEV_LIBINTS The Integral Package Library
 *  @brief
 *  Implementations of various one- and two-body
 *  integrals via McMurchie-Davidson recurrence scheme.
 *
 *  Here, we define the primitive Gaussian type functions (GTO's)
 *  \f{align*}{
 *    \phi_i({\bf r}) &\equiv x_A^{n_1} y_A^{l_1} z_A^{m_1} e^{-\alpha_1r_A^2} \\
 *    \phi_j({\bf r}) &\equiv x_B^{n_2} y_B^{l_2} z_B^{m_2} e^{-\alpha_2r_B^2} \\
 *    \phi_k({\bf r}) &\equiv x_C^{n_3} y_C^{l_3} z_C^{m_3} e^{-\alpha_3r_C^2} 
 *  \f}
 *  in which \f$ {\bf r}_A \equiv {\bf r} - {\bf A}\f$ and so on. \f$ {\bf A} \f$
 *  is the centre of the GTO, \f$\alpha_1\f$ its exponent, whereas \f$n_1,l_1,m_1\f$
 *  the Cartesian angular momenta, with the total angular momentum \f$\theta_1 = n_1+l_1+m_1\f$.
 *
 *  In OEPDev implementations, the following definition shall be in use:
 *  \f{align*}{
 *     {\bf P} &\equiv \frac{\alpha_1{\bf A} + \alpha_2{\bf B}}{\alpha_1+\alpha_2} \\
 *     {\bf Q} &\equiv \frac{\alpha_3{\bf C} + \alpha_4{\bf D}}{\alpha_3+\alpha_4} \\
 *     {\bf R} &\equiv \frac{\alpha_1{\bf A} + \alpha_2{\bf B} + \alpha_3{\bf C}}{\alpha_1+\alpha_2+\alpha_3} \\
 *     \alpha_P &\equiv \alpha_1+\alpha_2 \\
 *     \alpha_Q &\equiv \alpha_3+\alpha_4 \\
 *     \alpha_R &\equiv \alpha_1+\alpha_2+\alpha_3
 *  \f}
 *  The unnormalized products of primitive GTO's are denoted here as
 *  \f{align*}{
 *   [ij]  &\equiv \phi_i({\bf r}) \phi_j({\bf r}) \\
 *   [ijk] &\equiv \phi_i({\bf r}) \phi_j({\bf r}) \phi_k({\bf r})
 *  \f}
 *  \section libints_s1 Hermite Operators
 *  It is convenient to define 
 *  \f[
 *   \Lambda_j(x_P;\alpha_P) \equiv \left(\frac{\partial}{\partial P_x}\right)^j
 *      = \alpha_P^{j/2} H_j(\sqrt{\alpha_P}x_P)
 *  \f]
 *  where \f$ H_j(x)\f$ is the Hermite polynomial of order *j* evaluated at *x* .
 *  Introduction of the above Hermite operator can be used by invoking the
 *  recurrence relationship due to Hermite polynomial properties:
 *  \f[
 *    x_A \Lambda_j(x_P;\alpha_P) = j\Lambda_{j-1} + \vert {\bf P} - {\bf A}\vert_x \Lambda_j
 *                + \frac{1}{2\alpha_P} \Lambda_{j+1}
 *  \f]
 *  This can be directly used to derive very useful McMurchie-Davidson-Hermite coefficients
 *  as expansion coefficients of the polynomial expansions.
 *  \subsection libints_s2 Polynomial Expansions as Hermite Series
 *  By using the previous relation, it is possible to express the following expansions
 *  in Hermite series:
 *  \f{align*}{
 *    x_A^{n_1} &= \sum_{N=0}^{n_1} d_N^{n_1} \Lambda_N(x_A;\alpha_A) \\
 *    x_A^{n_1} x_B^{n_2} &= \sum_{N=0}^{n_1+n_2} d_N^{n_1n_2} \Lambda_N(x_P;\alpha_P) \\
 *    x_A^{n_1} x_B^{n_2} x_C^{n_3} &= \sum_{N=0}^{n_1+n_2+n_3} d_N^{n_1n_2n_3} \Lambda_N(x_R;\alpha_R) 
 *  \f}
 *  The recurrence relationships can be easily found and they read
 *  \f[
 *   d_N^{n_1+1} = \frac{1}{2\alpha_A}d_{N-1}^{n_1} + (N+1) d_{N+1}^{n_1}
 *  \f]
 *  as well as
 *  \f{align*}{
 *   d_N^{n_1+1,n_2} &= \frac{1}{2\alpha_P} d_{N-1}^{n_1n_2} 
 *                              + \vert {\bf P} - {\bf A}\vert_x d_N^{n_1n_2} + (N+1) d_{N+1}^{n_1n_2} \\
 *   d_N^{n_1,n_2+1} &= \frac{1}{2\alpha_P} d_{N-1}^{n_1n_2} 
 *                              + \vert {\bf P} - {\bf B}\vert_x d_N^{n_1n_2} + (N+1) d_{N+1}^{n_1n_2} 
 *  \f}
 *  and
 *  \f{align*}{
 *   d_N^{n_1+1,n_2,n_3} &= \frac{1}{2\alpha_R} d_{N-1}^{n_1n_2n_3} 
 *                              + \vert {\bf R} - {\bf A}\vert_x d_N^{n_1n_2n_3} + (N+1) d_{N+1}^{n_1n_2n_3} \\
 *   d_N^{n_1,n_2+1,n_3} &= \frac{1}{2\alpha_R} d_{N-1}^{n_1n_2n_3} 
 *                              + \vert {\bf R} - {\bf B}\vert_x d_N^{n_1n_2n_3} + (N+1) d_{N+1}^{n_1n_2n_3} \\
 *   d_N^{n_1,n_2,n_3+1} &= \frac{1}{2\alpha_R} d_{N-1}^{n_1n_2n_3} 
 *                              + \vert {\bf R} - {\bf C}\vert_x d_N^{n_1n_2n_3} + (N+1) d_{N+1}^{n_1n_2n_3} 
 *  \f}
 *  respectively.
 *  The first elements are given by
 *  \f{align*}{
 *   d_0^{0}   &= 1\\
 *   d_0^{00}  &= 1\\
 *   d_0^{000} &= 1
 *  \f}
 *  By using the above formalisms, it is strightforward to express the doublet of primitive GTO's
 *  as 
 *  \f[
 *   [ij] = E_{ij} \sum_{N=0}^{n_1+n_2} \sum_{L=0}^{l_1+l_2} \sum_{M=0}^{m_1+m_2}
 *          d_N^{n_1n_2} d_L^{l_1l_2} d_M^{m_1m_2}
 *          \Lambda_N(x_P)\Lambda_L(y_P)\Lambda_M(z_P)e^{-\alpha_Pr_P^2}
 *  \f]
 *  Analogously, the triplet of primitive GTO's is given by
 *  \f[
 *   [ijk] = E_{ijk} \sum_{N=0}^{n_1+n_2+n_3} \sum_{L=0}^{l_1+l_2+l+3} \sum_{M=0}^{m_1+m_2+m_3}
 *          d_N^{n_1n_2n_3} d_L^{l_1l_2l_3} d_M^{m_1m_2m_3}
 *          \Lambda_N(x_R)\Lambda_L(y_R)\Lambda_M(z_R)e^{-\alpha_Rr_R^2}
 *  \f]
 *  The multiplicative constants are given by
 *  \f{align*}{
 *  E_{ij}(\alpha_1,\alpha_2) &= \exp{\left[-\frac{\alpha_1\alpha_2}
 *                                       {\alpha_1+\alpha_2}\vert {\bf A}-{\bf B}\vert^2\right]} \\
 *  E_{ijk}(\alpha_1,\alpha_2,\alpha_3)  &= \exp{\left[-\frac{\alpha_1\alpha_2}
 *                                        {\alpha_1+\alpha_2}\vert {\bf A}-{\bf B}\vert^2\right]} 
 *                                         \exp{\left[-\frac{(\alpha_1+\alpha_2)\alpha_3}
 *                                        {\alpha_1+\alpha_2+\alpha_3}
 *                                         \vert {\bf P}-{\bf C}\vert^2\right]} 
 *  \f}
 *  \section libints_s3 One-Body Integrals over Hermite Functions
 *  The fundamental Hermite integrals that appear during computations
 *  of any kind of one-body integrals over GTO's are as follows
 *  \f[
 *    \left[ NLM \vert \Theta(1) \right] \equiv
 *    \int d{\bf r}_1 \Theta({\bf r}_1) 
 *    \Lambda_{N}(x_{1P};\alpha_P)
 *    \Lambda_{L}(y_{1P};\alpha_P)
 *    \Lambda_{M}(z_{1P};\alpha_P)
 *    e^{-\alpha_Pr^2_{1P}}
 *  \f]
 *  It immediately follows that the overlap, dipole, quadrupole and potential integrals
 *  are given as
 *  \f{align*}{
 *    \left[ NLM \vert 1 \right]   &= \delta_{N0}\delta_{L0}\delta_{M0} \left( \frac{\pi}{\alpha_P} \right)^{3/2} \\
 *    \left[ NLM \vert x_C \right] &= \left[ \delta_{N1} + \vert {\bf PC}\vert_x\delta_{N0} \right]\delta_{L0}\delta_{M0}
 *                                    \left( \frac{\pi}{\alpha_P} \right)^{3/2} \\
 *    \left[ NLM \vert x_C^2 \right] &= \left[ 
 *                                       2\delta_{N2} + 2\vert {\bf PC}\vert_x\delta_{N1} + 
 *                                           \left(
 *                                               \vert {\bf PC}\vert_x^2 + \frac{1}{2\alpha_P}
 *                                           \right) \delta_{N0} 
 *                                      \right] \delta_{L0}\delta_{M0} \left( \frac{\pi}{\alpha_P} \right)^{3/2} \\
 *    \left[ NLM \vert x_Cy_C \right] &= \left( 
 *                                         \delta_{N1} + \vert {\bf PC}\vert_x\delta_{N0}
 *                                        \right) \left(
 *                                         \delta_{L1} + \vert {\bf PC}\vert_y\delta_{L0}
 *                                        \right) \delta_{M0} \left( \frac{\pi}{\alpha_P} \right)^{3/2} \\
 *   \left[ NLM \vert r^{-1}_C \right] &= \frac{2\pi}{\alpha_P} R_{NLM}
 *  \f}
 *  The coefficients \f$ R_{NLM} \f$ are discussed in separate section below.
 *
 *  \section libints_s4 Two-Body Integrals over Hermite Functions
 *  The fundamental Hermite integrals that appear during computations
 *  of any kind of two-electron integrals over GTO's are as follows
 *  \f[
 *    \left[ N_1L_2M_2 \vert N_2L_2M_2\right] \equiv
 *    \iint d{\bf r}_1 d{\bf r}_2 
 *    \Lambda_{N_1}(x_{1P};\alpha_P)
 *    \Lambda_{L_1}(y_{1P};\alpha_P)
 *    \Lambda_{M_1}(z_{1P};\alpha_P)
 *    \Lambda_{N_2}(x_{2Q};\alpha_Q)
 *    \Lambda_{L_2}(y_{2Q};\alpha_Q)
 *    \Lambda_{M_2}(z_{2Q};\alpha_Q)
 *    e^{-\alpha_Pr^2_{1P}-\alpha_Qr^2_{2Q}}
 *  \f]
 *  The above formula dramatically reduces to the following
 *  \f[
 *    \left[ N_1L_2M_2 \vert N_2L_2M_2\right] = \lambda(-)^{N2+L2+M2} R_{N1+N2,L1+L2,M1+M2}
 *  \f]
 *  with
 *  \f[
 *   \lambda \equiv \frac{2\pi^{5/2}}{\alpha_P\alpha_Q\sqrt{\alpha_P+\alpha_Q}}
 *  \f]
 *  To compute the \f$R_{N1+N2,L1+L2,M1+M2}\f$ coefficients, the 
 *  parameter \f$T\f$ is given by
 *  \f[
 *    T = \frac{\alpha_P\alpha_Q}{\alpha_P+\alpha_Q} \vert {\bf P}-{\bf Q} \vert^2
 *  \f]
 *  \section libints_s5 The R(N,L,M) Coefficients
 *  The *R* coefficients are defined as
 *  \f[
 *    R_{NLM} \equiv 
 *          \left(\frac{\partial}{\partial a}\right)^N
 *          \left(\frac{\partial}{\partial b}\right)^L
 *          \left(\frac{\partial}{\partial c}\right)^M
 *          \int_0^1 e^{-Tu^2} \;du
 *  \f]
 *  with 
 *  \f[
 *    T \equiv \alpha \left( a^2+b^2+c^2 \right)
 *  \f]
 *  By extending the above definition to more general
 *  \f[
 *    R_{NLMj} \equiv \left(-\sqrt{\alpha}\right)^{N+L+M}\left(-2\alpha\right)^j
 *             \int_0^1 u^{N+L+M+2j}H_N(au\sqrt{\alpha})H_L(bu\sqrt{\alpha})H_M(cu\sqrt{\alpha})e^{-Tu^2} \;du
 *  \f]
 *  one can see that
 *  \f[
 *     R_{000j} = \left(-2\alpha\right)^j F_j(T)
 *  \f]
 *  The Boys function is here given by
 *  \f[
 *     F_j(T) \equiv \int_0^1 u^{2j}e^{-Tu^2} \;du
 *  \f]
 *  and its efficient implementation can be discussed elsewhere. In Psi4, `psi::Taylor_Fjt` class
 *  is used for this purpose.
 *
 *  Now, it is possible to show that the following recursion relationships 
 *  are true:
 *  \f{align*}{
 *   R_{0,0,M+1,j} &= cR_{0,0,M,j+1} + MR_{0,0,M-1,j+1} \\
 *   R_{0,L+1,M,j} &= bR_{0,L,M,j+1} + LR_{0,L-1,M,j+1} \\
 *   R_{N+1,L,M,j} &= aR_{N,L,M,j+1} + NR_{N-1,L,M,j+1} \\
 *  \f}
 *  This scheme is implemented in OEPDev.
 */

/**
 *  @defgroup OEPDEV_INTEGRAL_HELPERS The Integral Helper Library
 *  @brief
 *  You will find here various iterators to go through shells
 *  while computing ERI, or iterators over ERI itself.
 */

/**
 *  @defgroup OEPDEV_3DFIELDS The Three-Dimensional Scalar Fields Library
 *  @brief
 *  Handles all sorts of scalar distributions in 3D Euclidean space,
 *  such as potentials defined at particular collection of points.
 *  In this Module, you will also find handling both random and ordered
 *  points collections in a form of a G09 cube, as well as handling 
 *  G09 Cube files.
 */

/**
 *  @defgroup OEPDEV_ESP The Multipole Fitting Library
 *  @brief
 *  Implements methods to fit the generalized multipole moments
 *  of a generalized density distribution based on the input generalized
 *  potential scalar field. Among others, you will find here the 
 *  electrostatic potential (ESP) fitting method.
 */

/**
 *  @defgroup OEPDEV_DFT The Density Functional Theory Library
 *  @brief
 *  Implements the OEPDev ab initio DFT methods.
 */

/**
 *  @defgroup OEPDEV_UTILITIES The OEPDev Utilities
 *  @brief
 *  Contains utility functions such as printing OEPDev preambule 
 *  to the output file, class for wavefunction union etc.
 */

/**
 *  @defgroup OEPDEV_TESTS The OEPDev Testing Platform Library
 *  @brief
 *  Testing platform at C++ level of code. You should add more tests
 *  here when developing new functionalities, theories or models.
 */
