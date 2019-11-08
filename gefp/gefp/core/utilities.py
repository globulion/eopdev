#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 GEFP Utilities.
 Bartosz Błasiak, Gundelfingen, 18 Jul 2019
"""
import numpy
import psi4
import oepdev

__all__ = ["psi_molecule_from_file", 
           "psi_supermolecule_from_molecules",
           "wavefunction_union_from_dimer",
           "wavefunction_union_from_dfi_solver",
           "substitute_matrices_in_psi4_wavefunction"]

__all__+= ["UnitaryOptimizer", "UnitaryOptimizer_4_2"]



def psi_molecule_from_file(f, frm=None, no_com=True, no_reorient=True):
    "Construct psi4.core.Molecule object from structure file"
    if frm is None: frm = f.split('.')[-1].lower()
    #
    if   frm == 'xyz':
       qmol = psi4.qcdb.Molecule.init_with_xyz(f, no_com=no_com, no_reorient=no_reorient)  
       mol  = psi4.geometry(qmol.create_psi4_string_from_molecule())
    #
    elif frm == 'psi':
       log = open(f).read()
       mol = psi4.geometry('\n'+log)
    # 
    else: raise ValueError("Unrecognised format - %s -" % frm)
    #
    mol.update_geometry()
    return mol

def psi_supermolecule_from_molecules(*molecules, units='angs'):
    "Create one molecule object from multiple molecules. Assumes C1 symmetry"
    c = 1.0; log_units = "units bohr\n"
    if units.lower().startswith("angs"):
       c = psi4.constants.bohr2angstroms
       log_units = "units angstrom\n"
    log = "\n"
    for molecule in molecules:
        log+= "%d %d\n" % (molecule.molecular_charge(), molecule.multiplicity())
        for i in range(molecule.natom()):
            log += "%s"     %  molecule.symbol(i)
            log += "%16.6f" % (molecule.x(i) * c)
            log += "%16.6f" % (molecule.y(i) * c)
            log += "%16.6f" % (molecule.z(i) * c)
            log += "\n"
        log += log_units
        log += "symmetry c1\n"
        log += "no_reorient\n"
        log += "no_com\n"
                                                                                 
        log += "--\n"
    log = log[:-3]

    aggregate = psi4.geometry(log)
    aggregate.update_geometry()
    return aggregate

def _get_wavefunction_union_basis_sets(dimer):
    "Construct basis set objects necessary to build the oepdev.WavefunctionUnion object"
    dimer.update_geometry()
    pam = psi4.core.get_global_option("PUREAM")

    #
    molecule_A     = dimer.extract_subsets(1)                                                               
    molecule_B     = dimer.extract_subsets(2)

    # --- primary                                                                                                                 
    basis_A        = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("BASIS"), puream=pam)
    basis_B        = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("BASIS"), puream=pam)

    # --- auxiliary (DF-SCF)
    basis_df_scf_A = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=pam)
    basis_df_scf_B = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=pam)

    # --- auxiliary (OEP)
    opt_basis_df_oep_A = psi4.core.get_global_option("DF_BASIS_OEP_A")
    opt_basis_df_oep_B = psi4.core.get_global_option("DF_BASIS_OEP_B")
    opt_basis_df_oep   = psi4.core.get_global_option("DF_BASIS_OEP")
    if not opt_basis_df_oep_A: opt_basis_df_oep_A = opt_basis_df_oep
    if not opt_basis_df_oep_B: opt_basis_df_oep_B = opt_basis_df_oep
    basis_df_oep_A = psi4.core.BasisSet.build(molecule_A, "BASIS", opt_basis_df_oep_A, puream=pam)
    basis_df_oep_B = psi4.core.BasisSet.build(molecule_B, "BASIS", opt_basis_df_oep_B, puream=pam)

    # --- intermediate (OEP)
    if psi4.core.get_global_option("DF_BASIS_INT") == "":
       b_int = psi4.core.get_global_option("DF_BASIS_SCF")
    else:
       b_int = psi4.core.get_global_option("DF_BASIS_INT")

    basis_int_oep_A= psi4.core.BasisSet.build(molecule_A, "BASIS", b_int, puream=pam)
    basis_int_oep_B= psi4.core.BasisSet.build(molecule_B, "BASIS", b_int, puream=pam)

    basis        = psi4.core.BasisSet.build(dimer, "BASIS", psi4.core.get_global_option("BASIS"       ), puream=pam)
    basis_df_scf = psi4.core.BasisSet.build(dimer, "BASIS", psi4.core.get_global_option("DF_BASIS_SCF"), puream=pam)
    #
    return molecule_A, molecule_B, \
           basis, basis_df_scf, basis_A, basis_B, basis_df_oep_A, basis_df_oep_B, \
           basis_df_scf_A, basis_df_scf_B, basis_int_oep_A, basis_int_oep_B


def wavefunction_union_from_dimer(dimer, wfn_1=None, wfn_2=None,
                                         localize_orbitals=False, transform_integrals=False):
    "Create oepdev.WavefunctionUnion from the molecular dimer"

    molecule_A, molecule_B, \
    basis, basis_df_scf, \
    basis_A, basis_B,\
    basis_df_oep_A, basis_df_oep_B,\
    basis_df_scf_A, basis_df_scf_B,\
    basis_int_oep_A, basis_int_oep_B = _get_wavefunction_union_basis_sets(dimer)

    if wfn_1 is None:
       en_1, wfn_1 = psi4.energy('hf', molecule=molecule_A, return_wfn=True); del en_1
    if wfn_2 is None:
       en_2, wfn_2 = psi4.energy('hf', molecule=molecule_B, return_wfn=True); del en_2

    #psi4_io.set_default_path('/home/globulion/scr-d')
    psi4.core.clean()

    wfn_1.set_basisset("BASIS_DF_SCF", basis_df_scf_A)
    wfn_2.set_basisset("BASIS_DF_SCF", basis_df_scf_B)

    un = oepdev.WavefunctionUnion(dimer, basis, basis_df_scf,
                                  basis_A, basis_B,
                                  basis_df_oep_A, basis_df_oep_B,
                                  basis_df_scf_A, basis_df_scf_B,
                                  basis_int_oep_A, basis_int_oep_B,
                                  wfn_1, wfn_2,
                                  psi4.core.get_options())

    if localize_orbitals: 
       un.localize_orbitals()

    if transform_integrals: 
       un.transform_integrals()

    return un

def wavefunction_union_from_dfi_solver(dfi, **kwargs): 
    "Create oepdev.WavefunctionUnion from the oepdev.DFI solver"
    un = wavefunction_union_from_dimer(dfi.aggregate(), wfn_1=dfi.wfn(0), wfn_2=dfi.wfn(1), **kwargs)
    return un


def substitute_matrices_in_psi4_wavefunction(wfn, 
                                             Da, Ca, Fa, ea, 
                                             Db=None, Cb=None, Fb=None, eb=None):
    "Substitutes LCAO matrix, orbital energies and AO matrices in a given psi4.core.Wavefunction object"
    nbf = wfn.basisset().nbf()
    nmo = wfn.nmo()

    if Db is None: Db = Da
    if Cb is None: Cb = Ca
    if Fb is None: Fb = Fa
    if eb is None: eb = ea

    Ca_ = wfn.Ca(); Cb_ = wfn.Cb()
    Da_ = wfn.Da(); Db_ = wfn.Db()
    Fa_ = wfn.Fa(); Fb_ = wfn.Fb()
    ea_ = wfn.epsilon_a(); eb_ = wfn.epsilon_b()

    for a in range(nbf):
        for i in range(nmo):
            Ca_.set(a, i, Ca[a,i])
            Cb_.set(a, i, Cb[a,i])
        for b in range(nbf):
            Da_.set(a, b, Da[a,b])
            Db_.set(a, b, Db[a,b])
        for b in range(nbf):
            Fa_.set(a, b, Fa[a,b])
            Fb_.set(a, b, Fb[a,b])
    for i in range(nmo):
            ea_.set(i, ea[i])
            eb_.set(i, eb[i])
    return


class UnitaryOptimizer_4_2(object):
  """
 ---------------------------------------------------------------------------------

 Finds the unitary matrix X that optimizes the following function:
 
 Z(X) = \sum_{ijklmn} X_{ki} X_{lj} X_{mi} X_{nj} R_{ijklmn} 
      + \sum_{ijk}    X_{ji} X_{ki}               P_{ijk}
 
 where 
   * X is a square unitary matrix of size N x N
   * R is a general real 6-th rank tensor of size N^6
   * P is a general real 3-rd rank tensor of size N^3
 
 Usage:
   optimizer = UnitaryOptimizer(R, P, conv=1.0e-8, maxiter=100, verbose=True)
   optimizer.maximize() #or minimize()
   X = optimizer.X
   Z = optimizer.Z
 
 ---------------------------------------------------------------------------------
                                                       Last Revision: 07.04.2018
 """
  def __init__(self, R, P, conv=1.0e-8, maxiter=100, verbose=True):
      """Initialize with R and P matrix, as well as optimization options"""
      self._R   = R.copy()
      self._P   = P.copy()
      self._R0  = R.copy()
      self._P0  = P.copy()
      self.X    = None
      self._d   = len(P)
      # optimization options
      self.conv   = conv
      self.maxiter=maxiter
      self.verbose=verbose
      # initial Z value
      self._Zinit = self._eval_Z(numpy.identity(self._d), self._R0, self._P0)

  def maximize(self):
      """Maximize the Z function under unitary constraint for X""" 
      self.run('max')

  def minimize(self):
      """Minimize the Z function under unitary constraint for X"""
      self.run('min')

  def run(self, opt='minimize'):
      """Perform the optimization"""
      assert (opt.lower().startswith('min') or opt.lower().startswith('max')), 'Unrecognized optimization mode < %s >' % opt
      self._refresh()
      self._run(opt.lower())

  @property
  def Z(self):
      """Return the current value of objective function"""
      z = self._eval_Z(self.X, self._R0, self._P0)
      return z

  # -- protected
  def _run(self, opt):
      """Perform the optimization (protected interface)"""
      conv = 1e8
      #Xold = numpy.identity(self._d)
      Zold = self._Zinit
      Xacc = numpy.identity(self._d)
      success = False
      if self.verbose: 
         print    (" Start  : Z[1] = %15.6f" % Zold)
      niter = 0
      while (conv > self.conv):
         i, j, gamma = self._find_next(opt)
         #print(" Chosen x = %14.4f" % gamma)
         Xnew = self._form_X(i, j, gamma)
         self._update_RP(Xnew)
         Znew = self._eval_Z(numpy.identity(self._d), self._R, self._P)
         conv = abs(Znew-Zold)
         #Xold = Xnew.copy()
         Zold = Znew
         niter += 1
         Xacc = numpy.dot(Xacc, Xnew)
         if self.verbose:
            print (" Iter %2d: Z[X] = %15.6f  Conv= %15.6f" % (niter, Znew, conv))
         if niter > self.maxiter: 
            print(" Optimization unsuccesfull! Maximum iteration number exceeded!")
            success = False
            break
      success = True if niter <= self.maxiter else False
      self.X = Xacc.copy()
      if (self.verbose and success):
         print(" Optimization succesfull!\n")
         print(" Optimized Z[X] value: %15.6f" % self.Z)

  def _update_RP(self, X):
      """Update R and P tensors by transforming them by using X matrix"""
      P = numpy.tensordot(self._P, X, (1,0))
      P = P.transpose(0,2,1)
      P = numpy.tensordot(P, X, (2,0))

      R = numpy.tensordot(self._R, X, (2,0)) 
      R = R.transpose(0,1,5,2,3,4)
      R = numpy.tensordot(R  , X, (3,0))
      R = R.transpose(0,1,2,5,3,4)
      R = numpy.tensordot(R  , X, (4,0))
      R = R.transpose(0,1,2,3,5,4)
      R = numpy.tensordot(R  , X, (5,0))

      #P = self._P.copy()
      #R = self._R.copy()
      #P.fill(0); R.fill(0)
      #N = self._d

      #for i in range(N):
      #    for J in range(N): 
      #        for K in range(N):
      #            v = 0.0
      #            for j in range(N):
      #                for k in range(N):
      #                    v += X[j,J] * X[k,K] * self._P[i,j,k]
      #            P[i,J,K] = v


      #for i in range(N):
      #    for j in range(N):
      #        for K in range(N):
      #            for L in range(N):
      #                for M in range(N):
      #                    for NN in range(N):
      #                        v = 0.0
      #                        for k in range(N):
      #                            for l in range(N):
      #                                for m in range(N):
      #                                    for n in range(N):
      #                                        v += X[k,K] * X[l,L] * X[m,M] * X[n,NN] * self._R[i,j,k,l,m,n]
      #                        R[i,j,K,L,M,NN] = v
      self._P = P.copy()
      self._R = R.copy()

  def _refresh(self):
      """Restore the initial state of the optimizer"""
      self._R = self._R0.copy()
      self._P = self._P0.copy()
      self.X  = None

  def _find_next(self, opt):
      """Determine next pair of degrees of freedom for 2D rotation"""
      optfunc = operator.lt if opt.startswith('min') else operator.gt
      I, J = 0, 1
      Gamma = None
      dZold = 1e8 if opt.startswith('min') else -1e8
      for j in range(self._d):
          for i in range(j):
              a0, a1, a2, a3, a4, b1, b2, b3, b4 = self._get_Fourier_coeffs(i, j)
              gamma   = self._find_x(a0, a1, a2, a3, a4, b1, b2, b3, b4, i, j, opt)
              dZ = self._eval_dZ(gamma, self._P, self._R, i, j)
              if optfunc(dZ, dZold):
                 Gamma = gamma
                 I = i
                 J = j
                 dZold = dZ
      return I, J, Gamma

  def _find_x(self, a0, a1, a2, a3, a4, b1, b2, b3, b4, i, j, opt):
      """Find the optimal 2D rotation angle"""
      # Boyd's method in 4 dimensions: Boyd, J.P.; J. Eng. Math. (2006) 56:203-219
      d = a4 - 1.0j * b4
      K = numpy.zeros((8,8), numpy.complex64)
      K[0,1] = 1.0
      K[1,2] = 1.0
      K[2,3] = 1.0
      K[3,4] = 1.0
      K[4,5] = 1.0
      K[5,6] = 1.0
      K[6,7] = 1.0
      K[7,0] =-(a4 + 1.0j * b4) / d
      K[7,1] =-(a3 + 1.0j * b3) / d
      K[7,2] =-(a2 + 1.0j * b2) / d
      K[7,3] =-(a1 + 1.0j * b1) / d
      K[7,4] =- a0 * 2.0        / d
      K[7,5] =-(a1 - 1.0j * b1) / d
      K[7,6] =-(a2 - 1.0j * b2) / d
      K[7,7] =-(a3 - 1.0j * b3) / d
      #
      E, X = numpy.linalg.eig(K)
      X   = -1.0j * numpy.log(E)
      #print "Imaginary part of X: "
      #libbbg.utilities.PRINT(X.imag)
      X = X.real
      X[numpy.where(X<0.0)] += numpy.pi
      #print "Real      part of X: "
      #libbbg.utilities.PRINT(X)

      # Find optimal gamma 
      gamma = None
      if opt.startswith('min'):
         Zold = 1.0e8
         for x in X:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z < Zold:
                gamma = x
                Zold = Z
      else:
         Zold = -1e8
         for x in X:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z > Zold:
                gamma = x
                Zold = Z
      assert gamma is not None, "Error while searching for optimum!"
      return gamma

  def _get_Fourier_coeffs(self, I, J):
      """Retrieve ABCD parameters for root search"""
      a0, a1, a2, a3, a4, b1, b2, b3, b4 = 0, 0, 0, 0, 0, 0, 0, 0, 0
      d  = lambda i, j: 0.0 if (i!=j) else 1.0
      a_ = lambda i, k: (-d(I,i)*d(J,k)+d(I,k)*d(J,i))*(1.0-d(i,k))
      b_ = lambda i, k: d(i,k) * (d(I,k)+d(J,i))
      c_ = lambda i, k: d(i,k) * (1.0-d(I,k))*(1.0-d(J,i))
      N = self._d

      # P-contribution
      for i in range(N):
          for j in range(N):
              for k in range(N):
                  p   = self._P[i,j,k]
                  A =   (d(I,j)*d(J,i) - d(I,i)*d(J,j)) * (1-d(i,j))
                  B =   -d(i,j)*(d(I,j)+d(J,i))
                  C =    d(i,k)*(1-d(I,k))*(1-d(J,i))
                  D =    d(i,k)*(d(I,k)+d(J,i))
                  E =   -d(I,i)*d(J,k)*(1-d(i,k))
                  F =    d(I,k)*d(J,i)*(1-d(i,k))
                  G =    (1-d(i,k))*(d(I,k)*d(J,i) - d(I,i)*d(J,k))
                  H =   -d(i,k)*(d(I,k)+d(J,i))
                  II=    d(i,j)*(1-d(I,j))*(1-d(J,i))
                  JJ=    d(i,j)*(d(I,j)+d(J,i))
                  K =   -d(I,i)*d(J,j)*(1-d(i,j))
                  L =    d(I,j)*d(J,i)*(1-d(i,j))
                  a0 += p * (A*D+G*JJ+B*(E+F)+H*(K+L)) / 2.0
                  a1 += p * (A*C+G*II)
                  a2 += p * (A*D+G*JJ-B*(E+F)-H*(K+L)) / 2.0
                  b1 += p * (B*C+H*II)
                  b2 += p * (H*JJ+B*D+G*(K+L)+A*(E+F)) / 2.0

      # R-contribution
      for i in range(N):
          for j in range(N):
              for k in range(N):
                  for l in range(N):
                      for m in range(N):
                          for n in range(N):
                              r = self._R[i,j,k,l,m,n] 
                              # First R batch
                              A = a_(i,k)
                              B =-b_(i,k)
                              C = a_(i,m)
                              D = c_(i,m)
                              E = b_(i,m)
                              F = a_(j,l)
                              G = c_(j,l)
                              H = b_(j,l)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Second R batch
                              A = a_(i,m)
                              B =-b_(i,m)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(j,l)
                              G = c_(j,l)
                              H = b_(j,l)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Third R batch
                              A = a_(j,l)
                              B =-b_(j,l)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(i,m)
                              G = c_(i,m)
                              H = b_(i,m)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Fourth R batch
                              A = a_(j,n)
                              B =-b_(j,n)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(i,m)
                              G = c_(i,m)
                              H = b_(i,m)
                              II= a_(j,l)
                              JJ= c_(j,l)
                              K = b_(j,l)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
      
      return a0, a1, a2, a3, a4, b1, b2, b3, b4

  def _eval_Z(self, X, R, P):
      """Evaluate the objective Z function"""            
      Z = 0.0
      N = self._d
      p = numpy.tensordot(P, X, (1,0))
      p = numpy.tensordot(p, X, (1,0))
      Z = p.diagonal().diagonal().sum()
     
      r = numpy.tensordot(R, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      for i in range(N): 
          for j in range(N): 
              Z+=r[i,j,i,j,i,j]

      #for i in range(N):
      #    for j in range(N):
      #        for k in range(N):
      #            Z += P[i,j,k] * X[j,i] * X[k,i]
      #            for l in range(N):
      #                for m in range(N):
      #                    for n in range(N):
      #                        pass
      #                        Z += R[i,j,k,l,m,n] * X[k,i] * X[l,j] * X[m,i] * X[n,j]
      
      return Z

  def _eval_dZ(self, g, P, R, I, J):
      """Compute the change in Z"""
      X = self._form_X(I, J, g)
      E = numpy.identity(self._d)
      Z_new = self._eval_Z(X, R, P)
      Z_old = self._eval_Z(E, R, P)
      dZ = Z_new - Z_old
      return dZ

  def _form_X(self, i, j, gamma):
      """Form unitary matrix X"""                                        
      X = numpy.identity(self._d)
      T = numpy.cos(gamma)
      g = numpy.sin(gamma)
      X[i,i] = T
      X[j,j] = T
      X[i,j] = g
      X[j,i] =-g
      return X


class UnitaryOptimizer(object):
  """
 ---------------------------------------------------------------------------------

 Finds the unitary matrix X that optimizes the following function:
 
 Z(X) = \sum_{ijkl} X_{ij}X_{kl} R_{jl} - \sum_{ij} X_{ij}P_{j} 
 
 where 
   * X is a square unitary matrix of size N x N
   * R is a square, in general non-symmetric matrix of size N x N
   * P is a vector of length N   
 
 Usage:
   optimizer = UnitaryOptimizer(R, P, conv=1.0e-8, maxiter=100, verbose=True)
   optimizer.maximize() #or minimize()
   X = optimizer.X
   Z = optimizer.Z
 
 ---------------------------------------------------------------------------------
                                                       Last Revision: 25.03.2018
 """
  def __init__(self, R, P, conv=1.0e-8, maxiter=100, verbose=True):
      """Initialize with R and P matrix, as well as optimization options"""
      self._R   = R.copy()
      self._P   = P.copy()
      self._R0  = R.copy()
      self._P0  = P.copy()
      self.X    = None
      self._d   = P.size
      # optimization options
      self.conv   = conv
      self.maxiter=maxiter
      self.verbose=verbose
      # initial Z value
      self._Zinit = self._eval_Z(numpy.identity(self._d), self._R0, self._P0)
      # functions do find roots in 2D unitary optimization step
      self._f = lambda x, A, B, C, D:   A*numpy.sin(x)+B*numpy.cos(x)+  C*numpy.sin(2*x)+  D*numpy.cos(2*x)
      self._fg= lambda x, A, B, C, D:   A*numpy.cos(x)-B*numpy.sin(x)+2*C*numpy.cos(2*x)-2*D*numpy.sin(2*x)
      self._fh= lambda x, A, B, C, D: -(A*numpy.sin(x)+B*numpy.cos(x)+4*C*numpy.sin(2*x)+4*D*numpy.cos(2*x))

  def maximize(self):
      """Maximize the Z function under unitary constraint for X""" 
      self.run('max')

  def minimize(self):
      """Minimize the Z function under unitary constraint for X"""
      self.run('min')

  def run(self, opt='minimize'):
      """Perform the optimization"""
      assert (opt.lower().startswith('min') or opt.lower().startswith('max')), 'Unrecognized optimization mode < %s >' % opt
      self._refresh()
      self._run(opt.lower())

  @property
  def Z(self):
      """Return the current value of objective function"""
      z = self._eval_Z(self.X, self._R0, self._P0)
      return z

  # -- protected
  def _run(self, opt):
      """Perform the optimization (protected interface)"""
      conv = 1e8
      Xold = numpy.identity(self._d)
      Zold = self._Zinit
      Xacc = numpy.identity(self._d)
      success = False
      if self.verbose: 
         print    (" Start  : Z[1] = %15.6f" % Zold)
      niter = 0
      while (conv > self.conv):
         i, j, gamma = self._find_next(opt)
         Xnew = self._form_X(i, j, gamma)
         self._uptade_RP(Xnew)
         Znew = self._eval_Z(numpy.identity(self._d), self._R, self._P)
         conv = abs(Znew-Zold)
         Xold = Xnew.copy()
         Zold = Znew
         niter += 1
         Xacc = numpy.dot(Xnew, Xacc)
         if self.verbose:
            print (" Iter %2d: Z[X] = %15.6f  Conv= %15.6f" % (niter, Znew, conv))
         if niter > self.maxiter: 
            print(" Optimization unsuccesfull! Maximum iteration number exceeded!")
            success = False
            break
      success = True if niter <= self.maxiter else False
      self.X = Xacc
      if (self.verbose and success):
         print(" Optimization succesfull!\n")
         print(" Optimized Z[X] value: %15.6f" % self.Z)

  def _uptade_RP(self, X):
      self._P = numpy.dot(X, self._P)
      self._R = numpy.dot(X, numpy.dot(self._R, X.T))

  def _refresh(self):
      """Restore the initial state of the optimizer"""
      self._R = self._R0.copy()
      self._P = self._P0.copy()
      self.X  = None

  def _find_next(self, opt):
      """Determine next pair of degrees of freedom for 2D rotation"""
      optfunc = operator.lt if opt.startswith('min') else operator.gt
      I, J = 0, 1
      Gamma = 0.0
      dZold = 1e8 if opt.startswith('min') else -1e8
      for j in range(self._d):
          for i in range(j):
              A,B,C,D = self._get_ABCD(i, j)
              gamma   = self._find_x(A, B, C, D, i, j, opt)
              dZ = self._eval_dZ(gamma, self._P, self._R, i, j)
              if optfunc(dZ, dZold):
                 Gamma = gamma
                 I = i
                 J = j
                 dZold = dZ
      return I, J, Gamma

  def _find_x(self, A, B, C, D, i, j, opt):
      """Find the optimal 2D rotation angle"""
      #f = lambda x, A, B, C, D:   A*numpy.sin(x)+B*numpy.cos(x)+  C*numpy.sin(2*x)+  D*numpy.cos(2*x)
      #fg= lambda x, A, B, C, D:   A*numpy.cos(x)-B*numpy.sin(x)+2*C*numpy.cos(2*x)-2*D*numpy.sin(2*x)
      #fh= lambda x, A, B, C, D: -(A*numpy.sin(x)+B*numpy.cos(x)+4*C*numpy.sin(2*x)+4*D*numpy.cos(2*x))

      # Boyd's method in 4 dimensions: Boyd, J.P.; J. Eng. Math. (2006) 56:203-219
      d = D - 1.0j * C                        
      K = numpy.zeros((4,4), numpy.complex64)
      K[0,1] = 1.0
      K[1,2] = 1.0
      K[2,3] = 1.0
      K[3,0] = -(D + 1.0j * C)/d
      K[3,1] = -(B + 1.0j * A)/d
      K[3,2] = 0.0
      K[3,3] = -(B - 1.0j * A)/d
      E, X = numpy.linalg.eig(K)
      X   = -1.0j * numpy.log(E)
      X[numpy.where(X<0.0)] += 2.0 * numpy.pi
      X = X.real

      # discriminate between minima and maxima
      Xmin = list()
      Xmax = list()
      for x in X:
          F = self._f (x,A,B,C,D)
          g = self._fg(x,A,B,C,D)
          if   g> 0.0: Xmin.append(x)
          elif g< 0.0: Xmax.append(x)
          #else: raise ValueError("The Hessian of objective function is zero at X=%15.5f" % x)
      Xmin = numpy.array(Xmin)
      Xmax = numpy.array(Xmax)
    
      # Find optimal gamma 
      gamma = None
      if opt.startswith('min'):
         Zold = 1.0e8
         for x in Xmin:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z < Zold:
                gamma = x
             Zold = Z
      else:
         Zold = -1e8
         for x in Xmax:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z > Zold:
                gamma = x
             Zold = Z
      assert gamma is not None, "Error while searching for optimum!"
      return gamma

  def _get_ABCD(self, i, j):
      """Retrieve ABCD parameters for root search"""
      A = self._P[i]+self._P[j]
      B = self._P[i]-self._P[j] 
      C =-2.*(self._R[i,j]+self._R[j,i])
      D = 2.*(self._R[j,j]-self._R[i,i])
      ii = numpy.setdiff1d(numpy.arange(self._P.size), [i,j])
      A -= self._R[j,ii].sum() + self._R[i,ii].sum() + self._R[ii,j].sum() + self._R[ii,i].sum()
      B += self._R[j,ii].sum() - self._R[i,ii].sum() + self._R[ii,j].sum() - self._R[ii,i].sum()
      return A, B, C, D

  def _eval_Z(self, X, R, P):
      """Evaluate the objective Z function"""            
      z1 = numpy.dot(X, numpy.dot(R,X.T))
      z2 = numpy.dot(X, P)
      return z1.sum() - z2.sum()

  def _eval_dZ(self, g, P, R, i, j):
      """Compute the change in Z"""
      dZ = (1.0 - numpy.cos(g)) * (P[i] + P[j]) \
                + numpy.sin(g)  * (P[i] - P[j])       \
                + numpy.sin(2.*g) * (R[j,j] - R[i,i]) \
             - 2.*numpy.sin(g)**2 * (R[i,j] + R[j,i])
      ii = numpy.setdiff1d(numpy.arange(P.size), [i,j])
      dZ-= (1-numpy.cos(g))*(R[j,ii].sum() + R[i,ii].sum() + R[ii,j].sum() + R[ii,i].sum())
      dZ+=    numpy.sin(g) *(R[j,ii].sum() - R[i,ii].sum() + R[ii,j].sum() - R[ii,i].sum())
      return dZ

  def _form_X(self, i, j, gamma):
      """Form unitary matrix X"""                                        
      X = numpy.identity(self._d)
      T = numpy.cos(gamma)
      g = numpy.sin(gamma)
      X[i,i] = T
      X[j,j] = T
      X[i,j] = g
      X[j,i] =-g
      return X
 
  # private
  def __same(self, x, X, tol=0.0001):
      """Check if the solution was already obtained"""
      #pipi = 2*numpy.pi
      result = False
      for xp in X:
          dx = x-xp
          if (numpy.abs(dx) < tol) or (numpy.abs(dx+2*numpy.pi) < tol) or (numpy.abs(dx-2*numpy.pi) < tol): 
          #if not (numpy.abs(dx)%pipi):
             result = True
             break
      return result
