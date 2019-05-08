#!/usr/bin/python3
"""
 ---------------------------------------------

 The DMFT Universal Surface Module.

 a0 = f(x, y)

 at the FCI/STO-3G level of theory.

 ---------------------------------------------
       Created     : Gundelfingen, 5 May 2019
       Last Updated: Gundelgingen, 8 May 2019
"""
import gefp
import numpy, math, psi4
import scipy.optimize
import sys
sys.stdout.flush()

__all__ = ['gen_surface', 'gen_surface_mbb', 'Data', 'OUTFILE']

# defaults
OUTFILE='surface.dat'
ERRFILE='error.dat'
SYSTEM_ID=1

class Data:
    """\
Data structure with the universal surface information.

Contains:
    o number of alpha electrons
    o DMFT coefficient of interpolation
    o Total energies (REF, FUN, MBB)
    o Scalar and matrix electron correlation indicators (REF, FUN, MBB)
    o Dipole moments (REF, FUN, MBB)

 REF - Reference QM method
 FUN - DMFT interpolation functional
 MBB - Muller-Buijse-Baerends functional
"""
    nan = numpy.NaN

    def __init__(self, outfile=OUTFILE):
        "Namespace of all variables"
        self._outfile         = outfile
        self._data            = None  # All data
        self._data_n_alpha    = []  # N
        self._data_t1         = []  # Coefficient DMFT
        self._data_e_ref      = []  # Energy: Reference
        self._data_e_fun      = []  # Energy: DMFT
        self._data_e_mbb      = []  # Energy: MBB
        self._data_i_d_ref    = []  # I_DYN scalar: Reference
        self._data_i_n_ref    = []  # I_NON scalar: Reference
        self._data_i_d_fun    = []  # I_DYN scalar: DMFT
        self._data_i_n_fun    = []  # I_NON scalar: DMFT
        self._data_i_d_mbb    = []  # I_DYN scalar: MBB 
        self._data_i_n_mbb    = []  # I_NON scalar: MBB 
        self._data_I_d_ref    = []  # I_DYN matrix: Reference
        self._data_I_n_ref    = []  # I_NON matrix: Reference
        self._data_I_d_fun    = []  # I_DYN matrix: DMFT
        self._data_I_n_fun    = []  # I_NON matrix: DMFT
        self._data_I_d_mbb    = []  # I_DYN matrix: MBB 
        self._data_I_n_mbb    = []  # I_NON matrix: MBB 
        self._data_mom_dx_ref = []  # Dipole moment-X: Reference
        self._data_mom_dy_ref = []  # Dipole moment-Y: Reference
        self._data_mom_dz_ref = []  # Dipole moment-Z: Reference
        self._data_mom_dx_fun = []  # Dipole moment-X: DMFT
        self._data_mom_dy_fun = []  # Dipole moment-Y: DMFT
        self._data_mom_dz_fun = []  # Dipole moment-Z: DMFT
        self._data_mom_dx_mbb = []  # Dipole moment-X: MBB 
        self._data_mom_dy_mbb = []  # Dipole moment-Y: MBB 
        self._data_mom_dz_mbb = []  # Dipole moment-Z: MBB 
        self._create_data()
        #
        self._out  = open(self._outfile, 'w')
        pass

    @property
    def data(self):
        "Data structure in array format"
        return numpy.array(self._data).T

    @staticmethod
    def read(surface_file):
        "Read the data from the file"
        obj = Data(outfile='none')
        data = numpy.mafromtxt(surface_file).data
        for i in data:
            obj.add(\
            n_alpha=i[0], t1=i[1], 
            e_ref=i[2], e_fun=i[3], e_mbb=i[4],
            i_d_ref=i[5], i_n_ref=i[6],
            i_d_fun=i[7], i_n_fun=i[8],
            i_d_mbb=i[9], i_n_mbb=i[10],
            I_d_ref=i[11], I_n_ref=i[12],
            I_d_fun=i[13], I_n_fun=i[14],
            I_d_mbb=i[15], I_n_mbb=i[16],
            mom_d_ref=[i[17],i[18],i[19]], 
            mom_d_fun=[i[20],i[21],i[22]],
            mom_d_mbb=[i[23],i[24],i[25]])
        return obj

    def add(self, 
            n_alpha=nan, t1=nan, 
            e_ref=nan, e_fun=nan, e_mbb=nan,
            i_d_ref=nan, i_n_ref=nan,
            i_d_fun=nan, i_n_fun=nan,
            i_d_mbb=nan, i_n_mbb=nan,
            I_d_ref=nan, I_n_ref=nan,
            I_d_fun=nan, I_n_fun=nan,
            I_d_mbb=nan, I_n_mbb=nan,
            mom_d_ref=[nan,nan,nan], 
            mom_d_fun=[nan,nan,nan],
            mom_d_mbb=[nan,nan,nan],
            ):
        "Add new entry to the data"
        self._data_n_alpha.append(n_alpha)
        self._data_t1.append(t1)
        self._data_e_ref.append(e_ref)
        self._data_e_fun.append(e_fun)
        self._data_e_mbb.append(e_mbb)
        self._data_i_d_ref.append(i_d_ref)
        self._data_i_n_ref.append(i_n_ref)
        self._data_i_d_fun.append(i_d_fun)
        self._data_i_n_fun.append(i_n_fun)
        self._data_i_d_mbb.append(i_d_mbb)
        self._data_i_n_mbb.append(i_n_mbb)
        self._data_I_d_ref.append(I_d_ref)
        self._data_I_n_ref.append(I_n_ref)
        self._data_I_d_fun.append(I_d_fun)
        self._data_I_n_fun.append(I_n_fun)
        self._data_I_d_mbb.append(I_d_mbb)
        self._data_I_n_mbb.append(I_n_mbb)
        self._data_mom_dx_ref.append(mom_d_ref[0])
        self._data_mom_dy_ref.append(mom_d_ref[1])
        self._data_mom_dz_ref.append(mom_d_ref[2])
        self._data_mom_dx_fun.append(mom_d_fun[0])
        self._data_mom_dy_fun.append(mom_d_fun[1])
        self._data_mom_dz_fun.append(mom_d_fun[2])
        self._data_mom_dx_mbb.append(mom_d_mbb[0])
        self._data_mom_dy_mbb.append(mom_d_mbb[1])
        self._data_mom_dz_mbb.append(mom_d_mbb[2])
        return

    def dipole_norm(self):
        "Compute norm of dipole moments for REF, FUN and MBB"
        data = self.data
        dipole_ref = data[:,17:19+1]
        dipole_fun = data[:,20:22+1]
        dipole_mbb = data[:,23:25+1]
        d_ref = numpy.linalg.norm(dipole_ref, axis=1)
        d_fun = numpy.linalg.norm(dipole_fun, axis=1)
        d_mbb = numpy.linalg.norm(dipole_mbb, axis=1)
        return d_ref, d_fun, d_mbb

    def write_last(self):
        "Write last data to a file"
        data = self.data
        out = self._out
        assert(len(data)!=0)
        #
        if len(data) == 1:
           line = "%6s %14s" % ('N'.rjust(6), 'T'.rjust(14))                                               
           line+= "%18s %18s %18s" % ('E_ref'.rjust(18), 'E_fun'.rjust(18), 'E_mbb'.rjust(18)) 
           line+= 12*"%18s" \
                   % ('i_d_ref'.rjust(18), 'i_n_ref'.rjust(18), 
                      'i_d_fun'.rjust(18), 'i_n_fun'.rjust(18), 
                      'i_d_mbb'.rjust(18), 'i_n_mbb'.rjust(18),
                      'I_d_ref'.rjust(18), 'I_n_ref'.rjust(18), 
                      'I_d_fun'.rjust(18), 'I_n_fun'.rjust(18), 
                      'I_d_mbb'.rjust(18), 'I_n_mbb'.rjust(18))
           line+= 9*"%16s" \
                   % ('dx_ref'.rjust(16), 'dy_ref'.rjust(16), 'dz_ref'.rjust(16), 
                      'dx_fun'.rjust(16), 'dy_fun'.rjust(16), 'dz_fun'.rjust(16),
                      'dx_mbb'.rjust(16), 'dy_mbb'.rjust(16), 'dz_mbb'.rjust(16))
           out.write('#' + line + '\n')
           #
           line = "%6d %14d"  % (1, 2)
           line+= "%18d %18d %18d" % (3, 4, 5)
           line+=12*"%18d" % (6,7,8,9,10,11,12,13,14,15,16,17)
           line+= 9*"%16d" % (18,19,20,21,22,23,24,25,26)
           out.write('#' + line + '\n')
        #
        di = data[-1]                                                                  
        line = "%6.1f %14.6f"  % (di[0], di[1])
        line+= "%18.8f %18.8f %18.8f" % (di[2], di[3], di[4])
        line+=12*"%18.6f" % (di[5], di[6], di[7], di[8], di[9], di[10], di[11], di[12], di[13], di[14], di[15], di[16])
        line+= 9*"%16.4E" % (di[17], di[18], di[19], di[20], di[21], di[22], di[23], di[24], di[25])
        out.write(' ' + line + '\n')
        out.flush()
        #out.close()
        return

    def _create_data(self):
        "Create the data structure"
        self._data = [\
                      self._data_n_alpha          ,
                      self._data_t1               ,
                      self._data_e_ref            ,       
                      self._data_e_fun            ,
                      self._data_e_mbb            ,
                      self._data_i_d_ref          ,
                      self._data_i_n_ref          ,
                      self._data_i_d_fun          ,
                      self._data_i_n_fun          ,
                      self._data_i_d_mbb          ,
                      self._data_i_n_mbb          ,
                      self._data_I_d_ref          ,
                      self._data_I_n_ref          ,
                      self._data_I_d_fun          ,
                      self._data_I_n_fun          ,
                      self._data_I_d_mbb          ,
                      self._data_I_n_mbb          ,
                      self._data_mom_dx_ref       ,
                      self._data_mom_dy_ref       ,
                      self._data_mom_dz_ref       ,
                      self._data_mom_dx_fun       ,
                      self._data_mom_dy_fun       ,
                      self._data_mom_dz_fun       ,
                      self._data_mom_dx_mbb       ,
                      self._data_mom_dy_mbb       ,
                      self._data_mom_dz_mbb       ,]

# ---> Scan functions <--- #

def get_mol(xyzfile):
    "Create molecule object from XYZ file"
    qmol = psi4.qcdb.Molecule.init_with_xyz(xyzfile, no_com=True, no_reorient=True)
    lmol = psi4.geometry(qmol.create_psi4_string_from_molecule())
    lmol.fix_com(True)
    lmol.fix_orientation(True)
    lmol.update_geometry()
    return lmol

def scan_h2_stretch(x):
    global SYSTEM_ID
    mol = get_mol('xyz/h2.xyz')
    gefp.math.matrix.move_atom_along_bond(mol, 2, 1, x, units='ang')
    mol.save_xyz_file('xyz/ID_%03d_h2.s-sc.%03.3f.xyz' % (SYSTEM_ID, x), True)
    return mol

def scan_xhn_stretch_1h(x, molname):
    global SYSTEM_ID
    if   molname == 'h2o': fname = 'h2o.xyz'
    elif molname == 'h3o': fname = 'h3o.xyz'
    elif molname == 'nh3': fname = 'nh3.xyz'
    elif molname == 'ch4': fname = 'ch4.xyz'
    else: raise ValueError('This molecule is not included in surface yet!')
    mol = get_mol(moldir)
    gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])
    gefp.math.matrix.move_atom_along_bond(mol, 2, 1, x, units='ang')
    mol.save_xyz_file('xyz/ID_%03d_%s.s-1h.%03.3f.xyz' % (SYSTEM_ID, fname[:-4], x), True)
    return mol

def scan_xhn_stretch_nh(x, molname):
    global SYSTEM_ID
    if   molname == 'h2o': 
         fname  = 'h2o.xyz'
         atids  = [2,3]
    elif molname == 'h3o': 
         fname  = 'h3o.xyz'
         atids  = [2,3,4]
    elif molname == 'nh3': 
         fname  = 'nh3.xyz'
         atids  = [2,3,4]
    elif molname == 'ch4': 
         fname = 'ch4.xyz'
         atids  = [2,3,4,5]
    else: raise ValueError('This molecule is not included in surface yet!')
    mol = get_mol(moldir)
    gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])
    gefp.math.matrix.move_atom_symmetric_stretch(mol, atids, 1, x, units='ang')
    mol.save_xyz_file('xyz/ID_%03d_%s.s-nh.%03.3f.xyz' % (SYSTEM_ID, fname[:-4], x), True)
    return mol

def scan_h4_symmetric_stretch(x):
    global SYSTEM_ID
    mol = get_mol('xyz/h4.xyz')
    gefp.math.matrix.move_atom_scale_coordinates(mol, x)
    gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])
    mol.save_xyz_file('xyz/ID_%03d_h4.s-sy.%03.3f.xyz' % (x), True)
    return mol

def scan_h2h2_pull_1h(x):
    global SYSTEM_ID
    mol = get_mol('xyz/h2h2.xyz')
    gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])
    gefp.math.matrix.move_atom_along_bond(mol, 4, 3, x, units='ang')
    mol.save_xyz_file('xyz/ID_%03d_h2h2.p-1h.%03.3f.xyz' % (SYSTEM_ID, x), True)
    return mol

def scan_h2h2_pull_2h(x):
    global SYSTEM_ID
    mol = get_mol('xyz/h2h2.xyz')
    gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])
    gefp.math.matrix.move_atom_along_bond(mol, 3, 1, x, units='ang')
    gefp.math.matrix.move_atom_along_bond(mol, 4, 1, x, units='ang')
    mol.save_xyz_file('xyz/ID_%03d_h2h2.p-2h.%03.3f.xyz' % (SYSTEM_ID, x), True)
    return mol


# ---> Input data functions <--- #

def continue_if_runtime_error(func):
    "Do not stop the calculations when something goes wrong with full QM calculations"
    def wrapper(*args, **kwargs):
        errfile = open(ERRFILE, 'a')              
        try:
            func(*args,**kwargs)
        except RuntimeError:
            print("Reference method iterations not converged")
            log = "%s" % str(kwargs)
            errfile.write(log + '\n')
                                                  
        errfile.close()
        psi4.core.clean()
    return wrapper

@continue_if_runtime_error
def put_input(do_fci, mols, wfns, refs, dips, cors, scan_func, **kwargs):
    "Run the RHF and reference methods on each molecule"
    global SYSTEM_ID
    mol = scan_func(**kwargs)
    mols.append(mol)
    # RHF wavefunction
    e_hf, w_hf = psi4.driver.energy("hf", molecule=mol, return_wfn=True)
    wfns.append(w_hf)
    if do_fci:
       # FCI energy                                                                                                         
       e_ci, w_ci = psi4.driver.properties("fci", molecule=mol, return_wfn=True, properties=["NO_OCCUPATIONS"])
       refs.append(e_ci)
       # FCI dipole moment
       Da = w_ci.Da().to_array(dense=True)
       mints = psi4.core.MintsHelper(w_ci.basisset())
       T     = mints.ao_dipole()
       dip_x = 2.0 * numpy.dot(T[0].to_array(dense=True), Da).trace()
       dip_y = 2.0 * numpy.dot(T[1].to_array(dense=True), Da).trace()
       dip_z = 2.0 * numpy.dot(T[2].to_array(dense=True), Da).trace()
       dips.append([dip_x, dip_y, dip_z])
       # FCI correlation indicators
       S = w_ci.S().to_array(dense=True)
       Ca= w_ci.Ca_subset("AO","ALL").to_array(dense=True)
       n, c = gefp.density.partitioning.Density.natural_orbitals(Da, S, Ca, orthogonalize_mo=True,
                                          order='descending', no_cutoff=0.0, renormalize=False, return_ao_orthogonal=False,
                                          ignore_large_n=False)
       ns = n.copy(); ns[ns<0.0] = 0.0
       #S   = numpy.sqrt(ns).sum()
       #N   = ns.sum()
       i_n = (ns*(1.0 - ns)).sum()
       i_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - i_n
       I_n = (ns*(1.0 - ns))
       I_d = numpy.sqrt(abs(ns*(1.0 - ns))) / 2.0 - I_n
       I_n = numpy.linalg.multi_dot([c, numpy.diag(I_n), c.T])
       I_d = numpy.linalg.multi_dot([c, numpy.diag(I_d), c.T])
       I_n = 2.0 * numpy.dot(I_n, Da).trace()
       I_d = 2.0 * numpy.dot(I_d, Da).trace()
       cors.append([i_d, i_n, I_d, I_n])
       #
    SYSTEM_ID += 1
    return

def prepare(do_fci=True):
    "Prepare the input data: all the molecules, RHF wavefunctions and reference energies"
    mols = []
    wfns = []
    refs = []
    dips = []
    cors = []
    inp  = {'mol': mols, 'wfn': wfns, 'ref': refs, 'dip': dips, 'cor': cors}

    # 2-electron systems
    X = numpy.linspace(-0.8, 3.0, 2) 
    #X = numpy.array([-0.8, -0.6, -0.4, -0.2, -0. ,  0.2,  0.4,  0.6,  0.8,  1. ,  1.2,
    #    1.4,  1.6,  1.8,  2. ,  2.2,  2.4,  2.6,  2.8,  3. ])
    #X = [-0.4, -0.2, -0. ,  0.2,  0.4,  0.6, 0.8]
    for x in X:
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_h2_stretch, x=x)

    # 4-electron systems
    #S = numpy.linspace(2.3, 2.5, 3)
    #for s in S:
    #    put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_h4_symmetric_stretch, x=s)
    X = numpy.linspace(-0.3, 0.6, 10)
    for x in X:
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_h2h2_pull_1h, x=x)
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_h2h2_pull_2h, x=x)

    # 10-electron systems
    X = numpy.linspace(-0.3, 0.8, 10)
    #Y = X[::-1]
    for x in X:
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_1h, x=x, molname='h2o')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_1h, x=x, molname='h3o')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_1h, x=x, molname='nh3')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_1h, x=x, molname='ch4')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_nh, x=x, molname='h2o')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_nh, x=x, molname='h3o')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_nh, x=x, molname='nh3')
        put_input(do_fci, mols, wfns, refs, dips, cors, scan_func=scan_xhn_stretch_nh, x=x, molname='ch4')
    return inp

def run_dmft(wfn, f, kmax, t, g_0, verbose, conv=0.000005):
    "Run DMFT calculations for a particular wavefunction"
    if   f.lower() == 'v0': func = gefp.density.functional.XCFunctional.create('MBB'                              )
    elif f.lower() == 'v1': func = gefp.density.functional.XCFunctional.create('A1MEDI', coeff={'a0':t}, kmax=kmax)
    else:                   func = gefp.density.functional.XCFunctional.create('A2MEDI', coeff={'t' :t}, kmax=kmax)

    dmft = gefp.density.dmft.DMFT.create(wfn, func, 
                           algorithm='proj-p', 
                           guess='hcore', 
                           step_mode='search')
    if f.lower() == 'v0':
         dmft.run(maxit=300, conv=conv, g_0=g_0, verbose=verbose)
    elif f.lower() == 'v1':
         dmft.set_gradient_mode(num=True)
         dmft.run(maxit=300, conv=conv, g_0=g_0, verbose=verbose)
    elif f.lower() == 'v2':
         dmft.set_gradient_mode(approx=True)                   
         dmft.run(maxit=300, conv=0.000010, g_0=g_0, verbose=verbose)
         dmft.set_gradient_mode(exact=True)
         dmft.run(maxit=300, conv=conv, g_0=g_0, verbose=verbose, restart=True)
    else: raise ValueError("Only 'V0' - (MBB), 'V1' - (A1MEDI) and 'V2' - (A2MEDI) are available")
    return func, dmft

def evaluate_error(par, wfn, e_ref, f='v1', kmax=10, g_0=0.01, verbose=False):
    "Evaluate DMFT energy error wrt reference from input wavefunction"
    func, dmft = run_dmft(wfn, f, kmax, par[0], g_0, verbose)
    error = (dmft.E - e_ref)**2
    del dmft, func
    psi4.core.clean()
    return error

def obj(par, inp, f, kmax, g_0):
    "Objective function to minimize (to find the optimal DMFT parameter(s))"
    error = 0.0
    mols = inp['mol']
    for i, mol in enumerate(mols): 
        wfn   = inp['wfn'][i]
        ref   = inp['ref'][i]
        error += evaluate_error(par, wfn, ref, f, kmax, g_0)
    error /= len(mols)
    error = math.sqrt(error)
    error*= 627.5
    #print(" Current: ", coeffs)
    return error

def gen_surface(f, kmax, g_0, 
        tol=0.0000000000001, ftol=0.00000001, maxiter=100, disp=True, iprint=4, verbose=False,
        outfile=OUTFILE,
        errfile=ERRFILE):
    "Generate the universal surface"
    inp = prepare(do_fci=True)
    data= Data(outfile)
    err = open(errfile, 'a')

    if f.lower() == 'v1':
       coeffs = [0.50]
       bounds = [[0.3, 0.999],]
    elif f.lower() == 'v2':
       coeffs = [0.1]
       bounds = [[0.001, 0.99],]
    else: raise ValueError("Only 'V1' - (A1MEDI) and 'V2' - (A2MEDI) are available")

    options = {'disp': disp, 'iprint':iprint, 'ftol':ftol, 'maxiter':maxiter}

    for i in range(len(inp['mol'])):
        mol_i = inp['mol'][i]
        wfn_i = inp['wfn'][i]
        ref_i = inp['ref'][i]
        dip_i = inp['dip'][i]
        cor_i = inp['cor'][i]
        inp_i  = {'mol'  : [mol_i], 
                  'wfn'  : [wfn_i], 
                  'ref'  : [ref_i]}
        i_d_ref, i_n_ref, I_d_ref, I_n_ref = cor_i

        # find optimum DMFT parameter
        r_i = scipy.optimize.minimize(obj, coeffs, args=(inp_i, f, kmax, g_0), 
                   method='slsqp', 
                   tol=tol, 
                   options=options, 
                   bounds=bounds)

        t_i = r_i.x[0]
        print(" Fitted parameter: %14.6f" % t_i)

        # Run DMFT for the fitted parameter
        func_fun, dmft_fun \
                   = run_dmft(wfn_i, f, kmax, t_i, g_0, verbose, conv=1.0e-9)
        N          = int(dmft_fun.N.sum())
        E_fun      = dmft_fun.E
        D_fun      = dmft_fun.D
        d_fun      = dmft_fun.dipole
        i_d_fun, i_n_fun   = dmft_fun.scalar_correlation
        I_d_fun, I_n_fun   = dmft_fun.matrix_correlation
        I_d_fun        = numpy.dot(I_d_fun, D_fun).trace() * 2.0
        I_n_fun        = numpy.dot(I_n_fun, D_fun).trace() * 2.0
        del func_fun, dmft_fun
        psi4.core.clean()

        # Run MBB
        func_mbb, dmft_mbb \
                   = run_dmft(wfn_i, 'v0', None, None, g_0, verbose)
        N_mbb      = int(dmft_mbb.N.sum())
        if (N_mbb != N): 
            err.write("WARNING! %d-th system has invalid numbers of electrons!!! " % (i+1))
            err.write("N_fun = %d  N_mbb = %d\n" % (N, N_mbb))
        E_mbb      = dmft_mbb.E
        D_mbb      = dmft_mbb.D
        d_mbb      = dmft_mbb.dipole
        i_d_mbb, i_n_mbb   = dmft_mbb.scalar_correlation
        I_d_mbb, I_n_mbb   = dmft_mbb.matrix_correlation
        I_d_mbb        = numpy.dot(I_d_mbb, D_mbb).trace() * 2.0
        I_n_mbb        = numpy.dot(I_n_mbb, D_mbb).trace() * 2.0

        del func_mbb, dmft_mbb
        psi4.core.clean()

        # Append data
        data.add(n_alpha   = N, 
                 e_ref     = ref_i, 
                 e_fun     = E_fun, 
                 e_mbb     = E_mbb,
                 t1        = t_i, 
                 i_d_ref   = i_d_ref, i_n_ref = i_n_ref, I_d_ref = I_d_ref, I_n_ref = I_n_ref,
                 i_d_fun   = i_d_fun, i_n_fun = i_n_fun, I_d_fun = I_d_fun, I_n_fun = I_n_fun,
                 i_d_mbb   = i_d_mbb, i_n_mbb = i_n_mbb, I_d_mbb = I_d_mbb, I_n_mbb = I_n_mbb,
                 mom_d_ref = dip_i,
                 mom_d_fun = d_fun,
                 mom_d_mbb = d_mbb)

        # Save
        data.write_last()

        # Close all files
        err.close()
    return

def gen_surface_mbb(g_0, verbose, outfile=OUTFILE):
    "Generate the data for MBB functional only"
    inp = prepare(do_fci=False)
    data= Data(outfile)

    for i in range(len(inp['mol'])):
        wfn_i = inp['wfn'][i]

        func, dmft = run_dmft(wfn_i, 'v0', None, None, g_0, verbose)
        E          = dmft.E
        D          = dmft.D
        N          = int(wfn_i.nalpha())
        d          = dmft.dipole
        i_d, i_n   = dmft.scalar_correlation
        I_d, I_n   = dmft.matrix_correlation
        I_d        = numpy.dot(I_d, D).trace() * 2.0
        I_n        = numpy.dot(I_n, D).trace() * 2.0
        del func, dmft
        psi4.core.clean()

        data.add(n_alpha   = N, 
                 e_mbb     = E, 
                 i_d_mbb   = i_d, i_n_mbb = i_n, I_d_mbb = I_d, I_n_mbb = I_n,
                 mom_d_mbb = d)

        data.write_last()
    return
