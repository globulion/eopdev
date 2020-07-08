#!/usr/bin/python3
import gefp
mol = gefp.core.utilities.psi_molecule_from_file('reference.psi')
xyz_1 = mol.extract_subsets(1).geometry().to_array(dense=True)
xyz_2 = mol.extract_subsets(2).geometry().to_array(dense=True)

s = gefp.math.matrix.Superimposer()
s.set(xyz_2, xyz_1) # superimposition: mol_1 --> mol_2
s.run()
xyz = s.get_transformed()
rms = s.get_rms()
r,t = s.get_rotran()
print(rms)
print(xyz)
print(r)
print(t)
