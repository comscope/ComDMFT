from __future__ import print_function
try:
    from builtins import zip
except:
    pass

import sys, numpy
from pymatgen import Molecule


'''
Help funtions for extracting molecule from crystal and determine local rotational operations.
'''

# module global molecule index
imol = 0

class gmolecule(Molecule):
    '''adding equivalent indices to atoms.
    '''
    def __init__(self, species, coords, equivalent_indices=[], charge=0,
            fm_direction=None, spin_multiplicity=None,
            validate_proximity=False, site_properties=None):
        super(Molecule, self).__init__(species, coords, charge=charge,
                spin_multiplicity=spin_multiplicity,
                validate_proximity=validate_proximity,
                site_properties=site_properties)
        self.fm_direction = fm_direction
        self.equivalent_indices = equivalent_indices


def xtal_get_local_rot(symbols, scaled_positions, cell, iat, dist_cut_,
        equivalent_indices=[], locrot=None, fm_direction=None,
        Nmax=4, tol=1.e-5, log=sys.stdout):
    '''
    Get rotation operations of center atom iat in crystal.
    '''
    mol = xtal_extract_mol(symbols, scaled_positions, cell, iat,
            dist_cut_, equivalent_indices=equivalent_indices,
            locrot=locrot, fm_direction=fm_direction, Nmax=Nmax)
    return mol_get_rot_list(mol, tol=tol, log=log)


def xtal_extract_mol(symbols, scaled_positions, cell, iat, dist_cut_,
        locrot=None, fm_direction=None, Nmax=10, equivalent_indices=[]):
    '''
    Extracting molecule of center atom iat from crystal.
    atom iat is the molecule center
    '''
    center = numpy.copy(scaled_positions[iat])
    new_scaled_positions = numpy.copy(scaled_positions)
    new_scaled_positions = new_scaled_positions - center
    molecule_positions = []
    molecule_symbols = []
    pair_dist = []
    eq_indices = []
    # Get the (2*Nmax+1)x(2*Nmax+1)x(2*Nmax+1) block
    Nmax_list = numpy.arange(-Nmax, Nmax)
    for i in Nmax_list:
        for j in Nmax_list:
            for k in Nmax_list:
                for jat, sp in enumerate(new_scaled_positions):
                    n_sp = sp + numpy.array([i, j, k])
                    n_p = \
                            n_sp[0] * numpy.array(cell[0]) + \
                            n_sp[1] * numpy.array(cell[1]) + \
                            n_sp[2] * numpy.array(cell[2])
                    pair_dist.append(numpy.linalg.norm(n_p))
                    molecule_positions.append(n_p)
                    molecule_symbols.append(symbols[jat])
                    if len(equivalent_indices) > 0:
                        eq_indices.append(equivalent_indices[jat])
    # Get reasonable dist_cut
    dist = pair_dist[:]
    dist.sort()
    if dist_cut_ > dist[1]:
        dist_cut = dist_cut_
    else:
        dist_cut = dist[12] + 0.7
        print(" The default dist_cut for extracting a centered cluster" +\
                "\n for symmetry evaluation = {}".format(dist_cut))
    molecule_positions = [mpos for mpos, pdist in
            zip(molecule_positions, pair_dist) if pdist < dist_cut]
    molecule_symbols = [msym for msym, pdist in
            zip(molecule_symbols, pair_dist) if pdist < dist_cut]
    eq_indices = [eind for eind, pdist in
            zip(eq_indices, pair_dist) if pdist < dist_cut]

    from itertools import groupby
    mol_name = ''
    mol_symbols = sorted(molecule_symbols)
    for key, group in groupby(mol_symbols):
        mol_name += key + str(len(list(group)))

    # local rotation if required.
    if locrot is not None:
        for i, pos in enumerate(molecule_positions):
            molecule_positions[i] = locrot.T.dot(pos)

    print(' molecule extracted {}:'.format(mol_name))
    print(" atom   x      y       z    distance")
    for symbo, position in zip(molecule_symbols, molecule_positions):
        print(" {:3s} {:6.2f} {:6.2f} {:6.2f}  {:8.4f}".format(
                *([symbo] + position.tolist() +
                [numpy.sqrt(numpy.dot(position, position))])))

    global imol
    with open('{}_{}.xyz'.format(mol_name, imol), 'w') as f:
        print(' {}\n'.format(len(molecule_symbols)), file=f)
        for symbo, position in zip(molecule_symbols, molecule_positions):
            print(" {:3s} {:6.2f} {:6.2f} {:6.2f}".format(
                    *([symbo] + position.tolist())), file=f)
    imol += 1
    return gmolecule(molecule_symbols, molecule_positions,
            equivalent_indices=eq_indices, fm_direction=fm_direction)


def xyz_get_rot_list(symbols, positions, tol=1.e-5, log=sys.stdout):
    '''
    Get the rotation element list of a molecule,
    given symboal and position list.
    '''
    mol = gmolecule(symbols, positions)
    return mol_get_rot_list(mol, tol=tol, log=log)


def mol_get_rot_list0(mol, tol=1.e-5, log=sys.stdout):
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer,\
            generate_full_symmops
    from scipy.linalg import det

    analyzer = PointGroupAnalyzer(mol)
    print(" sch_symbol = {}".format(analyzer.sch_symbol), file=log)

    # Pick rotations only.
    symmops = [o for o in analyzer.symmops
            if abs(det(o.rotation_matrix)-1)<1.e-5]
    # symmops = generate_full_symmops(symmops, tol, max_recursion_depth=50)
    symmops = generate_full_symmops(symmops, tol)
    rot_list=[o.rotation_matrix for o in symmops]
    return rot_list


def chk_rot_keep(mol, rot, tol=1.e-5):
    if len(mol.equivalent_indices) > 1:
        coords = mol.cart_coords
        for i, coord in enumerate(coords):
            coordp = numpy.dot(rot, coord)
            diff = coords - coordp
            ind = numpy.where(numpy.all(numpy.abs(diff) < tol, axis=1))[0]
            if len(ind) != 1:
                print(' WARNING: identified rotation ids: {}'.format(ind))
                return False
            if mol.equivalent_indices[i] != mol.equivalent_indices[ind[0]]:
                return False
# not correct!
    if mol.fm_direction is not None:
        vec = numpy.dot(rot, mol.fm_direction)
# ok if in the same direction.
        if not numpy.allclose(vec, mol.fm_direction, rtol=tol):
            return False
    return True


def mol_get_rot_list_screening(mol, rot_list):
    '''Further screen rotations since atoms of the same symbol might be
    not equivalent, due to magnetism for instance.
    '''
    print(' number of rotations before screening: {}'.format(len(rot_list)))
    rot_list_tmp = rot_list
    rot_list = []
    numpy.set_printoptions(precision=3, suppress=True)
    for rot in rot_list_tmp:
        if chk_rot_keep(mol, rot):
            rot_list.append(rot)
        else:
            print(' rot \n{}\n removed!'.format(rot))
    print(' number of rotations after screening: {}'.format(len(rot_list)))
    return rot_list


def mol_get_rot_list(mol, tol=1.e-5, log=sys.stdout):
    rot_list = mol_get_rot_list0(mol, tol=tol, log=log)
    if len(mol.equivalent_indices) > 0 or mol.fm_direction is not None:
        rot_list = mol_get_rot_list_screening(mol, rot_list)
    return rot_list


if __name__ == "__main__":
    '''
    A test.
    '''
    # H cluster
    symbols = ['H', 'H', 'H', 'H', 'H', 'H']
    positions = [[0., 0., 0.], [2.0, 0., 0.], [0., 2.0, 0.0],
                 [0., -2.0, 0.], [0., 0., 2.0], [0., 0., -2.0]]
    rotations = xyz_get_rot_list(symbols, positions, log=sys.stdout)
    print(" H-dimer rotations:")
    print(rotations)
    # extract molecule
    symbols = ['Ce']
    scaled_positions = [[0, 0, 0]]
    cell = [[0.5, 0.5, 0.], [0.5, 0., 0.5], [0., 0.5, 0.5]]
    iat = 0
    dist_cut = -1.0
    molecule = xtal_extract_mol(symbols, scaled_positions, cell, iat, dist_cut)
    print(molecule.species)
    print(molecule.sites)
