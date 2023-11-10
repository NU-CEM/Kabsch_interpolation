#!/usr/bin/env python3

""" 
Linear interpolation of crystal structures along translation and rotation 
vectors. Designed primarily for application to hybrid organic-inorganic
materials. See https://github.com/NU-CEM/Kabsch_interpolation for more 
information. Permission to use software under the MIT license.

Author: Lucy Whalley (l.whalley@northumbria.ac.uk)
"""

import ase
from ase import io
from ase.geometry import analysis
from ase.build import molecule
from ase import Atoms
from ase import neighborlist

import numpy as np
from numpy.linalg import norm 
from scipy.spatial.transform import Rotation as R
from scipy.spatial.transform import Slerp
from scipy import sparse

import itertools
from collections import Counter
from copy import deepcopy
import argparse

def interpolate_structures(start_atoms, end_atoms, molecular_formulas=None, 
	number_intermediates=9, fformat="vasp", reverse=True, molecular_indices=None, 
	translation_species=None, mic_cutoff=0.5):
    """ Uses linear interpolation along translation and rotation vectors to 
    create intermediate structures that lie between those in start_filepath and 
    end_filepath.
    
    Args:

    start_atoms (str) - ASE Atoms object for the start structure.

    end_atoms (str) - ASE Atoms object for the end structure.

    molecular_formulas (list(str))(optional) - a list of molecular formulas to 
    specify which atoms will be rotated and translated. If not set then the 
    molecular_indices keyword argument must be set to explicitly identify atoms 
    for Kabsch interpolation. Defaults to None.

    number_intermediates (int)(optional) - number of intermediate structures 
    generated. Defaults to 9.

    fformat (string)(optional) - file format for writing intermediate structures. 
    For output options see https://wiki.fysik.dtu.dk/ase/ase/io/io.html. Defaults 
    to "vasp".

    reverse (bool)(optional) - create additional negative amplitude 
    interpolations along vectors. Defaults to True.

    molecular_indices (list(arrays))(optional) - list of numpy arrays. Each array 
    contains the atom indices for a molecule. If not set, molecular_indices will be 
    found automatically using the molecular_formulas keyword. Defaults to None.

    translation_species (list(str))(optional) - a list of elemental species which 
    will be translated (without rotation). If not set then all 
    non-molecular species will be translated. Defaults to None.
    
    mic_cutoff (float) - cutoff distance from edge of unit cell, below which the 
    minimum image convention is applied to any molecule."""
        
    # need to specify either molecular_formulas (most common use case) or 
    # molecular_indices (for awkward cases where ASE neighbour analysis doesn't work)
    if molecular_formulas is None and molecular_indices is None:
        raise ValueError("""either molecular_formulas or molecular_indices must 
        	be specified""")
        
    # if single entry for molecular formula, convert it to a list.
    if type(molecular_formulas) is str:
        molecular_formulas = [molecular_formulas]
        
    # need to ensure minimum image convention between start and end structures
    start_atoms, end_atoms = start_end_mic(start_atoms, end_atoms)
    
    # get index of every atom that belongs to a rotating molecule.
    if molecular_indices is not None:
        molecular_indices_list = molecular_indices
    else:
        molecular_indices_list = find_molecules(start_atoms, molecular_formulas)
    
    # get index of every atom that will be translated
    translation_indices = get_translation_indices(start_atoms, 
    	translation_species, molecular_indices_list)

    if reverse:
        iterator = range(-(number_intermediates+1),number_intermediates+2)
    else:
        iterator = range(0,number_intermediates+2)
    
    # generate structure at each interpolation step
    for step_index in iterator:

    	# deepcopies to avoid overwriting
    	# honestly, I don't understand how this works and atoms can't
    	# just be updated directly (without creating positions)
    	# it was a case of trial and error
        positions = deepcopy(start_atoms.get_positions())            
        atoms = deepcopy(start_atoms)   
        
        # "amplitude" of interpolation
        interval = step_index * (1/(number_intermediates+1))   
        
        # interpolate along translation vector 
        for atom_index in translation_indices:
            atom_translation = (end_atoms[atom_index].position - 
            start_atoms[atom_index].position)
            beta = atom_translation*interval
            atoms[atom_index].position = start_atoms[atom_index].position + beta
        
        # update positions with the interpolated positions
        positions = atoms.positions
        
        # interpolation along translation and rotation vectors
        for molecule_indices in molecular_indices_list: 
            
            # ASE atoms object to describe particular molecule
            start_molecule = start_atoms[molecule_indices] 
            end_molecule = end_atoms[molecule_indices] 
            
            # need to apply minimum image convention to any molecule that may 
            # move between neighbouring unit cells during relaxation
            if (((np.abs(start_molecule.positions) < mic_cutoff).any()) or 
            ((np.abs(end_molecule.positions) < mic_cutoff).any())):                           
                start_molecule.set_positions(ase.geometry.geometry.find_mic(
                	start_molecule.positions, start_atoms.cell)[0])
                end_molecule.set_positions(ase.geometry.geometry.find_mic(
                	end_molecule.positions, start_atoms.cell)[0])

            # get set of vectors that describe the molecule
            start_vectors = get_molecule_vectors(start_molecule)
            end_vectors = get_molecule_vectors(end_molecule)
            
            # get rotation and translation vectors 
            axis, angle = get_axis_angle(start_vectors, end_vectors)
            translation = get_translation(start_molecule, end_molecule)
            
            # scale vectors by the "amplitude" of the interpolation
            delta = translation * interval 
            alpha = angle * interval
        
            # apply the translation and rotation
            # deepcopy to avoid overwriting start_atoms
            molecule = deepcopy(start_molecule)
            molecule.rotate(alpha, axis, center="COM")
            molecule.translate(delta)
            
            # update atoms with the interpolated positions
            positions[molecule_indices] = molecule.positions
          
        # update atoms with interpolated positions
        atoms.set_positions(positions)
        # write out the interpolated structure      
        ase.io.write('POSCAR_'+str(step_index).zfill(3)+"."+fformat,atoms, 
        	format=fformat)

def find_molecules(atoms, molecular_formulas):
    """ Returns a list of arrays. 
    Each array contains the index of atoms in a molecule."""

    # create a matrix summarising the connectivity of the structure
    # this follows the example in the ASE documentation: 
    # https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html#ase.neighborlist.get_connectivity_matrix
    cutOff = ase.neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    # matrix is a scipy sparse matrix as most atoms are not connected to one another.
    # Also allows nice analysis of the connected components
    matrix = neighborList.get_connectivity_matrix()
    # calc number of components (connected atoms) and which atom belongs to which component
    n_components, component_list = sparse.csgraph.connected_components(
       matrix)

    # create list of lists. 
    # Each sub-list contains the indices for a single component (connected atoms).
    # will contain repeat components and components which are not our target molecules.
    molIdxs_list = []
    for idx in range(len(component_list)):
        molIdx = component_list[idx]
        molIdxs = [i for i in range(len(component_list))
               if component_list[i] == molIdx]
        molIdxs_list.append(molIdxs)

    # filter out repeat entries
    molIdxs_list.sort()
    molIdxs_list = list(molIdxs_list for molIdxs_list,_ 
        in itertools.groupby(molIdxs_list))

    # standardise the input molecular formulas so that formatted ASE-style.
    # this is for comparison against those we will find in the start_atoms object.
    molecular_formulas = [ase.Atoms(formula_string).get_chemical_formula() 
        for formula_string in molecular_formulas]

    # filter out molecules which are not the ones we want to rotate
    molIdxs_list = [molIdxs for molIdxs in molIdxs_list if 
        atoms[molIdxs].get_chemical_formula() in molecular_formulas]
    
    # report what has been found
    mol_counter = Counter(atoms[x].get_chemical_formula() for x in molIdxs_list)
    for key in mol_counter:
        print("{} {} molecules have been found".format(mol_counter[key],key))
    print("The molecules have the following indices:")
    for molIdxs in molIdxs_list:
        print(molIdxs)
        
    return [np.array(molIdxs) for molIdxs in molIdxs_list]
        
def get_axis_angle(start_vectors, end_vectors):
    """ Uses Kabsch algorithm to calculate the rotation axis and angle 
    between two sets of vectors"""
  
    transform = R.align_vectors(start_vectors, end_vectors)[0]
    transform.as_rotvec()
    angle = np.degrees(norm(transform.as_rotvec()))
    unit_axis = transform.as_rotvec()/norm(transform.as_rotvec())
    
    return unit_axis, -angle

def get_translation(start_molecule, end_molecule):
    """ Returns the displacement (in Angstrom) between the molecule 
    COM in the start position and end position"""
        
    start_COM = start_molecule.get_center_of_mass()
    end_COM = end_molecule.get_center_of_mass()
    translation = end_COM - start_COM
    
    return translation

def get_molecule_vectors(molecule_atoms):
    """Returns the distance vectors for all connected atoms in a molecule"""

    cutOff = ase.neighborlist.natural_cutoffs(molecule_atoms)
    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=False)
    neighborList.update(molecule_atoms)
    matrix = neighborList.get_connectivity_matrix()
    rows, columns = matrix.nonzero()
    pair_indices = np.column_stack((rows,columns))
    pair_indices = np.sort(pair_indices)
    pair_indices = pair_indices[np.lexsort([pair_indices[:, 1], pair_indices[:, 0]])]

    return [molecule_atoms.get_distance(i,j,mic=True,vector=True) for i,j in pair_indices]     

def get_translation_indices(atoms, translation_species, molecular_indices_list):
    """Returns the indices of all atoms that are to be translated (only, with 
    no rotational interpolation)."""
    
    # if translation_species is not specified then translate every atom not 
    # in molecular_indices_list
    if translation_species is None:
        flat_list = [item for sublist in molecular_indices_list for item in sublist]
        translation_indices = [atom.index for atom in atoms 
        	if atom.index not in flat_list]
    # if translation_species is specified then return the indices for every 
    # element in translation_species
    else:
        translation_indices = [atom.index for atom in atoms 
        	if atom.symbol in translation_species]
    print("{} translation-only atoms have been found"
    	.format(len(translation_indices)))
    return translation_indices

def start_end_mic(start_atoms,end_atoms):
    """In some cases an atom crosses a cell boundary during relaxation between 
    the start and end structures. Shift positions to ensure minimum image convention."""
    start_pos = start_atoms.get_scaled_positions()
    end_pos = end_atoms.get_scaled_positions()
    
    forward_cross = start_pos - end_pos > 0.5
    back_cross = start_pos - end_pos < -0.5
    
    new_end_pos = end_pos.copy()
    new_end_pos[forward_cross] += 1
    new_end_pos[back_cross] -= 1

    new_end_atoms = end_atoms.copy()
    new_end_atoms.set_scaled_positions(new_end_pos)
    return start_atoms, new_end_atoms

def main():
    parser = argparse.ArgumentParser(
                    prog = 'Kabsch_Interpolation',
                    description = """Linear interpolation of crystal structures 
                    along translation and rotation vectors""")

    parser.add_argument('start_POSCAR',
    	help='path to POSCAR for starting structure')
    parser.add_argument('end_POSCAR',
    	help='path to POSCAR for ending structure')
    parser.add_argument('-m', '--molecular_formulas',
		help='list of molecules to rotate and translate')
    parser.add_argument('-n', '--number_intermediates', default=9,
    	help='number of intermeditate structures to generate')
    parser.add_argument('-f', '--fformat', default="vasp",
    	help='''file format of output files. See ASE documentation for
		options: https://wiki.fysik.dtu.dk/ase/ase/io/io.html''')
    parser.add_argument('-r', '--reverse', default=True,
    	help='generate structures with negative amplitude interpolations')
    parser.add_argument('-i', '--molecular_indices',
		help='''indices of the atoms to rotate and translate.
		In most cases it will be easier to automatically find these by
		specifying molecular_formulas.''')
    parser.add_argument('-t', '--translation_species',
		help='''elements to translate without rotation.
		If not set, all atoms that do belong to a molecule are translated.''')
    parser.add_argument('-c', '--mic_cutoff', default=0.5,
		help='''cutoff for minimum image correction. If you generate unusual 
		structures it may be that this needs to be adjusted''')

    args = parser.parse_args()

    start_atoms = ase.io.read(args.start_POSCAR)
    end_atoms = ase.io.read(args.end_POSCAR)

    interpolate_structures(start_atoms, end_atoms, 
		molecular_formulas=args.molecular_formulas, 
		number_intermediates=args.number_intermediates, 
		fformat=args.fformat, 
		reverse=args.reverse, 
		molecular_indices=args.molecular_indices, 
		translation_species=args.translation_species, 
		mic_cutoff=args.mic_cutoff)

if __name__ == "__main__":
    main()