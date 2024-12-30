import copy
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, NeighborSearch, Selection


def protein_cut(structure, res_list_all, min_distance_res, protein_name):
    if protein_name == 'A1':
        del_list = [i for i in res_list_all if i > min_distance_res]
    elif protein_name == 'A2':
        del_list = [i for i in res_list_all if i < min_distance_res]
    elif protein_name == 'B1':
        del_list = [i for i in res_list_all if i >= min_distance_res]
    else:
        del_list = [i for i in res_list_all if i <= min_distance_res]

    print("Residues of protein "+protein_name+" removed:", del_list)
    for rid in del_list:
        res = (' ', rid, ' ')
        structure[0]['A'].detach_child(res)
    rest_res = [i for i in res_list_all if i not in del_list]

    return rest_res


def get_residues_list(structure):
    residues_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    residues_list.append(residue.get_id()[1])
    return residues_list


def get_ca_atoms(structure, residues_list):
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if (residue.get_id()[1] in residues_list) and ('CA' in residue):
                    ca_atoms.append(residue['CA'])
    return ca_atoms


def get_atom_matrix(structure, selected_residues):
    atom_coords = []
    alpha_carbon_coords = []

    # Iterate through each model in the PDB structure
    for model in structure:
        # Iterate through each chain in the model
        for chain in model:
            # Iterate through each residue in the chain
            for residue in chain:
                if residue.get_id()[1] in selected_residues:
                    # Find the alpha carbon atom of the selected residue
                    alpha_carbon = residue['CA']
                    alpha_carbon_coords.append(alpha_carbon.get_coord())
                # Iterate through each atom in the residue
                for atom in residue:
                    atom_coord = atom.get_coord()
                    atom_coords.append(atom_coord)
    return atom_coords, alpha_carbon_coords


def get_centered_matrix(atom_matrix):
    N = atom_matrix.shape[0]
    centered_coord = np.mean(atom_matrix, axis=0)
    # Create a center matrix by subtracting the mean from each data point
    centered_matrix = atom_matrix - np.tile(centered_coord, (N, 1))
    return centered_matrix, centered_coord


def get_seq(structure):
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    sequence = ""

    for residue in structure.get_residues():
        # Get the three-letter abbreviation of the residue
        residue_name = residue.get_resname()
        # Get the corresponding one-letter abbreviation, use 'X' if there is no corresponding one-letter code
        one_letter_code = three_to_one.get(residue_name, 'X')
        # Add the one-letter code to the sequence string
        sequence += one_letter_code

    return sequence


def get_clash_num(structure1, structure2, cutoff=1.5):
    clash_num = 0
    # Get lists of all atoms
    atoms_structure1 = Selection.unfold_entities(structure1, 'A')  # 'A' represents atoms
    atoms_structure2 = Selection.unfold_entities(structure2, 'A')

    # Create NeighborSearch object
    ns = NeighborSearch(atoms_structure1)

    # Check each atom in structure 2 to see if they have atoms in structure 1 within 1.5 Å
    for atom in atoms_structure2:
        close_atoms = ns.search(atom.get_coord(), cutoff, level='A')
        if close_atoms:
            clash_num += len(close_atoms)
    return clash_num


def get_intrastructure_clash_num(structure, cutoff=1.5):
    clash_num = 0
    # Get lists of all atoms
    atoms = Selection.unfold_entities(structure, 'A')  # 'A' represents atoms

    # Filter out hydrogen atoms and main chain atoms (N, CA, C, O)
    atoms = [atom for atom in atoms if atom.element != 'H' and atom.get_name() not in ['N', 'CA', 'C', 'O']]

    # Create NeighborSearch object
    ns = NeighborSearch(atoms)

    for atom in atoms:
        close_atoms = ns.search(atom.get_coord(), cutoff, level='A')
        # Remove self-comparison results
        close_atoms = [a for a in close_atoms if a != atom]

        # Further filter to exclude atoms from the same residue
        close_atoms = [a for a in close_atoms if a.get_parent() != atom.get_parent()]

        if close_atoms:
            clash_num += len(close_atoms)

    return clash_num // 2


def get_residue_clash_num(structure, cutoff=1.5, clash_threshold=0):
    # Only include side chain atoms, excluding main chain atoms N, CA, C, and O
    atoms = [atom for atom in structure.get_atoms() if not atom.element == 'H' and atom.get_name() not in ['N', 'CA', 'C', 'O']]
    ns = NeighborSearch(atoms)

    clash_res_id = []
    for model in structure:
        for chain in model:
            for residue in chain:
                conflict_atoms = 0
                # Only check side chain atoms for each residue
                sidechain_atoms = [atom for atom in residue if atom.get_name() not in ['N', 'CA', 'C', 'O']]
                for atom in sidechain_atoms:
                    neighbors = ns.search(atom.coord, cutoff)
                    # Exclude atoms within the same residue and main chain atoms
                    external_neighbors = [neigh for neigh in neighbors if not neigh.parent.id == atom.parent.id and neigh.get_name() not in ['N', 'CA', 'C', 'O']]
                    conflict_atoms += len(external_neighbors)
                # Use the clash_threshold to determine if there is a clash
                if conflict_atoms > clash_threshold:
                    clash_res_id.append(residue.get_id()[1])
    return clash_res_id


def residue_clash_check(pfile, key_residue_index_cat, cutoff=1.5, clash_threshold=0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pfile)
    clash_res_id = get_residue_clash_num(structure, cutoff=1.5, clash_threshold=0)
    print(f"Conflicting residues in {pfile}:", clash_res_id, '\n')
    key_clash_res_id = [i for i in clash_res_id if i in key_residue_index_cat]
    if len(key_clash_res_id) != 0:
        print(f"Conflicting specified key residues in {pfile}:", key_clash_res_id, '\n')
    return clash_res_id


def calculate_ca_parameters(structure1, structure2, residues_list1, residues_list2, residues_list_all1, residues_list_all2, ca_distance_cutoff=0.5):
    ca_atoms1 = get_ca_atoms(structure1, residues_list1)
    ca_atoms2 = get_ca_atoms(structure2, residues_list2)

    # Calculate the distance
    distances = []
    supatom_num = 0
    for i1, atom1 in enumerate(ca_atoms1):
        for i2, atom2 in enumerate(ca_atoms2):
            distance = atom1 - atom2
            distances.append(distance)
            if distance < ca_distance_cutoff:
                supatom_num += 1

    ca_atoms1_all = get_ca_atoms(structure1, residues_list_all1)
    ca_atoms2_all = get_ca_atoms(structure2, residues_list_all2)
    supatom_num_all = 0
    for i1, atom1 in enumerate(ca_atoms1_all):
        for i2, atom2 in enumerate(ca_atoms2_all):
            distance = atom1 - atom2
            if distance < ca_distance_cutoff:
                supatom_num_all += 1

    # Extract the minimum distance and its corresponding index
    min_distance = min(distances)
    min_distance_index = distances.index(min_distance)
    i1 = min_distance_index // len(ca_atoms2)
    i2 = min_distance_index % len(ca_atoms2)

    min_residue1 = ca_atoms1[i1].get_parent().get_resname()
    min_residue1_index = residues_list1[i1]
    min_residue2 = ca_atoms2[i2].get_parent().get_resname()
    min_residue2_index = residues_list2[i2]

    rmsd_forward = calculate_rmsd(ca_atoms1, ca_atoms2)
    rmsd_reverse = calculate_rmsd(ca_atoms1, ca_atoms2[::-1])
    ca_rmsd = min(rmsd_forward, rmsd_reverse)

    return supatom_num, supatom_num_all, min_distance, min_residue1, min_residue2, min_residue1_index, min_residue2_index, ca_rmsd


def calculate_rmsd(atoms1, atoms2):
    distances = np.array([atom1 - atom2 for atom1, atom2 in zip(atoms1, atoms2)])
    rmsd = np.sqrt(np.mean(np.square(distances)))
    return rmsd


def get_new_atom_matrix(pdb_file1, pdb_file2, residue_list1, residue_list2, align_num=5):
    new_structure, supatom_num_list, supatom_num_all_list, min_distance_list, min_residue1_list, \
        min_residue2_list, min_residue1_index_list, min_residue2_index_list, \
        ca_rmsd_list = [], [], [], [], [], [], [], [], []

    parser = PDBParser(QUIET=True)

    # Define the list of residue indices to select (in this example, residues 10 to 20 are selected)
    selected_residues1 = list(range(residue_list1[0], residue_list1[1]))
    selected_residues2 = list(range(residue_list2[0], residue_list2[1]))

    structure1 = parser.get_structure('protein', pdb_file1)
    structure2 = parser.get_structure('protein', pdb_file2)

    res_list_all1 = get_residues_list(structure1)
    res_list_all2 = get_residues_list(structure2)

    # Check if align_num is greater than the number of residues, if so, raise an error
    assert align_num <= len(selected_residues1), "align_num should be smaller than the number of residues"
    assert align_num <= len(selected_residues2), "align_num should be smaller than the number of residues"

    if align_num < len(selected_residues1):
        selected_residues1_new = []
        for i in range(len(selected_residues1) - align_num + 1):
            selected_residues1_new.append(selected_residues1[i:i + align_num])
        selected_residues1_use = selected_residues1_new
    else:
        selected_residues1_use = [selected_residues1]

    if align_num < len(selected_residues2):
        selected_residues2_new = []
        for i in range(len(selected_residues2) - align_num + 1):
            selected_residues2_new.append(selected_residues2[i:i + align_num])
        selected_residues2_use = selected_residues2_new
    else:
        selected_residues2_use = [selected_residues2]

    for i1 in range(len(selected_residues1_use)):
        for i2 in range(len(selected_residues2_use)):
            res_list1 = selected_residues1_use[i1]
            res_list2 = selected_residues2_use[i2]

            acl1, accl1 = get_atom_matrix(structure1, res_list1)
            acl2, accl2 = get_atom_matrix(structure2, res_list2)
            atom_matrix1, atom_matrix2 = np.array(acl1), np.array(acl2)
            alpha_carbon_matrix1, alpha_carbon_matrix2 = np.array(accl1), np.array(accl2)

            alpha_carbon_matrix1_centered, centered_coord1 = get_centered_matrix(alpha_carbon_matrix1)
            alpha_carbon_matrix2_centered, centered_coord2 = get_centered_matrix(alpha_carbon_matrix2)

            CM = np.dot(alpha_carbon_matrix1_centered.T, alpha_carbon_matrix2_centered)

            U, S, VT = np.linalg.svd(CM)

            R = np.dot(VT.T, U.T)
            T = centered_coord2 - np.dot(R, centered_coord1)

            atom_matrix_trans1 = np.dot(atom_matrix1, R.T) + T

            structure_copy = copy.deepcopy(structure1)
            for atom, new_coord in zip(structure_copy.get_atoms(), atom_matrix_trans1):
                atom.set_coord(new_coord)
            new_structure.append(structure_copy)

            supatom_num, supatom_num_all, min_distance, min_residue1, min_residue2, min_residue1_index, min_residue2_index, ca_rmsd = \
                calculate_ca_parameters(structure_copy, structure2, res_list1, res_list2, selected_residues1, selected_residues2)

            supatom_num_list.append(supatom_num)
            supatom_num_all_list.append(supatom_num_all)
            min_distance_list.append(min_distance)
            min_residue1_list.append(min_residue1)
            min_residue2_list.append(min_residue2)
            min_residue1_index_list.append(min_residue1_index)
            min_residue2_index_list.append(min_residue2_index)
            ca_rmsd_list.append(ca_rmsd)
    return new_structure, structure2, supatom_num_list, supatom_num_all_list, min_distance_list, min_residue1_list, \
        min_residue2_list, min_residue1_index_list, min_residue2_index_list, ca_rmsd_list, res_list_all1, res_list_all2