import copy
import re
from Bio.PDB import PDBParser, PDBIO


def ter_del(pdb_file):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
    filtered_lines = [line for line in lines if "TER" not in line]

    with open(pdb_file, 'w') as file:
        file.writelines(filtered_lines)


def protein_save(structure1, structure2, output_path, ds_list):
    parser = PDBParser(QUIET=True)
    s1 = parser.get_structure('p1', structure1)
    s2 = parser.get_structure('p2', structure2)

    ds_list.append(s2)

    io = PDBIO()
    io.set_structure(s1)
    io.save(f"{output_path}/original_structure.pdb", write_end=False)

    with open(f"{output_path}/original_structure.pdb", 'a') as file_handle:
        io.set_structure(s2)
        io.save(file_handle, write_end=True)


def renumber_residues(pdb_path, start_number=1):
    with open(pdb_path, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()

    # Define regular expression pattern to match residue numbers in ATOM or HETATM records
    pattern = r'(?<=\s)[0-9]+(?=\s+[\d\.-]+\s+[\d\.-]+\s+[\d\.-])'
    pattern_ter = r'TER\s+\d+\s+\w+\s+\w\s+(\d+)'

    # New residue number
    new_residue_number = start_number - 1

    res_now = 0
    # Iterate through each line of the PDB file
    for i, line in enumerate(pdb_lines):
        res_num = re.search(pattern, line)
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        if res_num is not None and res_num.group() != res_now:
            new_residue_number += 1

            line_res_num = re.search(pattern, line).group(0)
            pattern_use = pattern
            if len(line_res_num) > len(str(new_residue_number)):
                new_residue_number = ' ' * (len(line_res_num) - len(str(new_residue_number))) + str(new_residue_number)
            if len(line_res_num) < len(str(new_residue_number)):
                spaces = len(str(new_residue_number)) - len(line_res_num)
                pattern_use = rf'\s{{{spaces}}}[0-9]+(?=\s+[\d\.-]+\s+[\d\.-]+\s+[\d\.-])'

            pdb_lines[i] = re.sub(pattern_use, str(new_residue_number), line)
            new_residue_number = int(new_residue_number)
            res_now = res_num.group()
        elif res_num is not None and res_num.group() == res_now:
            line_res_num = re.search(pattern, line).group(0)
            pattern_use = pattern
            if len(line_res_num) > len(str(new_residue_number)):
                new_residue_number = ' ' * (len(line_res_num) - len(str(new_residue_number))) + str(new_residue_number)
            if len(line_res_num) < len(str(new_residue_number)):
                spaces = len(str(new_residue_number)) - len(line_res_num)
                pattern_use = rf'\s{{{spaces}}}[0-9]+(?=\s+[\d\.-]+\s+[\d\.-]+\s+[\d\.-])'

            pdb_lines[i] = re.sub(pattern_use, str(new_residue_number), line)
            new_residue_number = int(new_residue_number)

        else:
            continue

    # Write the modified content back to the PDB file
    with open(pdb_path, 'w') as modified_pdb_file:
        modified_pdb_file.writelines(pdb_lines)


def renumber_atoms(pdb_path, start_number=1):
    with open(pdb_path, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()

    # Define a regular expression pattern to match residue numbers in ATOM or HETATM records
    pattern = r'\b(ATOM|HETATM)\s+(\d+)(?=\s+[A-Z])'

    # New residue number
    new_residue_number = start_number - 1

    res_now = 0

    # Iterate through each line of the PDB file.
    for i, line in enumerate(pdb_lines):
        res_num = re.search(pattern, line)

        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        if res_num is not None and res_num.group() != res_now:
            new_residue_number += 1

            line_res_num = re.search(pattern, line).group(2)
            pattern_use = pattern
            if len(line_res_num) > len(str(new_residue_number)):
                new_residue_number = ' ' * (len(line_res_num) - len(str(new_residue_number))) + str(new_residue_number)
            if len(line_res_num) < len(str(new_residue_number)):
                spaces = len(str(new_residue_number)) - len(line_res_num)
                pattern_use = rf'\b(ATOM|HETATM)\s+(\s{{{spaces}}}\d+)(?=\s+[A-Z])'

            len_new_residue_number = len(str(new_residue_number))
            match = re.search(pattern_use, line).group(0)
            sub_match = match.replace(match[-len_new_residue_number:], str(new_residue_number))
            pdb_lines[i] = re.sub(pattern_use, sub_match, line)
            new_residue_number = int(new_residue_number)
            res_now = res_num.group()
        elif res_num is not None and res_num.group() == res_now:
            line_res_num = re.search(pattern, line).group(2)
            pattern_use = pattern
            if len(line_res_num) > len(str(new_residue_number)):
                new_residue_number = ' ' * (len(line_res_num) - len(str(new_residue_number))) + str(new_residue_number)
            if len(line_res_num) < len(str(new_residue_number)):
                spaces = len(str(new_residue_number)) - len(line_res_num)
                pattern_use = rf'\b(ATOM|HETATM)\s+(\s{{{spaces}}}\d+)(?=\s+[A-Z])'

            len_new_residue_number = len(str(new_residue_number))
            match = re.search(pattern_use, line).group(0)
            sub_match = match.replace(match[-len_new_residue_number:], str(new_residue_number))
            pdb_lines[i] = re.sub(pattern_use, sub_match, line)
            new_residue_number = int(new_residue_number)

        else:
            continue

    # Write the modified content back to the PDB file
    with open(pdb_path, 'w') as modified_pdb_file:
        modified_pdb_file.writelines(pdb_lines)


def res_ind_cat(krl1, krl2, rr1, rr2, print_mode):
    key_residue_index1, key_residue_index2 = [], []
    res_dict_cores = {}
    cat_residue_list = list(range(len(rr1) + len(rr2)))

    if print_mode == 'B1':
        for i in krl2:
            if i not in rr2:
                print(f"Note: Key residue {i} of protein A has been removed!")
            else:
                key_residue_index2.append(rr2.index(i))
        for i in krl1:
            if i not in rr1:
                print(f"Note: Key residue {i} of protein B has been removed!")
            else:
                key_residue_index1.append(rr1.index(i))
        res_dict_cores['A_ori'] = key_residue_index2
        res_dict_cores['B_ori'] = key_residue_index1
        res_dict_cores['A_new'] = [(cat_residue_list.index(i + len(rr1)) + 1) for i in key_residue_index2]
        res_dict_cores['B_new'] = [(cat_residue_list.index(i) + 1) for i in key_residue_index1]
    else:
        for i in krl1:
            if i not in rr1:
                print(f"Note: Key residue {i} of protein A has been removed!")
            else:
                key_residue_index1.append(rr1.index(i))
        for i in krl2:
            if i not in rr2:
                print(f"Note: Key residue {i} of protein B has been removed!")
            else:
                key_residue_index2.append(rr2.index(i))
        res_dict_cores['A_ori'] = key_residue_index1
        res_dict_cores['B_ori'] = key_residue_index2
        res_dict_cores['A_new'] = [(cat_residue_list.index(i) + 1) for i in key_residue_index1]
        res_dict_cores['B_new'] = [(cat_residue_list.index(i + len(rr1)) + 1) for i in key_residue_index2]

    key_residue_index_cat = [(cat_residue_list.index(i) + 1) for i in key_residue_index1] + \
                            [(cat_residue_list.index(i + len(rr1)) + 1) for i in key_residue_index2]
    return key_residue_index_cat, res_dict_cores