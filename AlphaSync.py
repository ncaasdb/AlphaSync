import os
import pandas as pd
import argparse
import datetime
from pdb_tools import *
from protein_process import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein_A", type=str, default="protein_example1.pdb",
                        help="path to the first protein")
    parser.add_argument("--protein_B", type=str, default="protein_example2.pdb",
                        help="path to the second protein")
    parser.add_argument('--align_ResNo_list', type=str, default="7, 8, 9",
                        help="list of the number of residues to align")
    parser.add_argument('--protein_align_range_A', type=str, default="5, 15",
                        help="range of residues to align in the first protein")
    parser.add_argument('--protein_align_range_B', type=str, default="40, 50",
                        help="range of residues to align in the second protein")
    parser.add_argument('--fixed_residue_list_A', type=str, default="20, 21, 22, 23, 24, 25",
                        help="list of key residues in the first protein")
    parser.add_argument('--fixed_residue_list_B', type=str, default="15, 16, 17, 18, 19, 20",
                        help="list of key residues in the second protein")
    parser.add_argument('--save_path', type=str, default="AlphaSync",
                        help="path to save the superimposition results")
    opt = parser.parse_args()

    print(opt)

    folder_name = opt.save_path + '_' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print(f"Creating a new HybridHelix Synthetica task, files will be saved in {folder_name}")

    sup_path = folder_name
    os.makedirs(sup_path, exist_ok=True)

    opt_align_num_list = [int(i) for i in opt.align_ResNo_list.split(',')]
    opt_protein_align_range_A = [int(i) for i in opt.protein_align_range_A.split(',')]
    opt_protein_align_range_B = [int(i) for i in opt.protein_align_range_B.split(',')]
    opt_fixed_residue_list_A = [int(i) for i in opt.fixed_residue_list_A.split(',')]
    opt_fixed_residue_list_B = [int(i) for i in opt.fixed_residue_list_B.split(',')]

    protein_file1 = opt.protein_A
    protein_file2 = opt.protein_B

    align_num_list = opt_align_num_list

    protein_align_range1 = [opt_protein_align_range_A[0], opt_protein_align_range_A[1] + 1]
    protein_align_range2 = [opt_protein_align_range_B[0], opt_protein_align_range_B[1] + 1]
    key_residue_list1 = opt_fixed_residue_list_A
    key_residue_list2 = opt_fixed_residue_list_B

    deep_structure_list = []
    protein_save(protein_file1, protein_file2, sup_path, deep_structure_list)

    file_index_output, supatom_num_list_output, clash_num_list_output, supatom_num_all_list_output, min_distance_list_output, min_residue1_list_output, \
        min_residue2_list_output, min_residue1_index_list_output, min_residue2_index_list_output, ca_rmsd_list_output, seq_list_output, seq_len_list_output, sup_clash_res_list_output = \
        [], [], [], [], [], [], [], [], [], [], [], [], []
    key_residue_dict, key_residue_dict_cores = {}, {}
    for an in align_num_list:
        new_structure_list, structure2, supatom_num_list, supatom_num_all_list, min_distance_list, min_residue1_list,\
        min_residue2_list, min_residue1_index_list, min_residue2_index_list, ca_rmsd_list, \
        res_list_all1, res_list_all2 = \
            get_new_atom_matrix(protein_file1, protein_file2, protein_align_range1, protein_align_range2, align_num=an)

        new_structure_list_copy = [copy.deepcopy(i) for i in new_structure_list]
        deep_structure_list += new_structure_list_copy

        min_distance_best_index = min_distance_list.index(min(min_distance_list))

        for i, new_structure in enumerate(new_structure_list):
            min_distance_res1 = min_residue1_index_list[i]
            min_distance_res2 = min_residue2_index_list[i]

            save_list = []

            print(f"Running: number of atoms for superimposition is {an}, {i + 1}th iteration")
            if min_distance_res1 - res_list_all1[0] > res_list_all1[-1] - min_distance_res1:
                rest_res1 = protein_cut(new_structure, res_list_all1, min_distance_res1, 'A1')
            else:
                rest_res1 = protein_cut(new_structure, res_list_all1, min_distance_res1, 'A2')

            structure2_copy = copy.deepcopy(structure2)
            if min_distance_res2 - res_list_all2[0] > res_list_all2[-1] - min_distance_res2:
                rest_res2 = protein_cut(structure2_copy, res_list_all2, min_distance_res2, 'B1')
                save_list.append(structure2_copy)
                save_list.append(new_structure)
                key_residue_index_cat, res_dict_cores = res_ind_cat(key_residue_list2, key_residue_list1, rest_res2, rest_res1, 'B1')
            else:
                rest_res2 = protein_cut(structure2_copy, res_list_all2, min_distance_res2, 'B2')
                save_list.append(new_structure)
                save_list.append(structure2_copy)
                key_residue_index_cat, res_dict_cores = res_ind_cat(key_residue_list1, key_residue_list2, rest_res1, rest_res2, 'B2')

            key_residue_dict[str(an) + '_' + str(i)] = key_residue_index_cat
            key_residue_dict_cores[str(an) + '_' + str(i)] = res_dict_cores

            io = PDBIO()
            io.set_structure(save_list[0])
            io.save(f"{sup_path}/sup_ca{an}_{i}.pdb", write_end=False)

            with open(f"{sup_path}/sup_ca{an}_{i}.pdb", 'a') as file_handle:
                io.set_structure(save_list[1])
                io.save(file_handle, write_end=True)
            ter_del(f"{sup_path}/sup_ca{an}_{i}.pdb")
            renumber_residues(f"{sup_path}/sup_ca{an}_{i}.pdb", start_number=1)
            renumber_atoms(f"{sup_path}/sup_ca{an}_{i}.pdb", start_number=1)

            parser_new = PDBParser(QUIET=True)
            new_sup_structure = parser_new.get_structure('protein', f"{sup_path}/sup_ca{an}_{i}.pdb")
            clash_num_now = get_intrastructure_clash_num(new_sup_structure, cutoff=1.5)
            clash_num_list_output.append(clash_num_now)

            seq_now = get_seq(new_sup_structure)
            seq_list_output.append(seq_now)

            seq_len_list_output.append(len(seq_now))

            print(f"New coordinates have been saved to {sup_path}/sup_ca{an}_{i}.pdb" + "\n")
            sup_clash_res_list = residue_clash_check(f"{sup_path}/sup_ca{an}_{i}.pdb", key_residue_index_cat, cutoff=1.5, clash_threshold=5)
            sup_clash_res_list_output.append(str(sup_clash_res_list))

            file_index_output.append(str(an) + '_' + str(i))

        supatom_num_list_output += supatom_num_list
        supatom_num_all_list_output += supatom_num_all_list
        min_distance_list_output += min_distance_list
        min_residue1_list_output += min_residue1_list
        min_residue2_list_output += min_residue2_list
        min_residue1_index_list_output += min_residue1_index_list
        min_residue2_index_list_output += min_residue2_index_list
        ca_rmsd_list_output += ca_rmsd_list

    io = PDBIO()
    io.set_structure(deep_structure_list[0])
    io.save(f"{sup_path}/saturated_superposition_structure.pdb", write_end=True)

    for structure in deep_structure_list[1:]:
        with open(f"{sup_path}/saturated_superposition_structure.pdb", 'a') as file_handle:
            io.set_structure(structure)
            io.save(file_handle, write_end=True)

    dict_output = {'file_index': file_index_output, 'clash_num': clash_num_list_output, 'supatom_num': supatom_num_list_output,
                   'supatom_num_all': supatom_num_all_list_output, 'min_distance': min_distance_list_output,
                   'min_residue1': min_residue1_list_output, 'min_residue2': min_residue2_list_output,
                   'min_residue1_index': min_residue1_index_list_output, 'min_residue2_index': min_residue2_index_list_output,
                   'ca_rmsd': ca_rmsd_list_output, 'seq': seq_list_output, 'seq_len': seq_len_list_output, 'sup_clash_res_list': sup_clash_res_list_output}
    df_output = pd.DataFrame(dict_output)
    df_output.to_csv(f'{sup_path}/output.csv', index=False)
    # with open(f'{sup_path}/key_residue_dict.pkl', 'wb') as file:
    #     pickle.dump(key_residue_dict, file)
    # with open(f'{sup_path}/key_residue_dict_cores.pkl', 'wb') as file:
    #     pickle.dump(key_residue_dict_cores, file)

    print("Superimposition completed. Please check the superimposition results in the sup_result folder and "
          "select the desired superimposition result for further optimization.")


if __name__ == '__main__':
    main()
