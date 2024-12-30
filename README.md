# Hybrid helix synthetica 

A tool for protein fusion.

## Requirements

* Python==3.7
* pandas==2.0.3
* numpy==1.25.1
* biopython==1.81

## Usage

Clone the repository locally, execute the following command.

```
python Hybrid_helix_synthetica.py \
    --protein_A "protein_example1.pdb" \
    --protein_B "protein_example2.pdb" \
    --align_num_list "7, 8, 9" \
    --protein_align_range_A "5, 15" \
    --protein_align_range_B "40, 50" \
    --save_path "HybridHelix_Synthetica"
```

## Flags

- **`--protein_A`** : Pathway of Protein A for superposition.
- **`--protein_B`** : Pathway of Protein B for superposition.
- **`--align_num_list`** : The number of atoms used in the superposition. (The default is 7, 8, 9, i.e., 7, 8, and 9 atoms are used for superposition respectively)
- **`--protein_align_range_A`** : The residue range of protein A to be superimposed. (The default parameter is 5, 15, i.e. residues 5-15 will be superimposed. Note that only two numbers are allowed here to define the start and end residues)
- **`--protein_align_range_B`** : The residue range of protein B to be superimposed.
- **`--save_path`** : Path for saving the superposition results.

## Output

Output sample is in AlphaSync_2024-12-30_20-02-25.

In outputs: <br>
`original_structure.pdb` is original structures of protein A and protein B; <br>
`saturated_superposition_structure.pdb` is a stack of fusion structures; <br>
`sup_ca7_0.pdb` is a fusion structure, ca7 means using 7 CA atoms for alignment, and 0 means using 7 CA atoms for alignment at the first position in the alignment range; <br>
`output.csv` records relevant indicators of the fusion structure.

In `output.csv`: <br>
`file_index` is the index of the fused structure; <br>
`clash_num` is the number of conflicts between side chains within the fusion structure, with a cutoff radius of 1.5 angstroms; <br>
`supatom_num` is the number of CA atoms/residues used for alignment; <br>
`supatom_num_all` is the number of residues successfully aligned within the selected alignment range (CA atom distance between two residues is less than 0.5 angstroms); <br>
`min_distance` is the minimum CA atomic distance between the aligned residues; <br>
`min_residue1` is the name of the residue that forms the smallest inter-residue distance in protein A; <br>
`min_residue2` is the name of the residue that forms the smallest inter-residue distance in protein B; <br>
`min_residue1_index` is the index of the residue that forms the smallest inter-residue distance in protein A; <br>
`min_residue2_index` is the index of the residue that forms the smallest inter-residue distance in protein B; <br>
`ca_rmsd` is the RMSD between the CA atoms used when fusion structures are superimposed; <br>
`seq` is the sequence of the fusion structure; <br>
`seq_len` is the sequence length of the fusion structure; <br>
`sup_clash_res_list` is a list of residue indices with atomic conflicts in the fusion structure.

```:test.csv
file_index,clash_num,supatom_num,supatom_num_all,min_distance,min_residue1,min_residue2,min_residue1_index,min_residue2_index,ca_rmsd,seq,seq_len,sup_clash_res_list
7_0,8,7,8,0.054805804,MET,LYS,10,45,0.13137935,TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEGVNALMRAAHEIRWLPNLTFDQRVAFIHKLEDDPSQSSELLSEAKKLNDSQAPK,93,"[14, 23, 27, 28, 31, 44, 52, 66, 67, 70]"
7_1,26,7,9,0.07894572,TRP,LYS,9,45,0.16165167,TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEGVNALWMRAAHEIRWLPNLTFDQRVAFIHKLEDDPSQSSELLSEAKKLNDSQAPK,94,"[14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 28, 29, 30, 31, 32, 33, 52, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 80, 81, 83]"
7_2,6,7,9,0.06367334,LYS,LEU,7,44,0.12121701,TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEGVNAKEWMRAAHEIRWLPNLTFDQRVAFIHKLEDDPSQSSELLSEAKKLNDSQAPK,95,"[14, 17, 21, 23, 28, 31, 44, 45, 48, 52, 79, 82]"
```

## Note

Confirm that the input protein pdb file format is consistent with the sample pdb file, otherwise it may cause some errors.

Ensure that the ranges of protein A and protein B used for alignment are located at different termini, i.e. (A-N, B-C) ​​or (A-C, B-N).

## References
