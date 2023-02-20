# Module 03 - Overlapping peptide sequence analysis tool

### The CLI - Command Line Interface
```
> python3 claudio_ops.py [-i <filepath>] [-o <directorypath>] [-v <int>] 

-i,   --input-filepath,         path to inputfile,
                                default="data/out/unique_protein_list/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs"
-o,   --output-directory,       output directory for produced csv-files, default="data/out/homo_signal/"
-v,    --verbose-level,         verbose level value, default=3:
                                    0: no outputs at all will be written to the commandline
                                    1: write tool inits and passed time
                                    2: write progressbars (where implemented)
                                    3: write alignments during data processing, and write extra information on process results
                                    4: write alignments during uniprot to pdb position translation
                                    5: write verfications during uniprot to pdb position translation
```

### Input
This tool requires a CSV-file containing multiple observed cross-linking interactions. Two columns have to contain 
UniProt IDs for each interacting residue, two columns the observed peptides for each interacting residue,
two columns the crosslinked residue's position within the full sequence (alternative: fill these with
Nans, but add two columns with the residue's positions in the respective peptides), two columns with the UniProt sequences,
and one column with the initial cross-link type estimation.\
It is recommended to start the analysis of your dataset with [module01](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module01),
which will generate said CSV-file with the extension '.sqcs'. If you wish to skip this step though, or if you customized
the output file used as input here, check this list to ensure that you can run this module:
* file has to have the extension '.sqcs'
* file is formatted like a normal CSV-file
* file has columns "unip_id", "pep_a", "pep_b", "seq", "pos_a", "pos_b", "res_pos_a", "res_pos_b", 
each filled with information (in that order) of UniProt IDs [strings], peptides [strings], UniProt sequences [strings], 
crosslinked residue positions in full UniProt sequences [integers], and (alternative to normal positions) positions of 
the crosslinked residue in the peptides [integers].

### Output
This tool will return a single CSV-file containing the results of the OPS-analysis tool and ending with the extension 
'.sqcs_ops.csv', as well as three histograms showing the results in detail, e.g. the one depicting the adjacencies of 
cross-linked residues, the relative overlaps of peptides, as well as a simple TRUE/FALSE-distribution of the dataset for
overlapping peptide sequences.

### Example
The project can be run like this:
```
python3 claudio_ops.py
```
This will result in a CSV-file in "data/out/homo_signal" containing the full dataset with the results of the OPS 
analysis tool, and three histograms depicting the analysis' results as distributions, all pertaining the project's 
default dataset ['liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'](https://github.com/KohlbacherLab/CLAUDIO/blob/main/data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv).
```
python3 claudio_ops.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id', 'peptide1': 'pep_a', 'peptide2': 'pep_b', 'position1': 'pos_a', 'position2': 'pos_b', 'k_pos1': 'res_pos_a', 'k_pos2': 'res_pos_b'}" -s True -o "c/user/documents/outs"
```
This will result in the respective outputs into the directory "c/user/documents/outs" for the placeholder 
'cross_links.csv' with a fresh uniprot search for the full protein sequences.

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
