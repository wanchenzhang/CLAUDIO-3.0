# Module 03 - Overlapping peptide sequence analysis tool

### The CLI - Command Line Interface
```
> python3 claudio_ops.py [-i <filepath>] [-p <"projection_dict">] [-s <True/False>] [-o <directorypath>] 

-i,   --input-filepath,         path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,   --projections,            string which can be parsed as dictionary, containing the column names for the uniprot 
                                entry columns (for naming convention see second example (Note: all values are mandatory,
                                only change the keys accordingly)),
                                default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant)
-s,   --uniprot-search,         if the tool has been run before with the same input a temporary file was saved, which
                                can be used to skip some of the steps, otherwise perform a full uniprot search again
                                (opposite of --read-temps), default=True
-o,   --output-directory,       output directory for produced csv-files, default="data/out/homo_signal/"
```

###Input
This tool requires a csv-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked lysin's position within the full sequence (alternative: fill these with
Nans, but add two columns with the lysin's positions in the respective peptides).

###Output
This tool will return a single csv-file containing the results of the ops analysis tool as well as three histograms
showing the results in detail.

### Example
The project can be run like this:
```
python3 claudio_ops.py
```
This will result in a csv-file in "data/out/homo_signal" containing the full dataset with the results of the OPS 
analysis tool, and three histograms depicting the analysis' results as distributions, all pertaining the project's 
default dataset 'liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'.
```
python3 claudio_ops.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id_a', 'my_entry2': 'unip_id_b', "peptide1": "pep_a", "peptide2": "pep_b", "position1": "pos_a", "position2": "pos_b", "k_pos1": "k_pos_a", "k_pos2": "k_pos_b", "gene1": "gene_a", "gene2": "gene_b", "Publication": "pub"}" -s True -o "c/user/documents/outs"
```
This will result in the respective outputs into the directory "c/user/documents/outs" for the placeholder 
'cross_links.csv' with a fresh uniprot search for the full protein sequences.

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
