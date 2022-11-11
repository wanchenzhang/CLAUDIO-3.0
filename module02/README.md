# Module 02 - Structural distance analysis tool

## Usage

### The CLI - Command Line Interface
**Input**

This tool requires a csv-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked lysin's position within the full sequence (alternative: fill these with
Nans, but add two columns with the lysin's positions in the respective peptides).

**Output**

This tool will return a directory full of pdb structure files, which were necessary for the computation, a csv-file 
ending with the extension ".sqcs.csv", which is the result of the uniprot (and needed structure) search including the 
structural distance computation's results, a ".log"-file, that lists eventual issues or errors in the 
inputfile that had to be changed/fixed for the computations, two histogram png-files depicting the computed 
structural distances and their differences, and two more which depict the experimental methods and their resolution for
each pdb structure.
```
> python3 claudio_structdi.py [-i <filepath>] [-p <"projection_dict">] [-rt <True/False>] [-t <"blastp"/"hhsearch">] [-pc <float>] [-e <float>] [-qi <float>] [-c <float>] [-r <float>] [-o <directorypath>]
-i, --input-filepath,   path to inputfile,
                        default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p, --projections,      string which can be parsed as dictionary, containing the column names for the uniprot entry
                        columns (for naming convention see second example (Note: all values are mandatory, only change 
                        the keys accordingly)),
                        default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant)
-rt,--read-temps,       if the tool has been run before with the same input a temporary file was saved, which can be 
                        used to retireve the previouse results of the uniprot and structure searches, default=False
-t, --search-tool,      can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb search,
                        default="blastp"
-pc,--plddt-cutoff,     float value used as cutoff for alphafold structure prediction confidences (plddt), default=70.0
-e, --e-value,          e-value used in structure search, default=1e-5
-qi,--query-id,         query identity used in structure search, default=90.0
-c, --coverage,         coverage used in structure search, default=50.0
-r, --res-cutoff,       float value used as cutoff in angstrom for resolution of structure files, default=6.5
-o, --output-directory, output directory, default="data/out/module02"
```

### Example
The project can be run like this:
```
python3 claudio_structdi.py
```
This will result in a directory called 'structures' in "data/out/module02" filled with pdb structures, two 
histograms describing their resolutions and experimental methods, and one csv-file with the associated results of the
structural distance analysis, all pertaining the project's default dataset 
'liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'.
```
python3 claudio_structdi.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id_a', 'my_entry2': 'unip_id_b', "peptide1": "pep_a", "peptide2": "pep_b", "position1": "pos_a", "position2": "pos_b", "k_pos1": "k_pos_a", "k_pos2": "k_pos_b", "gene1": "gene_a", "gene2": "gene_b", "Publication": "pub"}" -rt True -o "c/user/documents/outs"
```
This will result in the respective outputs into the directory "c/user/documents/outs" for the placeholder 
'cross_links.csv' using the temporary save files from the last execution (if this has not been run previously with the 
same input, this will return an error).

## Authors

* **Alexander Röhl**
