# Module 02 - Structural distance analysis tool

### The CLI - Command Line Interface
```
> python3 claudio_structdi.py [-i <filepath>] [-rt <True/False>] [-t <"blastp"/"hhsearch">] [-pc <float>] [-lmin <float>] [-lmax <float>] [-e <float>] [-qi <float>] [-c <float>] [-r <float>] [-o <directorypath>] [-bl <directorypath>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-tl <directorypath>] [-v <int>]

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-rt,   --read-temps,            if the tool has been run before with the same input a temporary file was saved, which 
                                can be used to retireve the previouse results of the uniprot and structure searches,
                                default=False
-t,    --search-tool,           can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-x,    --xl-residues,           comma-separated one-letter-code residues, optional: add two ':' after the 
                                one-letter-code symbol of the residue in order to specify full sequence position 
                                (either 1 for start, or -1 for end position) and/or the atom used for the distance
                                computation, default="K,M:N:1"
-pc,   --plddt-cutoff,          float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0
-lmin, --linker-minimum,        float value used as minimal crosslinker range in angstrom, default=5.0
-lmax, --linker-maximum,        float value used as maximal crosslinker range in angstrom, default=35.0
-e,    --e-value,               e-value used in structure search, default=1e-5
-qi,   --query-id,              query identity used in structure search, default=90.0
-c,    --coverage,              coverage used in structure search, default=50.0
-r,    --res-cutoff,            float value used as cutoff in angstrom for resolution of structure files, default=6.5
-o,    --output-directory,      output directory, default="data/out/module02"
-bl,   --blast-bin,             binary directory in blast installation, or None if binary directory has been added to 
                                PATH variable (e.g. if blast can be called from anywhere), default=None
-bldb, --blast-db,              database directory for blast installation, default="$BLASTDB"
-hh,   --hhsearch-bin,          binary directory in hh-suite installation, or None if binary directory has been added to
                                PATH variable (e.g. if hhsearch can be called from anywhere), default=None
-hhdb, --hhsearch-db,           database directory for hh-suite installation, default="$HHDB"
-tl,   --topolink-bin,          binary directory in topolink installation, or None if binary directory has been added to
                                PATH variable (e.g. if topolink can be called from anywhere), default=None
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
This tool will return a directory full of pdb structure files, which were necessary for the computation, a CSV-file 
ending with the extension '.sqcs_structdi.csv', which is the result of the structure search including the 
structural distance computation's results, two histogram PNG-files depicting the computed structural distances and their
differences, and two more which depict the experimental methods and their resolution for each pdb structure.

### Example
The project can be run like this:
```
python3 claudio_structdi.py
```
This will result in a directory called 'structures' in "data/out/module02" filled with pdb structures, two 
histograms describing their resolutions and experimental methods, and one CSV-file with the associated results of the
structural distance analysis, all pertaining the project's default dataset 
['liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'](https://github.com/KohlbacherLab/CLAUDIO/blob/main/data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv).
```
python3 claudio_structdi.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id', 'peptide1': 'pep_a', 'peptide2': 'pep_b', 'position1': 'pos_a', 'position2': 'pos_b', 'k_pos1': 'res_pos_a', 'k_pos2': 'res_pos_b'}" -rt True -o "c/user/documents/outs"
```
This will result in the respective outputs into the directory "c/user/documents/outs" for the placeholder 
'cross_links.csv' using the temporary save files from the last execution (if this has not been run previously with the 
same input, this will return an error).

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
