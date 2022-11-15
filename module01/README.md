# Module 01 - Unique protein (pair) listing tool

### The CLI - Command Line Interface
```
> python3 claudio_lists.py [-i <filepath>] [-p <"projection_dict">] [-t <"blastp"/"hhsearch">] [-o <directorypath>] [-bl <directorypath>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-hhout <directorypath>] 

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,    --projections,           string which can be parsed as dictionary, containing the column names for the uniprot 
                                entry columns (for naming convention see second example (Note: "unip_id_a" and 
                                "unip_id_b" are mandatory)), 
                                default=str(liu18_schweppe17_linked_residues_intra_homo_2370_nonredundant_unique)
-t,    --search-tool,           can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-o,    --output-directory,      output directory for produced csv-files, default="data/out/unique_protein_list"
-bl,   --blast-bin,             binary directory in blast installation, or None if binary directory has been added to 
                                PATH variable (e.g. if blast can be called from anywhere), default=None
-bldb, --blast-db,              database directory for blast installation, default="$BLASTDB"
-hh,   --hhsearch-bin,          binary directory in hh-suite installation, or None if binary directory has been added to
                                PATH variable (e.g. if hhsearch can be called from anywhere), default=None
-hhdb, --hhsearch-db,           database directory for hh-suite installation, default="$HHDB"
-hhout,--hhsearch-out,          output directory for hhsearch results, default="$HHOUT"
```
### Input
This tool requires a csv-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked lysin's position within the full sequence (alternative: fill these with
Nans, but add two columns with the lysin's positions in the respective peptides).

### Output
This tool will return two csv-files. One contains a list of all unique proteins, their full uniprot sequences,
additional information retrieved from uniprot, their pdb id and the number of occurences in the inputfile, while the 
other contains a list of all unique interaction pairs, their pdb ids and the number of occurences of this pair in the 
inputfile.

### Example
The project can be run like this:
```
python3 claudio_lists.py
```
This will result in two csv-files in "data/out/unique_protein_list" pertaining the project's default dataset 
'liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'.
```
python3 claudio_lists.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id_a', 'my_entry2': 'unip_id_b'}" -t "hhsearch" -o "c/user/documents/outs"
```
This will result in two csv-files in "c/user/documents/outs" for the placeholder 'cross_links.csv'.

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
