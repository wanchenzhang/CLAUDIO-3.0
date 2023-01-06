# Module 01 - Unique protein (pair) listing tool

### The CLI - Command Line Interface
```
> python3 claudio_lists.py [-i <filepath>] [-p <"projection_dict">] [-s <True/False>] [-t <"blastp"/"hhsearch">] [-o <directorypath>] [-bl <directorypath>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-hhout <directorypath>] 

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,    --projections,           comma-separated list used to map the input dataset column names to the ones used by the 
                                tool (for naming convention see second example (Note: all keys (values before the 
                                colons) are mandatory, and have to be mapped correctly onto the respective columns in 
                                the input dataset). The order does not matter here.),
                                default="pep_a:peptide1,pep_b:peptide2,pos_a:position1,pos_b:position2,res_pos_a:k_pos1,res_pos_b:k_pos2,unip_id_a:entry1,unip_id_b:entry2"
-s,    --uniprot_search,        if the tool has been run before with the same input a temporary file was saved, which 
                                can be used to retrieve the previouse results of the uniprot and structure searches. In 
                                this case you can set uniprot_search=False and it will try to retrieve said temporary 
                                save file, instead of rerunning the uniprot sequence search, default=False
-x,    --xl-residues,           commaseperated one-letter-code residues, optional: add ':' after the one-letter-code 
                                symbol of the residue in order to specify full sequence position (either 1 for start, or 
                                -1 for end position), default="K,M:1"
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
This tool requires a CSV-file containing multiple observed cross-linking interactions. One column has to contain the 
uniprot id associated with the crosslinked protein, two columns have to contain the observed peptides for each 
interacting residue and two columns have to contain the crosslinked residue's position within the full sequence 
(alternative: fill these with Nans, but add two columns with the crosslinked residue position in the respective 
peptide).\ It is recommended to start the analysis of your cross-linking dataset here, if you do not run the full 
pipeline. Most important for this first step is the customization of the parameter "-p / --projections". This parameter 
requires a python dictionary as input, which maps the column names of your dataset to the ones used in the tool. You may
see examples for this in [this module](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/src/dict/default_projections.py)
or in the [example configuration-file](https://github.com/KohlbacherLab/CLAUDIO/blob/main/config.txt). You can also add
your own projection to either [this file](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/src/dict/default_projections.py)
or to your own configuration file.

### Output
This tool will return two CSV-files and one LOG-file. One contains a list of all unique proteins, their full uniprot 
sequences, additional information retrieved from uniprot, their pdb id and the number of occurences in the inputfile. 
The second CSV-file will have the extension '.sqcs', which marks it as the input for module02 or module03. The LOG-file 
lists eventual issues or errors in cross-linked residue positions in the inputfile that had to be changed/fixed for the 
output.

### Example
The project can be run like this:
```
python3 claudio_lists.py
```
This will result in two CSV-files in "data/out/unique_protein_list" pertaining the project's default dataset 
['liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'](https://github.com/KohlbacherLab/CLAUDIO/blob/main/data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv).
```
python3 claudio_lists.py -i "c/user/documents/cross_links.csv -p "{'my_entry1': 'unip_id', 'peptide1': 'pep_a', 'peptide2': 'pep_b', 'position1': 'pos_a', 'position2': 'pos_b', 'k_pos1': 'res_pos_a', 'k_pos2': 'res_pos_b'}" -t "hhsearch" -o "c/user/documents/outs"
```
This will result in two csv-files in "c/user/documents/outs" for the placeholder 'cross_links.csv'.

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
