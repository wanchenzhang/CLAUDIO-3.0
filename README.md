# CLAUDIO 

*CLAUDIO*, the Cross-linking data analysis tool utilizing structural distances and overlapping peptide sequences, allows
for a methodical stepwise evaluation of cross-linking interaction types via in-depth analysis of structure and sequence 
information. It returns structural restraints, which can be applied in structure predictions. 

## Prerequisites

This tool is written in python and has the following python dependencies:

python 3.6

Packages:
* biopython
* click
* matplotlib
* pandas
* requests

This tool furthermore utilizes the following external softwares:
* *blast*
* *hhsearch*
* *topolink*

To ensure a flawless execution of *CLAUDIO* ensure that these tools can be executed directly from your terminal 
(e.g. add their respective `bin` directories to the `Path` variable of your OS).

For *blast* you have to download the `pdbaa` database, if it wasn't with the normal installation, and add 
the variable `$BLASTDB` (leading to the directory containing all database files) to the `Path` variable as a 
shortcut.

For *hhsearch* you have to download the `pdb70` database, and add both `$HHDB` (leading to the directory 
containing all database files) and `$HHOUT` (leading to the directory for hhsearch outputs) to the `Path` variable as a 
shortcut.

For *topolink* no additional steps have to be taken.

## Usage

### The CLI - Command Line Interface
#### CLAUDIO - Full pipeline
**Input**

This tool requires a csv-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked lysin's position within the full sequence (alternative: fill these with
Nans, but add two columns with the lysin's positions in the respective peptides).

All parameters can be given in a configuration file (see example: data/in/iiCheck_config.txt).

**Output**

This tool returns practically all the outputs listed in the other modules (see respective README.md-files).
Note: All csv-file outputs pertaining the input dataset are summarized into a single one.
```
> python3 claudio.py [-i <filepath>] [-p <"projection_dict">] [-rt <True/False>] [-t <"blastp"/"hhsearch">] [-e <float>] [-qi <float>] [-cv <float>] [-r <float>] [-pc <float>] [-lmin <float>] [-lmax <float>] [-es <float>] [-dm <float] [-ct <float>] [-o <directorypath>] [-c <filepath>] 

-i,   --input-filepath,         path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,   --projections,            string which can be parsed as dictionary, containing the column names for the uniprot 
                                entry columns (for naming convention see second example (Note: all values are mandatory,
                                only change the keys accordingly)),
                                default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant)
-rt,  --read-temps,             if the tool has been run before with the same input a temporary file was saved, which
                                can be used to skip some of the steps, default=False
-t,   --search-tool,            can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-e,   --e-value,                e-value used in structure search, default=1e-5
-qi,  --query-id,               query identity used in structure search, default=90.0
-cv,  --coverage,               coverage used in structure search, default=50.0
-r,   --res-cutoff,             float value used as cutoff in angstrom for resolution of structure files, default=6.5
-pc,  --plddt-cutoff,           float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0  
-lmin,--linker-minimum,         float value used as minimal crosslinker range in angstrom, default=0.0
-lmax,--linker-maximum,         float value used as maximal crosslinker range in angstrom, default=35.0
-es,  --euclidean-strictness,   float value substracted from the linker ranges for the euclidean distance scoring
                                (minimum will not go below 0), default=5.0
-dm,  --distance-maximum,       maximal distance value that seems realistic, if surpassed the distance will be set to 
                                this value during the confidence scoring, to ensure its consistency, default=50.0
-ct,  --cutoff,                 float value used as confidence score cutoff, if surpassed, the linker type will be set 
                                to inter, default=0.0  
-o,   --output-directory,       output directory for produced csv-files, default="data/out/full"
-c,   --config,                 filepath to configuration file containing all input parameters, if given all other 
                                parameters will be ignored (see example: data/in/iiCheck_config.txt), default=''
```

### Example
**CLAUDIO** (in full) can be run like this:
```
python3 claudio.py -c data/in/iiCheck_config.txt
```

## Authors

* **Alexander Röhl, Hadeer Elhabashy, Eugen Netz**
