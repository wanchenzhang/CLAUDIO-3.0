# CLAUDIO 

*CLAUDIO*, the tool for "**C**ross-**l**inking **a**nalysis **u**sing **di**stances and **o**verlaps", allows
for a methodical stepwise evaluation of cross-linking interaction types via in-depth analysis of structure and sequence 
information. It returns structural restraints, which can be applied in structure predictions, and the input dataset
extended by its analysis' results. 

## Prerequisites
### Python
This tool is written in python (v3.11) and has thus to be run with python 3.\
It has the following requirements:

Packages:
* biopython 1.79
* click 8.1.3
* matplotlib 3.6.3
* pandas 1.5.3
* requests 2.28.2

The packages may be installed all at once with the file [requirements.txt](https://github.com/KohlbacherLab/CLAUDIO/blob/main/requirements.txt):
```
pip install -r requirements.txt
```
or individually:
```
pip install biopython==1.79
pip install click==8.1.3
pip install matplotlib==3.6.3
pip install pandas==1.5.3
pip install requests==2.28.2
```
Note: Both approaches need to refer to the pip-installer associated to the python installation, that will be used to run
the tool.

### External Tools
In order to run *CLAUDIO* you need to install the following external tools:
* **Topolink** (for structural analysis)
* **BLASTP** and/or **HHsearch** (for finding suitable protein structures (recommended: BLASTP))

####Links to Installation Manuals
* *Blast* (see [Installation Manual](https://www.ncbi.nlm.nih.gov/books/NBK52640/))
  * For Blast, you also have to download the newest `pdbaa` database. You may do so by navigating into your 
  Blast installation directory and running the following command:
    ```
    perl bin/update_blastdb.pl --passive --decompress pdbaa
    ```
    Hint: If you followed the Installation Manual (linked up top) to-the-letter, you may have added the environmental 
    variable `$BLASTDB` to your paths. If so, you can use this variable instead of the full path in the input parameters
    of CLAUDIO.
* *HHsearch* (see [GitHub](https://github.com/soedinglab/hh-suite))
  * For HHsearch, you also have to download its newest `pdb70` database found 
  [here](https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/). 
* *TopoLink* (see [Installation Manual](http://leandro.iqm.unicamp.br/topolink/download.shtml))

To ease the execution of *CLAUDIO* you may want to ensure that these tools can be executed directly from your terminal 
(e.g. add their respective `bin` directories to the `Path` variable of your OS). Otherwise, you must specify the
location of the `bin` directories in the input parameters.

### Online connection
*CLAUDIO* calls upon the API of a number of bioinformatic online databases ([UniProt](https://www.uniprot.org/), 
[RCSB](https://www.rcsb.org/), [AlphaFold](https://alphafold.ebi.ac.uk/), and 
[SWISS-MODEL](https://swissmodel.expasy.org/)) during its computations. This means it cannot be run offline.\
It is furthermore recommended having a stable internet connection, as otherwise certain API calls may not be answered or
lead to empty results. This of course, also necessitates the database server's side to be running properly as well. If 
errors or suspicious inconsistencies in the results persist due to this, you may want to try again later.

## Usage
*CLAUDIO* consists of a total of 4 modules. Each module can be run independently as long as appropriate inputs are 
delivered. For details on how to run the modules individually see their respective README-files.
* [Module 01 - Unique protein (pair) listing tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/README.md)
* [Module 02 - Structural distance analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module02/README.md)
* [Module 03 - Overlapping peptide sequence analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module03/README.md)
* [Module 04 - XL-type evaluation](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module04/README.md)

For details on how to run the **full** pipeline continue below.

---
---

## CLAUDIO - Full pipeline
### The CLI - Command Line Interface
```
> python3 claudio.py [-i <filepath>] [-p <"projection_dict">] [-rt <True/False>] [-t <"blastp"/"hhsearch">] [-e <float>] [-qi <float>] [-cv <float>] [-r <float>] [-pc <float>] [-lmin <float>] [-lmax <float>] [-es <float>] [-dm <float] [-ct <float>] [-o <directorypath>] [-bl <directorypath>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-hhout <directorypath>] [-tl <directorypath>] [-c <filepath>] 

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,    --projections,           comma-separated position-sensitive list that names the column names of the users dataset
                                containing the necessary information for the tool. The column names should contain and 
                                should be given in the following order: crosslinked peptide_a, crosslinked peptide_b, 
                                crosslinked residue position_a, crosslinked residue position_b, position of crosslinked 
                                residue in peptide_a, position of crosslinked residue in peptide_b, UniProt ID of 
                                protein belonging to peptide_a, UniProt ID of protein belonging to peptide_b.
                                Note: The positions of the crosslinked residue in the peptides are information only 
                                accessed, if the given full sequence positions do not match into the retrieved UniProt 
                                sequence. If the positions are confirmed you may simply create two substitute columns 
                                for the positions in the peptides instead and leave them empty.
                                default="peptide1,peptide2,position1,position2,k_pos1,k_pos2,entry1,entry2"
-rt,   --read-temps,            if the tool has been run before with the same input a temporary file was saved, which
                                can be used to skip some of the steps, default=False
-x,    --xl-residues,           comma-separated one-letter-code residues, optional: add ':' after the one-letter-code 
                                symbol of the residue in order to specify full sequence position (either 1 for start, or 
                                -1 for end position), default="K,M:1"
-t,    --search-tool,           can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-e,    --e-value,               e-value used in structure search, default=1e-5
-qi,   --query-id,              query identity used in structure search, default=90.0
-cv,   --coverage,              coverage used in structure search, default=50.0
-r,    --res-cutoff,            float value used as cutoff in angstrom for resolution of structure files, default=6.5
-pc,   --plddt-cutoff,          float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0  
-lmin, --linker-minimum,        float value used as minimal crosslinker range in angstrom, default=0.0
-lmax, --linker-maximum,        float value used as maximal crosslinker range in angstrom, default=35.0
-es,   --euclidean-strictness,  float value substracted from the linker ranges for the euclidean distance scoring
                                (minimum will not go below 0), default=5.0
-dm,   --distance-maximum,      maximal distance value that seems realistic, if surpassed the distance will be set to 
                                this value during the confidence scoring, to ensure its consistency, default=50.0
-ct,   --cutoff,                float value used as confidence score cutoff, if surpassed, the linker type will be set 
                                to inter, default=0.0  
-o,    --output-directory,      output directory for produced csv-files, default="data/out/full"
-bl,   --blast-bin,             binary directory in blast installation, or None if binary directory has been added to 
                                PATH variable (e.g. if blast can be called from anywhere), default=None
-bldb, --blast-db,              database directory for blast installation, default="$BLASTDB"
-hh,   --hhsearch-bin,          binary directory in hh-suite installation, or None if binary directory has been added to
                                PATH variable (e.g. if hhsearch can be called from anywhere), default=None
-hhdb, --hhsearch-db,           database directory for hh-suite installation, default="$HHDB"
-hhout,--hhsearch-out,          output directory for hhsearch results, default="$HHOUT"
-tl,   --topolink-bin,          binary directory in topolink installation, or None if binary directory has been added to
                                PATH variable (e.g. if topolink can be called from anywhere), default=None)
-c,    --config,                filepath to configuration file containing all input parameters, if given all other 
                                parameters will be ignored (see example: config.txt), default=''
```
### Input
This tool requires a CSV-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked residue's position within the full sequence (alternative: fill these 
with Nans, but add two columns with the residue's positions in the respective peptides).\
First, it is important to customize the parameter "-p / --projections". This parameter requires a python dictionary as 
input, which maps the column names of your dataset to the ones used in the tool. You may see examples for this in 
[this module](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/src/dict/default_projections.py) or in the 
[example configuration-file](https://github.com/KohlbacherLab/CLAUDIO/blob/main/config.txt).

All parameters can be given in a configuration file (see example: [config.txt](https://github.com/KohlbacherLab/CLAUDIO/blob/main/config.txt)).

### Output
This tool returns all the outputs listed in the modules (see 
[module01](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module01),
[module02](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module02),
[module03](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module03),
[module04](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module04)).
Note: All CSV-file outputs pertaining the input dataset are summarized into a single one (marked with 
'_final.csv'-extension), e.g. the output CSV-file of module01 ending with '.sqcs', of module02 ending with 
'.sqcs_structdi.csv', and of module03 ending with '.sqcs_ops.csv' will be summarized here. 

### Example
**CLAUDIO** (in full) can be run like this:
```
python3 claudio.py -c config.txt
python3 claudio.py -i /home/user/docs/xl_dataset.csv -o /home/user/docs/claudio_outputs
```

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
